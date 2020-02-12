# DESCRIPTION =========================================================
# Author: Ashesh Rambachan
#
#  This script contains functions to implement the methods
#  described in Rambachan & Roth (2019) for robust inference
#  in difference-in-differences and event study designs.
#
#   First, this script contains function that are used to construct the
#   FLCI for a general choice of the vector l and Delta = Delta^{SD}(M).
#
#   Second, this script defines functions used to construct conditional
#   confidence intervals and hybrid confidence intervals for a general choice of the
#   vector l and Delta = Delta^{SD}, Delta^{SDPB}, Delta^{SDNB},
#   Delta^{SDI}, Delta^{SDD}. This implements functions to deal with
#   nuisance parameters as described in  Andrews, Roth, Pakes (2019).
#
# PRELIMINARIES =======================================================
library(tidyverse)
library(TruncatedNormal)
library(lpSolveAPI)
library(ROI)
library(Matrix)
library(pracma)
library(CVXR)
library(foreach)

# FLCI HELPER FUNCTIONS -----------------------------------------------
.createConstraints_AbsoluteValue <- function(sigma, numPrePeriods, UstackW){
  # This function creates linear constraints that help to minimize worst-case bias
  # over w using the optiSolve package, where we cast maximization of the absolute values
  # as a linear problem.
  # To do this, we assume that the optimization is over a vector (U, W) where
  # U_1 = |w_1|, U_2 = |w_1 + w_2|, etc.
  # This function creates a matrix of linear constraints that imposes these equalities
  # and can easily be passed to optimSolve.
  # To do this, we require that
  # -U_1 + w_1 <= 0
  # -U_1 - w_1 <= 0, and so on.
  # We return the values to be passed to optimSolve.
  K <- numPrePeriods
  lowerTriMat <- diag(K)
  lowerTriMat[lower.tri(lowerTriMat)] <- 1
  A_absolutevalue <-
    rbind(
      cbind( -diag(K), lowerTriMat  ),
      cbind( -diag(K), -lowerTriMat )
    )
  threshold_absolutevalue <- rep(0, NROW(A_absolutevalue))
  direction_absolutevalue <- rep("<=", NROW(A_absolutevalue))
  constraint = A_absolutevalue %*% UstackW <= threshold_absolutevalue
  return(constraint)
}

.createConstraints_SumWeights <- function(numPrePeriods, l_vec, UstackW){
  # Creates constraints on the sum of the weights.
  numPostPeriods = length(l_vec)
  A_sumweights <- c(rep(0,numPrePeriods), rep(1,numPrePeriods))
  threshold_sumweights <- t((1:numPostPeriods)) %*% l_vec
  direction_sumweights <- "=="

  constraint = t(A_sumweights) %*% UstackW == threshold_sumweights
  return(constraint)
}

.createObjectiveObjectForBias <- function(numPrePeriods, numPostPeriods, l_vec, UstackW){
  # Constructs the objective function for the worst-case bias.
  constant = sum(sapply(1:numPostPeriods, FUN = function(s) { abs(t(1:s) %*% l_vec[(numPostPeriods - s + 1):numPostPeriods]) })) - t((1:numPostPeriods)) %*% l_vec
  objective.UstackW = Minimize( constant +  t(c(rep(1,numPrePeriods), rep(0,numPrePeriods))) %*% UstackW )
  return(objective.UstackW)
}

.createMatricesForVarianceFromW <- function(sigma, numPrePeriods, l_vec, UstackW,
                                            prePeriodIndices = 1:numPrePeriods){
  # Constructs matrix to compute the variance of the affine estimator for a choice of weights W.
  SigmaPre <- sigma[prePeriodIndices, prePeriodIndices]
  SigmaPrePost <- sigma[prePeriodIndices, -prePeriodIndices]
  SigmaPost <- t(l_vec) %*% sigma[-prePeriodIndices, -prePeriodIndices] %*% l_vec
  WtoLPreMat <- diag(numPrePeriods)
  if (numPrePeriods == 1) {
    WtoLPreMat = 1
  } else {
    for(col in 1:(numPrePeriods-1) ){
      WtoLPreMat[col+1, col] <- -1
    }
  }

  UstackWtoLPreMat <- cbind(matrix(0, nrow = numPrePeriods, ncol = numPrePeriods), WtoLPreMat)
  A_quadratic_sd <-  t(UstackWtoLPreMat) %*% SigmaPre %*% UstackWtoLPreMat
  A_linear_sd <- 2 * t(UstackWtoLPreMat) %*% SigmaPrePost %*% l_vec
  A_constant_sd <- SigmaPost

  return(list(A_quadratic_sd = A_quadratic_sd,
              A_linear_sd = A_linear_sd,
              A_constant_sd = A_constant_sd))
}

.createConstraintsObject_SDLessThanH <- function(sigma, numPrePeriods, l_vec, UstackW, h, ...){
  # Creates constraint that the variance of the affine estimator must be less than
  # some level h.
  A_matrices <- .createMatricesForVarianceFromW(sigma, numPrePeriods, l_vec, ...)
  threshold_sd_constraint <- h^2
  constraint = CVXR::quad_form(UstackW, A_matrices$A_quadratic_sd) + t(A_matrices$A_linear_sd) %*% UstackW + A_matrices$A_constant_sd <= threshold_sd_constraint
  return(constraint)
}

.createObjectiveObject_MinimizeSD <- function(sigma, numPrePeriods, numPostPeriods, UstackW, l_vec, ...){
  # Create objective function to minimize the standard deviation.
  A_matrices <- .createMatricesForVarianceFromW(sigma, numPrePeriods, l_vec, ...)
  objective = Minimize(CVXR::quad_form(UstackW, A_matrices$A_quadratic_sd) + t(A_matrices$A_linear_sd) %*% UstackW + A_matrices$A_constant_sd)
  return(objective)
}

.qfoldednormal <- function(p, mu = 0, sd = 1, numSims = 10^6, seed = 0){
  # Computes the pth quantile of the folded normal distribution with mean mu and sd = sd
  # Vectorized over mu
  set.seed(seed)
  normDraws <- rnorm(n = numSims, sd = sd)
  pQuantiles <- purrr::map_dbl(.x = mu, .f = function(mu){quantile(abs(normDraws + mu), probs = p)})
  return(pQuantiles)
}

.wToLFn <- function(w){
  # Converts vector from w space to l space.
  numPrePeriods <- length(w)
  WtoLPreMat <- diag(numPrePeriods)
  if (numPrePeriods == 1) {
    WtoLPreMat = 1
  } else {
    for(col in 1:(numPrePeriods-1) ){
      WtoLPreMat[col+1, col] <- -1
    }
  }
  l <- c(WtoLPreMat %*% w)
  return(l)
}

.lToWFn <- function(l_vec) {
  # Converts vector from l space to w space
  numPostPeriods = length(l_vec)
  lToWPostMat <- diag(numPostPeriods)
  if (numPostPeriods == 1) {
    lToWPostMat = 1
  } else {
    for(col in 1:(numPostPeriods-1) ){
      lToWPostMat[col+1, col] <- 1
    }
  }
  w = c(lToWPostMat %*% l_vec)
  return(w)
}

# FLCI FUNCTIONS ------------------------------------------------------
.findWorstCaseBiasGivenH <- function(h, sigma, numPrePeriods, numPostPeriods, l_vec, M = 1, returnDF = F) {
  # This function minimizes worst-case bias over Delta^{SD}(M) subject to the constraint
  # that the SD of the estimator is <= h
  # Note: this function assumes M = 1 unless specified otherwise.
  UstackW <- CVXR::Variable(numPrePeriods + numPrePeriods)
  objectiveBias <- .createObjectiveObjectForBias(numPrePeriods = numPrePeriods, UstackW = UstackW, numPostPeriods = numPostPeriods, l_vec = l_vec)

  abs_constraint <- .createConstraints_AbsoluteValue(sigma = sigma, numPrePeriods = numPrePeriods, UstackW = UstackW)
  sum_constraint <- .createConstraints_SumWeights(numPrePeriods = numPrePeriods, l_vec = l_vec, UstackW = UstackW)
  quad_constraint <- .createConstraintsObject_SDLessThanH(sigma = sigma, numPrePeriods = numPrePeriods, l_vec = l_vec, UstackW = UstackW, h = h)

  biasProblem = CVXR::Problem(objectiveBias, constraints = list(abs_constraint, sum_constraint, quad_constraint))
  biasResult <- psolve(biasProblem)

  # Multiply objective by M (note that solution otherwise doesn't depend on M,
  # so no need to run many times with many different Ms)
  biasResult$value <- biasResult$value * M

  # Compute the implied w and l
  optimal.x <- biasResult$getValue(UstackW)
  optimal.w <- optimal.x[ (length(optimal.x)/2 + 1):length(optimal.x) ]
  optimal.l <- .wToLFn(optimal.w)

  if(returnDF == F){
    return(list(status = biasResult$status,
                value = biasResult$value,
                optimal.x = optimal.x,
                optimal.w = optimal.w,
                optimal.l = optimal.l))
  } else{
    temp = list(status = biasResult$status,
                value = biasResult$value,
                optimal.x = I(list(unlist(optimal.x))),
                optimal.w = I(list(unlist(optimal.w))),
                optimal.l = I(list(unlist(optimal.l)))
    )
    return(as.data.frame(temp))
  }
}

.findLowestH <- function(sigma, numPrePeriods, numPostPeriods, l_vec){
  # Finds the minimum variance affine estimator.
  UstackW <- Variable(numPrePeriods + numPrePeriods)
  abs_constraint <- .createConstraints_AbsoluteValue(sigma = sigma, numPrePeriods = numPrePeriods, UstackW = UstackW)
  sum_constraint <- .createConstraints_SumWeights(numPrePeriods = numPrePeriods, l_vec = l_vec, UstackW = UstackW)
  objectiveVariance = .createObjectiveObject_MinimizeSD(sigma = sigma, numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, UstackW = UstackW, l_vec = l_vec)
  varProblem = Problem(objectiveVariance, constraints = list(abs_constraint, sum_constraint))
  varResult <- psolve(varProblem)

  if(varResult$status != "optimal" & varResult$status != "optimal_inaccurate"){
    warning("Error in optimization for h0")
  }

  minimalVariance <- varResult$value
  minimalH <- sqrt(minimalVariance)
  return(minimalH)
}

.computeSigmaLFromW <- function(w, sigma, numPrePeriods, numPostPeriods, l_vec, ...){
  # Compute variance of affine estmiator with choice w
  A_matrices <- .createMatricesForVarianceFromW(sigma = sigma,
                                                numPrePeriods = numPrePeriods,
                                                l_vec = l_vec, ...)
  UstackW <- matrix(c( rep(0, length(w)), w ), ncol = 1)
  varL <- t(UstackW) %*% A_matrices$A_quadratic_sd %*% UstackW +
    t(A_matrices$A_linear_sd) %*% UstackW + A_matrices$A_constant_sd
  return(varL)
}

.findHForMinimumBias <- function(sigma, numPrePeriods, numPostPeriods, l_vec){
  hsquared <- .computeSigmaLFromW(c(rep(0,numPrePeriods-1), t(1:numPostPeriods) %*% l_vec),
                                  sigma = sigma,
                                  numPrePeriods = numPrePeriods,
                                  numPostPeriods = numPostPeriods,
                                  l_vec = l_vec)
  h <- sqrt(hsquared)
  return(h)
}

.findOptimalFLCI_helper <- function(sigma, M,
                                    numPrePeriods, numPostPeriods,
                                    l_vec, numPoints = 100, alpha){
  h0 <- .findHForMinimumBias(sigma = sigma, numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, l_vec = l_vec)
  hMin <- .findLowestH(sigma = sigma, numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, l_vec = l_vec)
  hGrid <- seq(from = hMin, to = h0, length.out = numPoints)
  biasDF <- purrr::map_dfr(.x = hGrid,
                           .f = function(h){ .findWorstCaseBiasGivenH(h, sigma = sigma, numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, l_vec = l_vec, returnDF = T) %>% mutate(h = h) } )
  biasDF <- biasDF %>% rename(bias = value)
  biasDF <-
    left_join(
      biasDF %>% mutate(id = 1),
      data.frame(m = M, id = 1),
      by = "id") %>% dplyr::select(-id)

  biasDF <- biasDF %>%
    rename(maxBias = bias) %>% filter(maxBias < Inf)
  biasDF <- biasDF %>%
    mutate(maxBias = maxBias * m) %>%
    mutate(CI.halflength = .qfoldednormal(p = 1-alpha, mu = maxBias/h) * h)

  optimalCIDF <- biasDF %>%
    group_by(m) %>%
    filter(status == "optimal" | status == "optimal_inaccurate") %>%
    filter(CI.halflength == min(CI.halflength))

  results = list(
    optimalVec = c(unlist(optimalCIDF$optimal.l), l_vec),
    optimalPrePeriodVec = unlist(optimalCIDF$optimal.l),
    optimalHalfLength = optimalCIDF$CI.halflength,
    M = optimalCIDF$m,
    status = optimalCIDF$status
  )
  return(results)
}

findOptimalFLCI <- function(betahat, sigma, M = 0,
                            numPrePeriods, numPostPeriods,
                            l_vec = .basisVector(index = 1, size = numPostPeriods),
                            numPoints = 100, alpha = 0.05) {
  FLCI_Results = .findOptimalFLCI_helper(sigma = sigma,
                                         M = M,
                                         numPrePeriods = numPrePeriods,
                                         numPostPeriods = numPostPeriods,
                                         l_vec = l_vec,
                                         numPoints = numPoints,
                                         alpha = alpha)
  FLCI = c(t(FLCI_Results$optimalVec) %*% betahat - FLCI_Results$optimalHalfLength,
           t(FLCI_Results$optimalVec) %*% betahat + FLCI_Results$optimalHalfLength)

  Results = list(
    FLCI = FLCI,
    optimalVec = FLCI_Results$optimalVec,
    optimalHalfLength = FLCI_Results$optimalHalfLength,
    M = FLCI_Results$M,
    status = FLCI_Results$status
  )
  return(Results)
}

# ARP HELPER FUNCTIONS ------------------------------------------------
.norminvp_generalized <- function(p, l, u, mu = 0, sd = 1){
  lnormalized <- (l-mu)/sd
  unormalized <- (u-mu)/sd
  qnormalized <- norminvp(p, lnormalized, unormalized)
  q <- mu + qnormalized * sd
  return(q)
}

.basisVector <- function(index = 1, size = 1){
  v <- matrix(0, nrow = size, ncol = 1)
  v[index] = 1
  return(v)
}

.max_program <- function(s_T, gamma_tilde, sigma, W_T, c) {
  # Define objective and constraints
  f = s_T + c(t(gamma_tilde) %*% sigma %*% gamma_tilde)^(-1)*(sigma %*% gamma_tilde)*c
  Aeq = t(W_T)
  beq = c(1, rep(0, dim(Aeq)[1]-1))
  # Set up linear program. Note: don't need to provide lower bound because default lb is zero in ROI package.
  linprog = Rglpk::Rglpk_solve_LP(obj = -c(f),
                                  mat = Aeq,
                                  dir = rep("==", NROW(Aeq)),
                                  rhs = beq)
  # Solve linear program and return negative of objective because we want the max.
  return(-linprog$optimum)
}

.check_if_solution_helper <- function(c, tol, s_T, gamma_tilde, sigma, W_T) {
  min_value = .max_program(s_T, gamma_tilde, sigma, W_T, c)
  solution = (abs(c - min_value) <= tol)
  return(solution)
}

.vlo_vup_dual_fn <- function(eta, s_T, gamma_tilde, sigma, W_T) {
  # This function computes vlo and vup for the dual linear program for the
  # conditional values using the bisection method described in Appendix
  # G in Algorithm 1.
  #
  # Inputs:
  #   eta         = solution to LP from test_delta_lp_fn
  #   s_T         =
  #   gamma_tilde = dual solution to LP from test_delta_lp_fn, given in output lambda.
  #   sigma       = covariance matrix of y_T
  #   W_T         =
  #
  # Output: list with elements
  #   vlo
  #   vup

  # Options for bisection algorithm
  tol_c = 10^(-6)
  tol_equality = 10^(-6)
  sigma_B = sqrt( t(gamma_tilde) %*% sigma %*% gamma_tilde )
  low_initial = min(-100, eta-20*sigma_B)
  high_initial = max(100, eta + 20*sigma_B)
  maxiters = 10000

  ### Compute vup ###
  dif = tol_c+1
  iters = 0
  low = eta
  high = high_initial

  if (is.na(.check_if_solution_helper(c = eta, tol = tol_equality,
                                      s_T = s_T, gamma_tilde = gamma_tilde,
                                      sigma = sigma, W_T = W_T))) {
    # warning('User-supplied eta is not a solution. Not rejecting automatically')
    return( list(vlo = eta, vup = Inf) )
  }
  if (!.check_if_solution_helper(c = eta, tol = tol_equality,
                                 s_T = s_T, gamma_tilde = gamma_tilde,
                                 sigma = sigma, W_T = W_T)) {
    # warning('User-supplied eta is not a solution. Not rejecting automatically')
    return( list(vlo = eta, vup = Inf) )
  }

  if (.check_if_solution_helper(c = high_initial,
                                tol = tol_equality,
                                s_T = s_T,
                                gamma_tilde = gamma_tilde,
                                sigma = sigma, W_T = W_T)) {
    vup = Inf
  } else {
    # Throughout the while loop, the high value is not a solution and the
    # low-value is a solution. If the midpoint between them is a solution,
    # then we set low to the midpoint; otherwise, we set high to the midpoint.
    while ( dif > tol_c & iters < maxiters ) {
      iters = iters+1
      mid = 0.5*(high + low)
      if (.check_if_solution_helper(c = mid, tol = tol_equality,
                                    s_T = s_T, gamma_tilde = gamma_tilde,
                                    sigma = sigma, W_T = W_T)) {
        low = mid
      } else {
        high = mid
      }
      dif = high-low
    }
    if (iters == maxiters) {
      # warning("vup: Reached max iterations without a solution!")
    }
    vup = mid
  }

  ### Compute vlo ###
  dif = tol_c+1
  iters = 0
  low = low_initial
  high = eta

  if (.check_if_solution_helper(low_initial, tol = tol_equality,
                                s_T = s_T, gamma_tilde = gamma_tilde,
                                sigma = sigma, W_T = W_T)) {
    vlo = -Inf
  } else {
    # Throughout the while loop, the low value is not a solution and the
    # high-value is a solution. If the midpoint between them is a solution,
    # then we set high to the midpoint; otherwise, we set low to the midpoint.
    while (dif > tol_c & iters < maxiters) {
      mid = 0.5*(low + high)
      iters = iters+1
      if (.check_if_solution_helper(mid, tol = tol_equality,
                                    s_T = s_T, gamma_tilde = gamma_tilde,
                                    sigma = sigma, W_T = W_T)) {
        high = mid
      } else {
        low = mid
      }
      dif = high-low
    }
    if (iters == maxiters) {
      # warning("vlo: Reached max iterations without a solution!")
    }
    vlo = mid
  }
  return(list(vlo = vlo, vup = vup))
}

.test_delta_lp_fn <- function(y_T, X_T, sigma) {
  # Returns the value of eta that solves
  #   min_{eta, delta} eta
  #     s.t. y_T - x_T delta <= eta * diag(sigma)
  #
  # Inputs:
  #   y_T = vector
  #   X_T = matrix
  #   sigma = covariance matrix of y_T
  # Outputs: list with elements
  #   etaStar   = minimum of linear program
  #   deltaStar = minimizer of linear program
  #   lambda    = lagrange multipliers on primal i.e. solution to dual
  #   error_flag = whether the linear program was solved.
  #
  # Note:
  #   To solve the lin. program, we need to transform the
  #   linear program such that it is in the form
  #     max f' X
  #     s.t. C X \leq b

  dimDelta = dim(X_T)[2] # dimension of delta
  sdVec = sqrt(diag(sigma)) # standard deviation of y_T elements

  # Define objective function
  f = c(1, rep(0, dimDelta))

  # Define linear constraint
  C = -cbind(sdVec, X_T)
  b = -y_T

  # Define linear program using lpSolveAPI
  linprog = make.lp(nrow = 0, ncol = length(f))
  set.objfn(linprog, f)
  for (r in 1:dim(C)[1]) {
    add.constraint(linprog, xt = c(C[r,]), type = "<=", rhs = b[r])
  }
  set.bounds(linprog, lower = rep(-Inf, length(f)), columns = 1:length(f))

  # Solve linear program using dual simplex
  lp.control(linprog, sense = "min", simplextype="dual", pivoting = "dantzig", verbose="neutral")
  error_flag = solve(linprog)
  eta = get.objective(linprog)
  primalSoln = get.primal.solution(linprog)
  delta = primalSoln[(length(primalSoln)-dimDelta+1):length(primalSoln)]
  dual = -get.sensitivity.rhs(linprog)$duals[1:dim(C)[1]]

  # Construct results  list
  results = list(eta_star = eta,
                 delta_star = delta,
                 lambda = dual,
                 error_flag = error_flag)
  return(results)
}

.lp_dual_fn <- function(y_T, X_T, eta, gamma_tilde, sigma) {
  # Wrapper function to calculate vlo, vlup using the bisection approach
  #
  # Inputs:
  #   y_T = vector
  #   X_T = matrix
  #   eta = solution to LP from test_delta_lp_fn
  #   gamma_tilde = vertex of the dual, the output lambda from test_delta_lp_fn
  #   sigma = covariance matrix of y_T.

  sdVec = sqrt(diag(sigma))
  W_T = cbind(sdVec, X_T)
  s_T = (diag(length(y_T)) - c(t(gamma_tilde) %*% sigma %*% gamma_tilde)^(-1)*(sigma %*% (gamma_tilde %*% t(gamma_tilde)))) %*% y_T
  vList = .vlo_vup_dual_fn(eta = eta, s_T = s_T, gamma_tilde = gamma_tilde, sigma = sigma, W_T = W_T)
  return(list(vlo = vList$vlo, vup = vList$vup,
              eta = eta, gamma_tilde = gamma_tilde))
}

# Functon to return column index of leading ones of matrix
.leading_one <- function(r, B) {
  i = 1
  while (B[r,i] == 0) {
    i = i+1
  }
  return(i)
}

.construct_Gamma <- function(l) {
  # This function constructs the invertible matrix Gamma that
  # has the 1 x \bar{T} vector l' as its first row. To do so, it uses
  # the following linear algebra trick:
  #
  # Construct matrix B = [l e_1 ... e_{\bar{T}}]. Convert it to reduced
  # row echelon form and find the columns containing the pivots of
  # each row. The colums in the original matrix associated with the pivots
  # form a basis of R^{\bar{T}} and we use this to construct \Gamma.
  #
  # This function is used in the change of basis to transform Delta into
  # a MI of the form that can be used by the ARP method.

  barT = length(l) # length of vector
  B = cbind(l, diag(barT))
  rrefB = pracma::rref(B) # compute reduced row echelon form of B.

  # Construct gamma and check if invertible
  leading_ones <- map_dbl(.x = 1:dim(rrefB)[1], .leading_one, B = rrefB)
  Gamma = t(B[, leading_ones])
  if (det(Gamma) == 0) {
    stop('Something went wrong in RREF algorithm.')
  } else{
    return(Gamma)
  }
}

# Compute maximum bias of a linear estimator
.maxBiasFN <- function(v, A, d){
  delta <- Variable(NCOL(A))
  objective <- Maximize(t(v) %*% delta)
  problem <- Problem(objective, constraints = list(A %*% delta - d<= 0) )
  soln <- solve(problem)
  return(soln)
}

# Compute minimum bias of a linear estimator
.minBiasFN <- function(v, A, d){
  delta <- CVXR::Variable(NCOL(A))
  objective <- CVXR::Minimize(t(v) %*% delta)
  problem <- Problem(objective, constraints = list(A %*% delta - d<= 0) )
  soln <- solve(problem)
  return(soln)
}

#Create a vector that extrapolates a linear trend to pre-period
.createV_Linear_NoIntercept <- function(lVec, tVec, referencePeriod = 0){
  relativeTVec <- tVec- referencePeriod
  tPre <- relativeTVec[ relativeTVec < 0 ]
  tPost <- relativeTVec[ relativeTVec > 0 ]
  slopePre <- solve( t(tPre) %*% (tPre)  ) %*%  tPre
  l_trend_NI <-
    c( -as.numeric(t(lVec) %*% tPost) * slopePre ,
       lVec )
  return(l_trend_NI)
}

#Create constraints related to max/min bias of linear estimator using pre-period
.createConstraintsLinearTrend <- function(A,d, lVec, tVec, referencePeriod = 0){
  v_Trend <- .createV_Linear_NoIntercept(lVec = lVec, tVec = tVec, referencePeriod = referencePeriod)
  maxBias <- .maxBiasFN(v = v_Trend, A = A, d= d)$value
  minBias <- .minBiasFN(v = v_Trend, A = A, d= d)$value
  A_trend <- rbind(v_Trend, -v_Trend)
  d_trend <- c(maxBias, -minBias)
  return( list(A_trend = A_trend, d_trend = d_trend))
}

# HYBRID HELPER FUNCTIONS ---------------------------------------------
.compute_least_favorable_cv <- function(X_T, sigma, hybrid_kappa, sims = 1000,
                                        rowsForARP = NULL) {
  # Computes the least favorable critical value following the algorithm in
  # Section 6.2 of Andrews, Roth, Pakes (2019).
  #
  # Inputs:
  #   X_T   = matrix
  #   sigma = covariance matrix of y_T
  #   hybrid_kappa = desired size
  #   sims  = number of simulations, default = 10000

  if(!is.null(rowsForARP)){
    X_T <- X_T[rowsForARP,]
    sigma <- sigma[rowsForARP, rowsForARP]
  }

  .compute_eta <- function(seed, f, C) {
    set.seed(seed)
    b = -mvtnorm::rmvnorm(n = 1, sigma = sigma) # Define linear constraint
    # Define linear program using lpSolveAPI
    linprog = make.lp(nrow = 0, ncol = length(f))
    set.objfn(linprog, f)
    for (r in 1:dim(C)[1]) {
      add.constraint(linprog, xt = c(C[r,]), type = "<=", rhs = b[r])
    }
    set.bounds(linprog, lower = rep(-Inf, length(f)), columns = 1:length(f))
    # Solve linear program using dual simplex
    lp.control(linprog, sense = "min", simplextype="dual", pivoting = "dantzig", verbose="neutral")
    error_flag = solve(linprog)
    # Return value of eta
    return(get.objective(linprog))
  }

  if (is.null(X_T)) { # no nuisance parameter case
    set.seed(0)
    xi.draws = mvtnorm::rmvnorm(n = sims, sigma = sigma)/sqrt(diag(sigma))
    eta_vec = matrixStats::rowMaxs(xi.draws)
    return(quantile(eta_vec, probs = 1 - hybrid_kappa, names = FALSE))
  } else { # Nuisance parameter case
    # Compute elements for LP that will be re-used
    sdVec = sqrt(diag(sigma))  # standard deviation of y_T elements
    dimDelta = dim(X_T)[2]     # dimension of delta
    f = c(1, rep(0, dimDelta)) # Define objective function
    C = -cbind(sdVec, X_T)     # Define linear constraint

    # Compute eta for each simulation
    eta_vec = purrr::map_dbl(.x = 1:sims, .f = .compute_eta, f = f, C = C)

    # We compute the 1-kappa quantile of eta_vec and return this value
    return(quantile(eta_vec, probs = 1-hybrid_kappa, names = FALSE))
  }
}

.FLCI_computeVloVup <- function(vbar, dbar, S, c) {
  # This function computes the values of Vlo, Vup modified for the FLCI hybrid.
  # Outputs:
  #   flci_vlo = value of vlo associated with FLCI
  #   flci_vup = value of vup associated with FLCI
  # Compute vbarMat
  VbarMat = rbind(t(vbar), -t(vbar))
  # Comptute max_min
  max_or_min = (dbar - (VbarMat %*% S))/(VbarMat %*% c)
  # Compute Vlo, Vup for the FLCI
  vlo = max(max_or_min[VbarMat %*% c < 0])
  vup = min(max_or_min[VbarMat %*% c > 0])
  # Return Vlo, Vup
  return(list(vlo = vlo, vup = vup))
}

# ARP FUNCTIONS -------------------------------------------------------
.lp_conditional_test_fn <- function(theta, y_T, X_T, sigma, alpha,
                                    hybrid_flag, hybrid_list, rowsForARP = 1:length(y_T)) {
  # Performs ARP test of moment inequality E[y_T - X_T \delta] <= 0.
  # It returns an indicator for whether the test rejects, as
  # as the maximum statistic and delta.
  #
  # Inputs:
  #   theta        = value of theta to test
  #   y_T          = vector
  #   X_T          = matrix
  #   sigma        = covariance matrix of y_T
  #   alpha        = desired size of test
  #   hybrid_flag  = Flag for whether implement hybrid test.
  #                   "LF"   = least-favorable hybrid,
  #                   "FLCI" = FLCI hybrid
  #                   "ARP"  = ARP
  #   hybrid_list   = list object with following elements:
  #                     hybrid_kappa    = size of first-stage test
  #                     lf_cv           = least-favorable cv (only if LF hybrid)
  #                     flci_vlo        = value of vlo for FLCI (only if FLCI hybrid)
  #                     flci_vup        = value of vup for FLCI (only if FLCI hybrid)
  #                     flci_pe         = FLCI point estimate (only if FLCI hybrid)
  #                     flci_halflength = half-length of FLCI (only if FLCI hybrid)
  # rowsForARP = an index of which rows of Y_T and X_T should be used in the ARP problem.
  # This allows for ARP to be based on a subset of moments used for the FLCI hybrid (e.g. excluding postperiods)
  #
  # Outputs: list with elements
  #   reject = indicator for wheter test rejects
  #   eta = eta value
  #   delta = delta value
  #   lambda = value of lagrange multipliers

  y_T_ARP <- y_T[rowsForARP]
  X_T_ARP <- X_T[rowsForARP,]
  sigma_ARP <- sigma[rowsForARP, rowsForARP]

  # Dimensions of objects
  M = dim(sigma_ARP)[1]
  k = dim(X_T_ARP)[2]

  # Compute eta, and the argmin delta
  linSoln = .test_delta_lp_fn(y_T = y_T_ARP, X_T = X_T_ARP, sigma = sigma_ARP)

  # Check to see if valid solution to LP obtained
  if (linSoln$error_flag > 0) {
    warning('LP for eta did not converge properly. Not rejecting');
    return(list(reject = 0,
                eta = linSoln$eta_star,
                delta = linSoln$delta_star))
  }

  # HYBRID: Implement first-stage test for hybrid.
  if (hybrid_flag == "LF") {
    # Least-favorable hybrid: check if estimated eta is larger
    # than least favorable critical value. If so, reject immediately.
    mod_size = (alpha - hybrid_list$hybrid_kappa)/(1 - hybrid_list$hybrid_kappa)
    if (linSoln$eta_star > hybrid_list$lf_cv) {
      return(list(reject = 1,
                  eta = linSoln$eta_star,
                  delta = linSoln$delta_star))
    }
  } else if (hybrid_flag == "FLCI") {
    # FLCI hybrid: Test whether parameter of interest falls within FLCI. If not, reject
    mod_size = (alpha - hybrid_list$hybrid_kappa)/(1 - hybrid_list$hybrid_kappa)
    VbarMat = rbind(t(hybrid_list$vbar), -t(hybrid_list$vbar))
    if ( (max(VbarMat %*% y_T - hybrid_list$dbar)) > 0 ) {
      return(list(reject = 1,
                  eta = linSoln$eta_star,
                  delta = linSoln$delta_star))
    }
  } else if (hybrid_flag == "ARP") {
    # No hybrid selected
    mod_size = alpha
  } else {
    # if hybrid flag does not equal LF, FLCI or ARP, return error
    stop("Hybrid flag must equal 'LF', 'FLCI' or 'ARP'")
  }

  # We now check for conditions under which the primal and dual solutions are equal to one another
  # If these conditions don't hold, we will switch to the dual.
  tol_lambda = 10^(-6)
  degenerate_flag = (sum(linSoln$lambda > tol_lambda) != (k+1))

  # Store which moments are binding: Do so using lambda rather than the moments.
  B_index = (linSoln$lambda > tol_lambda)
  Bc_index = (B_index == F)
  X_TB = matrix( X_T_ARP[B_index,], ncol = ncol(X_T_ARP) ) # select binding moments
  # Check whether binding moments have full rank.
  if (is.vector(X_TB)) {
    fullRank_flag = F
  } else {
    fullRank_flag = (rankMatrix(X_TB) == min(dim( X_TB ) ))
  }

  # If degenerate or binding moments don't have full rank, switch to dual
  if (!fullRank_flag | degenerate_flag) {
    ### Dual approach ###
    # warning('Using dual approach')

    # We calculate vlo and vup using the bisection approach that conditions on
    # having a gamma_tilde - a vertex of the dual. This is returned as lambda,
    # from the functioon test_delta_lp_fn.
    lpDualSoln = .lp_dual_fn(y_T = y_T_ARP, X_T = X_T_ARP, eta = linSoln$eta_star, gamma_tilde = linSoln$lambda, sigma = sigma_ARP)

    sigma_B_dual = sqrt( t(lpDualSoln$gamma_tilde) %*% sigma_ARP %*% lpDualSoln$gamma_tilde)

    #If sigma_B_dual is 0 to numerical precision, reject iff eta > 0
    if(sigma_B_dual < 10^(-10)){
      return(list(reject = ifelse(linSoln$eta_star > 0, 1, 0),
                  eta = linSoln$eta_star,
                  delta = linSoln$delta_star,
                  lambda = linSoln$lambda))
    }

    maxstat = lpDualSoln$eta/sigma_B_dual

    # HYBRID: Modify vlo, vup for the hybrid test
    if (hybrid_flag == "LF") {
      # Modify only vup using least-favorable CV for the least-favorable hybrid
      zlo_dual = lpDualSoln$vlo/sigma_B_dual
      zup_dual = min(lpDualSoln$vup, hybrid_list$lf_cv)/sigma_B_dual
    } else if (hybrid_flag == "FLCI") {
      # Compute vlo_FLCI, vup_FLCI.
      gamma_full <- matrix(0, nrow = length(y_T), ncol = 1)
      gamma_full[rowsForARP] <- lpDualSoln$gamma_tilde

      sigma_gamma = (sigma %*% gamma_full) %*% solve((t(gamma_full) %*% sigma %*% gamma_full))
      S = y_T - sigma_gamma %*% (t(gamma_full) %*% y_T)
      vFLCI = .FLCI_computeVloVup(vbar = hybrid_list$vbar, dbar = hybrid_list$dbar,
                                  S = S, c = sigma_gamma)

      # Modify vlo, vup using the FLCI vlo, vup values
      zlo_dual = max(lpDualSoln$vlo, vFLCI$vlo)/sigma_B_dual
      zup_dual = min(lpDualSoln$vup, vFLCI$vup)/sigma_B_dual

    } else if (hybrid_flag == "ARP") {
      # If ARP, don't modify
      zlo_dual = lpDualSoln$vlo/sigma_B_dual
      zup_dual = lpDualSoln$vup/sigma_B_dual
    } else {
      # if hybrid flag does not equal LF, FLCI or ARP, return error
      stop("Hybrid flag must equal 'LF', 'FLCI' or 'ARP'")
    }

    if (!(zlo_dual <= maxstat & maxstat <= zup_dual)) {
      # warning(sprintf("max stat (%3f) is not between z_lo (%3f) and z_up (%3f) in the dual approach", maxstat, zlo_dual, zup_dual))
      return(list(reject = 0,
                  eta = linSoln$eta_star,
                  delta = linSoln$delta_star,
                  lambda = linSoln$lambda))
    } else {
      cval = .norminvp_generalized(p = 1 - mod_size, l = zlo_dual, u = zup_dual)
      reject = as.numeric(c(maxstat) > cval)
      return(list(reject = reject,
                  eta = linSoln$eta_star,
                  delta = linSoln$delta_star,
                  lambda = linSoln$lambda))
    }
  } else {
    ### primal approach ###
    size_B = sum(B_index)
    size_Bc = sum(1 - B_index)

    sdVec = sqrt(diag(sigma_ARP))
    sdVec_B = sdVec[B_index]
    sdVec_Bc = sdVec[Bc_index]

    X_TBc = X_T_ARP[Bc_index, ]
    S_B = diag(M)
    S_B = S_B[B_index, ]
    S_Bc = diag(M)
    S_Bc = S_Bc[Bc_index, ]

    Gamma_B = cbind(sdVec_Bc, X_TBc) %*% solve(cbind(sdVec_B, X_TB)) %*% S_B - S_Bc
    e1 = c(1, rep(0, size_B-1))
    v_B = t( t(e1) %*% solve(cbind(sdVec_B, X_TB)) %*% S_B )
    sigma2_B = t(v_B) %*% sigma_ARP %*% v_B
    sigma_B = sqrt(sigma2_B)
    rho = Gamma_B %*% sigma_ARP %*% v_B %*% solve(sigma2_B)
    maximand_or_minimand = (-Gamma_B %*% y_T_ARP)/rho + c(t(v_B) %*% y_T_ARP)

    # Compute truncation values
    if ( sum(rho > 0) > 0) {
      vlo = max(maximand_or_minimand[rho > 0])
    } else {
      vlo = -Inf
    }
    if ( sum(rho < 0) > 0) {
      vup = min(maximand_or_minimand[rho < 0])
    } else {
      vup = Inf
    }

    # HYBRID: Modify vlo, vup for the hybrid
    if (hybrid_flag == "LF") {
      # if LF, modify vup.
      zlo = vlo/sigma_B
      zup = min(vup, hybrid_list$lf_cv)/sigma_B
    } else if (hybrid_flag == "FLCI") {
      # Compute vlo_FLCI, vup_FLCI.

      gamma_full <- matrix(0, nrow = length(y_T), ncol = 1)
      gamma_full[rowsForARP] <- v_B

      sigma_gamma = (sigma %*% gamma_full) %*% solve((t(gamma_full) %*% sigma %*% gamma_full))
      S = y_T - sigma_gamma %*% (t(gamma_full) %*% y_T)
      vFLCI = .FLCI_computeVloVup(vbar = hybrid_list$vbar, dbar = hybrid_list$dbar,
                                  S = S, c = sigma_gamma)

      # if FLCI, modify both vlo, vup
      zlo = max(vlo, vFLCI$vlo)/sigma_B
      zup = min(vup, vFLCI$vup)/sigma_B

      # if(vlo < vFLCI$vlo){
      #   ### vlo greater after conditioning
      #   x <- 1
      # }

    } else if (hybrid_flag == "ARP") {
      # If not hybrid, don't modify
      zlo = vlo/sigma_B
      zup = vup/sigma_B
    } else {
      # if hybrid flag does not equal LF, FLCI or ARP, return error
      stop("Hybrid flag must equal 'LF', 'FLCI' or 'ARP'")
    }

    # Compute test statistic
    maxstat = linSoln$eta_star/sigma_B

    if ( !(zlo <= maxstat & maxstat <= zup) ) {
      # warning(sprintf("max stat (%3f) is not between z_lo (%3f) and z_up (%3f) in the primal approach", maxstat, zlo, zup))
      return(list(reject = 0,
                  eta = linSoln$eta_star,
                  delta = linSoln$delta_star,
                  lambda = linSoln$lambda))
    } else {
      cval = .norminvp_generalized(p = 1 - mod_size, l = zlo, u = zup)
      reject = as.numeric(c(maxstat) > cval)
      return(list(reject = reject,
                  eta = linSoln$eta_star,
                  delta = linSoln$delta_star,
                  lambda = linSoln$lambda))
    }
  }
}

.test_delta_lp_fn_wrapper <- function(theta, y_T, X_T, sigma, alpha,
                                      hybrid_flag, hybrid_list, rowsForARP = NULL) {
  # This function providers a wrapper for the function
  # test_delta_lp_fn. it only returns an indicator for
  # whether the test rejects.
  # This is used by other functions to construct the ARP CIs.
  results = .lp_conditional_test_fn(theta = theta, y_T = y_T, X_T = X_T, sigma = sigma, alpha = alpha,
                                    hybrid_flag = hybrid_flag, hybrid_list = hybrid_list, rowsForARP = rowsForARP)
  return(results$reject)
}

.ARP_computeCI <- function(betahat, sigma, numPrePeriods, numPostPeriods,
                           A, d, l_vec, alpha,
                           hybrid_flag, hybrid_list,
                           returnLength, grid.lb, grid.ub, gridPoints, rowsForARP = NULL) {
  # This function computes the ARP confidence interval for delta \in Delta = {A delta <= d}.
  #
  # Inputs:
  #   betahat             = vector of estimated event study coefficients
  #   sigma               = covariance matrix of estimated event study coefficients
  #   l_vec               = vector that defines parameter of interest, l'beta.
  #   numPrePeriods       = number of pre-periods
  #   numPostPeriods      = number of post-periods
  #   A                   = matrix A that defines Delta.
  #   d                   = vector d that defines Delta
  #   alpha               = size of CI
  #   hybrid_flag         = flag for hybrid
  #   hybrid_list         = list of objects needed for hybrid
  #   returnLength        = True only returns the length of the CI, False returns the full grid of values
  #   grid.lb/ub          = upper and lower bound of the grid.
  #   gridPoints          = number of grid points to test over.
  #
  # Outputs:
  #   testResultsGrid

  # Construct grid of theta values to test over.
  thetaGrid <- seq(grid.lb, grid.ub, length.out = gridPoints)

  # Construct matrix Gamma and A*(Gamma)^{-1}
  Gamma = .construct_Gamma(l_vec)
  AGammaInv = A[,(numPrePeriods+1):(numPrePeriods+numPostPeriods)] %*% solve(Gamma)
  AGammaInv_one = AGammaInv[,1]
  AGammaInv_minusOne = AGammaInv[,-1] # This will be the matrix X_T in ARP functions

  # Compute Y = A betahat - d and its variance A*sigma*A'
  Y = c(A %*% betahat - d)
  sigmaY = A %*% sigma %*% t(A)

  # HYBRID: Check if hybrid is least favorable, then compute least-favorable CV and add it to hybrid list
  if (hybrid_flag == "LF") {
    hybrid_list$lf_cv = .compute_least_favorable_cv(X_T = AGammaInv_minusOne, sigma = sigmaY, hybrid_kappa = hybrid_list$hybrid_kappa, rowsForARP = rowsForARP)
  }

  # Define helper function to loop over theta grid
  testTheta <- function(theta) {
    # If FLCI Hybrid, compute dbar
    if (hybrid_flag == "FLCI") {
      hybrid_list$dbar = c(hybrid_list$flci_halflength - (t(hybrid_list$vbar) %*% d) + (1 - t(hybrid_list$vbar) %*% AGammaInv_one)*theta,
                           hybrid_list$flci_halflength + (t(hybrid_list$vbar) %*% d) - (1 - t(hybrid_list$vbar) %*% AGammaInv_one)*theta)
    }

    # For the FLCI, we compute the modified FLCI vlo, vup inside this function now.
    reject = .test_delta_lp_fn_wrapper(theta = theta, y_T = Y-AGammaInv_one*theta, X_T = AGammaInv_minusOne,
                                       sigma = sigmaY, alpha = alpha,
                                       hybrid_flag = hybrid_flag, hybrid_list = hybrid_list, rowsForARP = rowsForARP)
    accept = 1-reject
    return(accept)
  }

  # Loop over theta grid and store acceptances
  testResultsGrid <- purrr::map_dbl(.x = thetaGrid, .f = testTheta)
  testValsGrid <- thetaGrid
  resultsGrid <- cbind(testValsGrid, testResultsGrid)
  if( (resultsGrid[1, 2] == 1 | resultsGrid[NROW(resultsGrid), 2] == 1) & hybrid_flag != "FLCI"){
    warning("CI is open at one of the endpoints; CI length may not be accurate")
  }

  # Compute length, else return grid
  if (returnLength == T) {
    gridLength <- 0.5 * ( c(0, diff(thetaGrid)) + c(diff(thetaGrid), 0 ) )
    return(sum(resultsGrid[, 2]*gridLength))
  } else {
    return(tibble(grid = resultsGrid[, 1],
                  accept = resultsGrid[,2]))
  }
}

# DELTA^{SD}(M) FUNCTIONS ---------------------------------------------
# In this section, we implement helper functions to place testing with
# Delta^{SD}(M) into the form needed to use the ARP functions.

.create_A_SD <- function(numPrePeriods, numPostPeriods, postPeriodMomentsOnly = F) {
  # This function creates a matrix for the linear constraints that \delta \in Delta^SD(M).
  # It implements this using the general characterization of A, NOT the sharp
  # characterization of the identified set.
  #
  # Inputs:
  #   numPrePeriods = number of pre-periods. This is an element of resultsObjects.
  #   numPostPeriods = number of post-periods. This is an element of resultsObjects.
  #   postPeriodMomentsOnly = whether to exlude moments relating only to pre-period (which don't affect ID set)

  # First construct matrix Atilde -- (numPrePeriods+numPostPeriods-2) x (numPrePeriods+numPostPeriods+1)
  # Note Atilde is just the positive moments; is not related to Atilde, the rotate matrix, in the paper
  # Note: Atilde initially includes t = 0. We then drop it.
  Atilde = matrix(0, nrow = numPrePeriods+numPostPeriods-1, ncol = numPrePeriods+numPostPeriods+1)
  for (r in 1:(numPrePeriods+numPostPeriods-1)) {
    Atilde[r, r:(r+2)] = c(1, -2, 1)
  }
  Atilde = Atilde[, -(numPrePeriods+1)]

  # If postPeriodMomentsOnly == T, exclude moments that only involve pre-periods
  if(postPeriodMomentsOnly){
    postPeriodIndices <- (numPrePeriods +1):NCOL(Atilde)
    prePeriodOnlyRows <- which( rowSums( Atilde[ , postPeriodIndices] != 0 ) == 0 )
    Atilde <- Atilde[-prePeriodOnlyRows , ]

  }
  # Construct A = [Atilde; -Atilde]
  A = rbind(Atilde, -Atilde)
  return(A)
}

.create_d_SD <- function(numPrePeriods, numPostPeriods, M, postPeriodMomentsOnly = F) {
  # This function creates a vector for the linear constraints that \delta \in Delta^SD(M).
  # It implements this using the general characterization of d, NOT the sharp
  # characterization of the identified set.
  #
  # Inputs:
  #   numPrePeriods  = number of pre-periods. This is an element of resultsObjects.
  #   numPostPeriods = number of post-periods. This is an element of resultsObjects.
  #   M              = smoothness parameter of Delta^SD(M).
  #   postPeriodMomentsOnly = whether to exlude moments relating only to pre-period (which don't affect ID set)
  A_SD = .create_A_SD(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, postPeriodMomentsOnly = postPeriodMomentsOnly)
  d = rep(M, NROW(A_SD))
  return(d)
}

.compute_IDset_DeltaSD <- function(M, trueBeta, l_vec, numPrePeriods, numPostPeriods) {
  # This function computes the upper and lower bound of the identified set
  # given the event study coefficients, lvec and M.
  #
  # Note: lvec is assumed to be non-negative.
  #
  # Inputs:
  #   M              = smoothness param of Delta^SD
  #   trueBeta       = vector of population event study coefficients
  #   l_vec          = vector l defining parameter of interest
  #   numPrePeriods  = number of pre-periods
  #   numPostPeriods = number of post-periods
  #
  # Outputs:
  #   dataframe with columns
  #     id.ub = upper bound of ID set
  #     id.lb = lower bound of ID set
  #     M     = M value passed

  # Create objective function: Wish to min/max l'delta_post
  fDelta = c(rep(0, numPrePeriods), l_vec)

  # Create A, d that define Delta^SDPB(M)
  A_SD = .create_A_SD(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods)
  d_SD = .create_d_SD(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, M = M)
  dir_SD = rep("<=", NROW(A_SD))

  # Create equality constraint that delta_pre = beta_pre
  prePeriodEqualityMat = cbind(diag(numPrePeriods),
                               matrix(data = 0, nrow = numPrePeriods, ncol = numPostPeriods))
  A_SD = rbind(A_SD, prePeriodEqualityMat)
  d_SD = c(d_SD,
           trueBeta[1:numPrePeriods])
  dir_SD = c(dir_SD,
             rep("==", NROW(prePeriodEqualityMat)))
  bounds = list(lower = list(ind = 1:(numPrePeriods + numPostPeriods), val = rep(-Inf, numPrePeriods+numPostPeriods)),
                upper = list(ind = 1:(numPrePeriods + numPostPeriods), val = rep(Inf, numPrePeriods+numPostPeriods)))

  # Create and solve for max
  results.max = Rglpk::Rglpk_solve_LP(obj = fDelta,
                                      max = TRUE,
                                      mat = A_SD,
                                      dir = dir_SD,
                                      rhs = d_SD,
                                      bounds = bounds)

  # Create and solve for min
  results.min = Rglpk::Rglpk_solve_LP(obj = fDelta,
                                      max = FALSE,
                                      mat = A_SD,
                                      dir = dir_SD,
                                      rhs = d_SD,
                                      bounds = bounds)

  if (results.max$status != 0 & results.min$status != 0) {
    warning("Solver did not find an optimum")
    id.ub = (t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)])
    id.lb = (t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)])
  }
  else {
    # Construct upper/lower bound of identified set
    id.ub = (t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]) - results.min$optimum
    id.lb = (t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]) - results.max$optimum
  }
  # Return identified set
  return(tibble(
    id.lb = id.lb,
    id.ub = id.ub))
}

computeConditionalCS_DeltaSD <- function(betahat, sigma, numPrePeriods, numPostPeriods,
                                    l_vec = .basisVector(index = 1, size = numPostPeriods), M = 0,
                                    alpha = 0.05, hybrid_flag = "FLCI", hybrid_kappa = alpha/10,
                                    returnLength = F, postPeriodMomentsOnly = T,
                                    gridPoints =10^3, grid.midPoint = NA, grid.ub = NA, grid.lb = NA) {
  # This function computes the ARP CI that includes nuisance parameters
  # for Delta^{SD}(M). This functions uses ARP_computeCI for all
  # of its computations.
  #
  # Inputs:
  #   betahat             = vector of estimated event study coefficients
  #   sigma               = covariance matrix of estimated event study coefficients
  #   numPrePeriods       = number of pre-periods
  #   numPostPeriods      = number of post-periods
  #   l_vec               = vector that defines parameter of interest
  #   M                   = tuning parameter for Delta^SD(M), default M = 0.
  #   alpha               = desired size of CI, default alpha = 0.05
  #   hybrid_flag         = flag for hybrid, default = "ARP"
  #   hybrid_kappa        = desired size of first-stage hybrid test, default = NULL
  #   returnLength        = returns length of CI only. Otherwise, returns matrix with grid in col 1 and test result in col 2.
  #   numGridPoints       = number of gridpoints to test over, default = 1000
  #   postPeriodMomentsOnly = exclude moments for delta^SD that only include pre-period coefs
  #
  #  Outputs:
  #   data_frame containing upper and lower bounds of CI.

  # Construct A_SD, d_SD
  A_SD = .create_A_SD(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, postPeriodMomentsOnly = F)
  d_SD = .create_d_SD(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, M = M, postPeriodMomentsOnly = F)

  if(postPeriodMomentsOnly){
    postPeriodIndices <- (numPrePeriods +1):NCOL(A_SD)
    postPeriodRows <- which( rowSums( A_SD[ , postPeriodIndices] != 0 ) > 0 )
    rowsForARP <- postPeriodRows
  }else{
    rowsForARP <- 1:NROW(A_SD)
  }

  # Create hybrid_list object
  hybrid_list = list(hybrid_kappa = hybrid_kappa)

  # HYBRID: If hybrid, compute FLCI
  if (hybrid_flag == "FLCI") {

    flci = .findOptimalFLCI_helper(sigma = sigma, M = M,
                           numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                           l_vec = l_vec, alpha = hybrid_kappa)

    # Add objects to hybrid_list: flci l vector
    hybrid_list$flci_l = flci$optimalVec

    # Add vbar to flci l vector
    vbar = Variable(NROW(A_SD))
    obj <- Minimize( t(flci$optimalVec) %*% flci$optimalVec -
                       2 * t(flci$optimalVec) %*% t(A_SD) %*% vbar + quad_form(x = vbar, P = A_SD %*% t(A_SD)) )
    prob = Problem(obj)
    result = psolve(prob)
    hybrid_list$vbar = result$getValue(vbar)

    # Add objects to hybrid_list: flci half-length
    hybrid_list$flci_halflength = flci$optimalHalfLength

    # compute FLCI ub and FLCI lb
    if (is.na(grid.ub)){
      grid.ub = (t(flci$optimalVec) %*% betahat) + flci$optimalHalfLength
    }
    if (is.na(grid.lb)){
      grid.lb = (t(flci$optimalVec) %*% betahat) - flci$optimalHalfLength
    }

  } else {
    # Compute identified set under parallel trends
    IDset = .compute_IDset_DeltaSD(M = M, trueBeta = rep(0, numPrePeriods + numPostPeriods),
                                   l_vec = l_vec, numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods)
    sdTheta <- sqrt(t(l_vec) %*% sigma[(numPrePeriods+1):(numPrePeriods+numPostPeriods), (numPrePeriods+1):(numPrePeriods+numPostPeriods)] %*% l_vec)
    if(is.na(grid.ub)){
      grid.ub = IDset$id.ub + 20*sdTheta
    }
    if(is.na(grid.lb)){
      grid.lb = IDset$id.lb - 20*sdTheta
    }
  }

  # Compute ARP CI for l'beta using Delta^SD
  CI = .ARP_computeCI(betahat = betahat, sigma = sigma, numPrePeriods = numPrePeriods,
                      numPostPeriods = numPostPeriods, A = A_SD, d = d_SD,
                      l_vec = l_vec, alpha = alpha,
                      hybrid_flag = hybrid_flag, hybrid_list = hybrid_list,
                      returnLength = returnLength,
                      grid.lb = grid.lb, grid.ub = grid.ub,
                      gridPoints = gridPoints, rowsForARP = rowsForARP)
  # Returns CI
  return(CI)
}

# DELTA^{SDB}(M) FUNCTIONS --------------------------------------------
# In this section, we implement helper functions to place testing with
# Delta^{SDB}(M) into the form needed to use the ARP functions.
.create_A_SDB <- function(numPrePeriods, numPostPeriods,
                          biasDirection = "positive", postPeriodMomentsOnly = F) {
  # This function creates a matrix for the linear constraints that \delta \in Delta^SDPB(M).
  # It implements this using the general characterization of A.
  #
  # Inputs:
  #   numPrePeriods = number of pre-periods. This is an element of resultsObjects.
  #   numPostPeriods = number of post-periods. This is an element of resultsObjects.

  A_SD = .create_A_SD(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                      postPeriodMomentsOnly = postPeriodMomentsOnly)
  A_PB = -diag(numPrePeriods + numPostPeriods)
  A_PB = A_PB[(numPrePeriods+1):(numPrePeriods+numPostPeriods), ]

  if(biasDirection == "negative"){
    A_PB <- -A_PB
  }else if(biasDirection != "positive"){
    stop("Input biasDirection must equal either `positive' or `negative'")
  }

  A = rbind(A_SD, A_PB)
  return(A)
}

.create_d_SDB <- function(numPrePeriods, numPostPeriods, M, postPeriodMomentsOnly = F) {
  # This function creates a vector for the linear constraints that \delta \in Delta^SDPB(M).
  # It implements this using the general characterization of d.
  #
  # Inputs:
  #   numPrePeriods  = number of pre-periods. This is an element of resultsObjects.
  #   numPostPeriods = number of post-periods. This is an element of resultsObjects.
  #   M              = smoothness parameter of Delta^SD(M).

  d_SD = .create_d_SD(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, M = M, postPeriodMomentsOnly = postPeriodMomentsOnly)
  d_PB = rep(0, numPostPeriods)
  d = c(d_SD, d_PB)
  return(d)
}

.compute_IDset_DeltaSDB <- function(M, trueBeta, l_vec,
                                    numPrePeriods, numPostPeriods, biasDirection) {
  # This function computes the upper and lower bound of the identified set given the true
  # population event study coefficients.
  #
  # Inputs:
  #   M             = smoothness param
  #   trueBeta       = vec. of population event study coeff.
  #   l_vec          = vec. l defining param. of interest
  #   numPrePeriods  = number of pre-periods
  #   numPostPeriods = number of post-periods
  #   biasDirection  = direction of Bias

  # Create objective function: Wish to min/max l'delta_post
  fDelta = c(rep(0, numPrePeriods), l_vec)

  # Create A, d that define Delta^SDPB(M)
  A_SDB = .create_A_SDB(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                        biasDirection = biasDirection)
  d_SDB = .create_d_SDB(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, M = M)
  dir_SDB = rep("<=", NROW(A_SDPB))

  # Create equality constraint that delta_pre = beta_pre
  prePeriodEqualityMat = cbind(diag(numPrePeriods),
                               matrix(data = 0, nrow = numPrePeriods, ncol = numPostPeriods))
  A_SDB = rbind(A_SDB, prePeriodEqualityMat)
  d_SDB = c(d_SDB,
            trueBeta[1:numPrePeriods])
  dir_SDB = c(dir_SDB,
              rep("==", NROW(prePeriodEqualityMat)))
  bounds = list(lower = list(ind = 1:(numPrePeriods + numPostPeriods), val = rep(-Inf, numPrePeriods+numPostPeriods)),
                upper = list(ind = 1:(numPrePeriods + numPostPeriods), val = rep(Inf, numPrePeriods+numPostPeriods)))


  # Create and solve for max
  results.max = Rglpk::Rglpk_solve_LP(obj = fDelta,
                                      max = TRUE,
                                      mat = A_SDB,
                                      dir = dir_SDB,
                                      rhs = d_SDB,
                                      bounds = bounds)

  # Create and solve for min
  results.min = Rglpk::Rglpk_solve_LP(obj = fDelta,
                                      max = FALSE,
                                      mat = A_SDB,
                                      dir = dir_SDB,
                                      rhs = d_SDB,
                                      bounds = bounds)

  if (results.min$status != 0 & results.max$status != 0) {
    warning("Solver did not find an optimum")
    id.ub = (t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)])
    id.lb = (t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)])
  }
  else {
    # Construct upper/lower bound of identified set
    id.ub = (t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]) - results.min$optimum
    id.lb = (t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]) - results.max$optimum
  }
  # Return identified set
  return(tibble(
    id.lb = id.lb,
    id.ub = id.ub))
}

computeConditionalCS_DeltaSDB <- function(betahat, sigma, numPrePeriods, numPostPeriods, M = 0,
                                     l_vec = .basisVector(index = 1, size=numPostPeriods),
                                     alpha = 0.05,
                                     hybrid_flag = "FLCI", hybrid_kappa = alpha/10,
                                     returnLength = F, biasDirection = "positive",
                                     postPeriodMomentsOnly = T,
                                     gridPoints = 10^3,
                                     grid.lb = NA,
                                     grid.ub = NA) {
  # This function computes the ARP CI that includes nuisance parameters
  # for Delta^{SDPB}(M). This functions uses ARP_computeCI for all
  # of its computations.
  #
  # Inputs:
  #   betahat        = vector of estimated event study coefficients
  #   sigma          = covariance matrix of estimated event study coefficients
  #   numPrePeriods  = number of pre-periods
  #   l_vec          = vector that defines parameter theta
  #   postPeriod     = post-period of interest, set to first post-period by default
  #   M              = tuning parameter for Delta^SD(M), default M = 0.
  #   alpha          = desired size of CI, default alpha = 0.05
  #   hybrid_flag    = flag for hybrid default = "ARP".
  #   hybrid_kappa   = desired size of first-stage least favorable test, default = NULL
  #   returnLength   = returns length of CI only. Otherwise, returns matrix with grid in col 1 and test result in col 2.
  #   gridPoints     = number of gridpoints to test over, default = 1000
  #   postPeriodMomentsOnly = include the delta^SD moments corresponding with the first period only
  #
  #  Outputs:
  #   data_frame containing upper and lower bounds of CI.

  # Construct A_SDB, d_SDB
  A_SDB = .create_A_SDB(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                        biasDirection = biasDirection,
                        postPeriodMomentsOnly = F)
  d_SDB = .create_d_SDB(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, M = M, postPeriodMomentsOnly = F)

  if (postPeriodMomentsOnly) {
    postPeriodIndices <- (numPrePeriods +1):NCOL(A_SDB)
    postPeriodRows <- which( rowSums( A_SDB[ , postPeriodIndices] != 0 ) > 0 )
    rowsForARP <- postPeriodRows
  } else{
    rowsForARP <- 1:NROW(A_SDB)
  }

  # Create hybrid_list object
  hybrid_list = list(hybrid_kappa = hybrid_kappa)

  # HYBRID: If hybrid, compute FLCI
  if (hybrid_flag == "FLCI") {
    flci = .findOptimalFLCI_helper(sigma = sigma, M = M, numPrePeriods = numPrePeriods,
                           numPostPeriods = numPostPeriods, l_vec = l_vec, alpha = hybrid_kappa)

    # Add objects to hybrid_list: flci l vector
    hybrid_list$flci_l = flci$optimalVec

    # Add vbar to flci l vector
    vbar = Variable(NROW(A_SDB))
    obj <- Minimize( t(flci$optimalVec) %*% flci$optimalVec -
                       2 * t(flci$optimalVec) %*% t(A_SDB) %*% vbar + quad_form(x = vbar, P = A_SDB %*% t(A_SDB)) )
    prob = Problem(obj)
    result = psolve(prob)
    hybrid_list$vbar = result$getValue(vbar)

    # Add objects to hybrid_list: flci half-length
    hybrid_list$flci_halflength = flci$optimalHalfLength

    # compute FLCI ub and FLCI lb
    if (is.na(grid.ub)){
      grid.ub = (t(flci$optimalVec) %*% betahat) + flci$optimalHalfLength
    }
    if (is.na(grid.lb)){
      grid.lb = (t(flci$optimalVec) %*% betahat) - flci$optimalHalfLength
    }
  } else {
    # Compute ID set under parallel trends
    IDset = .compute_IDset_DeltaSDB(M = M, trueBeta = rep(0, numPrePeriods + numPostPeriods),
                                    l_vec = l_vec, numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods)
    # If direction is negative, flip the signs and the upper and lower bounds
    if(biasDirection == "negative"){
      new.lb <- - IDset$id.ub
      new.ub <- - IDset$id.lb
      IDset$id.lb <- new.lb
      IDset$id.ub <- new.ub
    }
    sdTheta <- sqrt(t(l_vec) %*% sigma[(numPrePeriods+1):(numPrePeriods+numPostPeriods), (numPrePeriods+1):(numPrePeriods+numPostPeriods)] %*% l_vec)

    if (is.na(grid.ub)) {
      grid.ub = IDset$id.ub + 20*sdTheta
    }
    if (is.na(grid.lb)) {
      grid.lb = IDset$id.lb - 20*sdTheta
    }
  }

  # Compute ARP CI for l'beta using Delta^SD
  CI = .ARP_computeCI(betahat = betahat, sigma = sigma, numPrePeriods = numPrePeriods,
                      numPostPeriods = numPostPeriods, A = A_SDB, d = d_SDB,
                      l_vec = l_vec, alpha = alpha,
                      hybrid_flag = hybrid_flag, hybrid_list = hybrid_list,
                      returnLength = returnLength,
                      grid.lb = grid.lb, grid.ub = grid.ub,
                      gridPoints = gridPoints, rowsForARP = rowsForARP)
  # Returns CI
  return(CI)
}

# DELTA^{SDM}(M) FUNCTIONS --------------------------------------------
# In this section, we implement helper functions to place testing with
# Delta^{SDM}(M) into the form needed to use the ARP functions.

.create_A_M <- function(numPrePeriods, numPostPeriods,
                        monotonicityDirection = "increasing", postPeriodMomentsOnly = F){
  #This function creates a matrix so that A \delta <= 0 implies delta is increasing/decreasing depending on what direction is specified
  A_M = matrix(0, nrow = numPrePeriods+numPostPeriods, ncol=numPrePeriods+numPostPeriods)
  for(r in 1:(numPrePeriods-1)){
    A_M[r, r:(r+1)] <- c(1,-1)
  }
  A_M[numPrePeriods, numPrePeriods] <- 1
  if(numPostPeriods > 0){
    A_M[numPrePeriods + 1, numPrePeriods + 1] <- -1
    for (r in (numPrePeriods + 2):(numPrePeriods+numPostPeriods)) {
      A_M[r, (r-1):r] <- c(1,-1)
    }
  }
  # If postPeriodMomentsOnly == T, exclude moments that only involve pre-periods
  if(postPeriodMomentsOnly){
    postPeriodIndices <- (numPrePeriods +1):NCOL(A_M)
    prePeriodOnlyRows <- which( rowSums( A_I[ , postPeriodIndices] != 0 ) == 0 )
    A_M <- A_M[-prePeriodOnlyRows , ]
  }
  if (monotonicityDirection == "decreasing") {
    A_M <- -A_M
  } else if(monotonicityDirection != "increasing") {
    stop("direction must be 'increasing' or 'decreasing'")
  }
  return(A_M)
}

.create_A_SDM <- function(numPrePeriods, numPostPeriods,
                          monotonicityDirection = "increasing", postPeriodMomentsOnly = F) {
  # This function creates a matrix for the linear constraints that \delta \in Delta^SDI(M).
  # It implements this using the general characterization of A.
  #
  # Inputs:
  #   numPrePeriods = number of pre-periods. This is an element of resultsObjects.
  #   numPostPeriods = number of post-periods. This is an element of resultsObjects.
  A_M <- .create_A_M(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                     monotonicityDirection = monotonicityDirection, postPeriodMomentsOnly = postPeriodMomentsOnly)
  A_SD = .create_A_SD(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, postPeriodMomentsOnly =  postPeriodMomentsOnly)

  A = rbind(A_SD, A_M)
  return(A)
}

.create_d_SDM <- function(numPrePeriods, numPostPeriods, M, postPeriodMomentsOnly = F) {
  # This function creates a vector for the linear constraints that \delta \in Delta^SDI(M).
  # It implements this using the general characterization of d.
  #
  # Inputs:
  #   numPrePeriods  = number of pre-periods. This is an element of resultsObjects.
  #   numPostPeriods = number of post-periods. This is an element of resultsObjects.
  #   M              = smoothness parameter of Delta^SD(M).

  d_SD = .create_d_SD(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, M = M, postPeriodMomentsOnly = postPeriodMomentsOnly)
  d_I = rep(0, ifelse(postPeriodMomentsOnly, numPostPeriods, numPrePeriods+numPostPeriods) )
  d = c(d_SD, d_I)
  return(d)
}

.compute_IDset_DeltaSDM <- function(M, trueBeta, l_vec, numPrePeriods,
                                    numPostPeriods, monotonicityDirection) {

  # Create objective function: Wish to min/max l'delta_post
  fDelta = c(rep(0, numPrePeriods), l_vec)

  # Create A, d that define Delta^SDM(M)
  A_SDM = .create_A_SDM(numPrePeriods = numPrePeriods,
                        numPostPeriods = numPostPeriods, monotonicityDirection = monotonicityDirection)
  d_SDM = .create_d_SDM(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, M = M)
  dir_SDM = rep("<=", NROW(A_SDM))

  # Create equality constraint that delta_pre = beta_pre
  prePeriodEqualityMat = cbind(diag(numPrePeriods),
                               matrix(data = 0, nrow = numPrePeriods, ncol = numPostPeriods))
  A_SDM = rbind(A_SDM, prePeriodEqualityMat)
  d_SDM = c(d_SDM,
            trueBeta[1:numPrePeriods])
  dir_SDM = c(dir_SDM,
              rep("==", NROW(prePeriodEqualityMat)))
  bounds = list(lower = list(ind = 1:(numPrePeriods + numPostPeriods),
                             val = rep(-Inf, numPrePeriods+numPostPeriods)),
                upper = list(ind = 1:(numPrePeriods + numPostPeriods),
                             val = rep(Inf, numPrePeriods+numPostPeriods)))

  # Create and solve for max
  max.results = Rglpk::Rglpk_solve_LP(obj = fDelta,
                                      max = TRUE,
                                      mat = A_SDM,
                                      dir = dir_SDM,
                                      rhs = d_SDM,
                                      bounds = bounds)

  # Create and solve for min
  min.results = Rglpk::Rglpk_solve_LP(obj = fDelta,
                                      max = FALSE,
                                      mat = A_SDM,
                                      dir = dir_SDM,
                                      rhs = d_SDM,
                                      bounds = bounds)

  if (max.results$status != 0 & min.results$status != 0) {
    warning("Solver did not find an optimum")
    id.ub = (t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)])
    id.lb = (t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)])
  }
  else {
    # Construct upper/lower bound of identified set
    id.ub = (t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]) - min.results$optimum
    id.lb = (t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]) - max.results$optimum
  }
  # Return identified set
  return(tibble(
    id.lb = id.lb,
    id.ub = id.ub))
}

computeConditionalCS_DeltaSDM <- function(betahat, sigma, numPrePeriods, numPostPeriods, M = 0,
                                     l_vec = .basisVector(index = 1, size = numPostPeriods),
                                     alpha = 0.05, monotonicityDirection = "increasing",
                                     hybrid_flag = "FLCI", hybrid_kappa = alpha/10,
                                     returnLength = F, postPeriodMomentsOnly = T,
                                     gridPoints=10^3,
                                     grid.lb = NA,
                                     grid.ub = NA) {
  # This function computes the ARP CI that includes nuisance parameters
  # for Delta^{SDI}(M). This functions uses ARP_computeCI for all
  # of its computations.
  #
  # Inputs:
  #   betahat        = vector of estimated event study coefficients
  #   sigma          = covariance matrix of estimated event study coefficients
  #   numPrePeriods  = number of pre-periods
  #   numPostPeriods = number of post-periods
  #   l_vec          = vector that defines parameter of interest theta
  #   M              = tuning parameter for Delta^SD(M), default M = 0.
  #   alpha          = desired size of CI, default alpha = 0.05
  #   hybrid_flag    = flag for hybrid default = "ARP".
  #   hybrid_kappa   = desired size of first-stage least favorable test, default = NULL
  #   returnLength   = returns length of CI only. Otherwise, returns matrix with grid in col 1 and test result in col 2.
  #   postPeriodMomentsOnly = exclude moments that relate only to the pre-period coeffs
  #   gridPoints     = number of gridpoints to test over, default = 1000
  #
  #  Outputs:
  #   data_frame containing upper and lower bounds of CI.

  # Construct A_SDM, d_SDM
  A_SDM = .create_A_SDM(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                        monotonicityDirection = monotonicityDirection, postPeriodMomentsOnly = F)
  d_SDM = .create_d_SDM(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, M = M, postPeriodMomentsOnly = F)

  if (postPeriodMomentsOnly) {
    postPeriodIndices <- (numPrePeriods +1):NCOL(A_SDM)
    postPeriodRows <- which( rowSums( A_SDM[ , postPeriodIndices] != 0 ) > 0 )
    rowsForARP <- postPeriodRows
  } else {
    rowsForARP <- 1:NROW(A_SDM)
  }

  # Create hybrid_list object
  hybrid_list = list(hybrid_kappa = hybrid_kappa)

  # HYBRID: If hybrid with FLCI, compute FLCI
  if (hybrid_flag == "FLCI") {
    flci = .findOptimalFLCI_helper(sigma = sigma, M = M, numPrePeriods = numPrePeriods,
                           numPostPeriods = numPostPeriods, l_vec = l_vec, alpha = hybrid_kappa)

    # Add objects to hybrid_list: flci l vector
    hybrid_list$flci_l = flci$optimalVec

    # Add vbar to flci l vector
    vbar = Variable(NROW(A_SDM))
    obj <- Minimize( t(flci$optimalVec) %*% flci$optimalVec -
                       2 * t(flci$optimalVec) %*% t(A_SDM) %*% vbar + quad_form(x = vbar, P = A_SDM %*% t(A_SDM)) )
    prob = Problem(obj)
    result = psolve(prob)
    hybrid_list$vbar = result$getValue(vbar)

    # Add objects to hybrid_list: flci half-length
    hybrid_list$flci_halflength = flci$optimalHalfLength

    # compute FLCI ub and FLCI lb
    grid.ub = (t(hybrid_list$flci_l) %*% betahat) + flci$optimalHalfLength
    grid.lb = (t(hybrid_list$flci_l) %*% betahat) - flci$optimalHalfLength
  } else {
    # Compute ID set
    IDset = .compute_IDset_DeltaSDM(M = M, trueBeta = rep(0, numPrePeriods + numPostPeriods),
                                    monotonicityDirection = monotonicityDirection, l_vec = l_vec,
                                    numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods)
    sdTheta <- sqrt(t(l_vec) %*% sigma[(numPrePeriods+1):(numPrePeriods+numPostPeriods), (numPrePeriods+1):(numPrePeriods+numPostPeriods)] %*% l_vec)
    if (is.na(grid.ub)) {
      grid.ub = IDset$id.ub + 20*sdTheta
    }
    if (is.na(grid.lb)) {
      grid.lb = IDset$id.lb - 20*sdTheta
    }
  }

  # Compute ARP CI for l'beta using Delta^SDI
  CI = .ARP_computeCI(betahat = betahat, sigma = sigma, numPrePeriods = numPrePeriods,
                      numPostPeriods = numPostPeriods, A = A_SDM, d = d_SDM,
                      l_vec = l_vec, alpha = alpha,
                      hybrid_flag = hybrid_flag, hybrid_list = hybrid_list,
                      returnLength = returnLength,
                      grid.lb = grid.lb, grid.ub = grid.ub,
                      gridPoints = gridPoints, rowsForARP = rowsForARP)

  # Returns CI
  return(CI)
}

# Delta^{RMI} functions -----------------------------------------------
.create_A_RMI <- function(numPrePeriods, numPostPeriods,
                          Mbar, monotonicityDirection, postPeriodMomentsOnly = F) {
  # This function creates a matrix for the linear constraints that \delta \in \Delta^{RM}(M) \intersect
  # with delta increasing or decreasing.
  #
  # Inputs:
  #   numPrePeriods = number of pre-periods.
  #   numPostPeriods = number of post-periods
  #   direction = "increasing" or "decreasing"
  #   postPeriodMomentsOnly = whether to exclude moments relating only to pre-period, which don't affect the ID set.
  #   Mbar = parameter that governs Delta^RM set.

  # Create matrix for the increasing/decreasing moments
  A_I = .create_A_M(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                    monotonicityDirection = "increasing", postPeriodMomentsOnly = postPeriodMomentsOnly)

  # Create matrix for the RM moments
  A_RMI_pos = matrix(0, nrow = numPrePeriods + numPostPeriods - 1, ncol = numPrePeriods + numPostPeriods + 1)
  for (r in 1:(numPrePeriods + numPostPeriods - 1)) {
    A_RMI_pos[r, r:(r+2)] = c(Mbar, -Mbar + 1, -1)
  }
  A_RMI_pos = A_RMI_pos[, -(numPrePeriods+1)]

  A_RMI_neg = matrix(0, nrow = numPrePeriods + numPostPeriods - 1, ncol = numPrePeriods + numPostPeriods + 1)
  for (r in 1:(numPrePeriods + numPostPeriods - 1)) {
    A_RMI_neg[r, r:(r+2)] = c(Mbar, -Mbar-1, 1)
  }
  A_RMI_neg = A_RMI_neg[, -(numPrePeriods+1)]
  A_RMI_tilde = rbind(A_RMI_pos, A_RMI_neg)

  # Impose postPeriodMomentsOnly on A_RMI
  if (postPeriodMomentsOnly) {
    postPeriodIndices <- (numPrePeriods +1):NCOL(A_RMI_tilde)
    prePeriodOnlyRows <- which( rowSums( A_RMI_tilde[ , postPeriodIndices] != 0 ) == 0 )
    A_RMI_tilde <- A_RMI_tilde[-prePeriodOnlyRows , ]
  }

  # Construct A_RMI matrix
  A_RMI = rbind(A_RMI_tilde, A_I)
  return(A_RMI)
}

.create_d_RMI <- function(numPrePeriods, numPostPeriods, postPeriodMomentsOnly = F) {
  # This function creates a vector for the linear constraints that \delta \in Delta^RM(M) \intersect increasing/decreasing.
  # It implements this using the general characterization of d.
  #
  # Inputs:
  #   numPrePeriods  = number of pre-periods. This is an element of resultsObjects.
  #   numPostPeriods = number of post-periods. This is an element of resultsObjects.
  d_RMI = rep(0, ifelse(postPeriodMomentsOnly, postPeriods, 2*(numPrePeriods + numPostPeriods-1)))
  d_I = rep(0, ifelse(postPeriodMomentsOnly, numPostPeriods, numPrePeriods+numPostPeriods))
  d = c(d_RMI, d_I)
  return(d)
}

.compute_IDset_DeltaRMI <- function(Mbar, trueBeta, l_vec, numPrePeriods, numPostPeriods) {

  # Create objective function: Wish to min/max l'delta_post
  fDelta = c(rep(0, numPrePeriods), l_vec)

  # Create A, d that define Delta^SDPB(M)
  A_RMI = .create_A_RMI(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, Mbar = Mbar)
  d_RMI = .create_d_RMI(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods)
  dir_RMI = rep("<=", NROW(A_RMI))

  # Create equality constraint that delta_pre = beta_pre
  prePeriodEqualityMat = cbind(diag(numPrePeriods),
                               matrix(data = 0, nrow = numPrePeriods, ncol = numPostPeriods))
  A_RMI = rbind(A_RMI, prePeriodEqualityMat)
  d_RMI = c(d_RMI,
            trueBeta[1:numPrePeriods])
  dir_RMI = c(dir_RMI,
              rep("==", NROW(prePeriodEqualityMat)))
  bounds = list(lower = list(ind = 1:(numPrePeriods + numPostPeriods), val = rep(-Inf, numPrePeriods+numPostPeriods)),
                upper = list(ind = 1:(numPrePeriods + numPostPeriods), val = rep(Inf, numPrePeriods+numPostPeriods)))

  # Create and solve for max
  max.results = Rglpk::Rglpk_solve_LP(obj = fDelta,
                                      max = TRUE,
                                      mat = A_RMI,
                                      dir = dir_RMI,
                                      rhs = d_RMI,
                                      bounds = bounds)

  # Create and solve for min
  min.results = Rglpk::Rglpk_solve_LP(obj = fDelta,
                                      max = FALSE,
                                      mat = A_RMI,
                                      dir = dir_RMI,
                                      rhs = d_RMI,
                                      bounds = bounds)

  if (max.results$status != 0 & min.results$status != 0) {
    warning("Solver did not find an optimum")
    id.ub = (t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)])
    id.lb = (t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)])
  }
  else {
    # Construct upper/lower bound of identified set
    id.ub = (t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]) - min.results$optimum
    id.lb = (t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]) - max.results$optimum
  }
  # Return identified set
  return(tibble(
    id.lb = id.lb,
    id.ub = id.ub))
}

computeConditionalCS_DeltaRMI <- function(betahat, sigma, numPrePeriods, numPostPeriods, Mbar = 0,
                                     l_vec = .basisVector(index = 1, size = numPostPeriods),
                                     alpha = 0.05,
                                     hybrid_flag = "LF", hybrid_kappa = alpha/10,
                                     returnLength = F, postPeriodMomentsOnly = T,
                                     gridPoints=10^3,
                                     grid.lb = NA, grid.ub = NA) {
  # This function computes the ARP CI that includes nuisance parameters
  # for Delta^{RM}(Mbar). This functions uses ARP_computeCI for all
  # of its computations.
  #
  # Inputs:
  #   betahat        = vector of estimated event study coefficients
  #   sigma          = covariance matrix of estimated event study coefficients
  #   numPrePeriods  = number of pre-periods
  #   numPostPeriods = number of post-periods
  #   l_vec          = vector that defines parameter of interest theta
  #   Mbar           = tuning parameter for Delta^RM(Mbar), default Mbar = 0.
  #   alpha          = desired size of CI, default alpha = 0.05
  #   hybrid_flag    = flag for hybrid default = "ARP".
  #   hybrid_kappa   = desired size of first-stage least favorable test, default = NULL
  #   returnLength   = returns length of CI only. Otherwise, returns matrix with grid in col 1 and test result in col 2.
  #   postPeriodMomentsOnly = exclude moments that relate only to the pre-period coeffs
  #   gridPoints     = number of gridpoints to test over, default = 1000
  #
  #  Outputs:
  #   data_frame containing upper and lower bounds of CI.

  # Construct A_SDI, d_SDI
  A_RMI = .create_A_RMI(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, Mbar = Mbar, postPeriodMomentsOnly = F)
  d_RMI = .create_d_RMI(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, postPeriodMomentsOnly = F)

  if(postPeriodMomentsOnly){
    postPeriodIndices <- (numPrePeriods +1):NCOL(A_RMI)
    postPeriodRows <- which( rowSums( A_RMI[ , postPeriodIndices] != 0 ) > 0 )
    rowsForARP <- postPeriodRows
  }else{
    rowsForARP <- 1:NROW(A_RMI)
  }

  # Create hybrid_list object
  hybrid_list = list(hybrid_kappa = hybrid_kappa)

  # HYBRID: If hybrid with FLCI, compute FLCI
  if (hybrid_flag == "FLCI") {
    stop("The FLCI hybrid is not an option for Delta^{RMI}!")
  } else {
    # Compute ID set
    IDset = .compute_IDset_DeltaRMI(Mbar = Mbar, trueBeta = rep(0, numPrePeriods + numPostPeriods), l_vec = l_vec,
                                    numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods)
    sdTheta <- sqrt(t(l_vec) %*% sigma[(numPrePeriods+1):(numPrePeriods+numPostPeriods), (numPrePeriods+1):(numPrePeriods+numPostPeriods)] %*% l_vec)
    if (is.na(grid.ub)) {
      grid.ub = IDset$id.ub + 20*sdTheta
    }
    if (is.na(grid.lb)) {
      grid.lb = IDset$id.lb - 20*sdTheta
    }
  }

  # Compute ARP CI for l'beta using Delta^SDI
  CI = .ARP_computeCI(betahat = betahat, sigma = sigma, numPrePeriods = numPrePeriods,
                      numPostPeriods = numPostPeriods, A = A_RMI, d = d_RMI,
                      l_vec = l_vec, alpha = alpha,
                      hybrid_flag = hybrid_flag, hybrid_list = hybrid_list,
                      returnLength = returnLength,
                      grid.lb = grid.lb, grid.ub = grid.ub,
                      gridPoints = gridPoints, rowsForARP = rowsForARP)

  # Returns CI
  return(CI)
}

# Lower and upper bounding M helper functions -------------------------
.selectionMat <- function(selection, size, select = "columns"){

  if(select == "rows"){
    m <- matrix(0, nrow = length(selection), ncol = size)

    m[ 1:length(selection), selection] <- diag(length(selection))
  }else{
    m <- matrix(0, ncol = length(selection), nrow = size)

    m[ selection , 1:length(selection)] <- diag(length(selection))
  }
  return(m)
}

.LeeCFN <- function(eta, Sigma){
  c = Sigma %*% eta / as.numeric( t(eta) %*% Sigma %*% eta )
  return(c)
}

.VLoVUpFN <- function(eta, Sigma, A, b, z){
  c = .LeeCFN(eta,Sigma)

  objective <- (b - A %*% z) / (A %*% c)

  ACNegativeIndex <- which( (A %*% c) < 0 )
  ACPositiveIndex <- which( (A %*% c) > 0 )

  if( length(ACNegativeIndex) == 0){
    VLo <- -Inf
  }else{
    VLo <- max(objective[ACNegativeIndex])
  }

  if( length(ACPositiveIndex) == 0){
    VUp <- Inf
  }else{
    VUp <- min(objective[ACPositiveIndex])
  }

  return(c(VLo,VUp))
}

.testInIdentifiedSet_Max <- function(M, y, sigma, A,alpha, d) {
  # Runs APR test of the moments E[AY - 1*M] <= 0, where Y ~ N(mu, sigma).
  # We construct this such that this tests whether the mean of the max moment equals thetabar.

  d_mod = d*M
  sigmaTilde <- as.vector( sqrt( diag(A %*% sigma %*% t(A)) ) )
  Atilde <- solve( diag(sigmaTilde) ) %*% A
  dtilde <- solve( diag(sigmaTilde) ) %*% d_mod

  normalizedMoments <- Atilde %*% y - dtilde
  maxLocation <- which.max(normalizedMoments)
  maxMoment <- normalizedMoments[maxLocation]

  T_B <- .selectionMat(maxLocation, size = NROW(Atilde), select = "rows")
  iota <- matrix(1, nrow = NROW(Atilde), ncol = 1)

  gamma <- t(T_B %*% Atilde)
  Abar <- Atilde - iota %*% T_B %*% Atilde
  dbar <- ( diag(NROW(dtilde)) - iota %*% T_B ) %*% dtilde

  sigmabar <- sqrt( t(gamma) %*% sigma %*% gamma )
  c <- sigma %*% gamma / as.numeric( t(gamma) %*% sigma %*% gamma  )
  z <- (diag(NROW(y)) - c %*% t(gamma)) %*% y
  VLoVUpVec <- .VLoVUpFN(eta = gamma, Sigma = sigma, A = Abar, b = dbar, z = z)
  criticalVal <- .norminvp_generalized(p = 1-alpha, l = VLoVUpVec[1], u = VLoVUpVec[2],
                                       mu = T_B %*% dtilde, sd = sigmabar)
  reject <- (maxMoment + T_B %*% dtilde > criticalVal)

  return(reject)
}

.create_A_and_D_SD_prePeriods <- function(numPrePeriods) {
  Atilde = matrix(0, nrow = numPrePeriods-1, ncol = numPrePeriods)
  if (numPrePeriods < 2) {
    stop("Can't estimate M in pre-period with < 2 pre-period coeffs")
  } else {
    Atilde[numPrePeriods-1, (numPrePeriods-1):(numPrePeriods)] <- c(1,-2)
    for (r in 1:(numPrePeriods-2)) {
      Atilde[r, r:(r+2)] = c(1, -2, 1)
    }
    A.pre <- rbind(Atilde, -Atilde)
    d = rep(1, NROW(A.pre))
    return(list(A = A.pre, d = d))
  }
}

.estimate_lowerBound_M_conditionalTest <- function(prePeriod.coef, prePeriod.covar,
                                                   grid.ub, alpha = 0.05, gridPoints) {
  # This function constructs a lower-bound for M using APR.
  numPrePeriods <- length(prePeriod.coef)
  # Helper function
  .APR_testOverMGrid <- function(prePeriod.coef, prePeriod.covar,
                                 mGrid, A, d, alpha) {
    # This function runs the APR test over a grid of possible values of M.
    .testMInSet <- function(maxVal) {
      reject <- .testInIdentifiedSet_Max(M = maxVal,
                                         y = prePeriod.coef, sigma =  prePeriod.covar,
                                         A = A, d = d, alpha = alpha)
      accept = 1 - reject
      return(accept)
    }
    ResultsGrid <-purrr::map_dbl(.x = mGrid, .testMInSet)
    ValsGrid <- mGrid
    return(cbind(ValsGrid, ResultsGrid))
  }
  Ad = .create_A_and_D_SD_prePeriods(numPrePeriods = numPrePeriods)
  # Construct grid of M values
  mGrid = seq(from = 0, to = grid.ub, length = gridPoints)
  # Compute APR test at each value of M in Grid
  resultsGrid = .APR_testOverMGrid(prePeriod.coef = prePeriod.coef, prePeriod.covar = prePeriod.covar,
                                   mGrid = mGrid, A = Ad$A, d = Ad$d, alpha = alpha)
  if (sum(resultsGrid[,2]) == 0) {
    warning("ARP conditional test rejects all values of M provided. User should increase upper bound of grid.")
    return(Inf)
  } else {
    return(min(resultsGrid[ resultsGrid[,2] == 1 ,1]))
  }
}

# Lower and upper bounding M functions --------------------------------
DeltaSD_upperBound_Mpre <- function(betahat, sigma, numPrePeriods, alpha = 0.05) {
  # This function constructs an upper-bound for M at the 1-alpha level
  # based on the observed pre-period coefficients.

  prePeriod.coef = betahat[1:numPrePeriods]
  prePeriod.sigma = sigma[1:numPrePeriods, 1:numPrePeriods]

  A_SD = .create_A_SD(numPrePeriods = numPrePeriods, numPostPeriods = 0)
  prePeriodCoefDiffs = A_SD %*% prePeriod.coef
  prePeriodSigmaDiffs = A_SD %*% prePeriod.sigma %*% t(A_SD)
  seDiffs = sqrt(diag(prePeriodSigmaDiffs))
  upperBoundVec = prePeriodCoefDiffs + qnorm(1-alpha)*seDiffs
  maxUpperBound = max(upperBoundVec)
  return(maxUpperBound)
}

DeltaSD_lowerBound_Mpre <- function(betahat, sigma, numPrePeriods, alpha = 0.05, grid.ub = NA, gridPoints = 1000) {
  # This function constructs a lower bound for M using the observed pre-period coefficients by
  # constructing a one-sided confidence interval on the maximal second difference of the observed
  # pre-period coefficients using the conditional test in Andrews, Roth, Pakes (2019)
  prePeriod.coef = betahat[1:numPrePeriods]
  prePeriod.sigma = sigma[1:numPrePeriods, 1:numPrePeriods]

  if (is.na(grid.ub)) {
    # If Mub not specified, use 3 times the max SE in the preperiod
    grid.ub <- 3 * max(sqrt(diag(prePeriod.sigma)))
  }
  results <- .estimate_lowerBound_M_conditionalTest(prePeriod.coef = prePeriod.coef,
                                                    prePeriod.covar = prePeriod.sigma,
                                                    alpha = alpha, grid.ub = grid.ub,
                                                    gridPoints = gridPoints)
  return(results)
}

# Construct Robust Results Function -----------------------------------
createSensitivityResults <- function(betahat, sigma,
                              numPrePeriods, numPostPeriods,
                              method = NULL,
                              Mvec = NULL,
                              l_vec = .basisVector(index = 1, size = numPostPeriods),
                              monotonicityDirection = NULL,
                              biasDirection = NULL,
                              alpha = 0.05,
                              parallel = FALSE) {

  # If Mvec is null, construct default Mvec
  if (is.null(Mvec)) {
    Mub = DeltaSD_upperBound_Mpre(betahat = betahat, sigma = sigma, numPrePeriods = numPrePeriods, alpha = 0.05)
    Mvec = seq(from = 0, to = Mub, length.out = 10)
  }

  # If Delta = Delta^{SD}, construct confidence intervals
  if (is.null(monotonicityDirection) & is.null(biasDirection)) {
    if (is.null(method)) {
      method = "FLCI"
    }

    if (method == "FLCI") { # If method = FLCI, construct FLCIs
      if (parallel == FALSE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = 0.05)
          tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                 method = "FLCI", Delta = "DeltaSD", M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = 0.05)
          tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                 method = "FLCI", Delta = "DeltaSD", M = Mvec[m])
        }
      }
    } else if (method == "Conditional") { # If method = Conditional, construct conditional confidence intervals
      if (parallel == FALSE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                              hybrid_flag = "ARP")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "Conditional",
                 Delta = "DeltaSD", M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                              hybrid_flag = "ARP")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "Conditional",
                 Delta = "DeltaSD", M = Mvec[m])
        }
       }
      } else if (method == "C-F") { # If method = C-F, construct conditional FLCI,
        if (parallel == FALSE) {
          Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
            temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                                numPrePeriods = numPrePeriods,
                                                numPostPeriods = numPostPeriods,
                                                l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                                hybrid_flag = "FLCI")
            tibble(lb = min(temp$grid[temp$accept == 1]),
                   ub = max(temp$grid[temp$accept == 1]),
                   method = "C-F",
                   Delta = "DeltaSD", M = Mvec[m])
          }
        } else {
          Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
            temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                                numPrePeriods = numPrePeriods,
                                                numPostPeriods = numPostPeriods,
                                                l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                                hybrid_flag = "FLCI")
            tibble(lb = min(temp$grid[temp$accept == 1]),
                   ub = max(temp$grid[temp$accept == 1]),
                   method = "C-F",
                   Delta = "DeltaSD", M = Mvec[m])
          }
        }
      } else if (method == "C-LF") {
        if (parallel == FALSE) {
          Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
            temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                                numPrePeriods = numPrePeriods,
                                                numPostPeriods = numPostPeriods,
                                                l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                                hybrid_flag = "LF")
            tibble(lb = min(temp$grid[temp$accept == 1]),
                   ub = max(temp$grid[temp$accept == 1]),
                   method = "C-LF",
                   Delta = "DeltaSD", M = Mvec[m])
          }
        } else {
          Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
            temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                                numPrePeriods = numPrePeriods,
                                                numPostPeriods = numPostPeriods,
                                                l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                                hybrid_flag = "LF")
            tibble(lb = min(temp$grid[temp$accept == 1]),
                   ub = max(temp$grid[temp$accept == 1]),
                   method = "C-LF",
                   Delta = "DeltaSD", M = Mvec[m])
          }
        }
      } else{
        stop("Method must equal one of: FLCI, Conditional, C-F or C-LF")
      }
  } else if (!is.null(biasDirection)) {
    if (is.null(method)) {
      method = "C-F"
    }

    if (biasDirection == "positive") {
      Delta = "DeltaSDPB"
     } else {
      Delta = "DeltaSDNB"
    }

    if (method == "FLCI") {
      warning("You specified a sign restriction but method = FLCI. The FLCI does not use the sign restriction!")
      if (parallel == FALSE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = 0.05)
          tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                 method = "FLCI", Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = 0.05)
          tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                 method = "FLCI", Delta = Delta, M = Mvec[m])
        }
      }
    } else if (method == "Conditional") { # If method = Conditional, construct conditional confidence intervals
      if (parallel == TRUE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "ARP")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "Conditional",
                 Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "ARP")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "Conditional",
                 Delta = Delta, M = Mvec[m])
        }
      }
    } else if (method == "C-F") { # If method = C-F, construct conditional FLCI,
      if (parallel == TRUE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "FLCI")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-F",
                 Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "FLCI")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-F",
                 Delta = Delta, M = Mvec[m])
        }
      }
    } else if (method == "C-LF") {
      if (parallel == TRUE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "LF")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-LF",
                 Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "LF")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-LF",
                 Delta = Delta, M = Mvec[m])
        }
      }
    } else{
      stop("Method must equal one of: FLCI, Conditional, C-F or C-LF")
    }
  } else {
    if (is.null(method)) {
      method = "C-F"
    }

    if (monotonicityDirection == "increasing") {
      Delta = "DeltaSDI"
    } else {
      Delta = "DeltaSDD"
    }

    if (method == "FLCI") {
      warning("You specified a shape restriction but method = FLCI. The FLCI does not use the shape restriction!")
      if (parallel == FALSE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = 0.05)
          tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                 method = "FLCI", Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = 0.05)
          tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                 method = "FLCI", Delta = Delta, M = Mvec[m])
        }
      }
    } else if (method == "Conditional") { # If method = Conditional, construct conditional confidence intervals
      if (parallel == FALSE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "ARP")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "Conditional",
                 Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "ARP")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "Conditional",
                 Delta = Delta, M = Mvec[m])
        }
      }
    } else if (method == "C-F") { # If method = C-F, construct conditional FLCI,
      if (parallel == FALSE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "FLCI")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-F",
                 Delta = Delta,
                 M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "FLCI")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-F",
                 Delta = Delta,
                 M = Mvec[m])
        }
      }
    } else if (method == "C-LF") {
      if (parallel == FALSE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "LF")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-LF",
                 Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "LF")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-LF",
                 Delta = Delta, M = Mvec[m])
        }
      }
    } else{
      stop("Method must equal one of: FLCI, Conditional, C-F or C-LF")
    }
  }

  return(Results)
}

constructOriginalCS <- function(betahat, sigma,
                           numPrePeriods, numPostPeriods,
                           l_vec = .basisVector(index = 1, size = numPostPeriods),
                           alpha = 0.05) {
  stdError = sqrt(t(l_vec) %*% sigma[(numPrePeriods+1):(numPrePeriods+numPostPeriods), (numPrePeriods+1):(numPrePeriods+numPostPeriods)] %*% l_vec)
  lb = t(l_vec) %*% betahat[(numPrePeriods+1):(numPrePeriods+numPostPeriods)] - qnorm(1-alpha)*stdError
  ub = t(l_vec) %*% betahat[(numPrePeriods+1):(numPrePeriods+numPostPeriods)] + qnorm(1-alpha)*stdError
  return(tibble(
    lb = lb,
    ub = ub,
    method = "Original",
    Delta = NA,
    M = 0
  ))
}

# Sensitivity plot functions ------------------------------------------
createEventStudyPlot <- function(betahat, stdErrors = NULL, sigma = NULL,
                                 numPrePeriods, numPostPeriods,
                                 timeVec, referencePeriod,
                                 useRelativeEventTime = F) {
  if (is.null(stdErrors) & is.null(sigma)) {
    stop("User must specify either vector of standard errors or vcv matrix!")
  } else if (is.null(stdErrors) & !is.null(sigma)) {
    stdErrors = sqrt(diag(sigma))
  }

  if (useRelativeEventTime == T) {
    timeVec = timeVec - referencePeriod
    referencePeriod = 0
  }

  EventStudyPlot <- ggplot(tibble(t = c(timeVec[1:numPrePeriods], referencePeriod, timeVec[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]),
                                  beta = c(betahat[1:numPrePeriods], 0, betahat[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]),
                                  se = c(stdErrors[1:numPrePeriods], NA, stdErrors[(numPrePeriods+1):(numPrePeriods+numPostPeriods)])),
                           aes(x = t)) +
    geom_point(aes(y = beta), color = "red") +
    geom_errorbar(aes(ymin = beta - qnorm(0.975)*se, ymax = beta + qnorm(0.975)*se), width = 0.5, colour = "#01a2d9") +
    theme(legend.position = "none") + labs(x = "Event time", y = "") +
    scale_x_continuous(breaks = seq(from = min(timeVec), to = max(timeVec), by = 1),
                       labels = as.character(seq(from = min(timeVec), to = max(timeVec), by = 1)))
  return(EventStudyPlot)
}

createSensitivityPlot <- function(robustResults, originalResults, rescaleFactor = 1, maxM = Inf, add_xAxis = TRUE) {
  # Set M for OLS to be the min M in robust results minus the gap between Ms in robust
  Mgap <- min( diff( sort( robustResults$M) ) )
  Mmin <- min( robustResults$M)
  originalResults$M <- Mmin - Mgap
  df <- bind_rows(originalResults, robustResults)

  # Rescale all the units by rescaleFactor
  df <- df %>% mutate_at( c("M", "ub", "lb"), ~ .x * rescaleFactor)

  # Filter out observations above maxM (after rescaling)
  df <- df %>% filter(M <= maxM)

  p <- ggplot(data = df, aes(x=M)) +
    geom_errorbar(aes(ymin = lb, ymax = ub, color = factor(method)),
                  width = Mgap * rescaleFactor / 2) +
    scale_color_manual(values = c("red", '#01a2d9')) +
    theme(legend.title=element_blank(), legend.position="bottom") +
    labs(x = "M", y = "")

  if (add_xAxis) {
    p <- p + geom_hline(yintercept = 0)
  }
  return(p)
}

# Other helper functions ----
basisVector <- function(index = 1, size = 1){
  v <- matrix(0, nrow = size, ncol = 1)
  v[index] = 1
  return(v)
}
