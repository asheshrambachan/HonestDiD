# DESCRIPTION =========================================================
# Author: Ashesh Rambachan <asheshr@g.harvard.edu>
#
#  This script contains functions to implement the methods
#  described in Rambachan & Roth (2019) for robust inference
#  in difference-in-differences and event study designs.
#
#  This script contains functions that are used to construct
#  the ARP test with nuisance parameters.

# PRELIMINARIES =======================================================
library(tidyverse)
library(TruncatedNormal)
library(lpSolveAPI)
library(ROI)
library(Matrix)
library(pracma)
library(CVXR)
library(foreach)

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
  return(linprog)
}

.roundeps <- function(x, eps = .Machine$double.eps^(3/4)) {
    if ( abs(x) < eps ) return (0) else return(x)
}

.check_if_solution_helper <- function(c, tol, s_T, gamma_tilde, sigma, W_T) {
  # Solve linear program and use negative of objective because we want the max.
  linprog = .max_program(s_T, gamma_tilde, sigma, W_T, c)
  linprog$honestsolution = (abs(c - (-linprog$optimum)) <= tol)
  return(linprog)
}

.vlo_vup_dual_fn <- function(eta, s_T, gamma_tilde, sigma, W_T) {
  # This function computes vlo and vup for the dual linear program for the
  # conditional values using the bisection method described in Appendix
  # D of ARP (2021) in Algorithm 1.
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
  tol_c        = 1e-6
  tol_equality = 1e-6
  sigma_B      = sqrt( t(gamma_tilde) %*% sigma %*% gamma_tilde )
  low_initial  = min(-100, eta-20*sigma_B)
  high_initial = max(100, eta + 20*sigma_B)
  maxiters     = 10000
  switchiters  = 10
  checksol     = .check_if_solution_helper(eta,
                                           tol_equality,
                                           s_T,
                                           gamma_tilde,
                                           sigma,
                                           W_T)$honestsolution

  if (is.na(checksol) || !checksol) {
    # warning('User-supplied eta is not a solution. Not rejecting automatically')
    return( list(vlo = eta, vup = Inf) )
  }

  ### Compute vup ###
  if ( (linprog = .check_if_solution_helper(high_initial,
                                            tol_equality,
                                            s_T,
                                            gamma_tilde,
                                            sigma,
                                            W_T))$honestsolution ) {
    vup = Inf
  } else {
    # Shortcut: We want to find c s.t. the maximized value is equal to c. This
    # takes the form a + b c = c, so for a given value we can find a, b and
    # compute the value of c s.t. a + (b - 1) c = 0. In practice iterating
    # over c in this way is faster than the bisection (it also appears to be
    # more precise); however, if it's taking too long we switch to the bisection.
    dif   = 0
    iters = 1
    b     = c(t(gamma_tilde) %*% sigma %*% gamma_tilde)^(-1)*(sigma %*% gamma_tilde)
    mid   = (.roundeps(linprog$solution %*% s_T) / (1 - linprog$solution %*% b))[1]
    while ( !(linprog = .check_if_solution_helper(mid,
                                                  tol_equality,
                                                  s_T,
                                                  gamma_tilde,
                                                  sigma,
                                                  W_T))$honestsolution && (iters < maxiters) ) {
      iters = iters + 1
      if ( iters >= switchiters ) {
        dif = tol_c + 1
        break
      }
      mid = (.roundeps(linprog$solution %*% s_T) / (1 - linprog$solution %*% b))[1]
    }

    # Throughout the while loop, the high value is not a solution and the
    # low-value is a solution. If the midpoint between them is a solution,
    # then we set low to the midpoint; otherwise, we set high to the midpoint.
    low   = eta
    high  = mid
    while ( dif > tol_c & iters < maxiters ) {
      iters = iters+1
      mid   = (high + low) / 2
      if ( .check_if_solution_helper(mid,
                                     tol_equality,
                                     s_T,
                                     gamma_tilde,
                                     sigma,
                                     W_T)$honestsolution ) {
        low = mid
      } else {
        high = mid
      }
      dif = high - low
    }
    if (iters == maxiters) {
      # warning("vup: Reached max iterations without a solution!")
    }
    vup = mid
  }

  ### Compute vlo ###
  if ( (linprog = .check_if_solution_helper(low_initial,
                                            tol_equality,
                                            s_T,
                                            gamma_tilde,
                                            sigma,
                                            W_T))$honestsolution ) {
    vlo = -Inf
  } else {
    # Shortcut: See shortcut notes above
    dif   = 0
    iters = 1
    b     = c(t(gamma_tilde) %*% sigma %*% gamma_tilde)^(-1)*(sigma %*% gamma_tilde)
    mid   = (.roundeps(linprog$solution %*% s_T) / (1 - linprog$solution %*% b))[1]
    while ( !(linprog = .check_if_solution_helper(mid,
                                                  tol_equality,
                                                  s_T,
                                                  gamma_tilde,
                                                  sigma,
                                                  W_T))$honestsolution && (iters < maxiters) ) {
      iters = iters + 1
      if ( iters >= switchiters ) {
        dif = tol_c + 1
        break
      }
      mid = (.roundeps(linprog$solution %*% s_T) / (1 - linprog$solution %*% b))[1]
    }

    # Throughout the while loop, the low value is not a solution and the
    # high-value is a solution. If the midpoint between them is a solution,
    # then we set high to the midpoint; otherwise, we set low to the midpoint.
    low  = mid
    high = eta
    while (dif > tol_c & iters < maxiters) {
      mid   = (low + high) / 2
      iters = iters + 1
      if (.check_if_solution_helper(mid,
                                    tol_equality,
                                    s_T,
                                    gamma_tilde,
                                    sigma,
                                    W_T)$honestsolution) {
        high = mid
      } else {
        low = mid
      }
      dif = high - low
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

  if (!is.null(rowsForARP)) {
    if (is.vector(X_T)) {
      X_T = matrix(X_T[rowsForARP], ncol = 1)
    } else {
      X_T <- X_T[rowsForARP,]
    }
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
    xi.draws = t(t(mvtnorm::rmvnorm(n = sims, sigma = sigma))/sqrt(diag(sigma)))
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
  if (is.vector(X_T)) {
    X_T = matrix(X_T, ncol = 1)
  }
  X_T_ARP <- X_T[rowsForARP,]
  if (is.vector(X_T_ARP)) {
    X_T_ARP = matrix(X_T_ARP, ncol = 1)
  }

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
      # Per ARP (2021), CV = max(0, c_{1-alpha}), where c_{1-alpha} is the 1-alpha
      # quantile of truncated normal.
      cval = max(0, .norminvp_generalized(p = 1 - mod_size, l = zlo_dual, u = zup_dual))
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
      # Per ARP (2021), CV = max(0, c_{1-alpha}), where c_{1-alpha} is the 1-alpha
      # quantile of truncated normal.
      cval = max(0, .norminvp_generalized(p = 1 - mod_size, l = zlo, u = zup))
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
