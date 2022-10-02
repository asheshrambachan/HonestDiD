# DESCRIPTION =========================================================
# Author: Ashesh Rambachan <asheshr@g.harvard.edu>
#
#  This script contains functions to implement the methods
#  described in Rambachan & Roth (2019) for robust inference
#  in difference-in-differences and event study designs.
#
#  This script contains functions that are used to construct
#  the FLCI for a general choice of vector l and Delta = Delta^{SD}(M)

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
    return(as.data.frame(temp, stringsAsFactors = FALSE))
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
  hstar <- .findOptimalCIDerivativeBisection(hMin, h0, M, numPoints, alpha, sigma, numPrePeriods, numPostPeriods, l_vec, T)

  if ( is.na(hstar) ) {
    # Numerical derivatives will occasionally fail; fall back into grid
    # search in that case
    hGrid <- seq(from = hMin, to = c(h0), length.out = numPoints)
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
  } else {
    optimalCIDF <- .findWorstCaseBiasGivenH(hstar, sigma, numPrePeriods, numPostPeriods, l_vec, T)
    optimalCIDF$m <- M
    optimalCIDF$CI.halflength <- .qfoldednormal(p = 1-alpha, mu = (M * optimalCIDF$value)/hstar) * hstar
  }

  results = list(
    optimalVec = c(unlist(optimalCIDF$optimal.l), l_vec),
    optimalPrePeriodVec = unlist(optimalCIDF$optimal.l),
    optimalHalfLength = optimalCIDF$CI.halflength,
    M = optimalCIDF$m,
    status = optimalCIDF$status
  )
  return(results)
}

.findOptimalCIDerivativeBisection <-  function(a, b, M, numPoints, alpha, sigma,
                                               numPrePeriods, numPostPeriods, l_vec, returnDF) {

  # Function of h, which is convex (returns CI half length)
  .f <- function(h, ...) {
    biasDF <- .findWorstCaseBiasGivenH(h, sigma, numPrePeriods, numPostPeriods, l_vec, returnDF)
    maxBias <- M * biasDF$value
    if (biasDF$value < Inf) {
      return(.qfoldednormal(p = 1-alpha, mu = maxBias/h) * h)
    } else {
      return(NaN)
    }
  }

  # Minimum search relying on convexity; check for failures at each
  # step; return failure and fall back on grid if failed
  hstar   <- NaN
  failtol <- (.Machine$double.eps)^(1/2)
  failed  <- FALSE
  dif     <- min((b - a) / numPoints, abs(b) * (.Machine$double.eps^(1/3)))
  fa      <- .f(a)
  fb      <- .f(b)
  fpa     <- (.f(a + dif) - fa) / dif  # Limit from the right for lb
  fpb     <- (.f(b - dif) - fb) / -dif # Limit from the left for ub
  iter    <- 1
  maxiter <- 10 * ceiling(log(abs(b - a) / dif) / log(2))

  if ( (fpa > fpb) | is.nan(fa) | is.nan(fb) ) {
    failed <-  TRUE
  } else if ( fpb < 0 ) {
    hstar <- b
  } else if ( fpa > 0 ) {
    hstar <- a
  } else {
    while ( !failed & abs(b - a) > dif ) {
      iter   <- iter + 1
      x      <- (a + b) / 2
      fpx    <- (.f(x + dif) - .f(x - dif)) / (2 * dif)
      failed <- (fpx > fpb + failtol) | (fpx + failtol < fpa) | iter > maxiter
      # fx     <- .f(x)
      if ( fpx > 0 ) {
        b <- x
      } else {
        a <- x
      }
    }
    hstar <- (a + b) / 2
  }

  # Only does 4 + 2 * iter function evaluations
  if (failed) {
    return(NaN)
  } else {
    return(hstar)
  }
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
