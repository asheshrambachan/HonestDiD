# DESCRIPTION =========================================================
# Author: Ashesh Rambachan <asheshr@g.havard.edu>
#
#  This script contains functions to implement the methods
#  described in Rambachan & Roth (2019) for robust inference
#  in difference-in-differences and event study designs.
#
#  Implements functions to upper and lower bound M.


.testInIdentifiedSet_Max <- function(M, y, sigma, A,alpha, d) {
  # Runs APR test of the moments E[AY - 1*M] <= 0, where Y ~ N(mu, sigma).
  # We construct this such that this tests whether the mean of the max moment equals thetabar.

  d_mod = d*M
  sigmaTilde <- base::as.vector( base::sqrt( base::diag(A %*% sigma %*% base::t(A)) ) )
  Atilde <- base::solve( base::diag(sigmaTilde) ) %*% A
  dtilde <- base::solve( base::diag(sigmaTilde) ) %*% d_mod

  normalizedMoments <- Atilde %*% y - dtilde
  maxLocation <- base::which.max(normalizedMoments)
  maxMoment <- normalizedMoments[maxLocation]

  T_B <- .selectionMat(maxLocation, size = base::NROW(Atilde), select = "rows")
  iota <- base::matrix(1, nrow = base::NROW(Atilde), ncol = 1)

  gamma <- base::t(T_B %*% Atilde)
  Abar <- Atilde - iota %*% T_B %*% Atilde
  dbar <- ( base::diag(base::NROW(dtilde)) - iota %*% T_B ) %*% dtilde

  sigmabar <- base::sqrt( base::t(gamma) %*% sigma %*% gamma )
  c <- sigma %*% gamma / base::as.numeric( base::t(gamma) %*% sigma %*% gamma  )
  z <- (base::diag(base::NROW(y)) - c %*% base::t(gamma)) %*% y
  VLoVUpVec <- .VLoVUpFN(eta = gamma, Sigma = sigma, A = Abar, b = dbar, z = z)
  criticalVal <- .norminvp_generalized(p = 1-alpha, l = VLoVUpVec[1], u = VLoVUpVec[2],
                                       mu = T_B %*% dtilde, sd = sigmabar)
  reject <- (maxMoment + T_B %*% dtilde > criticalVal)

  base::return(reject)
}

.create_A_and_D_SD_prePeriods <- function(numPrePeriods) {
  Atilde = base::matrix(0, nrow = numPrePeriods-1, ncol = numPrePeriods)
  if (numPrePeriods < 2) {
    base::stop("Can't estimate M in pre-period with < 2 pre-period coeffs")
  } else {
    Atilde[numPrePeriods-1, (numPrePeriods-1):(numPrePeriods)] <- base::c(1,-2)
    for (r in 1:(numPrePeriods-2)) {
      Atilde[r, r:(r+2)] = base::c(1, -2, 1)
    }
    A.pre <- base::rbind(Atilde, -Atilde)
    d = base::rep(1, base::NROW(A.pre))
    base::return(list(A = A.pre, d = d))
  }
}

.estimate_lowerBound_M_conditionalTest <- function(prePeriod.coef, prePeriod.covar,
                                                   grid.ub, alpha = 0.05, gridPoints) {
  # This function constructs a lower-bound for M using APR.
  numPrePeriods <- base::length(prePeriod.coef)
  # Helper function
  .APR_testOverMGrid <- function(prePeriod.coef, prePeriod.covar,
                                 mGrid, A, d, alpha) {
    # This function runs the APR test over a grid of possible values of M.
    .testMInSet <- function(maxVal) {
      reject <- .testInIdentifiedSet_Max(M = maxVal,
                                         y = prePeriod.coef, sigma =  prePeriod.covar,
                                         A = A, d = d, alpha = alpha)
      accept = 1 - reject
      base::return(accept)
    }
    ResultsGrid <-purrr::map_dbl(.x = mGrid, .testMInSet)
    ValsGrid <- mGrid
    base::return(base::cbind(ValsGrid, ResultsGrid))
  }
  Ad = .create_A_and_D_SD_prePeriods(numPrePeriods = numPrePeriods)
  # Construct grid of M values
  mGrid = base::seq(from = 0, to = grid.ub, length = gridPoints)
  # Compute APR test at each value of M in Grid
  resultsGrid = .APR_testOverMGrid(prePeriod.coef = prePeriod.coef, prePeriod.covar = prePeriod.covar,
                                   mGrid = mGrid, A = Ad$A, d = Ad$d, alpha = alpha)
  if (base::sum(resultsGrid[,2]) == 0) {
    base::warning("ARP conditional test rejects all values of M provided. User should increase upper bound of grid.")
    base::return(Inf)
  } else {
    base::return(base::min(resultsGrid[ resultsGrid[,2] == 1 ,1]))
  }
}

# Lower and upper bounding M functions --------------------------------
DeltaSD_upperBound_Mpre <- function(betahat, sigma, numPrePeriods, alpha = 0.05) {
  # This function constructs an upper-bound for M at the 1-alpha level
  # based on the observed pre-period coefficients.

  base::stopifnot(numPrePeriods > 1)

  prePeriod.coef = betahat[1:numPrePeriods]
  prePeriod.sigma = sigma[1:numPrePeriods, 1:numPrePeriods]

  A_SD = .create_A_SD(numPrePeriods = numPrePeriods, numPostPeriods = 0)
  prePeriodCoefDiffs = A_SD %*% prePeriod.coef
  prePeriodSigmaDiffs = A_SD %*% prePeriod.sigma %*% base::t(A_SD)
  seDiffs = base::sqrt(base::diag(prePeriodSigmaDiffs))
  upperBoundVec = prePeriodCoefDiffs + stats::qnorm(1-alpha)*seDiffs
  maxUpperBound = base::max(upperBoundVec)
  base::return(maxUpperBound)
}

DeltaSD_lowerBound_Mpre <- function(betahat, sigma, numPrePeriods, alpha = 0.05, grid.ub = NA, gridPoints = 1000) {
  # This function constructs a lower bound for M using the observed pre-period coefficients by
  # constructing a one-sided confidence interval on the maximal second difference of the observed
  # pre-period coefficients using the conditional test in Andrews, Roth, Pakes (2019)

  base::stopifnot(numPrePeriods > 1)

  prePeriod.coef = betahat[1:numPrePeriods]
  prePeriod.sigma = sigma[1:numPrePeriods, 1:numPrePeriods]

  if (base::is.na(grid.ub)) {
    # If Mub not specified, use 3 times the max SE in the preperiod
    grid.ub <- base::c(3 * base::max(base::sqrt(base::diag(prePeriod.sigma))))
  }
  results <- .estimate_lowerBound_M_conditionalTest(prePeriod.coef = prePeriod.coef,
                                                    prePeriod.covar = prePeriod.sigma,
                                                    alpha = alpha, grid.ub = grid.ub,
                                                    gridPoints = gridPoints)
  base::return(results)
}
