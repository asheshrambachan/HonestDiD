# DESCRIPTION =========================================================
# Author: Ashesh Rambachan <asheshr@g.harvard.edu>
#
#  This script contains functions to implement the methods
#  described in Rambachan & Roth (2019) for robust inference
#  in difference-in-differences and event study designs.
#
#  This script contains functions that are used to construct
#  the ARP test without nuisance parameters.

.testInIdentifiedSet <- function(y, sigma, A, d,
                                    Abar_additional = NULL, dbar_additional = NULL, alpha) {
  # Runs APR test of the moments E[AY] - d <= 0, where Y ~ N(mu, sigma) and mu <= 0 under the null.
  # The APR test conditions on the location of the binding moment, which we write as Abar Y <= dbar.
  # This can be used, e.g. for hybrid values with the FLCI

  sigmaTilde <- as.vector( sqrt( diag(A %*% sigma %*% t(A)) ) )
  Atilde <- solve( diag(sigmaTilde) ) %*% A
  dtilde <- solve( diag(sigmaTilde) ) %*% d

  normalizedMoments <- Atilde %*% y - dtilde
  maxLocation <- which.max(normalizedMoments)
  maxMoment <- normalizedMoments[maxLocation]

  T_B <- .selectionMat(maxLocation, size = NROW(Atilde), select = "rows")
  iota <- matrix(1, nrow = NROW(Atilde), ncol = 1)

  gamma <- t(T_B %*% Atilde)
  Abar <- Atilde - iota %*% T_B %*% Atilde
  dbar <- ( diag(NROW(dtilde)) - iota %*% T_B ) %*% dtilde

  # If statement, modifies Abar for the FLCI hybrid
  if (!is.null(Abar_additional)) {
    Abar <- rbind(Abar, Abar_additional)
    dbar <- c(dbar, dbar_additional)
  }

  sigmabar <- sqrt( t(gamma) %*% sigma %*% gamma )
  c <- sigma %*% gamma / as.numeric( t(gamma) %*% sigma %*% gamma  )
  z <- (diag(NROW(y)) - c %*% t(gamma)) %*% y
  VLoVUpVec <- .VLoVUpFN(eta = gamma, Sigma = sigma, A = Abar, b = dbar, z = z)
  criticalVal <- .norminvp_generalized(p = 1-alpha, l = VLoVUpVec[1], u = VLoVUpVec[2],
                                       mu = T_B %*% dtilde, sd = sigmabar)
  reject <- (maxMoment + T_B %*% dtilde > criticalVal)

  return(reject)
}

.testInIdentifiedSet_FLCI_Hybrid <- function(y, sigma, A, d, alpha, hybrid_list) {
  # This function does a hybrid test where we first test if |l'y| > halflength
  # where l and halflength are for the FLCI of size beta
  # If so, we reject in the first stage
  # If not, we add the event that |l'y| <= halflength to the conditioning event,
  # and we adjust second stage size accordingly

  # Note that if y = (betahat - tau), then can derive that $tau \in l'\betahat \pm \chi$ iff |l'y| \leq \chi.
  # Note that this assume l places weight of 1 on \tau

  # Create A_firststage and d_firststage to capture the absolutevalue constraints
  # We reject if A_firstsage %*%y - d_firstage has any positive elements
  # otherwise we add these to the constraints

  A_firststage <- rbind(hybrid_list$flci_l,
                        -hybrid_list$flci_l)
  d_firststage <- c(hybrid_list$flci_halflength,
                    hybrid_list$flci_halflength)

  # Run the first-stage test
  if (max(A_firststage %*% y - d_firststage) > 0) {
    reject <- T
  } else {
    alphatilde <- (alpha - hybrid_list$hybrid_kappa) / (1 - hybrid_list$hybrid_kappa)
    reject <- .testInIdentifiedSet(y = y, sigma = sigma,
                                      A = A, d = d,
                                      Abar_additional = A_firststage,
                                      dbar_additional = d_firststage,
                                      alpha = alphatilde)
  }
  return(reject)
}

.testInIdentifiedSet_LF_Hybrid <- function(y, sigma, A, d, alpha, hybrid_list) {

  sigmaTilde <- as.vector( sqrt( diag(A %*% sigma %*% t(A)) ) )
  Atilde <- solve( diag(sigmaTilde) ) %*% A
  dtilde <- solve( diag(sigmaTilde) ) %*% d

  normalizedMoments <- Atilde %*% y - dtilde
  maxLocation <- which.max(normalizedMoments)
  maxMoment <- normalizedMoments[maxLocation]

  if (maxMoment > hybrid_list$lf_cv) {
    reject = 1
    return(reject)
  } else {
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
}

.testOverThetaGrid <- function(betahat, sigma, A, d, thetaGrid,
                               numPrePeriods, alpha, testFn = .testInIdentifiedSet, ...) {
  # Tests whether values in a grid lie in the identified set.
  testTauInSetFn <- function(theta) {
    reject <- testFn(y = betahat - basisVector(index = numPrePeriods + 1, size = length(betahat))*theta,
                     sigma =  sigma, A = A, d = d, alpha = alpha, ...)
    inSet <- !reject
    return(inSet)
  }
  testResultsGrid <- purrr::map_dbl(.x = thetaGrid, .f = testTauInSetFn)
  testValsGrid <- thetaGrid
  resultsGrid <- cbind(testValsGrid, testResultsGrid)
  return(resultsGrid)
}

.APR_computeCI_NoNuis <- function(betahat, sigma, A, d,
                                  numPrePeriods, numPostPeriods, l_vec,
                                  alpha, returnLength,
                                  hybrid_flag, hybrid_list, grid.ub, grid.lb,
                                  gridPoints, postPeriodMomentsOnly) {
  # This function computes the APR confidence interval for Delta^SD(M) for a given event study.
  # It takes the following inputs:
  #   betahat = vector of estimated event study coefficients
  #   sigma   = covariance matrix of estimated event study coefficients
  #   A       = matrix defining the set Delta
  #   d       = vector defining the set Delta
  #   numPrePeriods = number of pre-periods
  #   numPostPeriods = number of post-periods
  #   M = tuning parameter of Delta^SD(M), default M = 0
  #   postPeriod = post-period of interest
  #   alpha = size of CI, default 0.05.

  # Construct grid of tau values to test and test which values of tau lie in ID set.
  thetaGrid <- seq(grid.lb, grid.ub, length.out = gridPoints)

  if (hybrid_flag == "ARP") {
    resultsGrid <- .testOverThetaGrid(betahat = betahat, sigma = sigma,
                                      numPrePeriods = numPrePeriods,
                                      A = A, d = d, thetaGrid = thetaGrid, alpha = alpha)
  } else if (hybrid_flag == "FLCI") {
    resultsGrid <- .testOverThetaGrid(betahat = betahat, sigma = sigma, thetaGrid = thetaGrid,
                                      A = A, d = d, alpha = alpha,
                                      numPrePeriods = numPrePeriods,
                                      testFn = .testInIdentifiedSet_FLCI_Hybrid,
                                      hybrid_list = hybrid_list)
  } else if (hybrid_flag == "LF") {
    resultsGrid <- .testOverThetaGrid(betahat = betahat, sigma = sigma, thetaGrid = thetaGrid,
                                      A = A, d = d, alpha = alpha,
                                      numPrePeriods = numPrePeriods,
                                      testFn = .testInIdentifiedSet_LF_Hybrid,
                                      hybrid_list = hybrid_list)
  } else {
    stop("hybrid_flag must equal 'APR' or 'FLCI' or 'LF'")
  }

  if( resultsGrid[1] == 1 | resultsGrid[length(resultsGrid)] == 1 ){
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
