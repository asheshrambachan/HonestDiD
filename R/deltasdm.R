# DESCRIPTION =========================================================
# Author: Ashesh Rambachan <asheshr@g.harvard.edu>
#
#  This script contains functions to implement the methods
#  described in Rambachan & Roth (2021) for robust inference
#  in difference-in-differences and event study designs.
#
#  This script contains functions that are used to construct
#  the confidence sets for Delta^{SDM}(M).

# DELTA^{SDM}(M) FUNCTIONS --------------------------------------------
# In this section, we implement helper functions to place testing with
# Delta^{SDM}(M) into the form needed to use the ARP functions.

.create_A_SDM <- function(numPrePeriods, numPostPeriods,
                          monotonicityDirection = "increasing", postPeriodMomentsOnly = FALSE) {
  # This function creates a matrix for the linear constraints that \delta \in Delta^SDI(M).
  # It implements this using the general characterization of A.
  #
  # Inputs:
  #   numPrePeriods = number of pre-periods. This is an element of resultsObjects.
  #   numPostPeriods = number of post-periods. This is an element of resultsObjects.
  A_M <- .create_A_M(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                     monotonicityDirection = monotonicityDirection, postPeriodMomentsOnly = postPeriodMomentsOnly)
  A_SD = .create_A_SD(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, postPeriodMomentsOnly =  postPeriodMomentsOnly)

  A = base::rbind(A_SD, A_M)
  base::return(A)
}

.create_d_SDM <- function(numPrePeriods, numPostPeriods, M, postPeriodMomentsOnly = FALSE) {
  # This function creates a vector for the linear constraints that \delta \in Delta^SDI(M).
  # It implements this using the general characterization of d.
  #
  # Inputs:
  #   numPrePeriods  = number of pre-periods. This is an element of resultsObjects.
  #   numPostPeriods = number of post-periods. This is an element of resultsObjects.
  #   M              = smoothness parameter of Delta^SD(M).

  d_SD = .create_d_SD(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, M = M, postPeriodMomentsOnly = postPeriodMomentsOnly)
  d_M = base::rep(0, base::ifelse(postPeriodMomentsOnly, numPostPeriods, numPrePeriods+numPostPeriods) )
  d = c(d_SD, d_M)
  base::return(d)
}

.compute_IDset_DeltaSDM <- function(M, trueBeta, l_vec, numPrePeriods,
                                    numPostPeriods, monotonicityDirection) {

  # Create objective function: Wish to min/max l'delta_post
  fDelta = base::c(base::rep(0, numPrePeriods), l_vec)

  # Create A, d that define Delta^SDM(M)
  A_SDM = .create_A_SDM(numPrePeriods = numPrePeriods,
                        numPostPeriods = numPostPeriods, monotonicityDirection = monotonicityDirection)
  d_SDM = .create_d_SDM(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, M = M)
  dir_SDM = base::rep("<=", base::NROW(A_SDM))

  # Create equality constraint that delta_pre = beta_pre
  prePeriodEqualityMat = base::cbind(base::diag(numPrePeriods),
                                     base::matrix(data = 0, nrow = numPrePeriods, ncol = numPostPeriods))
  A_SDM = base::rbind(A_SDM, prePeriodEqualityMat)
  d_SDM = base::c(d_SDM, trueBeta[1:numPrePeriods])
  dir_SDM = base::c(dir_SDM, base::rep("==", base::NROW(prePeriodEqualityMat)))
  bounds = base::list(lower = base::list(ind = 1:(numPrePeriods + numPostPeriods),
                                         val = base::rep(-Inf, numPrePeriods+numPostPeriods)),
                      upper = base::list(ind = 1:(numPrePeriods + numPostPeriods),
                                         val = base::rep(Inf, numPrePeriods+numPostPeriods)))

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
    base::warning("Solver did not find an optimum")
    id.ub = (base::t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)])
    id.lb = (base::t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)])
  }
  else {
    # Construct upper/lower bound of identified set
    id.ub = (base::t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]) - min.results$optimum
    id.lb = (base::t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]) - max.results$optimum
  }
  # Return identified set
  base::return(tibble::tibble(
    id.lb = id.lb,
    id.ub = id.ub))
}

computeConditionalCS_DeltaSDM <- function(betahat, sigma, numPrePeriods, numPostPeriods,
                                          M = 0, l_vec = .basisVector(index = 1, size = numPostPeriods),
                                          alpha = 0.05, monotonicityDirection = "increasing",
                                          hybrid_flag = "FLCI", hybrid_kappa = alpha/10,
                                          returnLength = FALSE, postPeriodMomentsOnly = TRUE,
                                          gridPoints=10^3, grid.lb = NA, grid.ub = NA, seed = 0) {
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
                        monotonicityDirection = monotonicityDirection, postPeriodMomentsOnly = FALSE)
  d_SDM = .create_d_SDM(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, M = M, postPeriodMomentsOnly = FALSE)

  if (postPeriodMomentsOnly & numPostPeriods > 1) {
    postPeriodIndices <- (numPrePeriods +1):base::NCOL(A_SDM)
    postPeriodRows <- base::which( base::rowSums( A_SDM[ , postPeriodIndices] != 0 ) > 0 )
    rowsForARP <- postPeriodRows
  } else {
    rowsForARP <- 1:base::NROW(A_SDM)
  }

  # Create hybrid_list object
  hybrid_list = base::list(hybrid_kappa = hybrid_kappa)

  # if there is only one post-period, we use the no-nuisance parameter functions
  if (numPostPeriods == 1) {
    if (hybrid_flag == "FLCI") {
      # Compute FLCI
      flci = .findOptimalFLCI_helper(sigma = sigma, M = M,
                                     numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                     l_vec = l_vec, alpha = hybrid_kappa)

      # Add objects to hybrid_list: flci l vector
      hybrid_list$flci_l = flci$optimalVec

      # Add objects to hybrid_list: flci half-length
      hybrid_list$flci_halflength = flci$optimalHalfLength

      # compute FLCI ub and FLCI lb
      if (base::is.na(grid.ub)){
        grid.ub = base::c(base::t(flci$optimalVec) %*% betahat) + flci$optimalHalfLength
      }
      if (base::is.na(grid.lb)){
        grid.lb = base::c(base::t(flci$optimalVec) %*% betahat) - flci$optimalHalfLength
      }
    } else if (hybrid_flag == "LF") {
      # Compute LF CV
      lf_cv = .compute_least_favorable_cv(X_T = NULL, sigma = A_SDM %*% sigma %*% t(A_SDM), hybrid_kappa = hybrid_kappa, seed = seed)

      # Store lf cv
      hybrid_list$lf_cv = lf_cv

      # construct theta grid
      if (base::is.na(grid.ub) & base::is.na(grid.lb)) {
        # Compute identified set under parallel trends
        IDset = .compute_IDset_DeltaSDM(M = M, trueBeta = base::rep(0, numPrePeriods + numPostPeriods),
                                        monotonicityDirection = monotonicityDirection,
                                        l_vec = l_vec, numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods)
        sdTheta <- base::c(base::sqrt(base::t(l_vec) %*% sigma[(numPrePeriods+1):(numPrePeriods+numPostPeriods), (numPrePeriods+1):(numPrePeriods+numPostPeriods)] %*% l_vec))
        if(base::is.na(grid.ub)){
          grid.ub = IDset$id.ub + 20*sdTheta
        }
        if(base::is.na(grid.lb)){
          grid.lb = IDset$id.lb - 20*sdTheta
        }
      }
    } else if (hybrid_flag == "ARP") {
      # construct theta grid
      if (base::is.na(grid.ub) & base::is.na(grid.lb)) {
        # Compute identified set under parallel trends
        IDset = .compute_IDset_DeltaSDM(M = M, trueBeta = rep(0, numPrePeriods + numPostPeriods),
                                        monotonicityDirection = monotonicityDirection,
                                        l_vec = l_vec, numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods)
        sdTheta <- base::c(base::sqrt(base::t(l_vec) %*% sigma[(numPrePeriods+1):(numPrePeriods+numPostPeriods), (numPrePeriods+1):(numPrePeriods+numPostPeriods)] %*% l_vec))
        if(base::is.na(grid.ub)){
          grid.ub = IDset$id.ub + 20*sdTheta
        }
        if(base::is.na(grid.lb)){
          grid.lb = IDset$id.lb - 20*sdTheta
        }
      }
    } else {
      base::stop("hybrid_flag must equal 'APR' or 'FLCI' or 'LF'")
    }

    # Compute confidence set
    CI <- .APR_computeCI_NoNuis(betahat = betahat, sigma = sigma,
                                A = A_SDM, d = d_SDM, numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                l_vec = l_vec, alpha = alpha, returnLength = returnLength,
                                hybrid_flag = hybrid_flag, hybrid_list = hybrid_list,
                                grid.ub = grid.ub, grid.lb = grid.lb,
                                gridPoints = gridPoints)
    base::return(CI)
  } else { # if there are multiple post-periods, we use the nuisance parameter functions
    # HYBRID: If hybrid with FLCI, compute FLCI
    if (hybrid_flag == "FLCI") {
      flci = .findOptimalFLCI_helper(sigma = sigma, M = M, numPrePeriods = numPrePeriods,
                                     numPostPeriods = numPostPeriods, l_vec = l_vec, alpha = hybrid_kappa)

      # Add objects to hybrid_list: flci l vector
      hybrid_list$flci_l = flci$optimalVec

      # Add vbar to flci l vector
      vbar = CVXR::Variable(base::NROW(A_SDM))
      obj <- CVXR::Minimize( base::t(flci$optimalVec) %*% flci$optimalVec -
                         2 * base::t(flci$optimalVec) %*% base::t(A_SDM) %*% vbar + CVXR::quad_form(x = vbar, P = A_SDM %*% base::t(A_SDM)) )
      prob = CVXR::Problem(obj)
      result = CVXR::psolve(prob)
      hybrid_list$vbar = result$getValue(vbar)

      # Add objects to hybrid_list: flci half-length
      hybrid_list$flci_halflength = flci$optimalHalfLength

      # compute FLCI ub and FLCI lb
      grid.ub = base::c(base::t(hybrid_list$flci_l) %*% betahat) + flci$optimalHalfLength
      grid.lb = base::c(base::t(hybrid_list$flci_l) %*% betahat) - flci$optimalHalfLength
    } else {
      # Compute ID set
      IDset = .compute_IDset_DeltaSDM(M = M, trueBeta = base::rep(0, numPrePeriods + numPostPeriods),
                                      monotonicityDirection = monotonicityDirection, l_vec = l_vec,
                                      numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods)
      sdTheta <- base::c(base::sqrt(base::t(l_vec) %*% sigma[(numPrePeriods+1):(numPrePeriods+numPostPeriods), (numPrePeriods+1):(numPrePeriods+numPostPeriods)] %*% l_vec))
      if (base::is.na(grid.ub)) {
        grid.ub = IDset$id.ub + 20*sdTheta
      }
      if (base::is.na(grid.lb)) {
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
    base::return(CI)
    }
}
