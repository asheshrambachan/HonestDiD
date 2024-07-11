# DESCRIPTION =========================================================
# Author: Ashesh Rambachan <asheshr@g.harvard.edu>
#
#  This script contains functions to implement the methods
#  described in Rambachan & Roth (2019) for robust inference
#  in difference-in-differences and event study designs.
#
#  This script contains functions that are used to construct
#  the confidence sets for Delta^{SDB}(M).

# DELTA^{SDB}(M) FUNCTIONS --------------------------------------------
# In this section, we implement helper functions to place testing with
# Delta^{SDB}(M) into the form needed to use the ARP functions.
.create_A_SDB <- function(numPrePeriods, numPostPeriods,
                          biasDirection = "positive", postPeriodMomentsOnly = FALSE) {
  # This function creates a matrix for the linear constraints that \delta \in Delta^SDB(M).
  # It implements this using the general characterization of A.
  #
  # Inputs:
  #   numPrePeriods = number of pre-periods. This is an element of resultsObjects.
  #   numPostPeriods = number of post-periods. This is an element of resultsObjects.

  A_SD = .create_A_SD(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                      postPeriodMomentsOnly = postPeriodMomentsOnly)
  A_B = .create_A_B(numPrePeriods = numPrePeriods,
                    numPostPeriods = numPostPeriods, biasDirection = biasDirection)

  A = base::rbind(A_SD, A_B)
  base::return(A)
}

.create_d_SDB <- function(numPrePeriods, numPostPeriods, M, postPeriodMomentsOnly = FALSE) {
  # This function creates a vector for the linear constraints that \delta \in Delta^SDB(M).
  # It implements this using the general characterization of d.
  #
  # Inputs:
  #   numPrePeriods  = number of pre-periods. This is an element of resultsObjects.
  #   numPostPeriods = number of post-periods. This is an element of resultsObjects.
  #   M              = smoothness parameter of Delta^SD(M).

  d_SD = .create_d_SD(numPrePeriods = numPrePeriods,
                      numPostPeriods = numPostPeriods, M = M, postPeriodMomentsOnly = postPeriodMomentsOnly)
  d_B = base::rep(0, numPostPeriods)
  d = base::c(d_SD, d_B)
  base::return(d)
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
  fDelta = base::c(base::rep(0, numPrePeriods), l_vec)

  # Create A, d that define Delta^SDPB(M)
  A_SDB = .create_A_SDB(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                        biasDirection = biasDirection)
  d_SDB = .create_d_SDB(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, M = M)
  dir_SDB = base::rep("<=", base::NROW(A_SDB))

  # Create equality constraint that delta_pre = beta_pre
  prePeriodEqualityMat = base::cbind(base::diag(numPrePeriods),
                                     base::matrix(data = 0, nrow = numPrePeriods, ncol = numPostPeriods))
  A_SDB = base::rbind(A_SDB, prePeriodEqualityMat)
  d_SDB = base::c(d_SDB, trueBeta[1:numPrePeriods])
  dir_SDB = base::c(dir_SDB, base::rep("==", base::NROW(prePeriodEqualityMat)))
  bounds = base::list(lower = base::list(ind = 1:(numPrePeriods + numPostPeriods), val = base::rep(-Inf, numPrePeriods+numPostPeriods)),
                upper = base::list(ind = 1:(numPrePeriods + numPostPeriods), val = base::rep(Inf, numPrePeriods+numPostPeriods)))


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
    base::warning("Solver did not find an optimum")
    id.ub = (base::t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)])
    id.lb = (base::t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)])
  }
  else {
    # Construct upper/lower bound of identified set
    id.ub = (base::t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]) - results.min$optimum
    id.lb = (base::t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]) - results.max$optimum
  }
  # Return identified set
  base::return(tibble::tibble(
    id.lb = id.lb,
    id.ub = id.ub))
}

computeConditionalCS_DeltaSDB <- function(betahat, sigma, numPrePeriods, numPostPeriods,
                                          M = 0, l_vec = .basisVector(index = 1, size=numPostPeriods),
                                          alpha = 0.05, hybrid_flag = "FLCI", hybrid_kappa = alpha/10,
                                          returnLength = FALSE, biasDirection = "positive",
                                          postPeriodMomentsOnly = TRUE,
                                          gridPoints = 10^3, grid.lb = NA, grid.ub = NA, seed = 0) {
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
  #   biasDirection       = "positive" or "negative"; specifies direction of sign restriction.
  #   postPeriodMomentsOnly = include the delta^RM moments corresponding with the first period only
  #
  #  Outputs:
  #   data_frame containing upper and lower bounds of CI.

  # Construct A_SDB, d_SDB
  A_SDB = .create_A_SDB(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                        biasDirection = biasDirection,
                        postPeriodMomentsOnly = FALSE)
  d_SDB = .create_d_SDB(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods, M = M, postPeriodMomentsOnly = FALSE)

  if (postPeriodMomentsOnly  & numPostPeriods > 1) {
    postPeriodIndices <- (numPrePeriods +1):base::NCOL(A_SDB)
    postPeriodRows <- base::which( base::rowSums( A_SDB[ , postPeriodIndices] != 0 ) > 0 )
    rowsForARP <- postPeriodRows
  } else{
    rowsForARP <- 1:base::NROW(A_SDB)
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
      lf_cv = .compute_least_favorable_cv(X_T = NULL, sigma = A_SDB %*% sigma %*% base::t(A_SDB), hybrid_kappa = hybrid_kappa, seed = seed)

      # Store lf cv
      hybrid_list$lf_cv = lf_cv

      # construct theta grid
      if (base::is.na(grid.ub) & base::is.na(grid.lb)) {
        # Compute identified set under parallel trends
        IDset = .compute_IDset_DeltaSDB(M = M, trueBeta = base::rep(0, numPrePeriods + numPostPeriods),
                                        biasDirection = biasDirection,
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
        IDset = .compute_IDset_DeltaSDB(M = M, trueBeta = base::rep(0, numPrePeriods + numPostPeriods),
                                        biasDirection = biasDirection,
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
                                A = A_SDB, d = d_SDB, numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                l_vec = l_vec, alpha = alpha, returnLength = returnLength,
                                hybrid_flag = hybrid_flag, hybrid_list = hybrid_list,
                                grid.ub = grid.ub, grid.lb = grid.lb,
                                gridPoints = gridPoints)
    base::return(CI)
  } else { # If multiple post-periods, we use the nuisance parameter functions
    # HYBRID: If hybrid, compute FLCI
    if (hybrid_flag == "FLCI") {
      flci = .findOptimalFLCI_helper(sigma = sigma, M = M, numPrePeriods = numPrePeriods,
                                     numPostPeriods = numPostPeriods, l_vec = l_vec, alpha = hybrid_kappa)

      # Add objects to hybrid_list: flci l vector
      hybrid_list$flci_l = flci$optimalVec

      # Add vbar to flci l vector
      vbar = CVXR::Variable(base::NROW(A_SDB))
      obj <- CVXR::Minimize( base::t(flci$optimalVec) %*% flci$optimalVec -
                         2 * base::t(flci$optimalVec) %*% base::t(A_SDB) %*% vbar + CVXR::quad_form(x = vbar, P = A_SDB %*% base::t(A_SDB)) )
      prob = CVXR::Problem(obj)
      result = CVXR::psolve(prob)
      hybrid_list$vbar = result$getValue(vbar)

      # Add objects to hybrid_list: flci half-length
      hybrid_list$flci_halflength = flci$optimalHalfLength

      # compute FLCI ub and FLCI lb
      if (base::is.na(grid.ub)){
        grid.ub = base::c(base::t(flci$optimalVec) %*% betahat) + flci$optimalHalfLength
      }
      if (base::is.na(grid.lb)){
        grid.lb = base::c(base::t(flci$optimalVec) %*% betahat) - flci$optimalHalfLength
      }
    } else {
      # Compute ID set under parallel trends
      IDset = .compute_IDset_DeltaSDB(M = M, trueBeta = base::rep(0, numPrePeriods + numPostPeriods), biasDirection = biasDirection,
                                      l_vec = l_vec, numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods)
      # If direction is negative, flip the signs and the upper and lower bounds
      if(biasDirection == "negative"){
        new.lb <- - IDset$id.ub
        new.ub <- - IDset$id.lb
        IDset$id.lb <- new.lb
        IDset$id.ub <- new.ub
      }
      sdTheta <- base::c(base::sqrt(base::t(l_vec) %*% sigma[(numPrePeriods+1):(numPrePeriods+numPostPeriods), (numPrePeriods+1):(numPrePeriods+numPostPeriods)] %*% l_vec))

      if (base::is.na(grid.ub)) {
        grid.ub = IDset$id.ub + 20*sdTheta
      }
      if (base::is.na(grid.lb)) {
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
    base::return(CI)
  }
}
