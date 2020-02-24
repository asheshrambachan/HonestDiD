# DESCRIPTION =========================================================
# Author: Ashesh Rambachan <asheshr@g.harvard.edu>
#
#  This script contains functions to implement the methods
#  described in Rambachan & Roth (2019) for robust inference
#  in difference-in-differences and event study designs.
#
#  This script contains functions that are used to construct
#  the confidence sets for Delta^{SD}(M).

# PRELIMINARIES =======================================================
library(tidyverse)
library(TruncatedNormal)
library(lpSolveAPI)
library(ROI)
library(Matrix)
library(pracma)
library(CVXR)
library(foreach)

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
  #   hybrid_flag         = flag for hybrid, default = "FLCI"
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

  if (postPeriodMomentsOnly & numPostPeriods > 1){
    postPeriodIndices <- (numPrePeriods +1):NCOL(A_SD)
    postPeriodRows <- which( rowSums( A_SD[ , postPeriodIndices] != 0 ) > 0 )
    rowsForARP <- postPeriodRows
  } else{
    rowsForARP <- 1:NROW(A_SD)
  }

  # Create hybrid_list object
  hybrid_list = list(hybrid_kappa = hybrid_kappa)

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
      if (is.na(grid.ub)){
        grid.ub = (t(flci$optimalVec) %*% betahat) + flci$optimalHalfLength
      }
      if (is.na(grid.lb)){
        grid.lb = (t(flci$optimalVec) %*% betahat) - flci$optimalHalfLength
      }
    } else if (hybrid_flag == "LF") {
      # Compute LF CV
      lf_cv = .compute_least_favorable_cv(X_T = NULL, sigma = A_SD %*% sigma %*% t(A_SD), hybrid_kappa = hybrid_kappa)

      # Store lf cv
      hybrid_list$lf_cv = lf_cv

      # construct theta grid
      if (is.na(grid.ub) & is.na(grid.lb)) {
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
    } else if (hybrid_flag == "ARP") {
      # construct theta grid
      if (is.na(grid.ub) & is.na(grid.lb)) {
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
    } else {
      stop("hybrid_flag must equal 'APR' or 'FLCI' or 'LF'")
    }

    # Compute confidence set
    CI <- .APR_computeCI_NoNuis(betahat = betahat, sigma = sigma,
                                A = A_SD, d = d_SD, numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                l_vec = l_vec, alpha = alpha, returnLength = returnLength,
                                hybrid_flag = hybrid_flag, hybrid_list = hybrid_list,
                                grid.ub = grid.ub, grid.lb = grid.lb,
                                gridPoints = gridPoints)
    return(CI)
  } else { # CASE: NumPostPeriods > 1
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

}
