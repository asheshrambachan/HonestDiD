# DESCRIPTION =========================================================
# Author: Ashesh Rambachan <asheshr@g.harvard.edu>
#
#  This script contains functions to implement the methods
#  described in Rambachan & Roth (2019) for robust inference
#  in difference-in-differences and event study designs.
#
#  This script contains functions that are used to construct
#  the confidence sets for Delta^{RMI}(M).

# PRELIMINARIES =======================================================
library(tidyverse)
library(TruncatedNormal)
library(lpSolveAPI)
library(ROI)
library(Matrix)
library(pracma)
library(CVXR)
library(foreach)

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

  if(postPeriodMomentsOnly & numPostPeriods > 1){
    postPeriodIndices <- (numPrePeriods +1):NCOL(A_RMI)
    postPeriodRows <- which( rowSums( A_RMI[ , postPeriodIndices] != 0 ) > 0 )
    rowsForARP <- postPeriodRows
  } else{
    rowsForARP <- 1:NROW(A_RMI)
  }

  # Create hybrid_list object
  hybrid_list = list(hybrid_kappa = hybrid_kappa)

  # If only one post-period, we use the no-nuisance parameter functions
  if (numPostPeriods == 1) {
    if (hybrid_flag == "FLCI") {
      stop("The FLCI hybrid is not an option for Delta^{RMI}!")
    } else if (hybrid_flag == "LF") {
      # Compute LF CV
      lf_cv = .compute_least_favorable_cv(X_T = NULL, sigma = A_RMI %*% sigma %*% t(A_RMI), hybrid_kappa = hybrid_kappa)

      # Store lf cv
      hybrid_list$lf_cv = lf_cv

      # construct theta grid
      if (is.na(grid.ub) & is.na(grid.lb)) {
        # Compute identified set under parallel trends
        IDset = .compute_IDset_DeltaRMI(Mbar = Mbar, trueBeta = rep(0, numPrePeriods + numPostPeriods),
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
        IDset = .compute_IDset_DeltaRMI(Mbar = Mbar, trueBeta = rep(0, numPrePeriods + numPostPeriods),
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
                                A = A_RMI, d = d_RMI, numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                l_vec = l_vec, alpha = alpha, returnLength = returnLength,
                                hybrid_flag = hybrid_flag, hybrid_list = hybrid_list,
                                grid.ub = grid.ub, grid.lb = grid.lb,
                                gridPoints = gridPoints)
    return(CI)
  } else { # If there are multiple post-periods, we use the implementation with nuisance parameters
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
}
