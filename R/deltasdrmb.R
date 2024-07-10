# DESCRIPTION =========================================================
# Author: Ashesh Rambachan <asheshr@g.harvard.edu>
#
#  This script contains functions to implement the methods
#  described in Rambachan & Roth (2021) for robust inference
#  in difference-in-differences and event study designs.
#
#  This script contains functions that are used to construct
#  the confidence sets for Delta^{SDRM}(Mbar), which intersects Delta^{SDRM}(Mbar)
#  with a sign restriction (i.e., Delta^{PB} or Delta^{NB}). See the discussion
#  Section 2.3.3 of Rambachan and Roth (2021).

# DELTA^{SDRMB}(Mbar) FUNCTIONS ---------------------------------------------
# In this section, we implement helper functions to place testing with
# Delta^{SDRMB}(Mbar) into the form needed to use the ARP functions.
.create_A_SDRMB <- function(numPrePeriods, numPostPeriods,
                           Mbar = 1, s, max_positive = TRUE,
                           dropZero = TRUE, biasDirection) {
  # This function creates a matrix for the linear constraints that \delta \in Delta^SDRMB_{s,.}(Mbar), where . is + if max_positve = TRUE and - if max_positive = FALSE.
  #
  # Inputs:
  #   numPrePeriods = number of pre-periods. This is an element of resultsObjects.
  #   numPostPeriods = number of post-periods. This is an element of resultsObjects.
  #   biasDirection = must be either positive or negative.

  # First construct matrix Atilde used in DeltaSD -- (numPrePeriods+numPostPeriods-1) x (numPrePeriods+numPostPeriods+1)
  # Note Atilde is just the positive moments; is not related to Atilde, the rotate matrix, in the paper
  # Note: Atilde initially includes t = 0. We then drop it.
  Atilde = base::matrix(0, nrow = numPrePeriods+numPostPeriods-1, ncol = numPrePeriods+numPostPeriods+1)
  for (r in 1:(numPrePeriods+numPostPeriods-1)) {
    Atilde[r, r:(r+2)] = c(1, -2, 1)
  }

  # Create a vector to extract the max second dif, which corresponds with the second dif for period s, or minus this if max_positive == FALSE
  v_max_dif <- base::matrix(0, nrow = 1, ncol = numPrePeriods + numPostPeriods + 1)
  v_max_dif[(numPrePeriods+1+s-2):(numPrePeriods+1+s)] <- c(1,-2, 1)

  if(max_positive == FALSE){
    v_max_dif <- -v_max_dif
  }

  # The bounds for the 2nd dif starting with period t are 1*v_max_dif if t<=0 and M*v_max_dif if t>0
  A_UB <- base::rbind( pracma::repmat(v_max_dif, n=numPrePeriods-1, m = 1),
                       pracma::repmat(Mbar*v_max_dif, n=numPostPeriods, m = 1))

  # Construct A that imposes |Atilde * delta | <= A_UB * delta
  A = base::rbind(Atilde - A_UB, -Atilde - A_UB)

  # Remove all-zero rows of the matrix Atilde, corresponding with the constraint (delta_s - delta_s-1) - (delta_s - delta_s-1) <= (delta_s - delta_s-1) - (delta_s - delta_s-1)
  zerorows <- base::apply(X = A, MARGIN = 1, FUN = function(x) t(x) %*% x) <= 10^-10
  A <- A[!zerorows, ]
  # Construct the sign restriction matrix
  A_B <- .create_A_B(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                     biasDirection = biasDirection)

  # Remove the period corresponding with t=0
  if (dropZero) {
    A = A[, -(numPrePeriods+1)]
    # Bind rows of A for Delta^{RM}_s with A_B and return
    base::return(base::rbind(A, A_B))
  } else {
    # Bind rows of A for Delta^{RM}_s with A_B and return
    base::return(base::rbind(A, A_B))
  }
}

.create_d_SDRMB <- function(numPrePeriods, numPostPeriods, dropZero = TRUE){
  # This function creates a vector for the linear constraints that delt is in Delta^SDRM_{s,.}(Mbar), where . is + if max_positve = TRUE and - if max_positive = FALSE.
  # It implements this using the general characterization of d, NOT the sharp
  # characterization of the identified set.
  #
  # Inputs:
  #   numPrePeriods  = number of pre-periods. This is an element of resultsObjects.
  #   numPostPeriods = number of post-periods. This is an element of resultsObjects.

  A_SDRM = .create_A_SDRM(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                         Mbar = 0, s = 0, dropZero = dropZero) # d doesn't depend on Mbar or s
  d_SDRM = base::rep(0, base::NROW(A_SDRM))
  d_B = base::rep(0, numPostPeriods)
  d = base::c(d_SDRM, d_B)
  base::return(d)
}

# DELTA^{SDRMB}(Mbar) Identified Set Helper Functions --------------------
.compute_IDset_DeltaSDRMB_fixedS <- function(s, Mbar, max_positive,
                                             trueBeta, l_vec, numPrePeriods, numPostPeriods, biasDirection) {
  # This helper function computes the upper and lower bound of the identified set
  # given the event study coefficients, lvec and Mbar. It computes the identified
  # set for a user-specified choice of s and (+), (-). This is used by
  # the function compute_IDset_DeltaSDRMB below.

  # Create objective function: Wish to min/max l'delta_post
  fDelta = base::c(base::rep(0, numPrePeriods), l_vec)

  # Create A_SDRMB, d_SDRMB for this choice of s, max_positive
  A_SDRMB_s = .create_A_SDRMB(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                              Mbar = Mbar, s = s, max_positive = max_positive, biasDirection = biasDirection)
  d_SDRMB = .create_d_SDRMB(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods)

  # Create vector for direction of inequalities associated with MB
  dir_SDRMB = base::rep("<=", base::length(d_SDRMB))

  # Add equality constraint for pre-period coefficients
  prePeriodEqualityMat = base::cbind(base::diag(numPrePeriods),
                                     base::matrix(data = 0, nrow = numPrePeriods, ncol = numPostPeriods))
  A_SDRMB_s = base::rbind(A_SDRMB_s, prePeriodEqualityMat)
  d_SDRMB = base::c(d_SDRMB, trueBeta[1:numPrePeriods])
  dir_SDRMB = base::c(dir_SDRMB, base::rep("==", base::NROW(prePeriodEqualityMat)))

  # Specify variables between (-inf, inf)
  bounds = base::list(lower = base::list(ind = 1:(numPrePeriods + numPostPeriods), val = base::rep(-Inf, numPrePeriods+numPostPeriods)),
                      upper = base::list(ind = 1:(numPrePeriods + numPostPeriods), val = base::rep(Inf, numPrePeriods+numPostPeriods)))

  # Create and solve for max
  results.max = Rglpk::Rglpk_solve_LP(obj = fDelta,
                                      max = TRUE,
                                      mat = A_SDRMB_s,
                                      dir = dir_SDRMB,
                                      rhs = d_SDRMB,
                                      bounds = bounds)

  # Create and solve for min
  results.min = Rglpk::Rglpk_solve_LP(obj = fDelta,
                                      max = FALSE,
                                      mat = A_SDRMB_s,
                                      dir = dir_SDRMB,
                                      rhs = d_SDRMB,
                                      bounds = bounds)

  if (results.max$status != 0 & results.min$status != 0) {
      # If the solver does not return solution, we just return the l_vec'trueBeta.
      id.ub = (base::t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)])
      id.lb = (base::t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)])
  }
  else {
    # Construct upper/lower bound of identified set
    id.ub = (base::t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]) - results.min$optimum
    id.lb = (base::t(l_vec) %*% trueBeta[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]) - results.max$optimum
  }
  base::return(
    tibble::tibble(id.lb = id.lb, id.ub = id.ub)
  )
}

.compute_IDset_DeltaSDRMB <- function(Mbar, trueBeta, l_vec, numPrePeriods, numPostPeriods, biasDirection) {
  # This function computes the upper and lower bound of the identified set
  # given the event study coefficients, lvec and Mbar.
  #
  # To do so, we construct the identified set at each choice of s, +, -. We
  # then take the union of these intervals.
  #
  # Note: lvec is assumed to be non-negative.
  #
  # Inputs:
  #   Mbar          = smoothness param of Delta^SDRM
  #   trueBeta       = vector of population event study coefficients
  #   l_vec          = vector l defining parameter of interest
  #   numPrePeriods  = number of pre-periods
  #   numPostPeriods = number of post-periods
  #
  # Outputs:
  #   dataframe with columns
  #     id.ub = upper bound of ID set
  #     id.lb = lower bound of ID set

  # Construct identified sets for (+) at each value of s
  min_s = -(numPrePeriods - 2)
  id_bounds_plus = purrr::map_dfr(
    .x = min_s:0,
    .f = ~.compute_IDset_DeltaSDRMB_fixedS(s = .x, Mbar = Mbar, max_positive = TRUE,
                                           trueBeta = trueBeta, l_vec = l_vec,
                                           numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                           biasDirection = biasDirection)
  )
  id_bounds_minus = purrr::map_dfr(
    .x = min_s:0,
    .f = ~.compute_IDset_DeltaSDRMB_fixedS(s = .x, Mbar = Mbar, max_positive = FALSE,
                                           trueBeta = trueBeta, l_vec = l_vec,
                                           numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                           biasDirection = biasDirection)
  )

  # Construct the identified set by taking the max of the upper bound and the min of the lower bound
  id.lb = base::min(base::min(id_bounds_plus$id.lb), base::min(id_bounds_minus$id.lb))
  id.ub = base::max(base::max(id_bounds_plus$id.ub), base::max(id_bounds_minus$id.ub))

  # Return identified set
  base::return(tibble::tibble(
    id.lb = id.lb,
    id.ub = id.ub))
}

# Delta^{SDRMB}(Mbar) Inference Helper Functions -------------------------
.computeConditionalCS_DeltaSDRMB_fixedS <- function(s, max_positive, Mbar,
                                                    betahat, sigma, numPrePeriods, numPostPeriods, l_vec,
                                                    alpha, hybrid_flag, hybrid_kappa, biasDirection,
                                                    postPeriodMomentsOnly, gridPoints, grid.ub, grid.lb, seed = 0) {
  # This function computes the ARP CI that includes nuisance parameters
  # for Delta^{SDRMB}(Mbar) for a fixed s and (+),(-). This functions uses ARP_computeCI for all
  # of its computations. It is used as a helper function in computeConditionalCS_DeltaSDRMB below.

  # Check that hybrid_flag equals LF or ARP
  if (hybrid_flag != "LF" & hybrid_flag != "ARP") {
    base::stop("hybrid_flag must equal 'ARP' or 'FLCI' or 'LF'")
  }

  # Create hybrid_list object
  hybrid_list = base::list(hybrid_kappa = hybrid_kappa)

  # Create matrix A_SDRM_s, and vector d_SDRM
  A_SDRMB_s = .create_A_SDRMB(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                             Mbar = Mbar, s = s, max_positive = max_positive, biasDirection = biasDirection)
  d_SDRMB = .create_d_SDRMB(numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods)

  # If only use post period moments, construct indices for the post period moments only.
  if (postPeriodMomentsOnly){
    if(numPostPeriods > 1){
      postPeriodIndices <- (numPrePeriods +1):base::NCOL(A_SDRMB_s)
      postPeriodRows <- base::which( base::rowSums( A_SDRMB_s[ , postPeriodIndices] != 0 ) > 0 )
      rowsForARP <- postPeriodRows
    }else{
      #If only one post-period, then it is the last column
      postPeriodRows <- base::which(A_SDRMB_s[ ,NCOL(A_SDRMB_s)] != 0 )
      A_SDRMB_s <- A_SDRMB_s[postPeriodRows, ]
      d_SDRMB <- d_SDRMB[postPeriodRows]

    }
  } else{
    rowsForARP <- 1:base::NROW(A_SDRMB_s)
  }

  # if there is only one post-period, we use the no-nuisance parameter functions
  if (numPostPeriods == 1) {
    if (hybrid_flag == "LF") {
      # Compute LF CV and store it in hybrid_list
      lf_cv = .compute_least_favorable_cv(X_T = NULL, sigma = A_SDRMB_s %*% sigma %*% base::t(A_SDRMB_s), hybrid_kappa = hybrid_kappa, seed = seed)
      hybrid_list$lf_cv = lf_cv
    }
    # Compute confidence set
    CI <- .APR_computeCI_NoNuis(betahat = betahat, sigma = sigma,
                                A = A_SDRMB_s, d = d_SDRMB,
                                numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                l_vec = l_vec, alpha = alpha, returnLength = FALSE,
                                hybrid_flag = hybrid_flag, hybrid_list = hybrid_list,
                                grid.ub = grid.ub, grid.lb = grid.lb,
                                gridPoints = gridPoints)
  } else { # CASE: NumPostPeriods > 1
    # Compute ARP CI for l'beta using Delta^SDRMB
    CI = .ARP_computeCI(betahat = betahat, sigma = sigma, numPrePeriods = numPrePeriods,
                        numPostPeriods = numPostPeriods, A = A_SDRMB_s, d = d_SDRMB,
                        l_vec = l_vec, alpha = alpha,
                        hybrid_flag = hybrid_flag, hybrid_list = hybrid_list,
                        returnLength = FALSE,
                        grid.lb = grid.lb, grid.ub = grid.ub,
                        gridPoints = gridPoints, rowsForARP = rowsForARP)
  }
  base::return(CI)
}

computeConditionalCS_DeltaSDRMB <- function(betahat, sigma, numPrePeriods, numPostPeriods,
                                            l_vec = .basisVector(index = 1, size = numPostPeriods),
                                            Mbar = 0, alpha = 0.05, hybrid_flag = "LF",
                                            hybrid_kappa = alpha/10, returnLength = FALSE,
                                            postPeriodMomentsOnly = TRUE, biasDirection = "positive",
                                            gridPoints = 10^3, grid.ub = NA, grid.lb = NA, seed = 0) {
  # This function computes the ARP CI that includes nuisance parameters
  # for Delta^{SDRMB}(Mbar). This functions uses ARP_computeCI for all
  # of its computations.
  #
  # Inputs:
  #   betahat             = vector of estimated event study coefficients
  #   sigma               = covariance matrix of estimated event study coefficients
  #   numPrePeriods       = number of pre-periods
  #   numPostPeriods      = number of post-periods
  #   l_vec               = vector that defines parameter of interest
  #   Mbar                = tuning parameter for Delta^{SDRM}(Mbar), default Mbar = 0.
  #   alpha               = desired size of CI, default alpha = 0.05
  #   hybrid_flag         = flag for hybrid, default = "LF". Must be either "LF" or "ARP"
  #   hybrid_kappa        = desired size of first-stage hybrid test, default = NULL
  #   returnLength        = returns length of CI only. Otherwise, returns matrix with grid in col 1 and test result in col 2.
  #   numGridPoints       = number of gridpoints to test over, default = 1000
  #   postPeriodMomentsOnly = exclude moments for delta^SDRM that only include pre-period coefs
  #
  #  Outputs:
  #   data_frame containing upper and lower bounds of CI.

  # Note: Since this choice of Delta^{SDRM} bounds the variation in post-treatment trends,
  # based on observed variation in the pre-treatment trends, we provide an error
  # if the user tries to provide data with only one pre-treatment period.
  if (numPrePeriods == 1) {
    base::stop("Error: not enough pre-periods (Delta^{SDRM} as base choice)! Plese see documentation.")
  }

  # Create minimal s index for looping.
  min_s = -(numPrePeriods - 2)
  s_indices = min_s:0

  # Construct theta grid by computing id set under parallel trends.
  # The default sets the grid to be equal to [-20*sdTheta, 20*sdTheta]
  sdTheta <- base::c(base::sqrt(base::t(l_vec) %*% sigma[(numPrePeriods+1):(numPrePeriods+numPostPeriods), (numPrePeriods+1):(numPrePeriods+numPostPeriods)] %*% l_vec))
  if (base::is.na(grid.ub)) { grid.ub = 20*sdTheta }
  if (base::is.na(grid.lb)) { grid.lb = -20*sdTheta }

  # Loop over s values for (+), (-), left join the resulting CIs based on the grid
  CIs_SDRMB_plus_allS = base::matrix(0, nrow = gridPoints, ncol = base::length(s_indices))
  CIs_SDRMB_minus_allS = base::matrix(0, nrow = gridPoints, ncol = base::length(s_indices))
  for (s_i in 1:base::length(s_indices)) {
    # Compute CI for s, (+) and bind it to all CI's for (+)
    CI_s_plus = .computeConditionalCS_DeltaSDRMB_fixedS(s = s_indices[s_i], max_positive = TRUE, Mbar = Mbar,
                                                        betahat = betahat, sigma = sigma, numPrePeriods = numPrePeriods,
                                                        numPostPeriods = numPostPeriods, l_vec = l_vec,
                                                        alpha = alpha, hybrid_flag = hybrid_flag, hybrid_kappa = hybrid_kappa,
                                                        postPeriodMomentsOnly = postPeriodMomentsOnly, biasDirection = biasDirection,
                                                        gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb, seed = seed)
    CIs_SDRMB_plus_allS[,s_i] = CI_s_plus$accept

    # Compute CI for s, (-) and bind it to all CI's for (-)
    CI_s_minus = .computeConditionalCS_DeltaSDRMB_fixedS(s = s_indices[s_i], max_positive = FALSE, Mbar = Mbar,
                                                        betahat = betahat, sigma = sigma, numPrePeriods = numPrePeriods,
                                                        numPostPeriods = numPostPeriods, l_vec = l_vec,
                                                        alpha = alpha, hybrid_flag = hybrid_flag, hybrid_kappa = hybrid_kappa,
                                                        postPeriodMomentsOnly = postPeriodMomentsOnly, biasDirection = biasDirection,
                                                        gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb, seed = seed)
    CIs_SDRMB_minus_allS[,s_i] = CI_s_minus$accept
  }
  CIs_SDRMB_plus_maxS = base::apply(CIs_SDRMB_plus_allS, MARGIN = 1, FUN = base::max)
  CIs_SDRMB_minus_maxS = base::apply(CIs_SDRMB_minus_allS, MARGIN = 1, FUN = base::max)

  # Take the max between (+), (-) and Construct grid containing theta points and whether any CI accepted
  CI_SDRMB = tibble::tibble(grid   = base::seq(grid.lb, grid.ub, length.out = gridPoints),
                            accept = base::pmax(CIs_SDRMB_plus_maxS, CIs_SDRMB_minus_maxS))

  # Compute length, else return grid
  if (returnLength == TRUE) {
    gridLength <- 0.5 * ( base::c(0, base::diff(CI_SDRMB$grid)) + base::c(base::diff(CI_SDRMB$grid), 0 ) )
    base::return(base::sum(CI_SDRMB$accept*gridLength))
  } else {
    base::return(CI_SDRMB)
  }
}
