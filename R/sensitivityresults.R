# DESCRIPTION =========================================================
# Author: Ashesh Rambachan <asheshr@g.havard.edu>
#
#  This script contains functions to implement the methods
#  described in Rambachan & Roth (2019) for robust inference
#  in difference-in-differences and event study designs.
#
#  Implements functions to perform sensitivity analysis on event study coefficients


# Construct Robust Results Function for Smoothness Restrictions -----------------------------------
createSensitivityResults <- function(betahat, sigma,
                                     numPrePeriods, numPostPeriods,
                                     method = NULL,
                                     Mvec = NULL,
                                     l_vec = .basisVector(index = 1, size = numPostPeriods),
                                     monotonicityDirection = NULL,
                                     biasDirection = NULL,
                                     alpha = 0.05,
                                     parallel = FALSE,
                                     seed = 0) {

  .CustomErrorHandling({
  .stopIfNotConformable(betahat, sigma, numPrePeriods, numPostPeriods, l_vec)
  .warnIfNotSymmPSD(sigma)
  m <- NULL

  # If Mvec is null, construct default Mvec
  if (base::is.null(Mvec)) {
    if (numPrePeriods == 1) {
      Mvec = base::seq(from = 0, to = base::c(base::sqrt(sigma[1, 1])), length.out = 10)
    } else {
      Mub = DeltaSD_upperBound_Mpre(betahat = betahat, sigma = sigma, numPrePeriods = numPrePeriods, alpha = 0.05)
      Mvec = base::seq(from = 0, to = Mub, length.out = 10)
    }
  }

  # If Delta = Delta^{SD}, construct confidence intervals
  if (base::is.null(monotonicityDirection) & base::is.null(biasDirection)) {
    if (base::is.null(method)) {
      method = "FLCI"
    }

    if (method == "FLCI") { # If method = FLCI, construct FLCIs
      if (parallel == FALSE) {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %do% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = alpha, seed = seed)
          tibble::tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                         method = "FLCI", Delta = "DeltaSD", M = Mvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %dopar% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = alpha, seed = seed)
          tibble::tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                         method = "FLCI", Delta = "DeltaSD", M = Mvec[m])
        }
      }
    } else if (method == "Conditional") { # If method = Conditional, construct conditional confidence intervals
      if (parallel == FALSE) {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                              hybrid_flag = "ARP", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "Conditional",
                         Delta = "DeltaSD", M = Mvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                              hybrid_flag = "ARP", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "Conditional",
                         Delta = "DeltaSD", M = Mvec[m])
        }
      }
    } else if (method == "C-F") { # If method = C-F, construct conditional FLCI,
      if (parallel == FALSE) {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                              hybrid_flag = "FLCI", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "C-F",
                         Delta = "DeltaSD", M = Mvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                              hybrid_flag = "FLCI", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "C-F",
                         Delta = "DeltaSD", M = Mvec[m])
        }
      }
    } else if (method == "C-LF") {
      if (parallel == FALSE) {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                              hybrid_flag = "LF", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "C-LF",
                         Delta = "DeltaSD", M = Mvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                              hybrid_flag = "LF", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "C-LF",
                         Delta = "DeltaSD", M = Mvec[m])
        }
      }
    } else{
      base::stop("Method must equal one of: FLCI, Conditional, C-F or C-LF")
    }
  } else if (!base::is.null(biasDirection)) {
    if (base::is.null(method)) {
      method = "C-F"
    }

    if (biasDirection == "positive") {
      Delta = "DeltaSDPB"
    } else {
      Delta = "DeltaSDNB"
    }

    if (method == "FLCI") {
      base::warning("You specified a sign restriction but method = FLCI. The FLCI does not use the sign restriction!")
      if (parallel == FALSE) {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %do% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = alpha, seed = seed)
          tibble::tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                         method = "FLCI", Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %dopar% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = alpha, seed = seed)
          tibble::tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                         method = "FLCI", Delta = Delta, M = Mvec[m])
        }
      }
    } else if (method == "Conditional") { # If method = Conditional, construct conditional confidence intervals
      if (parallel == FALSE) {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "ARP", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "Conditional",
                         Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "ARP", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "Conditional",
                         Delta = Delta, M = Mvec[m])
        }
      }
    } else if (method == "C-F") { # If method = C-F, construct conditional FLCI,
      if (parallel == FALSE) {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "FLCI", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "C-F",
                         Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "FLCI", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "C-F",
                         Delta = Delta, M = Mvec[m])
        }
      }
    } else if (method == "C-LF") {
      if (parallel == FALSE) {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "LF", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "C-LF",
                         Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "LF", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "C-LF",
                         Delta = Delta, M = Mvec[m])
        }
      }
    } else{
      base::stop("Method must equal one of: FLCI, Conditional, C-F or C-LF")
    }
  } else {
    if (base::is.null(method)) {
      method = "C-F"
    }

    if (monotonicityDirection == "increasing") {
      Delta = "DeltaSDI"
    } else {
      Delta = "DeltaSDD"
    }

    if (method == "FLCI") {
      base::warning("You specified a shape restriction but method = FLCI. The FLCI does not use the shape restriction!")
      if (parallel == FALSE) {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %do% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = alpha, seed = seed)
          tibble::tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                         method = "FLCI", Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %dopar% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = alpha, seed = seed)
          tibble::tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                         method = "FLCI", Delta = Delta, M = Mvec[m])
        }
      }
    } else if (method == "Conditional") { # If method = Conditional, construct conditional confidence intervals
      if (parallel == FALSE) {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "ARP", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "Conditional",
                         Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "ARP", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "Conditional",
                         Delta = Delta, M = Mvec[m])
        }
      }
    } else if (method == "C-F") { # If method = C-F, construct conditional FLCI,
      if (parallel == FALSE) {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "FLCI", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "C-F",
                         Delta = Delta,
                         M = Mvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "FLCI", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "C-F",
                         Delta = Delta,
                         M = Mvec[m])
        }
      }
    } else if (method == "C-LF") {
      if (parallel == FALSE) {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "LF", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "C-LF",
                         Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "LF", seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = "C-LF",
                         Delta = Delta, M = Mvec[m])
        }
      }
    } else{
      base::stop("Method must equal one of: FLCI, Conditional, C-F or C-LF")
    }
  }
  base::return(Results)
  }, "Error in computing main results")
}

createSensitivityPlot <- function(robustResults, originalResults, rescaleFactor = 1, maxM = Inf, add_xAxis = TRUE) {
  # Set M for OLS to be the min M in robust results minus the gap between Ms in robust
  Mgap <- base::min( base::diff( base::sort( robustResults$M) ) )
  Mmin <- base::min( robustResults$M)

  originalResults$M <- Mmin - Mgap
  df <- dplyr::bind_rows(originalResults, robustResults)

  # Rescale all the units by rescaleFactor
  df <- df %>% dplyr::mutate_at( c("M", "ub", "lb"), ~ .x * rescaleFactor)

  # Filter out observations above maxM (after rescaling)
  df <- df %>% dplyr::filter(.data$M <= maxM)

  p <- ggplot2::ggplot(data = df, ggplot2::aes(x=.data$M)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lb, ymax = .data$ub, color = base::factor(.data$method)),
                           width = Mgap * rescaleFactor / 2) +
    ggplot2::scale_color_manual(values = base::c("red", '#01a2d9')) +
    ggplot2::theme(legend.title=ggplot2::element_blank(), legend.position="bottom") +
    ggplot2::labs(x = "M", y = "")

  if (add_xAxis) {
    p <- p + ggplot2::geom_hline(yintercept = 0)
  }
  base::return(p)
}

# Construct Robust Results Function for Relative Magnitude Bounds -----------------------------------
createSensitivityResults_relativeMagnitudes <- function(betahat, sigma,
                                                        numPrePeriods, numPostPeriods,
                                                        bound = "deviation from parallel trends",
                                                        method = "C-LF",
                                                        Mbarvec = NULL,
                                                        l_vec = .basisVector(index = 1, size = numPostPeriods),
                                                        monotonicityDirection = NULL,
                                                        biasDirection = NULL,
                                                        alpha = 0.05,
                                                        gridPoints = 10^3,
                                                        grid.ub = NA,
                                                        grid.lb = NA,
                                                        parallel = FALSE,
                                                        seed = 0) {

  .CustomErrorHandling({
  .stopIfNotConformable(betahat, sigma, numPrePeriods, numPostPeriods, l_vec)
  .warnIfNotSymmPSD(sigma)
  m <- NULL

  # If Mbarvec is null, construct default Mbarvec to be 10 values on [0,2].
  if (base::is.null(Mbarvec)) {
    Mbarvec = base::seq(from = 0, to = 2, length.out = 10)
  }

  # Check if bound is specified correctly
  if (bound != "deviation from parallel trends" & bound != "deviation from linear trend") {
    base::stop("bound must equal either 'deviation from parallel trends' or 'deviation from linear trend'.")
  }

  # Check if both monotonicity direction and bias direction are specified:
  if (!base::is.null(monotonicityDirection) & !base::is.null(biasDirection)) {
    base::stop("Please select either a shape restriction or sign restriction (not both).")
  }

  # Use method to specify the choice of confidence set
  if (method == "C-LF") {
    hybrid_flag = "LF"
    method_named = "C-LF"
  } else if (method == "Conditional") {
    hybrid_flag = "ARP"
    method_named = "Conditional"
  } else {
    base::stop("method must be either NULL, Conditional or C-LF.")
  }

  # If bound = "parallel trends violation", we select Delta^{RM} and its variants.
  if (bound == "deviation from parallel trends") {
    # use monotonicity direction and biasDirection to select the choice of Delta
    if (base::is.null(monotonicityDirection) & base::is.null(biasDirection)) { # if both null, Delta^{RM}
      Delta = "DeltaRM"

      if (parallel) {
        Results = foreach::foreach(m = 1:base::length(Mbarvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaRM(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                              hybrid_flag = hybrid_flag,
                                              gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb, seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = method_named,
                         Delta = Delta, Mbar = Mbarvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mbarvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaRM(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                              hybrid_flag = hybrid_flag,
                                              gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb, seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = method_named,
                         Delta = Delta, Mbar = Mbarvec[m])
        }
      }
    } else if (!base::is.null(monotonicityDirection)) { # if monotonicityDirection specified, Delta^{RMM}.
      if (monotonicityDirection == "increasing") {
        Delta = "DeltaRMI"
      } else if (monotonicityDirection == "decreasing") {
        Delta = "DeltaRMD"
      } else {
        base::stop("monotonicityDirection must equal either increasing or decreasing.")
      }

      if (parallel) {
        Results = foreach::foreach(m = 1:base::length(Mbarvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaRMM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = hybrid_flag,
                                               gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb, seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = method_named,
                         Delta = Delta, Mbar = Mbarvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mbarvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaRMM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = hybrid_flag,
                                               gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb, seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = method_named,
                         Delta = Delta, Mbar = Mbarvec[m])
        }
      }
    } else { # if biasDirection is specified, DeltaRMB.
      if (biasDirection == "positive") {
        Delta = "DeltaRMPB"
      } else if (biasDirection == "negative") {
        Delta = "DeltaRMNB"
      } else {
        base::stop("monotonicityDirection must equal either positive or negative.")
      }

      if (parallel) {
        Results = foreach::foreach(m = 1:base::length(Mbarvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaRMB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = hybrid_flag,
                                               gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb, seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = method_named,
                         Delta = Delta, Mbar = Mbarvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mbarvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaRMB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = hybrid_flag,
                                               gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb, seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = method_named,
                         Delta = Delta, Mbar = Mbarvec[m])
        }
      }
    }
  } else { # if bound = "deviation from linear trend", we select Delta^{SDRM} and its variants.

    # Note: Since this choice of Delta^{SDRM} bounds the variation in post-treatment trends,
    # based on observed variation in the pre-treatment trends, we provide an error
    # if the user tries to provide data with only one pre-treatment period.
    if (numPrePeriods == 1) {
      base::stop("Error: not enough pre-periods for 'deviation from linear trend' (Delta^{SDRM} as base choice). Plese see documentation.")
    }

    # use monotonicity direction and biasDirection to select the choice of Delta
    if (base::is.null(monotonicityDirection) & base::is.null(biasDirection)) { # if both null, Delta^{SDRM}
      Delta = "DeltaSDRM"

      if (parallel) {
        Results = foreach::foreach(m = 1:base::length(Mbarvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDRM(betahat = betahat, sigma = sigma,
                                                numPrePeriods = numPrePeriods,
                                                numPostPeriods = numPostPeriods,
                                                l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                                hybrid_flag = hybrid_flag,
                                                gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb, seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = method_named,
                         Delta = Delta, Mbar = Mbarvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mbarvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDRM(betahat = betahat, sigma = sigma,
                                                numPrePeriods = numPrePeriods,
                                                numPostPeriods = numPostPeriods,
                                                l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                                hybrid_flag = hybrid_flag,
                                                gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb, seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = method_named,
                         Delta = Delta, Mbar = Mbarvec[m])
        }
      }
    } else if (!base::is.null(monotonicityDirection)) { # if monotonicityDirection specified, Delta^{SDRMM}.
      if (monotonicityDirection == "increasing") {
        Delta = "DeltaSDRMI"
      } else if (monotonicityDirection == "decreasing") {
        Delta = "DeltaSDRMD"
      } else {
        base::stop("monotonicityDirection must equal either increasing or decreasing.")
      }

      if (parallel) {
        Results = foreach::foreach(m = 1:base::length(Mbarvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDRMM(betahat = betahat, sigma = sigma,
                                                 numPrePeriods = numPrePeriods,
                                                 numPostPeriods = numPostPeriods,
                                                 l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                                 monotonicityDirection = monotonicityDirection,
                                                 hybrid_flag = hybrid_flag,
                                                 gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb, seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = method_named,
                         Delta = Delta, Mbar = Mbarvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mbarvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDRMM(betahat = betahat, sigma = sigma,
                                                 numPrePeriods = numPrePeriods,
                                                 numPostPeriods = numPostPeriods,
                                                 l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                                 monotonicityDirection = monotonicityDirection,
                                                 hybrid_flag = hybrid_flag,
                                                 gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb, seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = method_named,
                         Delta = Delta, Mbar = Mbarvec[m])
        }
      }
    } else { # if biasDirection is specified, DeltaRMB.
      if (biasDirection == "positive") {
        Delta = "DeltaSDRMPB"
      } else if (biasDirection == "negative") {
        Delta = "DeltaSDRMNB"
      } else {
        base::stop("monotonicityDirection must equal either positive or negative.")
      }

      if (parallel) {
        Results = foreach::foreach(m = 1:base::length(Mbarvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDRMB(betahat = betahat, sigma = sigma,
                                                 numPrePeriods = numPrePeriods,
                                                 numPostPeriods = numPostPeriods,
                                                 l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                                 biasDirection = biasDirection,
                                                 hybrid_flag = hybrid_flag,
                                                 gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb, seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = method_named,
                         Delta = Delta, Mbar = Mbarvec[m])
        }
      } else {
        Results = foreach::foreach(m = 1:base::length(Mbarvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDRMB(betahat = betahat, sigma = sigma,
                                                 numPrePeriods = numPrePeriods,
                                                 numPostPeriods = numPostPeriods,
                                                 l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                                 biasDirection = biasDirection,
                                                 hybrid_flag = hybrid_flag,
                                                 gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb, seed = seed)
          tibble::tibble(lb = base::min(temp$grid[temp$accept == 1]),
                         ub = base::max(temp$grid[temp$accept == 1]),
                         method = method_named,
                         Delta = Delta, Mbar = Mbarvec[m])
        }
      }
    }
  }
  # Return tibble of results
  base::return(Results)
  }, "Error in computing main results")
}

createSensitivityPlot_relativeMagnitudes <- function(robustResults, originalResults,
                                                     rescaleFactor = 1, maxMbar = Inf,
                                                     add_xAxis = TRUE) {
  # Set Mbar for OLS to be the min Mbar in robust results minus the gap between Mbars in robust
  Mbargap <- base::min( base::diff( base::sort( robustResults$Mbar) ) )
  Mbarmin <- base::min( robustResults$Mbar)

  originalResults$Mbar <- Mbarmin - Mbargap
  df <- dplyr::bind_rows(originalResults, robustResults)

  # Rescale all the units by rescaleFactor
  df <- df %>% dplyr::mutate_at( base::c("Mbar", "ub", "lb"), ~ .x * rescaleFactor)

  # Filter out observations above maxM (after rescaling)
  df <- df %>% dplyr::filter(.data$Mbar <= maxMbar)

  p <- ggplot2::ggplot(data = df, ggplot2::aes(x=.data$Mbar)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lb, ymax = .data$ub, color = base::factor(.data$method)),
                           width = Mbargap * rescaleFactor / 2) +
    ggplot2::scale_color_manual(values = base::c("red", '#01a2d9')) +
    ggplot2::theme(legend.title=ggplot2::element_blank(), legend.position="bottom") +
    ggplot2::labs(x = latex2exp::TeX("$\\bar{M}$"), y = "")

  if (add_xAxis) {
    p <- p + ggplot2::geom_hline(yintercept = 0)
  }
  base::return(p)
}

# Construct Original CS -----------------------------------------------
constructOriginalCS <- function(betahat, sigma,
                                numPrePeriods, numPostPeriods,
                                l_vec = .basisVector(index = 1, size = numPostPeriods),
                                alpha = 0.05) {
  .stopIfNotConformable(betahat, sigma, numPrePeriods, numPostPeriods, l_vec)
  .warnIfNotSymmPSD(sigma)
  stdError = base::sqrt(base::t(l_vec) %*% sigma[(numPrePeriods+1):(numPrePeriods+numPostPeriods), (numPrePeriods+1):(numPrePeriods+numPostPeriods)] %*% l_vec)
  lb = base::t(l_vec) %*% betahat[(numPrePeriods+1):(numPrePeriods+numPostPeriods)] - stats::qnorm(1-alpha/2)*stdError
  ub = base::t(l_vec) %*% betahat[(numPrePeriods+1):(numPrePeriods+numPostPeriods)] + stats::qnorm(1-alpha/2)*stdError
  base::return(tibble::tibble(
    lb = lb,
    ub = ub,
    method = "Original",
    Delta = NA
  ))
}

createEventStudyPlot <- function(betahat, stdErrors = NULL, sigma = NULL,
                                 numPrePeriods, numPostPeriods, alpha = 0.05,
                                 timeVec, referencePeriod,
                                 useRelativeEventTime = FALSE) {
  if (base::is.null(stdErrors) & base::is.null(sigma)) {
    base::stop("User must specify either vector of standard errors or vcv matrix!")
  } else if (base::is.null(stdErrors) & !is.null(sigma)) {
    stdErrors = base::sqrt(base::diag(sigma))
  }

  if (useRelativeEventTime == TRUE) {
    timeVec = timeVec - referencePeriod
    referencePeriod = 0
  }

  EventStudyPlot <- ggplot2::ggplot(tibble::tibble(t    = base::c(timeVec[1:numPrePeriods], referencePeriod, timeVec[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]),
                                                   beta = base::c(betahat[1:numPrePeriods], 0, betahat[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]),
                                                   se   = base::c(stdErrors[1:numPrePeriods], NA, stdErrors[(numPrePeriods+1):(numPrePeriods+numPostPeriods)])),
                                    ggplot2::aes(x = t)) +
    ggplot2::geom_point(ggplot2::aes(y = beta), color = "red") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = beta - stats::qnorm(1-alpha/2)*.data$se, ymax = beta + stats::qnorm(1-alpha/2)*.data$se), width = 0.5, colour = "#01a2d9") +
    ggplot2::theme(legend.position = "none") + 
    ggplot2::labs(x = "Event time", y = "") +
    ggplot2::scale_x_continuous(breaks = base::seq(from = base::min(timeVec), to = base::max(timeVec), by = 1),
                                labels = base::as.character(base::seq(from = base::min(timeVec), to = base::max(timeVec), by = 1)))
  base::return(EventStudyPlot)
}
