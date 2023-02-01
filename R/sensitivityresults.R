# DESCRIPTION =========================================================
# Author: Ashesh Rambachan <asheshr@g.havard.edu>
#
#  This script contains functions to implement the methods
#  described in Rambachan & Roth (2019) for robust inference
#  in difference-in-differences and event study designs.
#
#  Implements functions to perform sensitivity analysis on event study coefficients


# PRELIMINARIES =======================================================
library(tidyverse)
library(TruncatedNormal)
library(lpSolveAPI)
library(ROI)
library(Matrix)
library(pracma)
library(CVXR)
library(foreach)

# Construct Robust Results Function for Smoothness Restrictions -----------------------------------
createSensitivityResults <- function(betahat, sigma,
                                     numPrePeriods, numPostPeriods,
                                     method = NULL,
                                     Mvec = NULL,
                                     l_vec = .basisVector(index = 1, size = numPostPeriods),
                                     monotonicityDirection = NULL,
                                     biasDirection = NULL,
                                     alpha = 0.05,
                                     parallel = FALSE) {

  # If Mvec is null, construct default Mvec
  if (is.null(Mvec)) {
    if (numPrePeriods == 1) {
      Mvec = seq(from = 0, to = sqrt(sigma[1, 1]), length.out = 10)
    } else {
      Mub = DeltaSD_upperBound_Mpre(betahat = betahat, sigma = sigma, numPrePeriods = numPrePeriods, alpha = 0.05)
      Mvec = seq(from = 0, to = Mub, length.out = 10)
    }
  }

  # If Delta = Delta^{SD}, construct confidence intervals
  if (is.null(monotonicityDirection) & is.null(biasDirection)) {
    if (is.null(method)) {
      method = "FLCI"
    }

    if (method == "FLCI") { # If method = FLCI, construct FLCIs
      if (parallel == FALSE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = alpha)
          tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                 method = "FLCI", Delta = "DeltaSD", M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = alpha)
          tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                 method = "FLCI", Delta = "DeltaSD", M = Mvec[m])
        }
      }
    } else if (method == "Conditional") { # If method = Conditional, construct conditional confidence intervals
      if (parallel == FALSE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                              hybrid_flag = "ARP")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "Conditional",
                 Delta = "DeltaSD", M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                              hybrid_flag = "ARP")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "Conditional",
                 Delta = "DeltaSD", M = Mvec[m])
        }
      }
    } else if (method == "C-F") { # If method = C-F, construct conditional FLCI,
      if (parallel == FALSE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                              hybrid_flag = "FLCI")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-F",
                 Delta = "DeltaSD", M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                              hybrid_flag = "FLCI")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-F",
                 Delta = "DeltaSD", M = Mvec[m])
        }
      }
    } else if (method == "C-LF") {
      if (parallel == FALSE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                              hybrid_flag = "LF")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-LF",
                 Delta = "DeltaSD", M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSD(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                              hybrid_flag = "LF")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-LF",
                 Delta = "DeltaSD", M = Mvec[m])
        }
      }
    } else{
      stop("Method must equal one of: FLCI, Conditional, C-F or C-LF")
    }
  } else if (!is.null(biasDirection)) {
    if (is.null(method)) {
      method = "C-F"
    }

    if (biasDirection == "positive") {
      Delta = "DeltaSDPB"
    } else {
      Delta = "DeltaSDNB"
    }

    if (method == "FLCI") {
      warning("You specified a sign restriction but method = FLCI. The FLCI does not use the sign restriction!")
      if (parallel == FALSE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = alpha)
          tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                 method = "FLCI", Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = alpha)
          tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                 method = "FLCI", Delta = Delta, M = Mvec[m])
        }
      }
    } else if (method == "Conditional") { # If method = Conditional, construct conditional confidence intervals
      if (parallel == TRUE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "ARP")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "Conditional",
                 Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "ARP")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "Conditional",
                 Delta = Delta, M = Mvec[m])
        }
      }
    } else if (method == "C-F") { # If method = C-F, construct conditional FLCI,
      if (parallel == TRUE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "FLCI")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-F",
                 Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "FLCI")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-F",
                 Delta = Delta, M = Mvec[m])
        }
      }
    } else if (method == "C-LF") {
      if (parallel == TRUE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "LF")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-LF",
                 Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = "LF")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-LF",
                 Delta = Delta, M = Mvec[m])
        }
      }
    } else{
      stop("Method must equal one of: FLCI, Conditional, C-F or C-LF")
    }
  } else {
    if (is.null(method)) {
      method = "C-F"
    }

    if (monotonicityDirection == "increasing") {
      Delta = "DeltaSDI"
    } else {
      Delta = "DeltaSDD"
    }

    if (method == "FLCI") {
      warning("You specified a shape restriction but method = FLCI. The FLCI does not use the shape restriction!")
      if (parallel == FALSE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = alpha)
          tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                 method = "FLCI", Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = alpha)
          tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                 method = "FLCI", Delta = Delta, M = Mvec[m])
        }
      }
    } else if (method == "Conditional") { # If method = Conditional, construct conditional confidence intervals
      if (parallel == FALSE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "ARP")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "Conditional",
                 Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "ARP")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "Conditional",
                 Delta = Delta, M = Mvec[m])
        }
      }
    } else if (method == "C-F") { # If method = C-F, construct conditional FLCI,
      if (parallel == FALSE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "FLCI")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-F",
                 Delta = Delta,
                 M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "FLCI")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-F",
                 Delta = Delta,
                 M = Mvec[m])
        }
      }
    } else if (method == "C-LF") {
      if (parallel == FALSE) {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "LF")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-LF",
                 Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, M = Mvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = "LF")
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = "C-LF",
                 Delta = Delta, M = Mvec[m])
        }
      }
    } else{
      stop("Method must equal one of: FLCI, Conditional, C-F or C-LF")
    }
  }
  return(Results)
}

createSensitivityPlot <- function(robustResults, originalResults, rescaleFactor = 1, maxM = Inf, add_xAxis = TRUE) {
  # Set M for OLS to be the min M in robust results minus the gap between Ms in robust
  Mgap <- min( diff( sort( robustResults$M) ) )
  Mmin <- min( robustResults$M)

  originalResults$M <- Mmin - Mgap
  df <- bind_rows(originalResults, robustResults)

  # Rescale all the units by rescaleFactor
  df <- df %>% mutate_at( c("M", "ub", "lb"), ~ .x * rescaleFactor)

  # Filter out observations above maxM (after rescaling)
  df <- df %>% dplyr::filter(M <= maxM)

  p <- ggplot(data = df, aes(x=M)) +
    geom_errorbar(aes(ymin = lb, ymax = ub, color = factor(method)),
                  width = Mgap * rescaleFactor / 2) +
    scale_color_manual(values = c("red", '#01a2d9')) +
    theme(legend.title=element_blank(), legend.position="bottom") +
    labs(x = "M", y = "")

  if (add_xAxis) {
    p <- p + geom_hline(yintercept = 0)
  }
  return(p)
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
                                                        gridPoints = 10^3, grid.ub = NA, grid.lb = NA,
                                                        parallel = FALSE) {

  # If Mbarvec is null, construct default Mbarvec to be 10 values on [0,2].
  if (is.null(Mbarvec)) {
    Mbarvec = seq(from = 0, to = 2, length.out = 10)
  }

  # Check if bound is specified correctly
  if (bound != "deviation from parallel trends" & bound != "deviation from linear trend") {
    stop("bound must equal either 'deviation from parallel trends' or 'deviation from linear trend'.")
  }

  # Check if both monotonicity direction and bias direction are specified:
  if (!is.null(monotonicityDirection) & !is.null(biasDirection)) {
    stop("Please select either a shape restriction or sign restriction (not both).")
  }

  # Use method to specify the choice of confidence set
  if (method == "C-LF") {
    hybrid_flag = "LF"
    method_named = "C-LF"
  } else if (method == "Conditional") {
    hybrid_flag = "ARP"
    method_named = "Conditional"
  } else {
    stop("method must be either NULL, Conditional or C-LF.")
  }

  # If bound = "parallel trends violation", we select Delta^{RM} and its variants.
  if (bound == "deviation from parallel trends") {
    # use monotonicity direction and biasDirection to select the choice of Delta
    if (is.null(monotonicityDirection) & is.null(biasDirection)) { # if both null, Delta^{RM}
      Delta = "DeltaRM"

      if (parallel) {
        Results = foreach(m = 1:length(Mbarvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaRM(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                              hybrid_flag = hybrid_flag,
                                              gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb)
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = method_named,
                 Delta = Delta, Mbar = Mbarvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mbarvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaRM(betahat = betahat, sigma = sigma,
                                              numPrePeriods = numPrePeriods,
                                              numPostPeriods = numPostPeriods,
                                              l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                              hybrid_flag = hybrid_flag,
                                              gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb)
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = method_named,
                 Delta = Delta, Mbar = Mbarvec[m])
        }
      }
    } else if (!is.null(monotonicityDirection)) { # if monotonicityDirection specified, Delta^{RMM}.
      if (monotonicityDirection == "increasing") {
        Delta = "DeltaRMI"
      } else if (monotonicityDirection == "decreasing") {
        Delta = "DeltaRMD"
      } else {
        stop("monotonicityDirection must equal either increasing or decreasing.")
      }

      if (parallel) {
        Results = foreach(m = 1:length(Mbarvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaRMM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = hybrid_flag,
                                               gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb)
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = method_named,
                 Delta = Delta, Mbar = Mbarvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mbarvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaRMM(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                               monotonicityDirection = monotonicityDirection,
                                               hybrid_flag = hybrid_flag,
                                               gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb)
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
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
        stop("monotonicityDirection must equal either positive or negative.")
      }

      if (parallel) {
        Results = foreach(m = 1:length(Mbarvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaRMB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = hybrid_flag,
                                               gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb)
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = method_named,
                 Delta = Delta, Mbar = Mbarvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mbarvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaRMB(betahat = betahat, sigma = sigma,
                                               numPrePeriods = numPrePeriods,
                                               numPostPeriods = numPostPeriods,
                                               l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                               biasDirection = biasDirection,
                                               hybrid_flag = hybrid_flag,
                                               gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb)
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
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
      stop("Error: not enough pre-periods for 'deviation from linear trend' (Delta^{SDRM} as base choice). Plese see documentation.")
    }

    # use monotonicity direction and biasDirection to select the choice of Delta
    if (is.null(monotonicityDirection) & is.null(biasDirection)) { # if both null, Delta^{SDRM}
      Delta = "DeltaSDRM"

      if (parallel) {
        Results = foreach(m = 1:length(Mbarvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDRM(betahat = betahat, sigma = sigma,
                                                numPrePeriods = numPrePeriods,
                                                numPostPeriods = numPostPeriods,
                                                l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                                hybrid_flag = hybrid_flag,
                                                gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb)
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = method_named,
                 Delta = Delta, Mbar = Mbarvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mbarvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDRM(betahat = betahat, sigma = sigma,
                                                numPrePeriods = numPrePeriods,
                                                numPostPeriods = numPostPeriods,
                                                l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                                hybrid_flag = hybrid_flag,
                                                gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb)
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = method_named,
                 Delta = Delta, Mbar = Mbarvec[m])
        }
      }
    } else if (!is.null(monotonicityDirection)) { # if monotonicityDirection specified, Delta^{SDRMM}.
      if (monotonicityDirection == "increasing") {
        Delta = "DeltaSDRMI"
      } else if (monotonicityDirection == "decreasing") {
        Delta = "DeltaSDRMD"
      } else {
        stop("monotonicityDirection must equal either increasing or decreasing.")
      }

      if (parallel) {
        Results = foreach(m = 1:length(Mbarvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDRMM(betahat = betahat, sigma = sigma,
                                                 numPrePeriods = numPrePeriods,
                                                 numPostPeriods = numPostPeriods,
                                                 l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                                 monotonicityDirection = monotonicityDirection,
                                                 hybrid_flag = hybrid_flag,
                                                 gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb)
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = method_named,
                 Delta = Delta, Mbar = Mbarvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mbarvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDRMM(betahat = betahat, sigma = sigma,
                                                 numPrePeriods = numPrePeriods,
                                                 numPostPeriods = numPostPeriods,
                                                 l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                                 monotonicityDirection = monotonicityDirection,
                                                 hybrid_flag = hybrid_flag,
                                                 gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb)
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
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
        stop("monotonicityDirection must equal either positive or negative.")
      }

      if (parallel) {
        Results = foreach(m = 1:length(Mbarvec), .combine = 'rbind') %do% {
          temp = computeConditionalCS_DeltaSDRMB(betahat = betahat, sigma = sigma,
                                                 numPrePeriods = numPrePeriods,
                                                 numPostPeriods = numPostPeriods,
                                                 l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                                 biasDirection = biasDirection,
                                                 hybrid_flag = hybrid_flag,
                                                 gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb)
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = method_named,
                 Delta = Delta, Mbar = Mbarvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mbarvec), .combine = 'rbind') %dopar% {
          temp = computeConditionalCS_DeltaSDRMB(betahat = betahat, sigma = sigma,
                                                 numPrePeriods = numPrePeriods,
                                                 numPostPeriods = numPostPeriods,
                                                 l_vec = l_vec, alpha = alpha, Mbar = Mbarvec[m],
                                                 biasDirection = biasDirection,
                                                 hybrid_flag = hybrid_flag,
                                                 gridPoints = gridPoints, grid.ub = grid.ub, grid.lb = grid.lb)
          tibble(lb = min(temp$grid[temp$accept == 1]),
                 ub = max(temp$grid[temp$accept == 1]),
                 method = method_named,
                 Delta = Delta, Mbar = Mbarvec[m])
        }
      }
    }
  }
  # Return tibble of results
  return(Results)
}

createSensitivityPlot_relativeMagnitudes <- function(robustResults, originalResults, rescaleFactor = 1, maxMbar = Inf, add_xAxis = TRUE) {
  # Set Mbar for OLS to be the min Mbar in robust results minus the gap between Mbars in robust
  Mbargap <- min( diff( sort( robustResults$Mbar) ) )
  Mbarmin <- min( robustResults$Mbar)

  originalResults$Mbar <- Mbarmin - Mbargap
  df <- bind_rows(originalResults, robustResults)

  # Rescale all the units by rescaleFactor
  df <- df %>% mutate_at( c("Mbar", "ub", "lb"), ~ .x * rescaleFactor)

  # Filter out observations above maxM (after rescaling)
  df <- df %>% dplyr::filter(Mbar <= maxMbar)

  p <- ggplot(data = df, aes(x=Mbar)) +
    geom_errorbar(aes(ymin = lb, ymax = ub, color = factor(method)),
                  width = Mbargap * rescaleFactor / 2) +
    scale_color_manual(values = c("red", '#01a2d9')) +
    theme(legend.title=element_blank(), legend.position="bottom") +
    labs(x = latex2exp::TeX("$\\bar{M}$"), y = "")

  if (add_xAxis) {
    p <- p + geom_hline(yintercept = 0)
  }
  return(p)
}

# Construct Original CS -----------------------------------------------
constructOriginalCS <- function(betahat, sigma,
                                numPrePeriods, numPostPeriods,
                                l_vec = .basisVector(index = 1, size = numPostPeriods),
                                alpha = 0.05) {
  stdError = sqrt(t(l_vec) %*% sigma[(numPrePeriods+1):(numPrePeriods+numPostPeriods), (numPrePeriods+1):(numPrePeriods+numPostPeriods)] %*% l_vec)
  lb = t(l_vec) %*% betahat[(numPrePeriods+1):(numPrePeriods+numPostPeriods)] - qnorm(1-alpha/2)*stdError
  ub = t(l_vec) %*% betahat[(numPrePeriods+1):(numPrePeriods+numPostPeriods)] + qnorm(1-alpha/2)*stdError
  return(tibble(
    lb = lb,
    ub = ub,
    method = "Original",
    Delta = NA
  ))
}

createEventStudyPlot <- function(betahat, stdErrors = NULL, sigma = NULL,
                                 numPrePeriods, numPostPeriods, alpha = 0.05,
                                 timeVec, referencePeriod,
                                 useRelativeEventTime = F) {
  if (is.null(stdErrors) & is.null(sigma)) {
    stop("User must specify either vector of standard errors or vcv matrix!")
  } else if (is.null(stdErrors) & !is.null(sigma)) {
    stdErrors = sqrt(diag(sigma))
  }

  if (useRelativeEventTime == T) {
    timeVec = timeVec - referencePeriod
    referencePeriod = 0
  }

  EventStudyPlot <- ggplot(tibble(t = c(timeVec[1:numPrePeriods], referencePeriod, timeVec[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]),
                                  beta = c(betahat[1:numPrePeriods], 0, betahat[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]),
                                  se = c(stdErrors[1:numPrePeriods], NA, stdErrors[(numPrePeriods+1):(numPrePeriods+numPostPeriods)])),
                           aes(x = t)) +
    geom_point(aes(y = beta), color = "red") +
    geom_errorbar(aes(ymin = beta - qnorm(1-alpha/2)*se, ymax = beta + qnorm(1-alpha/2)*se), width = 0.5, colour = "#01a2d9") +
    theme(legend.position = "none") + labs(x = "Event time", y = "") +
    scale_x_continuous(breaks = seq(from = min(timeVec), to = max(timeVec), by = 1),
                       labels = as.character(seq(from = min(timeVec), to = max(timeVec), by = 1)))
  return(EventStudyPlot)
}
