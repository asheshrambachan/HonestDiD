# DESCRIPTION =========================================================
# Author: Ashesh Rambachan
#
#  This script contains functions to implement the methods
#  described in Rambachan & Roth (2019) for robust inference
#  in difference-in-differences and event study designs.
#
#   First, this script contains function that are used to construct the
#   FLCI for a general choice of the vector l and Delta = Delta^{SD}(M).
#
#   Second, this script defines functions used to construct conditional
#   confidence intervals and hybrid confidence intervals for a general choice of the
#   vector l and Delta = Delta^{SD}, Delta^{SDPB}, Delta^{SDNB},
#   Delta^{SDI}, Delta^{SDD}. This implements functions to deal with
#   nuisance parameters as described in  Andrews, Roth, Pakes (2019).

# PRELIMINARIES =======================================================
library(TruncatedNormal)
library(lpSolveAPI)
library(Matrix)
library(pracma)
library(CVXR)
library(foreach)

# Construct Robust Results Function -----------------------------------
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
      Mvec = seq(from = 0, to = c(sigma[1, 1]), length.out = 10)
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
                                 l_vec = l_vec, M = Mvec[m], alpha = 0.05)
          tibble::tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                         method = "FLCI", Delta = "DeltaSD", M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = 0.05)
          tibble::tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
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
          tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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
          tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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
          tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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
          tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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
          tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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
          tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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
                                 l_vec = l_vec, M = Mvec[m], alpha = 0.05)
         tibble::tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                        method = "FLCI", Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = 0.05)
         tibble::tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
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
          tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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
          tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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
          tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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
          tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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
          tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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
          tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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
                                 l_vec = l_vec, M = Mvec[m], alpha = 0.05)
         tibble::tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
                        method = "FLCI", Delta = Delta, M = Mvec[m])
        }
      } else {
        Results = foreach(m = 1:length(Mvec), .combine = 'rbind') %dopar% {
          temp = findOptimalFLCI(betahat = betahat, sigma = sigma,
                                 numPrePeriods = numPrePeriods, numPostPeriods = numPostPeriods,
                                 l_vec = l_vec, M = Mvec[m], alpha = 0.05)
         tibble::tibble(lb = temp$FLCI[1], ub = temp$FLCI[2],
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
          tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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
          tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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
         tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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
          tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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
          tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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
          tibble::tibble(lb = min(temp$grid[temp$accept == 1]),
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

constructOriginalCS <- function(betahat, sigma,
                                numPrePeriods, numPostPeriods,
                                l_vec = .basisVector(index = 1, size = numPostPeriods),
                                alpha = 0.05) {
  stdError = sqrt(t(l_vec) %*% sigma[(numPrePeriods+1):(numPrePeriods+numPostPeriods), (numPrePeriods+1):(numPrePeriods+numPostPeriods)] %*% l_vec)
  lb = t(l_vec) %*% betahat[(numPrePeriods+1):(numPrePeriods+numPostPeriods)] - qnorm(1-alpha/2)*stdError
  ub = t(l_vec) %*% betahat[(numPrePeriods+1):(numPrePeriods+numPostPeriods)] + qnorm(1-alpha/2)*stdError
  return(tibble::tibble(
    lb = lb,
    ub = ub,
    method = "Original",
    Delta = NA,
    M = 0
  ))
}

# Sensitivity plot functions ------------------------------------------
createEventStudyPlot <- function(betahat, stdErrors = NULL, sigma = NULL,
                                 numPrePeriods, numPostPeriods, alpha = 0.05,
                                 timeVec, referencePeriod,
                                 useRelativeEventTime = F) {
  if (is.null(stdErrors) & is.null(sigma)) {
    stop("User must specify either vector of standard errors or vcv matrix!")
  } else if (is.null(stdErrors) & !is.null(sigma)) {
    stdErrors = sqrt(diag(sigma))
  }

  if (useRelativeEventTime == TRUE) {
    timeVec = timeVec - referencePeriod
    referencePeriod = 0
  }

  EventStudyPlot <- ggplot(tibble::tibble(t = c(timeVec[1:numPrePeriods], referencePeriod, timeVec[(numPrePeriods+1):(numPrePeriods+numPostPeriods)]),
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

createSensitivityPlot <- function(robustResults, originalResults, rescaleFactor = 1, maxM = Inf, add_xAxis = TRUE) {
  # Set M for OLS to be the min M in robust results minus the gap between Ms in robust
  Mgap <- min( diff( sort( robustResults$M) ) )
  Mmin <- min( robustResults$M)
  originalResults$M <- Mmin - Mgap
  df <- bind_rows(originalResults, robustResults)

  # Rescale all the units by rescaleFactor
  df <- df %>% mutate_at( c("M", "ub", "lb"), ~ .x * rescaleFactor)

  # Filter out observations above maxM (after rescaling)
  df <- df %>% filter(M <= maxM)

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

# Other helper functions ----
basisVector <- function(index = 1, size = 1){
  v <- matrix(0, nrow = size, ncol = 1)
  v[index] = 1
  return(v)
}
