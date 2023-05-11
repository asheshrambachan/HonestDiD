#' @title honest_did
#'
#' @description a function to compute a sensitivity analysis
#'  using the approach of Rambachan and Roth (2021)
honest_did <- function(...) UseMethod("honest_did")

#' @title honest_did.AGGTEobj
#'
#' @description a function to compute a sensitivity analysis
#'  using the approach of Rambachan and Roth (2021) when
#'  the event study is estimating using the `did` package
#'
#' @param e event time to compute the sensitivity analysis for.
#'  The default value is `e=0` corresponding to the "on impact"
#'  effect of participating in the treatment.
#' @param type Options are "smoothness" (which conducts a
#'  sensitivity analysis allowing for violations of linear trends
#'  in pre-treatment periods) or "relative_magnitude" (which
#'  conducts a sensitivity analysis based on the relative magnitudes
#'  of deviations from parallel trends in pre-treatment periods).
#' @inheritParams HonestDiD::createSensitivityResults
#' @inheritParams HonestDid::createSensitivityResults_relativeMagnitudes
honest_did.AGGTEobj <- function(es,
                                e          = 0,
                                type       = c("smoothness", "relative_magnitude"),
                                gridPoints = 100,
                                ...) {

  type <- match.arg(type)

  # Make sure that user is passing in an event study
  if (es$type != "dynamic") {
    stop("need to pass in an event study")
  }

  # Check if used universal base period and warn otherwise
  if (es$DIDparams$base_period != "universal") {
    stop("Use a universal base period for honest_did")
  }

  # Recover influence function for event study estimates
  es_inf_func <- es$inf.function$dynamic.inf.func.e

  # Recover variance-covariance matrix
  n <- nrow(es_inf_func)
  V <- t(es_inf_func) %*% es_inf_func / n / n

  # Remove the coefficient normalized to zero
  referencePeriodIndex <- which(es$egt == -1)
  V    <- V[-referencePeriodIndex,-referencePeriodIndex]
  beta <- es$att.egt[-referencePeriodIndex]

  nperiods <- nrow(V)
  npre     <- sum(1*(es$egt < -1))
  npost    <- nperiods - npre
  baseVec1 <- basisVector(index=(e+1),size=npost)
  orig_ci  <- constructOriginalCS(betahat        = beta,
                                  sigma          = V,
                                  numPrePeriods  = npre,
                                  numPostPeriods = npost,
                                  l_vec          = baseVec1)

  if (type=="relative_magnitude") {
    robust_ci <- createSensitivityResults_relativeMagnitudes(betahat        = beta,
                                                             sigma          = V,
                                                             numPrePeriods  = npre,
                                                             numPostPeriods = npost,
                                                             l_vec          = baseVec1,
                                                             gridPoints     = gridPoints,
                                                             ...)

  } else if (type == "smoothness") {
    robust_ci <- createSensitivityResults(betahat        = beta,
                                          sigma          = V,
                                          numPrePeriods  = npre,
                                          numPostPeriods = npost,
                                          l_vec          = baseVec1,
                                          ...)
  }

  return(list(robust_ci=robust_ci, orig_ci=orig_ci, type=type))
}
