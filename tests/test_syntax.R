# remove.packages("HonestDiD")
# install.packages(".", repos=NULL, type="source")
# testthat::test_dir("tests")
#
# Syntax smoke tests: exercise as many exported entry points and option
# combinations as possible with tiny grids so the full suite runs fast.
# The goal is breadth (catching regressions in call-site plumbing, e.g. the
# CVXR API change behind issue #68) rather than numerical correctness.

library(testthat)
library(HonestDiD)
data(BCdata_EventStudy)

BC_betahat        <- BCdata_EventStudy$betahat
BC_sigma          <- BCdata_EventStudy$sigma
BC_numPrePeriods  <- length(BCdata_EventStudy$prePeriodIndices)
BC_numPostPeriods <- length(BCdata_EventStudy$postPeriodIndices)
BC_l_vec          <- basisVector(index = 1, size = BC_numPostPeriods)

# Tiny knobs to keep the suite fast. The tests verify that each entry point
# runs and returns the expected object type, not that results are accurate,
# so single-value M/Mbar vectors and a 10-point grid are enough.
SYN_Mvec       <- 0.1
SYN_Mbarvec    <- 0.5
SYN_gridPoints <- 10
SYN_grid.ub    <- 1
SYN_grid.lb    <- -1

# Suppress benign noise: several low-level helpers print progress to stdout
# and the tiny grid occasionally produces empty accept-sets, which trip
# min()/max() warnings. Wrap runners so the tests stay focused on errors.
.quiet <- function(expr) base::suppressWarnings(base::suppressMessages(
  utils::capture.output(val <- expr)))
.run <- function(expr) {
  .quiet(val <- expr)
  val
}
expect_runs <- function(expr) {
  val <- .run(expr)
  expect_false(base::is.null(val))
  base::invisible(val)
}
expect_tibble <- function(expr, expected_rows = NA) {
  val <- .run(expr)
  expect_true(base::inherits(val, "data.frame"))
  if ( !base::is.na(expected_rows) ) expect_equal(base::nrow(val), expected_rows)
  base::invisible(val)
}

if ( Sys.getenv("HONESTDID_RUN_TESTS") == "1" ) {
  test_that("HonestDiD utilities and low-level helpers run with no errors", {
    # basisVector: both forms used in examples
    expect_runs(basisVector(index = 1, size = BC_numPostPeriods))
    expect_runs(basisVector(index = BC_numPostPeriods, size = BC_numPostPeriods))

    # Pre-period-based M bounds
    expect_runs(DeltaSD_upperBound_Mpre(betahat       = BC_betahat,
                                        sigma         = BC_sigma,
                                        numPrePeriods = BC_numPrePeriods,
                                        alpha         = 0.05))
    expect_runs(DeltaSD_lowerBound_Mpre(betahat       = BC_betahat,
                                        sigma         = BC_sigma,
                                        numPrePeriods = BC_numPrePeriods,
                                        alpha         = 0.05,
                                        grid.ub       = SYN_grid.ub,
                                        gridPoints    = SYN_gridPoints))

    # Direct FLCI entry point
    flci <- expect_runs(findOptimalFLCI(betahat        = BC_betahat,
                                        sigma          = BC_sigma,
                                        numPrePeriods  = BC_numPrePeriods,
                                        numPostPeriods = BC_numPostPeriods,
                                        l_vec          = BC_l_vec,
                                        M              = SYN_Mvec,
                                        alpha          = 0.05,
                                        seed           = 0))
    expect_true(base::all(c("FLCI", "optimalVec", "optimalHalfLength", "status") %in% base::names(flci)))

    # Original (non-robust) CS
    expect_tibble(constructOriginalCS(betahat        = BC_betahat,
                                      sigma          = BC_sigma,
                                      numPrePeriods  = BC_numPrePeriods,
                                      numPostPeriods = BC_numPostPeriods,
                                      l_vec          = BC_l_vec))
  })
} else {
  print("HonestDiD utilities run was skipped")
}

if ( Sys.getenv("HONESTDID_RUN_TESTS") == "1" ) {
  test_that("createSensitivityResults covers all method x restriction combos", {
    methods <- c("FLCI", "Conditional", "C-F", "C-LF")

    # No shape / sign restriction
    for ( m in methods ) {
      expect_tibble(
        createSensitivityResults(betahat        = BC_betahat,
                                 sigma          = BC_sigma,
                                 numPrePeriods  = BC_numPrePeriods,
                                 numPostPeriods = BC_numPostPeriods,
                                 l_vec          = BC_l_vec,
                                 method         = m,
                                 Mvec           = SYN_Mvec,
                                 alpha          = 0.05),
        expected_rows = base::length(SYN_Mvec))
    }

    # Bias direction restriction. FLCI ignores sign restrictions (it
    # delegates to findOptimalFLCI and only warns about the extra argument),
    # so its call path is already covered in (a); skip here to save time.
    for ( m in base::setdiff(methods, "FLCI") ) {
      for ( bd in c("positive", "negative") ) {
        expect_tibble(
          createSensitivityResults(betahat        = BC_betahat,
                                   sigma          = BC_sigma,
                                   numPrePeriods  = BC_numPrePeriods,
                                   numPostPeriods = BC_numPostPeriods,
                                   l_vec          = BC_l_vec,
                                   method         = m,
                                   biasDirection  = bd,
                                   Mvec           = SYN_Mvec,
                                   alpha          = 0.05),
          expected_rows = base::length(SYN_Mvec))
      }
    }

    # Monotonicity (shape) restriction. Same caveat as (b): FLCI ignores
    # shape restrictions, so we skip it here.
    for ( m in base::setdiff(methods, "FLCI") ) {
      for ( md in c("increasing", "decreasing") ) {
        expect_tibble(
          createSensitivityResults(betahat               = BC_betahat,
                                   sigma                 = BC_sigma,
                                   numPrePeriods         = BC_numPrePeriods,
                                   numPostPeriods        = BC_numPostPeriods,
                                   l_vec                 = BC_l_vec,
                                   method                = m,
                                   monotonicityDirection = md,
                                   Mvec                  = SYN_Mvec,
                                   alpha                 = 0.05),
          expected_rows = base::length(SYN_Mvec))
      }
    }
  })
} else {
  print("HonestDiD createSensitivityResults run was skipped")
}

if ( Sys.getenv("HONESTDID_RUN_TESTS") == "1" ) {
  test_that("createSensitivityResults_relativeMagnitudes covers all bound x method x restriction combos", {
    methods <- c("C-LF", "Conditional")
    bounds  <- c("deviation from parallel trends", "deviation from linear trend")

    # No shape / sign restriction
    for ( b in bounds ) {
      for ( m in methods ) {
        expect_tibble(
          createSensitivityResults_relativeMagnitudes(
            betahat        = BC_betahat,
            sigma          = BC_sigma,
            numPrePeriods  = BC_numPrePeriods,
            numPostPeriods = BC_numPostPeriods,
            l_vec          = BC_l_vec,
            bound          = b,
            method         = m,
            Mbarvec        = SYN_Mbarvec,
            gridPoints     = SYN_gridPoints,
            grid.ub        = SYN_grid.ub,
            grid.lb        = SYN_grid.lb,
            alpha          = 0.05),
          expected_rows = base::length(SYN_Mbarvec))
      }
    }

    # Bias direction restriction
    for ( b in bounds ) {
      for ( m in methods ) {
        for ( bd in c("positive", "negative") ) {
          expect_tibble(
            createSensitivityResults_relativeMagnitudes(
              betahat        = BC_betahat,
              sigma          = BC_sigma,
              numPrePeriods  = BC_numPrePeriods,
              numPostPeriods = BC_numPostPeriods,
              l_vec          = BC_l_vec,
              bound          = b,
              method         = m,
              biasDirection  = bd,
              Mbarvec        = SYN_Mbarvec,
              gridPoints     = SYN_gridPoints,
              grid.ub        = SYN_grid.ub,
              grid.lb        = SYN_grid.lb,
              alpha          = 0.05),
            expected_rows = base::length(SYN_Mbarvec))
        }
      }
    }

    # Monotonicity (shape) restriction
    for ( b in bounds ) {
      for ( m in methods ) {
        for ( md in c("increasing", "decreasing") ) {
          expect_tibble(
            createSensitivityResults_relativeMagnitudes(
              betahat               = BC_betahat,
              sigma                 = BC_sigma,
              numPrePeriods         = BC_numPrePeriods,
              numPostPeriods        = BC_numPostPeriods,
              l_vec                 = BC_l_vec,
              bound                 = b,
              method                = m,
              monotonicityDirection = md,
              Mbarvec               = SYN_Mbarvec,
              gridPoints            = SYN_gridPoints,
              grid.ub               = SYN_grid.ub,
              grid.lb               = SYN_grid.lb,
              alpha                 = 0.05),
            expected_rows = base::length(SYN_Mbarvec))
        }
      }
    }
  })
} else {
  print("HonestDiD createSensitivityResults_relativeMagnitudes run was skipped")
}

if ( Sys.getenv("HONESTDID_RUN_TESTS") == "1" ) {
  test_that("Low-level computeConditionalCS_* entry points all run", {
    # Baseline calls hit every hybrid_flag to catch solver-plumbing regressions
    # (this is the path the CVXR API change in issue #68 broke). Restricted
    # variants only need to exercise one hybrid to cover the call-site syntax.
    for ( hf in c("ARP", "FLCI", "LF") ) {
      expect_runs(computeConditionalCS_DeltaSD(
        betahat        = BC_betahat,
        sigma          = BC_sigma,
        numPrePeriods  = BC_numPrePeriods,
        numPostPeriods = BC_numPostPeriods,
        l_vec          = BC_l_vec,
        M              = SYN_Mvec,
        alpha          = 0.05,
        hybrid_flag    = hf,
        gridPoints     = SYN_gridPoints,
        grid.ub        = SYN_grid.ub,
        grid.lb        = SYN_grid.lb,
        seed           = 0))
    }
    for ( hf in c("ARP", "LF") ) {
      expect_runs(computeConditionalCS_DeltaRM(
        betahat        = BC_betahat,
        sigma          = BC_sigma,
        numPrePeriods  = BC_numPrePeriods,
        numPostPeriods = BC_numPostPeriods,
        l_vec          = BC_l_vec,
        Mbar           = SYN_Mbarvec, alpha = 0.05,
        hybrid_flag    = hf,
        gridPoints     = SYN_gridPoints,
        grid.ub        = SYN_grid.ub,
        grid.lb        = SYN_grid.lb, seed  = 0))
      expect_runs(computeConditionalCS_DeltaSDRM(
        betahat        = BC_betahat,
        sigma          = BC_sigma,
        numPrePeriods  = BC_numPrePeriods,
        numPostPeriods = BC_numPostPeriods,
        l_vec          = BC_l_vec,
        Mbar           = SYN_Mbarvec, alpha = 0.05,
        hybrid_flag    = hf,
        gridPoints     = SYN_gridPoints,
        grid.ub        = SYN_grid.ub,
        grid.lb        = SYN_grid.lb, seed  = 0))
    }

    # Bias-direction restricted variants (one hybrid flag each)
    for ( bd in c("positive", "negative") ) {
      expect_runs(computeConditionalCS_DeltaSDB(
        betahat        = BC_betahat,
        sigma          = BC_sigma,
        numPrePeriods  = BC_numPrePeriods,
        numPostPeriods = BC_numPostPeriods,
        l_vec          = BC_l_vec,
        M              = SYN_Mvec,
        alpha          = 0.05,
        biasDirection  = bd,
        hybrid_flag    = "ARP",
        gridPoints     = SYN_gridPoints,
        grid.ub        = SYN_grid.ub,
        grid.lb        = SYN_grid.lb,
        seed           = 0))
      expect_runs(computeConditionalCS_DeltaRMB(
        betahat        = BC_betahat,
        sigma          = BC_sigma,
        numPrePeriods  = BC_numPrePeriods,
        numPostPeriods = BC_numPostPeriods,
        l_vec          = BC_l_vec,
        Mbar           = SYN_Mbarvec,
        alpha          = 0.05,
        biasDirection  = bd,
        hybrid_flag    = "LF",
        gridPoints     = SYN_gridPoints,
        grid.ub        = SYN_grid.ub,
        grid.lb        = SYN_grid.lb,
        seed           = 0))
      expect_runs(computeConditionalCS_DeltaSDRMB(
        betahat        = BC_betahat,
        sigma          = BC_sigma,
        numPrePeriods  = BC_numPrePeriods,
        numPostPeriods = BC_numPostPeriods,
        l_vec          = BC_l_vec,
        Mbar           = SYN_Mbarvec,
        alpha          = 0.05,
        biasDirection  = bd,
        hybrid_flag    = "LF",
        gridPoints     = SYN_gridPoints,
        grid.ub        = SYN_grid.ub,
        grid.lb        = SYN_grid.lb,
        seed           = 0))
    }

    # Monotonicity-restricted variants (one hybrid flag each)
    for ( md in c("increasing", "decreasing") ) {
      expect_runs(computeConditionalCS_DeltaSDM(
        betahat               = BC_betahat,
        sigma                 = BC_sigma,
        numPrePeriods         = BC_numPrePeriods,
        numPostPeriods        = BC_numPostPeriods,
        l_vec                 = BC_l_vec,
        M                     = SYN_Mvec,
        alpha                 = 0.05,
        monotonicityDirection = md,
        hybrid_flag           = "ARP",
        gridPoints            = SYN_gridPoints,
        grid.ub               = SYN_grid.ub,
        grid.lb               = SYN_grid.lb,
        seed                  = 0))
      expect_runs(computeConditionalCS_DeltaRMM(
        betahat               = BC_betahat,
        sigma                 = BC_sigma,
        numPrePeriods         = BC_numPrePeriods,
        numPostPeriods        = BC_numPostPeriods,
        l_vec                 = BC_l_vec,
        Mbar                  = SYN_Mbarvec,
        alpha                 = 0.05,
        monotonicityDirection = md,
        hybrid_flag           = "LF",
        gridPoints            = SYN_gridPoints,
        grid.ub               = SYN_grid.ub,
        grid.lb               = SYN_grid.lb,
        seed                  = 0))
      expect_runs(computeConditionalCS_DeltaSDRMM(
        betahat               = BC_betahat,
        sigma                 = BC_sigma,
        numPrePeriods         = BC_numPrePeriods,
        numPostPeriods        = BC_numPostPeriods,
        l_vec                 = BC_l_vec,
        Mbar                  = SYN_Mbarvec,
        alpha                 = 0.05,
        monotonicityDirection = md,
        hybrid_flag           = "LF",
        gridPoints            = SYN_gridPoints,
        grid.ub               = SYN_grid.ub,
        grid.lb               = SYN_grid.lb,
        seed                  = 0))
    }
  })
} else {
  print("HonestDiD computeConditionalCS run was skipped")
}

if ( Sys.getenv("HONESTDID_RUN_TESTS") == "1" ) {
  test_that("Plotting wrappers run with no errors", {
    origCS <- .run(constructOriginalCS(betahat        = BC_betahat,
                                       sigma          = BC_sigma,
                                       numPrePeriods  = BC_numPrePeriods,
                                       numPostPeriods = BC_numPostPeriods,
                                       l_vec          = BC_l_vec))

    # Need >= 2 M values so createSensitivityPlot can compute the between-M gap
    robustSD <- .run(createSensitivityResults(betahat        = BC_betahat,
                                              sigma          = BC_sigma,
                                              numPrePeriods  = BC_numPrePeriods,
                                              numPostPeriods = BC_numPostPeriods,
                                              l_vec          = BC_l_vec,
                                              method         = "FLCI",
                                              Mvec           = c(0, 0.1)))

    robustRM <- .run(createSensitivityResults_relativeMagnitudes(
      betahat        = BC_betahat,
      sigma          = BC_sigma,
      numPrePeriods  = BC_numPrePeriods,
      numPostPeriods = BC_numPostPeriods,
      l_vec          = BC_l_vec,
      bound          = "deviation from parallel trends",
      method         = "C-LF",
      Mbarvec        = c(0, 0.5),
      gridPoints     = SYN_gridPoints,
      grid.ub        = SYN_grid.ub,
      grid.lb        = SYN_grid.lb))

    expect_runs(createSensitivityPlot(robustResults   = robustSD,
                                      originalResults = origCS))
    expect_runs(createSensitivityPlot_relativeMagnitudes(robustResults   = robustRM,
                                                         originalResults = origCS))

    # Event-study plot: one with stdErrors, one with full sigma
    timeVec <- base::c(base::seq(from = -BC_numPrePeriods, to = -1),
                       base::seq(from = 1, to = BC_numPostPeriods))
    stdErrors <- base::sqrt(base::diag(BC_sigma))

    expect_runs(createEventStudyPlot(betahat         = BC_betahat,
                                     stdErrors       = stdErrors,
                                     numPrePeriods   = BC_numPrePeriods,
                                     numPostPeriods  = BC_numPostPeriods,
                                     timeVec         = timeVec,
                                     referencePeriod = 0))

    expect_runs(createEventStudyPlot(betahat         = BC_betahat,
                                     sigma           = BC_sigma,
                                     numPrePeriods   = BC_numPrePeriods,
                                     numPostPeriods  = BC_numPostPeriods,
                                     timeVec         = timeVec,
                                     referencePeriod = 0))
  })
} else {
  print("HonestDiD plotting run was skipped")
}
