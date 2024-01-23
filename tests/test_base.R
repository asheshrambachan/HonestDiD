# remove.packages("HonestDiD")
# install.packages(".", repos=NULL, type="source")
# testthat::test_dir("tests")

library(lfe)
library(haven)
library(testthat)
library(HonestDiD)
data(BCdata_EventStudy)
data(LWdata_EventStudy)

BC_numPrePeriods  <- length(BCdata_EventStudy$prePeriodIndices)
BC_numPostPeriods <- length(BCdata_EventStudy$postPeriodIndices)
BC_l_vec          <- basisVector(index = 1, size = BC_numPostPeriods)
BC_l_vec          <- cbind(c(1, 0, 0, 0))

if ( Sys.getenv("HONESTDID_RUN_TESTS") == "1" ) {
  test_that("HonestDiD base ran with no errors", {
    BC_DeltaSDNB_RobustResults <-
        createSensitivityResults(betahat        = BCdata_EventStudy$betahat,
                                 sigma          = BCdata_EventStudy$sigma,
                                 numPrePeriods  = BC_numPrePeriods,
                                 numPostPeriods = BC_numPostPeriods,
                                 l_vec          = BC_l_vec,
                                 method         = "FLCI",
                                 Mvec           = seq(from=0, to=0.3, by=0.1))

    BC_DeltaSDNB_RobustResultsConditional <-
        createSensitivityResults(betahat        = BCdata_EventStudy$betahat,
                                 sigma          = BCdata_EventStudy$sigma,
                                 numPrePeriods  = BC_numPrePeriods,
                                 numPostPeriods = BC_numPostPeriods,
                                 l_vec          = BC_l_vec,
                                 method         = "Conditional",
                                 Mvec           = seq(from=0, to=0.3, by=0.1))

    BC_DeltaSDNB_RobustResultsCF <-
        createSensitivityResults(betahat        = BCdata_EventStudy$betahat,
                                 sigma          = BCdata_EventStudy$sigma,
                                 numPrePeriods  = BC_numPrePeriods,
                                 numPostPeriods = BC_numPostPeriods,
                                 l_vec          = BC_l_vec,
                                 method         = "C-F",
                                 Mvec           = seq(from=0, to=0.3, by=0.1))

    BC_DeltaSDNB_RobustResultsCLF <-
        createSensitivityResults(betahat        = BCdata_EventStudy$betahat,
                                 sigma          = BCdata_EventStudy$sigma,
                                 numPrePeriods  = BC_numPrePeriods,
                                 numPostPeriods = BC_numPostPeriods,
                                 l_vec          = BC_l_vec,
                                 method         = "C-LF",
                                 Mvec           = seq(from=0, to=0.3, by=0.1))

    BC_OriginalResults <-
        constructOriginalCS(betahat        = BCdata_EventStudy$betahat,
                            sigma          = BCdata_EventStudy$sigma,
                            numPrePeriods  = BC_numPrePeriods,
                            numPostPeriods = BC_numPostPeriods,
                            l_vec          = BC_l_vec)

    BC_DeltaSDNB_SensitivityPlot <-
        createSensitivityPlot(robustResults   = BC_DeltaSDNB_RobustResults,
                              originalResults = BC_OriginalResults)

    expect_silent(BC_DeltaSDNB_RobustResults)
    expect_silent(BC_DeltaSDNB_RobustResultsConditional)
    expect_silent(BC_DeltaSDNB_RobustResultsCF)
    expect_silent(BC_DeltaSDNB_RobustResultsCLF)
    expect_silent(BC_OriginalResults)
    expect_silent(BC_DeltaSDNB_SensitivityPlot)
  })
} else {
  print("HonestDiD base run was skipped")
}

if ( Sys.getenv("HONESTDID_RUN_TESTS") == "1" ) {
  test_that("HonestDiD options ran with no errors", {
    LWdata_RawData = haven::read_dta(system.file("extdata", "LWdata_RawData.dta", package = "HonestDiD"))
    sum(LWdata_RawData$nobs)

    # Estimate event study using lfe package
    EmpFemale.EventStudy = lfe::felm(emp ~
                                     rtESV13 + rtESV14 + rtESV15 +
                                     rtESV16 + rtESV17 + rtESV18 +
                                     rtESV19 + rtESV110 + rtESV111 + # End Pre-periods
                                     rtESV113 + rtESV114 + rtESV115 +
                                     rtESV116 + rtESV117 + rtESV118 +
                                     rtESV119 + rtESV120 + rtESV121 +
                                     rtESV122 + rtESV123 + rtESV124 +
                                     rtESV125 + rtESV126 + rtESV127 +
                                     rtESV128 + rtESV129 + rtESV130 +
                                     rtESV131 + rtESV132 + rtESV133 +
                                     rtESV134 + rtESV135 + # End post-periods
                                     yearsfcor + yearsflr + aveitc + fscontrol +
                                     asian + black + hispanic + other |
                                     factor(PUS_SURVEY_YEAR)*factor(BIRTHYEAR) +
                                     factor(PUS_SURVEY_YEAR) + factor(BIRTHSTATE) |
                                     0 | BIRTHSTATE,
                                     data = LWdata_RawData,
                                     weights = LWdata_RawData$nobs)
    summary(EmpFemale.EventStudy)

    coefIndex = which(grepl(x = dimnames(EmpFemale.EventStudy$coefficients)[[1]], pattern = "rtESV"))
    betahat = EmpFemale.EventStudy$beta[coefIndex, ]

    # Extract estimated variance-covariance matrix of event study coefficients
    sigma = EmpFemale.EventStudy$clustervcv[coefIndex, coefIndex]

    # Construct vector of event times and the scalar reference period
    timeVec = c(seq(from = -11, to = -3, by = 1), seq(from = -1, to = 21, by = 1))
    referencePeriod   <- -2
    postPeriodIndices <- which(timeVec > -2)
    prePeriodIndices  <- which(timeVec < -2)
    LW_numPrePeriods  <- length(prePeriodIndices)
    LW_numPostPeriods <- length(postPeriodIndices)
    LW_l_vec          <- basisVector(index = 1, size = LW_numPostPeriods)

    for( method in c("C-F", "C-LF", "Conditional", "FLCI") ) {
      for( monotonicityDirection in c("increasing", "decreasing") ) {
        for ( biasDirection in c("positive", "negative") ) {
          LW_DeltaSDNB_RobustResults <-
              createSensitivityResults(betahat               = betahat,
                                       sigma                 = sigma,
                                       numPrePeriods         = LW_numPrePeriods,
                                       numPostPeriods        = LW_numPostPeriods,
                                       l_vec                 = LW_l_vec,
                                       method                = method,
                                       monotonicityDirection = monotonicityDirection,
                                       biasDirection         = biasDirection,
                                       Mvec                  = seq(from=0, to=0.3, by=0.1))
          print(c(method, monotonicityDirection, biasDirection, LW_DeltaSDNB_RobustResults))
          expect_silent(LW_DeltaSDNB_RobustResults)
        }
      }
    }

    for ( method in c(NULL, "C-LF", "Conditional") ) {
      for ( monotonicityDirection in c("increasing", "decreasing", NULL) ) {
        for ( bound in c("deviation from parallel trends", "deviation from linear trend") ) {
          BC_DeltaRM_RobustResults <-
            createSensitivityResults_relativeMagnitudes(betahat               = BCdata_EventStudy$betahat,
                                                        sigma                 = BCdata_EventStudy$sigma,
                                                        numPrePeriods         = BC_numPrePeriods,
                                                        numPostPeriods        = BC_numPostPeriods,
                                                        l_vec                 = BC_l_vec,
                                                        gridPoints            = 100,
                                                        grid.ub               = 1,
                                                        grid.lb               = -1,
                                                        bound                 = bound,
                                                        method                = method,
                                                        monotonicityDirection = monotonicityDirection,
                                                        Mbarvec               = seq(from=0, to=1, by=0.5))
          print(c(method, monotonicityDirection, biasDirection, bound, BC_DeltaRM_RobustResults))
          expect_silent(BC_DeltaRM_RobustResults)
        }
      }
    }

    for ( method in c(NULL, "C-LF", "Conditional") ) {
      for ( biasDirection in c("positive", "negative", NULL) ) {
        for ( bound in c("deviation from parallel trends", "deviation from linear trend") ) {
          BC_DeltaRM_RobustResults <-
            createSensitivityResults_relativeMagnitudes(betahat               = BCdata_EventStudy$betahat,
                                                        sigma                 = BCdata_EventStudy$sigma,
                                                        numPrePeriods         = BC_numPrePeriods,
                                                        numPostPeriods        = BC_numPostPeriods,
                                                        l_vec                 = BC_l_vec,
                                                        gridPoints            = 100,
                                                        grid.ub               = 1,
                                                        grid.lb               = -1,
                                                        bound                 = bound,
                                                        method                = method,
                                                        biasDirection         = biasDirection,
                                                        Mbarvec               = seq(from=0, to=1, by=0.5))
          print(c(method, monotonicityDirection, biasDirection, bound, BC_DeltaRM_RobustResults))
          expect_silent(BC_DeltaRM_RobustResults)
        }
      }
    }
  })
} else {
  print("HonestDiD options run was skipped")
}
