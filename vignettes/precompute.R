library(HonestDiD)

# data('VignetteResults', package="HonestDiD")
data('BCdata_EventStudy', package="HonestDiD")

# Number of pre-periods
BC_numPrePeriods = length(BCdata_EventStudy$prePeriodIndices)
BC_numPostPeriods = length(BCdata_EventStudy$postPeriodIndices)

# Create l_vec to define the parameter of interest, the first post-treatment period.
BC_l_vec = basisVector(index = 1, size = BC_numPostPeriods)

# Construct robust confidence intervals for Delta^{SDRM}(Mbar) for first post-treatment period.
# We specify 100 gridPoints over [-1, 1] for the underlying test inversion to construct the robust confidence set.
# Users may wish to leave this at the default values.
BC_DeltaSDRM_RobustResults =
  createSensitivityResults_relativeMagnitudes(betahat = BCdata_EventStudy$betahat,
                                              sigma = BCdata_EventStudy$sigma,
                                              bound = "deviation from linear trend",
                                              numPrePeriods = BC_numPrePeriods,
                                              numPostPeriods = BC_numPostPeriods,
                                              l_vec = BC_l_vec,
                                              Mbarvec = seq(from = 0, to = 2, by = 0.5),
                                              gridPoints = 100, grid.lb = -1, grid.ub = 1)

head(BC_DeltaSDRM_RobustResults)

# Construct dataframe with OLS confidence interval for theta.
BC_OriginalResults = constructOriginalCS(betahat = BCdata_EventStudy$betahat,
                                         sigma = BCdata_EventStudy$sigma,
                                         numPrePeriods = BC_numPrePeriods,
                                         numPostPeriods = BC_numPostPeriods,
                                         l_vec = BC_l_vec )

# Construct sensitivity plot.
BC_DeltaSDRM_SensitivityPlot =
  createSensitivityPlot_relativeMagnitudes(robustResults = BC_DeltaSDRM_RobustResults,
                                           originalResults = BC_OriginalResults)

BC_DeltaSDRM_SensitivityPlot

# Construct robust confidence intervals for Delta^{SDNB}(M) for first post-treatment period
BC_DeltaSDNB_RobustResults = createSensitivityResults(betahat = BCdata_EventStudy$betahat,
                                                      sigma = BCdata_EventStudy$sigma,
                                                      numPrePeriods = BC_numPrePeriods,
                                                      numPostPeriods = BC_numPostPeriods,
                                                      l_vec = BC_l_vec,
                                                      Mvec = seq(from = 0, to = 0.3, by = 0.1),
                                                      biasDirection = "negative")
BC_DeltaSDNB_SensitivityPlot = createSensitivityPlot(robustResults = BC_DeltaSDNB_RobustResults,
                                                     originalResults = BC_OriginalResults)
BC_DeltaSDNB_SensitivityPlot

data('LWdata_EventStudy', package = "HonestDiD")

# Number of pre-periods
LW_numPrePeriods = length(LWdata_EventStudy$prePeriodIndices)
LW_numPostPeriods = length(LWdata_EventStudy$postPeriodIndices)

# Create l_vec corresponding with 15 years of exposure
# Reference is -2 years of exposure, so want effect 17 pds later
LW_l_vec = basisVector(15 - (-2), LW_numPostPeriods)

# Construct robust confidence intervals for Delta^{SD}(M) for 15 years of exposure
LW_DeltaSD_RobustResults = createSensitivityResults(betahat = LWdata_EventStudy$betahat,
                                                    sigma = LWdata_EventStudy$sigma,
                                                    numPrePeriods = LW_numPrePeriods,
                                                    numPostPeriods = LW_numPostPeriods,
                                                    l_vec = LW_l_vec,
                                                    Mvec = seq(from = 0, to = 0.04, by = 0.005))
head(LW_DeltaSD_RobustResults)

# Construct dataframe with OLS confidence interval for theta
LW_OriginalResults = constructOriginalCS(betahat = LWdata_EventStudy$betahat,
                                         sigma = LWdata_EventStudy$sigma,
                                         numPrePeriods = LW_numPrePeriods,
                                         numPostPeriods = LW_numPostPeriods,
                                         l_vec = LW_l_vec )

# Construct sensitivity plot
LW_DeltaSD_SensitivityPlot = createSensitivityPlot(robustResults = LW_DeltaSD_RobustResults,
                                                               originalResults = LW_OriginalResults)
LW_DeltaSD_SensitivityPlot

# Construct robust confidence intervals for Delta^{SDD}(M)
LW_DeltaSDD_RobustResults = createSensitivityResults(betahat = LWdata_EventStudy$betahat,
                                                     sigma = LWdata_EventStudy$sigma,
                                                     monotonicityDirection = "decreasing",
                                                     numPrePeriods = LW_numPrePeriods,
                                                     numPostPeriods = LW_numPostPeriods,
                                                     l_vec = LW_l_vec,
                                                     Mvec = seq(from = 0, to = 0.04, by = 0.005))

# Construct sensitivity plot
LW_DeltaSDD_SensitivityPlot = createSensitivityPlot(robustResults = LW_DeltaSDD_RobustResults,
                                                    originalResults = LW_OriginalResults)
LW_DeltaSDD_SensitivityPlot

LW_lowerBound_M = DeltaSD_lowerBound_Mpre(betahat = LWdata_EventStudy$betahat,
                                          sigma = LWdata_EventStudy$sigma,
                                          numPrePeriods = LW_numPrePeriods)
LW_upperBound_M = DeltaSD_upperBound_Mpre(betahat = LWdata_EventStudy$betahat,
                                          sigma = LWdata_EventStudy$sigma,
                                          numPrePeriods = LW_numPrePeriods)

# Load in LWdata_RawData.dta
LWdata_RawData = haven::read_dta(system.file("extdata", "LWdata_RawData.dta",
                                             package = "HonestDiD"))

# Estimate event study using lfe package
EmpFemale.EventStudy = lfe::felm(emp ~ rtESV13 + rtESV14 + rtESV15 +
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

# Extract coefficients of regression associated with event study coefficients
coefIndex = which(grepl(x = dimnames(EmpFemale.EventStudy$coefficients)[[1]],
                        pattern = "rtESV"))
betahat = EmpFemale.EventStudy$beta[coefIndex, ]

# Extract estimated variance-covariance matrix of event study coefficients
sigma = EmpFemale.EventStudy$clustervcv[coefIndex, coefIndex]

#Rescale by 100 so that results will be in units of percentage points
betahat = 100 * betahat
sigma = 100^2 * sigma

# Construct vector of event times and the scalar reference period
timeVec = c(seq(from = -11, to = -3, by = 1), seq(from = -1, to = 21, by = 1))
referencePeriod = -2
postPeriodIndices = which(timeVec > -2)
prePeriodIndices = which(timeVec < -2)

# Construct standard errors associated with event study coefficients
stdErrors = summary(EmpFemale.EventStudy)$coefficients[coefIndex,2]

# Create list containing objects produced by the event study
LWdata_EventStudy = list(
  betahat = betahat,
  sigma = sigma,
  timeVec = timeVec,
  referencePeriod = referencePeriod,
  prePeriodIndices = prePeriodIndices,
  postPeriodIndices = postPeriodIndices,
  stdErrors = stdErrors
)

VignetteResults <- list(BC_DeltaSDRM_RobustResults = BC_DeltaSDRM_RobustResults,
                        BC_OriginalResults         = BC_OriginalResults,
                        BC_DeltaSDNB_RobustResults = BC_DeltaSDNB_RobustResults,
                        LW_DeltaSD_RobustResults   = LW_DeltaSD_RobustResults,
                        LW_DeltaSDD_RobustResults  = LW_DeltaSDD_RobustResults)
save(VignetteResults, file="data/VignetteResults.rda")
