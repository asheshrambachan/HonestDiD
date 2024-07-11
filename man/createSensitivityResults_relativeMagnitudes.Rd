\name{createSensitivityResults_relativeMagnitudes}
\alias{createSensitivityResults_relativeMagnitudes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Constructs robust confidence intervals for \eqn{\Delta = \Delta^{RM}(Mbar)}, \eqn{\Delta^{SDRM}(Mbar)} and their variants that incorporate shape or sign restrictions for a vector of possible Mbar values.
}
\description{
  Constructs robust confidence intervals for \eqn{\Delta = \Delta^{RM}(Mbar)}, \eqn{\Delta^{SDRM}(Mbar)} and their variants that incorporate shape or sign restrictions for a vector of possible Mbar values. By default, the function constructs sensitivity results for \eqn{\Delta^{RM}(Mbar)} and its variants. The confidence sets are constructed through test inversion.
}
\usage{
createSensitivityResults_relativeMagnitudes(betahat, sigma,
                                            numPrePeriods, numPostPeriods,
                                            bound = "deviation from parallel trends",
                                            method = "C-LF",
                                            Mbarvec = NULL,
                                            l_vec = .basisVector(index = 1,
                                                                 size = numPostPeriods),
                                            monotonicityDirection = NULL,
                                            biasDirection = NULL,
                                            alpha = 0.05,
                                            gridPoints = 10^3,
                                            grid.ub = NA,
                                            grid.lb = NA,
                                            parallel = FALSE,
                                            seed = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{betahat}{
  Vector of estimated event study coefficients.
  }
  \item{sigma}{
  Covariance matrix of event study coefficients.
  }
  \item{numPrePeriods}{
  Number of pre-periods. If user selects bound = "deviation from linear trends" (\eqn{Delta^{SDRM}} as base choice of Delta), then numPrePeriods must be greater than one. See details for further explanation.
  }
  \item{numPostPeriods}{
  Number of post-periods.
  }
  \item{bound}{
  String that specifies the base choice of Delta (to which additional sign and shape restrictions will be incorporated if specified by the user). This must be either "deviation from parallel trends" or "deviation from linear trend". If bound equals "deviation from parallel trends", then the function will select \eqn{\Delta^{RM}(Mbar)} as the base choice of \eqn{\Delta}. If bound equals "deviation from linear trends", then the function will select \eqn{\Delta^{SDRM}} as the base choice of \eqn{\Delta}. By default, this is set to "deviation from parallel trends". See Section 2.3.1 and 2.3.2 of Rambachan & Roth (2021) for a discussion of these choices of \eqn{\Delta}.
  }
  \item{method}{
  String that specifies the choice of method for constructing robust confidence intervals. This must be either "Conditional", or "C-LF" (conditional least-favorable hybrid). Default equals "C-LF" and the function automatically sets method to be "C-LF" based on the recommendations in Rambachan & Roth (2021).
  }
  \item{Mbarvec}{
  Vector of Mbar values for which the user wishes to construct robust confidence intervals. If NULL, the function constructs a grid of length 10 that starts at Mbar = 0 and ends at Mbar = 2. Default equals null.
  }
  \item{l_vec}{
  Vector of length numPostPeriods that describes the scalar parameter of interest, theta = l_vec'tau. Default equals to first basis vector, (1, 0, ..., 0)
  }
  \item{biasDirection}{
  This must be specified if the user wishes to add an additional bias restriction to base choice of Delta. If "positive", bias is restricted to be positive, delta >= 0. If "negative", bias is restricted to be negative, delta <= 0. Default equals NULL.
  }
  \item{monotonicityDirection}{
  This must be specified if the user wishes to add an additional monotonicity restriction to base choice of Delta. If "increasing", underlying trend specified to be increasing, \eqn{\delta_t \ge \delta_{t-1}}. If "decreasing", underlying trend specified to be decreasing \eqn{\delta_t \le \delta_{t-1}}. Default equals NULL.
  }
  \item{alpha}{
  Desired size of the robust confidence sets. Default equals 0.05 (corresponding to 95\% confidence interval)
  }
  \item{parallel}{
  Logical to indicate whether the user would like to construct the robust confidence intervals in parallel. This uses the Foreach package and doParallel package. Default equals FALSE.
  }
  \item{gridPoints}{
  Number of grid points used for the underlying test inversion. Default equals 1000. User may wish to change the number of grid points for computational reasons.
  }
  \item{grid.ub}{
  Upper bound of grid used for underlying test inversion. Default sets grid.ub to be equal to twenty times the standard deviation of the estimated target parameter, l_vec * betahat. User may wish to change the upper bound of the grid to suit their application.
  }
  \item{grid.lb}{
  Lower bound of grid used for underlying test inversion. Default sets grid.lb to be equal to negative twenty times the standard deviation of the estimated target parameter, l_vec * betahat. User may wish to change the lower bound of the grid to suit their application.
  }
  \item{seed}{
  Random seed for internal computations; included for reproducibility.
  }
}
\details{
Note: If the user specifies bound = "deviation from linear trends", then numPrePeriods must be greater than one. By specifying bound = "deviation from linear trends", then the function selects \eqn{\Delta^{SDRM}} as the base choice of \eqn{\Delta}. As discussed in Section 2.3.2 of Rambachan & Roth (2021), \eqn{\Delta^{SDRM}} uses observed non-linearities in the pre-treatment difference in trends to bound the possible non-linearities in the post-treatment difference in trends. This is only possible if there are multiple pre-treatment periods (i.e., numPrePeriods > 1).
}
\value{
Returns a dataframe with columns
  \item{lb}{
  Lower bound of robust confidence sets.
  }
  \item{ub}{
  Upper bound of robust confidence sets.
  }
  \item{method}{
  Method for constructing robust confidence sets
  }
  \item{Delta}{
  The set Delta that was specified.
  }
  \item{M}{
  Values of M associated with each robust confidence set.
  }
}
\references{
Rambachan, Ashesh and Jonathan Roth. "An Honest Approach to Parallel Trends." 2021.
}
\author{
Ashesh Rambachan
}
\examples{
\donttest{
  # Simple use case. For more detailed examples,
  # see <https://github.com/asheshrambachan/HonestDiD#honestdid>
  kwargs <- list(betahat        = BCdata_EventStudy$betahat,
                 sigma          = BCdata_EventStudy$sigma,
                 numPrePeriods  = length(BCdata_EventStudy$prePeriodIndices),
                 numPostPeriods = length(BCdata_EventStudy$postPeriodIndices),
                 alpha          = 0.05)
  do.call(createSensitivityResults_relativeMagnitudes, kwargs)
}
}
