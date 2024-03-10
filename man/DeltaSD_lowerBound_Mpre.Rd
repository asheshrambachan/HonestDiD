\name{DeltaSD_lowerBound_Mpre}
\alias{DeltaSD_lowerBound_Mpre}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct lower bound for M for \eqn{\Delta = \Delta^{SD}(M)} based on observed pre-period coefficients.
}
\description{
  Constructs a lower bound for M using the observed pre-period coefficients. It constructs a one-sided confidence interval for the maximal second difference of the observed pre-period using the conditional test developed in Andrews, Roth & Pakes (2019). The number of pre-periods (not including the reference period) must be larger than or equal to two.
}
\usage{
DeltaSD_lowerBound_Mpre(betahat, sigma, numPrePeriods,
                        alpha = 0.05, grid.ub = NA, gridPoints = 1000)
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
  Number of pre-periods. Must be larger than or equal to two.
  }
  \item{alpha}{
  Desired size of the one-sided confidence set. Default equals 0.05 (corresponding to 95\% confidence interval)
  }
  \item{grid.ub}{
  Upper bound of grid of values of M that is used to construct the confidence interval by test inversion. Default equals NA and the upper bound of the grid is set equal to three times the maximum standard error of the observed pre-period event-study coefficients.
  }
  \item{gridPoints}{
  Number of points to include in the grid that is used to construct the confidence interval by test inversion. Default equals 1000 points.
  }
}
\value{
Returns a scalar that equals the lower bound of a one-sided confidence interval for the maximal second difference of the observed pre-period coefficients.
}
\references{
Andrews, Isaiah, Jonathan Roth and Ariel Pakes. "Inference for Linear Conditional Moment Inequalities." 2019.
Rambachan, Ashesh and Jonathan Roth. "An Honest Approach to Parallel Trends." 2021.
}
\author{
Ashesh Rambachan
}
