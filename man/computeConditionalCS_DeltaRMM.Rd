\name{computeConditionalCS_DeltaRMM}
\alias{computeConditionalCS_DeltaRMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes conditional and hybridized confidence set for \eqn{\Delta = \Delta^{RMM}(Mbar)}.
}
\description{
Computes the conditional confidence set and hybridized confidence set for \eqn{\Delta = \Delta^{RMM}(Mbar)}. The set \eqn{\Delta^{RMM}(Mbar)} adds an additional shape restriction to \eqn{\Delta^{RM}(Mbar)} that restricts the underlying trend to be monotone. It may either be increasing (\eqn{\delta_t \ge \delta_{t-1}}) or decreasing (\eqn{\delta_t \le \delta_{t-1}}).
}
\usage{
computeConditionalCS_DeltaRMM(betahat, sigma, numPrePeriods, numPostPeriods,
                         l_vec = .basisVector(index = 1,
                                              size = numPostPeriods),
                                              Mbar = 0,
                         alpha = 0.05, hybrid_flag = "LF", hybrid_kappa = alpha/10,
                         returnLength = FALSE, postPeriodMomentsOnly = TRUE,
                         monotonicityDirection = "increasing",
                         gridPoints=10^3, grid.ub = NA, grid.lb = NA, seed = 0)
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
  Number of pre-periods.
  }
  \item{numPostPeriods}{
  Number of post-periods.
  }
  \item{l_vec}{
  Vector of length numPostPeriods that describes the scalar parameter of interest, theta = l_vec'tau. Default equals to first basis vector, (1, 0, ..., 0)
  }
  \item{Mbar}{
  Tuning parameter Mbar for \eqn{\Delta^{RM}(Mbar)} that governs how different the maximal pre-period violation of parallel trends may be from the post-period differential trend. Default sets Mbar = 0. See Section 2.3.2 of Rambachan & Roth (2021) for more details.
  }
  \item{alpha}{
  Desired level of the confidence set. Default equals 0.05 (corresponding to 95\% confidence interval)
  }
  \item{hybrid_flag}{
  Flag for whether user wishes to compute a hybridized confidence set. "ARP" specifies the conditional confidence set "LF" specifies the conditional least-favorable confidence set. The conditional FLCI hybrid confidence set is not available for \eqn{\Delta^{RM}(Mbar)} since the FLCI is infinite length for this choice of \eqn{\Delta}. See Section 3.3 and Section 5.3 of Rambachan & Roth (2021) for details. Default equals "LF".
  }
  \item{hybrid_kappa}{
  Desired first-stage size of hybridized confidence set. Only specify this value if the user wishes to compute a hybridized confidence set. Default equals alpha/10. If user specifies hybrid_flag = "ARP", set this value to NULL.
  }
  \item{returnLength}{
Logical value. If TRUE, function only returns the length of the robust confidence. If FALSE, function returns dataframe that contains a grid of possible parameter values and a vector of zeros and ones associated with each value in the grid (one denotes that the grid value lies in the confidence set and zero denotes that the grid value does not fall within the confidence set.) Default equals FALSE.
  }
  \item{postPeriodMomentsOnly}{
  Logical value. If TRUE, function excludes moments for \eqn{\Delta^{RM}(Mbar)} that only include pre-period coefficients. Default equals TRUE.
  }
  \item{monotonicityDirection}{
Specifies direction of monotonicity restriction. If "increasing", underlying trend specified to be increasing, \eqn{\delta_t \ge \delta_{t-1}}. If "decreasing", underlying trend specified to be decreasing \eqn{\delta_t \le \delta_{t-1}}. Default equals "increasing."
  }
  \item{gridPoints}{
  Number of grid points used in test inversion step. Default equals 1000.
  }
  \item{grid.ub}{
  Upper bound of grid for test inversion. The user should only specify this if she wishes to manually specify the upper bound of the grid. Default equals NA and sets grid upper bound to equal the upper bound of the identified set under parallel trends plus 20*standard deviation of the point estimate, l_vec'betahat.
  }
  \item{grid.lb}{
  Lower bound of grid for test inversion. The user should only specify this if she wishes to manually specify the upper bound of the grid. Default equals NA sets grid lower bound to equal the lower bound of the identified set under parallel trends minus 20*standard deviation of the point estimate, l_vec'betahat.
  }
  \item{seed}{
  Random seed for internal computations; included for reproducibility.
  }
}
\value{
If returnLength equals TRUE, function returns a scalar that equals the length of the confidence interval. If returnLength equals FALSE, function returns a dataframe with columns
  \item{grid}{
    Vector of grid values used to construct the confidence interval by test inversion.
  }
  \item{accept}{
    Vector of zeros-ones associated with grid values, where one denotes a grid value that falls within the confidence interval and zero denotes a grid value that falls outside the confidence interval.
  }
}
\references{
Rambachan, Ashesh and Jonathan Roth. "An Honest Approach to Parallel Trends." 2021.
}
\author{
 Ashesh Rambachan
}
