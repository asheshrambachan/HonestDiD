\name{createSensitivityPlot_relativeMagnitudes}
\alias{createSensitivityPlot_relativeMagnitudes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constructs sensitivity plot for \eqn{\Delta = \Delta^{RM}(Mbar)}, \eqn{\Delta^{SDRM}{Mbar}} and their variants that incorporate additional shape or sign restrictions.
}
\description{
  This function constructs sensitivity plots that examine how the robust confidence sets change as the parameter Mbar varies for \eqn{\Delta = \Delta^{RM}(Mbar)}, \eqn{\Delta^{SDRM}(Mbar)} and their variants that incorporate additional shape or sign restrictions. Similar plots are constructed in Section 6 of Rambachan & Roth (2021).
}
\usage{
createSensitivityPlot_relativeMagnitudes(robustResults, originalResults,
                                         rescaleFactor = 1, maxMbar = Inf,
                                         add_xAxis = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{robustResults}{
  Dataframe that contains the upper/lower bounds of robust confidence sets for each choice of Mbar. Contains columns: method --  Method of constructing robust confidence set, lb -- Lower bound of robust confidence set, ub -- Upper bound of robust confidence set, Mbar -- M values associated with each robust confidence set.
  }
  \item{originalResults}{
  Dataframe that contains the original confidence set for the parameter of interest. Contains columns: method --  Method of constructing confidence set (e.g., "Original"), lb -- Lower bound of confidence set, ub -- Upper bound of confidence set.
  }
  \item{rescaleFactor}{
  Scalar that is used to rescale the user specified choices of M and the upper/lower bounds of the confidence sets. Default equals one.
  }
  \item{maxMbar}{
  Scalar that specifies the maximum Mbar value to plot in the sensitivity plot. Default equals infinity (no truncation).
  }
  \item{add_xAxis}{
  Logical specifying whether to plot the x-axis in the sensitivity plot. Default equals TRUE.
  }
}
\value{
Returns ggplot object of the sensitivity plot.
}
\references{
  Rambachan, Ashesh and Jonathan Roth. "An Honest Approach to Parallel Trends." 2019.
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
  robustResults   <- do.call(createSensitivityResults_relativeMagnitudes, kwargs)
  originalResults <- do.call(constructOriginalCS, kwargs)
  createSensitivityPlot_relativeMagnitudes(robustResults, originalResults)
}
}
