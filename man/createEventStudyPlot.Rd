\name{createEventStudyPlot}
\alias{createEventStudyPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constructs event study plot
}
\description{
Constructs event study plot using the estimated event study coefficients and standard errors.
}
\usage{
createEventStudyPlot(betahat, stdErrors = NULL, sigma = NULL,
                     numPrePeriods, numPostPeriods, alpha = 0.05, timeVec,
                     referencePeriod, useRelativeEventTime = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{betahat}{
   Vector of estimated event study coefficients.
  }
  \item{stdErrors}{
  Vector of standard errors associated with the estimated event study coefficients. Default equals NULL. Either stdErrors or sigma must be specified by the user. If stdErrors is not specified but sigma is, the stdErrors are set to equal the square root of the diagonal elements of sigma.
  }
  \item{sigma}{
   Covariance matrix of event study coefficients. Default equals NULL. Either stdErrors or sigma must be specified by the user.
  }
  \item{numPrePeriods}{
  Number of pre-periods.
  }
  \item{numPostPeriods}{
  Number of post-periods.
  }
  \item{alpha}{
  Desired size of confidence intervals. Default = 0.05.
  }
  \item{timeVec}{
  Vector that contains the time periods associated with the event study coefficients. This vector should not include the reference period that is normalized to zero.
  }
  \item{referencePeriod}{
  Scalar that contains the time period associated with the reference period.
  }
  \item{useRelativeEventTime}{
  Logical that specifies whether user would like the plot to be in relative event time (normalizes the reference period to be zero). Default equals FALSE.
  }
}
\value{
Returns ggplot object of the event study plot.
}
\references{
Rambachan, Ashesh and Jonathan Roth. "An Honest Approach to Parallel Trends." 2021.
}
\author{
Ashesh Rambachan
}
\examples{
\donttest{
  # Simple use case; for more detailed examples,
  # see <https://github.com/asheshrambachan/HonestDiD#honestdid>
  createEventStudyPlot(betahat         = BCdata_EventStudy$betahat,
                       sigma           = BCdata_EventStudy$sigma,
                       numPrePeriods   = length(BCdata_EventStudy$prePeriodIndices),
                       numPostPeriods  = length(BCdata_EventStudy$postPeriodIndices),
                       alpha           = 0.05,
                       timeVec         = BCdata_EventStudy$timeVec,
                       referencePeriod = BCdata_EventStudy$referencePeriod)
}
}
