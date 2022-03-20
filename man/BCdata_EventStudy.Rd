\docType{data}
\name{BCdata_EventStudy}
\alias{BCdata_EventStudy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Event study estimates from baseline event study specification on profits in Benzarti & Carloni (2019). See discussion in Section 6.1 of Rambachan & Roth (2021).
}
\description{
This list contains the event study estimates from baseline event study specification on profits in Benzarti & Carloni (2019). See discussion in Section 6.1 of Rambachan & Roth (2021).
}
\format{
  \describe{A list, containing 7 objects:
    \item{betahat}{Vector of estimated event study coefficients.}
    \item{sigma}{Estimated variance-covariance matrix.}
    \item{timeVec}{Vector that contains the time periods associated with the event study coefficients.}
    \item{referencePeriod}{Reference period that is normalized to zero.}
    \item{prePeriodIndices}{Vector containing elements of timeVec that correspond to the pre-periods.}
    \item{postPeriodIndices}{Vector containing elements of timeVec that correspond to the post-periods.}
  }
}
