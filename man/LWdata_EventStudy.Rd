\docType{data}
\name{LWdata_EventStudy}
\alias{LWdata_EventStudy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Event study estimates from baseline female specification on employment in Lovenheim & Willen (2019). See discussion in Section 6.2 of Rambachan & Roth (2021).
}
\description{
  This list contains the event study estimates from the baseline female specification on employment in Lovenheim & Willen (2019). See discussion in Section 6.2 of Rambachan & Roth (2021).
}
\format{
  \describe{A list, containing 7 objects:
    \item{betahat}{Vector of estimated event study coefficients.}
    \item{sigma}{Estimated variance-covariance matrix.}
    \item{timeVec}{Vector that contains the time periods associated with the event study coefficients.}
    \item{referencePeriod}{Reference period that is normalized to zero.}
    \item{prePeriodIndices}{Vector containing elements of timeVec that correspond to the pre-periods.}
    \item{postPeriodIndices}{Vector containing elements of timeVec that correspond to the post-periods.}
    \item{stdErrors}{Vector of standard errors associated with estimated event study coefficients}
  }
}
