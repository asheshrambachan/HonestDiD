\name{constructOriginalCS}
\alias{constructOriginalCS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Constructs original confidence interval for parameter of interest, theta = l_vec'tau.
}
\description{
  Constructs original confidence interval for parameter of interest, theta = l_vec'tau using the user-specified estimated event study coefficients and variance-covariance matrix.
}
\usage{
constructOriginalCS(betahat, sigma,
                    numPrePeriods, numPostPeriods,
                    l_vec = .basisVector(index = 1, size = numPostPeriods),
                    alpha = 0.05)
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
  \item{alpha}{
  Desired size of the robust confidence sets. Default equals 0.05 (corresponding to 95\% confidence interval)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
