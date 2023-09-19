\name{basisVector}
\alias{basisVector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a standard basis vector.
}
\description{
Creates a basis vector of length \code{size} with a 1 in the \code{index} position.
}
\usage{
# Create the index basis vector in R^size
basisVector(index, size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{index}{
    The index at which there should be a one. Default equals one.
  }
  \item{size}{
    The length of the vector. Default equals one.
  }
}
\value{
Returns a basis vector of length \code{size} with a 1 in the \code{index} position.
}
\examples{
  # Create the first basis vector in R^2
  basisVector(index = 1, size = 2)
  # Create the third basis vector in R^6
  basisVector(index = 3, size = 6)
}
