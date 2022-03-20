\name{basisVector}
\alias{basisVector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a standard basis vector.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Creates a basis vector of length \code{size} with a 1 in the \code{index} position.
}
\usage{
# Create the third basis vector in R^6
basisVector(index = 1, size = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
%%     ~~Describe \code{index} here~~
  \item{index}{ The index at which there should be a one. Default equals one.

}
  \item{size}{
%%     ~~Describe \code{size} here~~
The length of the vector. Default equals one.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a basis vector of length \code{size} with a 1 in the \code{index} position.
}

