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
basisVector(index = 3, size = 6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
%%     ~~Describe \code{index} here~~
  \item{index}{ The index at which there should be a 1.

}
  \item{size}{
%%     ~~Describe \code{size} here~~
The length of the vector.

}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a basis vector of length \code{size} with a 1 in the \code{index} position.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
