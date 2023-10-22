# DESCRIPTION =========================================================
# Author: Ashesh Rambachan <asheshr@g.harvard.edu>
#
#  This script contains functions to implement the methods
#  described in Rambachan & Roth (2021) for robust inference
#  in difference-in-differences and event study designs.
#
#  This script contains utility functions to construct
#  shape and sign restrictions on Delta.


.create_A_M <- function(numPrePeriods, numPostPeriods,
                        monotonicityDirection, postPeriodMomentsOnly) {
  # This function creates a matrix so that A \delta <= 0 implies
  # delta is increasing/decreasing depending on what direction is specified.
  # This is used to impose monotonicity restrictions on the user's choice of Delta.

  A_M = base::matrix(0, nrow = numPrePeriods+numPostPeriods, ncol=numPrePeriods+numPostPeriods)
  for(r in 1:(numPrePeriods-1)){
    A_M[r, r:(r+1)] <- base::c(1,-1)
  }
  A_M[numPrePeriods, numPrePeriods] <- 1
  if(numPostPeriods > 0){
    A_M[numPrePeriods + 1, numPrePeriods + 1] <- -1
    if (numPostPeriods > 1) {
      for (r in (numPrePeriods + 2):(numPrePeriods+numPostPeriods)) {
        A_M[r, (r-1):r] <- c(1,-1)
      }
    }
  }

  # If postPeriodMomentsOnly == TRUE, exclude moments that only involve pre-periods
  if(postPeriodMomentsOnly){
    postPeriodIndices <- (numPrePeriods +1):base::NCOL(A_M)
    prePeriodOnlyRows <- base::which( base::rowSums( A_M[ , postPeriodIndices] != 0 ) == 0 )
    A_M <- A_M[-prePeriodOnlyRows , ]
  }
  if (monotonicityDirection == "decreasing") {
    A_M <- -A_M
  } else if(monotonicityDirection != "increasing") {
    base::stop("direction must be 'increasing' or 'decreasing'")
  }
  base::return(A_M)
}

.create_A_B <- function(numPrePeriods, numPostPeriods, biasDirection) {
  # This function creates a matrix for the linear constraints that \delta satisfy a sign restriction.
  #
  # Inputs:
  #   numPrePeriods = number of pre-periods. This is an element of resultsObjects.
  #   numPostPeriods = number of post-periods. This is an element of resultsObjects.
  #   biasDirection = must be "positive" or "negative".
  A_B = -base::diag(numPrePeriods + numPostPeriods)
  A_B = A_B[(numPrePeriods+1):(numPrePeriods+numPostPeriods), ]

  if (biasDirection == "negative"){
    A_B <- -A_B
  } else if(biasDirection != "positive"){
    base::stop("Input biasDirection must equal either `positive' or `negative'")
  }
  base::return(A_B)
}
