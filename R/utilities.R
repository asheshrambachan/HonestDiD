# Other helper functions ----
.selectionMat <- function(selection, size, select = "columns"){
  if(select == "rows"){
    m <- base::matrix(0, nrow = base::length(selection), ncol = size)

    m[ 1:base::length(selection), selection] <- base::diag(base::length(selection))
  }else{
    m <- base::matrix(0, ncol = base::length(selection), nrow = size)

    m[ selection , 1:base::length(selection)] <- base::diag(base::length(selection))
  }
  base::return(m)
}

.LeeCFN <- function(eta, Sigma){
  c = Sigma %*% eta / base::as.numeric( base::t(eta) %*% Sigma %*% eta )
  base::return(c)
}

.VLoVUpFN <- function(eta, Sigma, A, b, z){
  c = .LeeCFN(eta,Sigma)

  objective <- (b - A %*% z) / (A %*% c)

  ACNegativeIndex <- base::which( (A %*% c) < 0 )
  ACPositiveIndex <- base::which( (A %*% c) > 0 )

  if( base::length(ACNegativeIndex) == 0){
    VLo <- -Inf
  }else{
    VLo <- base::max(objective[ACNegativeIndex])
  }

  if( base::length(ACPositiveIndex) == 0){
    VUp <- Inf
  }else{
    VUp <- base::min(objective[ACPositiveIndex])
  }

  base::return(base::c(VLo,VUp))
}

basisVector <- function(index = 1, size = 1){
  v <- base::matrix(0, nrow = size, ncol = 1)
  v[index] = 1
  base::return(v)
}

.warnIfNotSymmPSD <- function (sigma) {
  # Check that sigma is PSD
  if ( !base::isSymmetric(sigma, check.attributes=FALSE) ) {
    # sigma <- (sigma + t(sigma)) / 2
    base::warning(base::sprintf("matrix sigma not exactly symmetric (largest asymmetry was %g)",
                                base::max(base::abs(sigma - t(sigma)))))
  }

  lambda <- base::eigen(sigma, TRUE, only.values=TRUE)$values
  if ( base::any(lambda < 0) ) {
    base::warning(base::sprintf("matrix sigma not numerically positive semi-definite (smallest eigenvalue was %g)",
                                base::min(lambda)))
  }

  # return(sigma)
}

.stopIfNotConformable <- function (betahat, sigma, numPrePeriods, numPostPeriods, l_vec) {
  betaDim <- c(base::NROW(betahat), base::NCOL(betahat))

  betaL <- base::max(betaDim)
  beta1 <- base::min(betaDim)
  if ( beta1 > 1 ) {
      base::stop(base::sprintf("expected a vector but betahat was %d by %d", betaDim[1], betaDim[2]))
  }

  sigmaR <- base::NROW(sigma)
  sigmaC <- base::NCOL(sigma)
  if ( sigmaR != sigmaC ) {
      base::stop(base::sprintf("expected a square matrix but sigma was %d by %d", sigmaR, sigmaC))
  }

  if ( sigmaR != betaL ) {
      base::stop(base::sprintf("betahat (%d by %d) and sigma (%d by %d) were non-conformable",
                 betaDim[1], betaDim[2], sigmaR, sigmaC))
  }

  numPeriods <- numPrePeriods + numPostPeriods
  if ( numPeriods != betaL ) {
      base::stop(base::sprintf("betahat (%d by %d) and pre + post periods (%d + %d) were non-conformable",
                 betaDim[1], betaDim[2], numPrePeriods, numPostPeriods))
  }

  if ( base::length(l_vec) != numPostPeriods ) {
      base::stop(base::sprintf("l_vec (length %d) and post periods (%d) were non-conformable",
                 base::length(l_vec), numPostPeriods))
  }
}

# Custom error messages
.CustomErrorHandling <- function(code, custom) {
  base::tryCatch({
    code
  }, error = function(e) {
    base::message(custom)
    base::message("\tcall:  ", base::paste0(rlang::expr_deparse(e$call), collapse='\n\t'))
    base::message("\terror: ", e$message)
    base::message("Traceback:"); traceback()
  })
}
