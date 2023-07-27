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

.stopIfNotSymmPSD <- function (sigma) {
  # Check that sigma is PSD
  if ( base::isSymmetric(sigma) ) {
    if ( base::any(base::eigen(sigma, TRUE, only.values=TRUE)$values < 0) ) {
      base::stop("sigma must be a positive semi-definite symmetric matrix")
    }
  } else {
    base::stop("sigma must be a symmetric matrix")
  }
}
