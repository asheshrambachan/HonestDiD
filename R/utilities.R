# Other helper functions ----
.selectionMat <- function(selection, size, select = "columns"){
  if(select == "rows"){
    m <- matrix(0, nrow = length(selection), ncol = size)

    m[ 1:length(selection), selection] <- diag(length(selection))
  }else{
    m <- matrix(0, ncol = length(selection), nrow = size)

    m[ selection , 1:length(selection)] <- diag(length(selection))
  }
  return(m)
}

.LeeCFN <- function(eta, Sigma){
  c = Sigma %*% eta / as.numeric( t(eta) %*% Sigma %*% eta )
  return(c)
}

.VLoVUpFN <- function(eta, Sigma, A, b, z){
  c = .LeeCFN(eta,Sigma)

  objective <- (b - A %*% z) / (A %*% c)

  ACNegativeIndex <- which( (A %*% c) < 0 )
  ACPositiveIndex <- which( (A %*% c) > 0 )

  if( length(ACNegativeIndex) == 0){
    VLo <- -Inf
  }else{
    VLo <- max(objective[ACNegativeIndex])
  }

  if( length(ACPositiveIndex) == 0){
    VUp <- Inf
  }else{
    VUp <- min(objective[ACPositiveIndex])
  }

  return(c(VLo,VUp))
}

basisVector <- function(index = 1, size = 1){
  v <- matrix(0, nrow = size, ncol = 1)
  v[index] = 1
  return(v)
}
