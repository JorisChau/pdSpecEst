#' Zonoid and Geodesic Distance Depth for HPD matrices
#'
#' @importFrom ddalpha depth.zonoid
#' @export
pdDepth <- function(y = NULL, X, method = c('zonoid', 'gdd')){

  method <- match.arg(method, c('zonoid', 'gdd'))
  d <- dim(X)[1]
  if(length(dim(X)) == 3){
    n <- dim(X)[3]
  } else if(length(dim(X)) == 4){
    n <- dim(X)[4]
    N <- dim(X)[3]
  }
  E <- E_basis(d)

  ## Zonoid depth
  if(method == 'zonoid'){
    if(length(dim(X)) == 3){
      X.vec <- sapply(1:n, function(i) E_coeff(Logm(y, X[,,i]), E))
      depth <- depth.zonoid(t(as.matrix(rep(0, d^2))), t(X.vec))
    } else if(length(dim(X)) == 4){
      depth.t <- rep(NA, N)
      for(t in 1:N){
        X.vec <- sapply(1:n, function(i) E_coeff(Logm(y[,,t], X[,,t,i]), E))
        depth.t[t] <- depth.zonoid(t(as.matrix(rep(0, d^2))), t(X.vec))
      }
      depth <- mean(depth.t)
    }
  }

  ## Geodesic distance depth
  if(method == 'gdd'){
    if(!is.null(y)){
      if(length(dim(X)) == 3){
        depth <- exp(-mean(sapply(1:n, function(i) RiemmDist(y, X[,,i]))))
      } else if(length(dim(X)) == 4){
        depth <- exp(-mean(sapply(1:n, function(i) mean(sapply(1:N, function(t)
          RiemmDist(y[,,t], X[,,t,i]))))))
      }
    } else if(is.null(y)){
      if(length(dim(X)) == 3){
        dist <- matrix(0, nrow = n, ncol = n)
        for(i in 1:n){
          for(j in 1:i){
            if(j < i){
              dist[i,j] <- dist[j,i] <- RiemmDist(X[,,i], X[,,j])
            }
          }
        }
        depth <- exp(-colMeans(dist))
      } else if(length(dim(X)) == 4){
        depth <- NaN
      }
    }
  }
  return(depth)
}
