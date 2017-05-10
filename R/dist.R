#' Vairous distance measures between HPD matrices
#'
#' @export
pdDist <- function(A, B, method = 'Riemannian'){

  method <- match.arg(method, c('Riemannian', 'logEuclidean', 'Cholesky', 'Euclidean', 'Procrustes'))
  d <- nrow(A)

  if(method == 'Riemannian'){
    dd <- RiemmDist(A, B)
  }
  if(method == 'logEuclidean'){
    dd <- NormF(Logm(diag(d), A) - Logm(diag(d), B))
  }
  if(method == 'Cholesky'){
    dd <- NormF(Chol(A) - Chol(B))
  }
  if(method == 'Euclidean'){
    dd <- NormF(A - B)
  }
  if(method == 'Procrustes'){
    l1 <- Sqrt(A)
    l2 <- Sqrt(B)
    dd <- sqrt(NormF(l1)^2 + NormF(l2)^2 - 2 * sum(svd(t(Conj(l2)) %*% l1)$d))
  }
  dd
}

