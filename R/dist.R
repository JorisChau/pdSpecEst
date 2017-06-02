#' Various distance measures between HPD matrices
#'
#' @export
pdDist <- function(A, B, method = 'Riemannian'){

  if (!(is.TRUE(all.equal(dim(Ad), dim(B)) & (length(dim(A)) == 3)))) {
    stop("Incorrect input dimensions for arguments: 'A' and/or 'B',
             consult the function documentation for the requested inputs.")
  }
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

