#' Compute distance between two HPD matrices
#'
#' \code{pdDist} calculates a distance between two Hermitian PD matrices.
#'
#' Available distance measures between two Hermitian PD matrices are (i) Riemannian distance (default) as in
#' (Bhatia, 2009, Chapter 6), (ii) log-Euclidean distance, the Euclidean distance between matrix logarithms,
#' (iii) Cholesky distance, the Euclidean distance between Cholesky decompositions, (iv) Euclidean distance,
#' (v) root-Euclidean distance and (vi) Procrustes distance as in (Dryden et al., 2009). In particular, \code{pdDist} generalizes the function
#' \code{shapes::distcov} to compute the distance between two symmetric positive definite matrices to the
#' distance between two Hermitian positive definite matrices.
#'
#' @param A,B Hermitian positive definite matrices (of equal dimension).
#' @param method the distance measure, one of \code{'Riemannian'},
#' \code{'logEuclidean'}, \code{'Cholesky'}, \code{'Euclidean'}, \code{'squareRoot'} or \code{'Procrustes'}. Defaults to \code{'Riemannian'}.
#'
#' @examples
#'  a <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
#'  A <- t(Conj(a)) %*% a
#'  b <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
#'  B <- t(Conj(b)) %*% b
#'  pdDist(A, B) ## Riemannian distance
#'
#' @references Bhatia, R. (2009). \emph{Positive Definite Matrices}. New Jersey: Princeton University Press.
#' @references Dryden, I.L., Koloydenko, A., Zhou, D. (2009). Non-Euclidean statistics for covariance matrices,
#' with applications to diffusion tensor imaging. \emph{Annals of Applied Statistics}, 3(3), 1102-1123.
#'
#' @export
pdDist <- function(A, B, method = "Riemannian") {

  if (!(isTRUE(all.equal(dim(A), dim(B)) & (dim(A)[1] == dim(A)[2]) & (length(dim(A)) == 2)))) {
    stop("Incorrect input dimensions for arguments: 'A' and/or 'B',
         consult the function documentation for the requested inputs.")
  }
  method <- match.arg(method, c("Riemannian", "logEuclidean", "Cholesky", "Euclidean", "rootEuclidean", "Procrustes"))
  d <- nrow(A)

  if (method == "Riemannian") {
    dd <- RiemmDist(A, B)
  }
  if (method == "logEuclidean") {
    dd <- NormF(Logm(diag(d), A) - Logm(diag(d), B))
  }
  if (method == "Cholesky") {
    dd <- NormF(Chol(A) - Chol(B))
  }
  if (method == "Euclidean") {
    dd <- NormF(A - B)
  }
  if (method == "rootEuclidean") {
    dd <- NormF(Sqrt(A) - Sqrt(B))
  }
  if (method == "Procrustes") {
    l1 <- Sqrt(A)
    l2 <- Sqrt(B)
    dd <- sqrt(NormF(l1)^2 + NormF(l2)^2 - 2 * sum(svd(t(Conj(l2)) %*% l1)$d))
  }
  return(dd)
}

#' Weighted geometric mean of HPD matrices
#'
#' \code{pdMean} calculates an (approximate) weighted geometric mean of \eqn{S} different
#' \eqn{(d \times d)}-dimensional Hermitian PD matrices based on the Riemannian metric by
#' the fast recursive algorithm in (Chau and von Sachs, 2017) or the slower but more accurate
#' gradient descent algorithm in (Pennec, 2006). By default, the unweighted geometric mean is computed.
#'
#' @param M a \eqn{(d,d,S)}-dimensional array of Hermitian PD matrices.
#' @param w an \eqn{S}-dimensional nonnegative weight vector, such that \code{sum(w) = 1}.
#' @param grad_desc a logical value deciding if the gradient descent algorithm be used, defaults to
#' \code{FALSE}.
#' @param max_iter maximum number of iterations in gradient descent algorithm, only used if
#' \code{isTRUE(grad_desc)}.
#' @param tol optional tolerance parameter in gradient descent algorithm, only used if
#' \code{isTRUE(grad_desc)}, defaults to \code{.Machine$double.eps}.
#' @param ... additional arguments for internal usage.
#'
#' @examples
#' m <- function(){
#'  X <- matrix(complex(real=rnorm(9), imaginary=rnorm(9)), nrow=3)
#'  t(Conj(X)) %*% X
#' }
#' M <- replicate(100, m())
#' z <- rnorm(100)
#' w <- abs(z)/sum(abs(z))
#' Ave <- pdMean(M, w)
#'
#' @references Chau, J. and von Sachs, R. (2017). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#' @references Pennec, X. (2006). Intrinsic statistics on Riemannian manifolds: Basic tools for geometric
#' measurements. \emph{Journal of Mathematical Imaging and Vision} 25(1), 127-154.
#'
#' @seealso \code{\link{Mid}}
#'
#' @export
pdMean <- function(M, w, grad_desc = F, max_iter = 1000, tol, ...) {

  ## Set parameters
  dots = list(...)
  metric = (if(is.null(dots$metric)) "Riemannian" else dots$metric)
  tol = (if(missing(tol)) NA else tol)
  w = (if(missing(w)) rep(1/dim(M)[3], dim(M)[3]) else w)
  d = dim(M)[1]

  if(dim(M)[3] == 1){
    Mean <- M[, , 1]
  } else{

    if(metric == "Riemannian"){

      ## Recursive algorithm
      Mean <- kMean(do.call(rbind, lapply(1:dim(M)[3], function(s) M[, , s])), w)

      ## Gradient descent algorithm
      if(grad_desc){
        tol <- (if(is.na(tol)) .Machine$double.eps else tol)
        if(!isTRUE(is.numeric(tol))){
          stop("'tol' should be NA (default) or a numeric value.")
        }
        Mean_new <- Mean
        i <- 0
        while((pdDist(Mean_new, Mean) > tol) & (i < max_iter)){
          Mean <- Mean_new
          Mean_new <- Expm(Mean, apply(sapply(1:dim(M)[3], function(i) w[i] *
                                                Logm(Mean, M[, , i]), simplify = "array"), c(1, 2), sum))
          i <- i+1
        }
        Mean <- Mean_new
      }
    } else if(metric == "Euclidean"){

      ## Euclidean weighted mean
      Mean <- apply(array(rep(w, each = d^2), dim = c(d, d, dim(M)[3])) * M, c(1, 2), sum)
    }
  }
  return(Mean)
}

