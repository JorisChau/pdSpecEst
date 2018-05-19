#' Compute distance between two HPD matrices
#'
#' \code{pdDist} calculates a distance between two Hermitian PD matrices.
#'
#' Available distance measures between two Hermitian PD matrices are (i) affine-invariant Riemannian distance (default) as in
#' (Bhatia, 2009, Chapter 6), (ii) log-Euclidean distance, the Euclidean distance between matrix logarithms,
#' (iii) Cholesky distance, the Euclidean distance between Cholesky decompositions, (iv) Euclidean distance,
#' (v) root-Euclidean distance and (vi) Procrustes distance as in (Dryden et al., 2009). In particular, \code{pdDist} generalizes the function
#' \code{shapes::distcov} to compute the distance between two symmetric positive definite matrices to the
#' distance between two Hermitian positive definite matrices.
#'
#' @param A,B Hermitian positive definite matrices (of equal dimension).
#' @param metric the distance measure, one of \code{'Riemannian'},
#' \code{'logEuclidean'}, \code{'Cholesky'}, \code{'Euclidean'}, \code{'rootEuclidean'} or \code{'Procrustes'}.
#' Defaults to \code{'Riemannian'}.
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
pdDist <- function(A, B, metric = "Riemannian") {

  if (!(isTRUE(all.equal(dim(A), dim(B)) & (dim(A)[1] == dim(A)[2]) & (length(dim(A)) == 2)))) {
    stop("Incorrect input dimensions for arguments: 'A' and/or 'B',
         consult the function documentation for the requested inputs.")
  }
  metric <- match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "Euclidean", "rootEuclidean", "Procrustes"))

  if (metric != "Procrustes") {
    dd <- pdDist_C(A, B, metric)
  } else {
    l1 <- Sqrt(A)
    l2 <- Sqrt(B)
    dd <- sqrt(NormF(l1)^2 + NormF(l2)^2 - 2 * sum(svd(t(Conj(l2)) %*% l1)$d))
  }
  return(dd)
}

#' Weighted Karcher mean of HPD matrices
#'
#' \code{pdMean} calculates an (approximate) weighted Karcher or Frechet mean of \eqn{S} different
#' \eqn{(d \times d)}-dimensional Hermitian PD matrices intrinsic to a specific metric. For the
#' affine-invariant Riemannian metric, the weighted Karcher mean is either approximated via
#' the fast recursive algorithm in (Ho et al., 2016) or computed via the slower, but more accurate,
#' gradient descent algorithm in (Pennec, 2006). By default, the unweighted Karcher mean is computed.
#'
#' @param M a \eqn{(d,d,S)}-dimensional array of Hermitian PD matrices.
#' @param w an \eqn{S}-dimensional nonnegative weight vector, such that \code{sum(w) = 1}.
#' @param metric the distance measure, one of \code{'Riemannian'}, \code{'logEuclidean'},
#' \code{'Cholesky'}, \code{'Euclidean'} or \code{'rootEuclidean'}. Defaults to \code{'Riemannian'}.
#' @param grad_desc if \code{metric == "Riemannian"}, a logical value indicating if the
#' gradient descent algorithm should be used, defaults to \code{FALSE}.
#' @param maxit maximum number of iterations in gradient descent algorithm, only used if
#' \code{isTRUE(grad_desc & metric == "Riemannian")}. Defaults to \code{1000}
#' @param reltol optional tolerance parameter in gradient descent algorithm, only used if
#' \code{isTRUE(grad_desc & metric == "Riemannian")}. Defaults to \code{1E-10}.
#'
#' @examples
#' m <- function(){
#'  X <- matrix(complex(real=rnorm(9), imaginary=rnorm(9)), nrow=3)
#'  t(Conj(X)) %*% X
#' }
#' M <- replicate(100, m())
#' z <- rnorm(100)
#' w <- abs(z)/sum(abs(z))
#' pdMean(M, w)
#'
#' @references Chau, J. and von Sachs, R. (2017). \emph{Intrinsic wavelet regression for curves of
#' Hermitian positive definite matrices}. Available at \url{http://arxiv.org/abs/1701.03314}.
#' @references Pennec, X. (2006). Intrinsic statistics on Riemannian manifolds: Basic tools for geometric
#' measurements. \emph{Journal of Mathematical Imaging and Vision} 25(1), 127-154.
#'
#' @seealso \code{\link{Mid}}, \code{\link{pdMedian}}
#'
#' @export
pdMean <- function(M, w, metric = "Riemannian", grad_desc = F, maxit = 1000, reltol) {

  if (!(isTRUE(is.array(M) & (dim(M)[1] == dim(M)[2]) & (length(dim(M)) == 3)))) {
    stop("Incorrect input dimensions for arguments: 'M',
         consult the function documentation for the requested inputs.")
  }
  ## Set parameters
  metric = match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "Euclidean", "rootEuclidean"))
  w = (if(missing(w)) rep(1/dim(M)[3], dim(M)[3]) else w)
  d = dim(M)[1]

  if(dim(M)[3] == 1){
    Mean <- M[, , 1]
  } else{

    if(metric == "Riemannian"){
      ## Recursive algorithm
      Mean <- pdMean_C_approx(M, w)
      ## Gradient descent algorithm
      if(grad_desc){
        reltol <- (if(missing(reltol)) 1E-10 else reltol)
        if(!isTRUE(is.numeric(reltol) & is.numeric(maxit) & maxit > 0)){
          stop("Incorrect input for arguments: 'reltol' or 'maxit'.")
        }
        Mean <- pdMean_C(Mean, M, w, round(maxit), reltol)
      }
    } else {
      ## Transform
      M1 <- switch(metric,
                   logEuclidean = sapply(1:dim(M)[3], function(i) Logm(diag(d), M[, , i]), simplify = "array"),
                   Cholesky = sapply(1:dim(M)[3], function(i) Chol_C(M[, , i]), simplify = "array"),
                   rootEuclidean = sapply(1:dim(M)[3], function(i) Sqrt(M[, , i]), simplify = "array"),
                   Euclidean = M)
      ## Euclidean mean
      Mean0 <- apply(sweep(M1, 3, w, "*"), c(1, 2), sum)
      ## Transform back
      Mean <- switch(metric,
                     logEuclidean = Expm(diag(d), Mean0),
                     Cholesky = t(Conj(Mean0)) %*% Mean0,
                     rootEuclidean = t(Conj(Mean0)) %*% Mean0,
                     Euclidean = Mean0)
    }
  }
  return(Mean)
}

#' Weighted intrinsic median of HPD matrices
#'
#' \code{pdMedian} calculates a weighted intrinsic median of \eqn{S} different
#' \eqn{(d \times d)}-dimensional Hermitian PD matrices based on a Weiszfeld algorithm intrinsic
#' to the chosen metric. For the affine-invariant Riemannian metric, the intrinsic Weisfeld algorithm in
#' (Fletcher et al., 2009) is used. By default, the unweighted intrinsic median is computed.
#'
#' @param M a \eqn{(d,d,S)}-dimensional array of Hermitian PD matrices.
#' @param w an \eqn{S}-dimensional nonnegative weight vector, such that \code{sum(w) = 1}.
#' @param metric the distance measure, one of \code{'Riemannian'}, \code{'logEuclidean'},
#' \code{'Cholesky'}, \code{'Euclidean'} or \code{'rootEuclidean'}. Defaults to \code{'Riemannian'}.
#' @param maxit maximum number of iterations in gradient descent algorithm. Defaults to \code{1000}
#' @param reltol optional tolerance parameter in gradient descent algorithm. Defaults to \code{1E-10}.
#'
#' @examples
#' m <- function(){
#'  X <- matrix(complex(real=rnorm(9), imaginary=rnorm(9)), nrow=3)
#'  t(Conj(X)) %*% X
#' }
#' M <- replicate(100, m())
#' z <- rnorm(100)
#' w <- abs(z)/sum(abs(z))
#' pdMedian(M, w)
#'
#' @references Fletcher, P.T., Venkatasubramanian, S. and Joshi, S. (2009). The geometric median
#' on Riemannian manifolds with application to robust atlas estimation. \emph{NeuroImage} 45(1),
#' S143-S152.
#'
#' @seealso \code{\link{pdMean}}
#'
#' @export
pdMedian <- function(M, w, metric = "Riemannian", maxit = 1000, reltol) {

  if (!(isTRUE(is.array(M) & (dim(M)[1] == dim(M)[2]) & (length(dim(M)) == 3)))) {
    stop("Incorrect input dimensions for arguments: 'M',
         consult the function documentation for the requested inputs.")
  }
  ## Set parameters
  metric = match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "Euclidean", "rootEuclidean"))
  w = (if(missing(w)) rep(1/dim(M)[3], dim(M)[3]) else w)
  d = dim(M)[1]
  reltol = (if(missing(reltol)) 1E-10 else reltol)
  if(!isTRUE(is.numeric(reltol) & is.numeric(maxit) & maxit > 0)){
    stop("Incorrect input for arguments: 'reltol' or 'maxit'.")
  }

  if(dim(M)[3] == 1){
    Med <- M[, , 1]
  } else{

    if(metric == "Riemannian"){
      ## Initial estimate
      Med0 <- pdMean_C_approx(M, w)
      ## Weiszfeld algorithm
      Med <- pdMedian_C(Med0, M, w, round(maxit), reltol)
    } else {
      ## Transform
      M1 <- switch(metric,
                   logEuclidean = sapply(1:dim(M)[3], function(i) Logm(diag(d), M[, , i]), simplify = "array"),
                   Cholesky = sapply(1:dim(M)[3], function(i) Chol_C(M[, , i]), simplify = "array"),
                   rootEuclidean = sapply(1:dim(M)[3], function(i) Sqrt(M[, , i]), simplify = "array"),
                   Euclidean = M)
      ## Initial estimate
      Med0 <- apply(sweep(M1, 3, w, "*"), c(1, 2), sum)
      ## Euclidean Weiszfeld algorithm
      Med1 <- Euclid_Median_C(Med0, M1, w, round(maxit), reltol)
      ## Transform back
      Med <- switch(metric,
                    logEuclidean = Expm(diag(d), Med1),
                    Cholesky = t(Conj(Med1)) %*% Med1,
                    rootEuclidean = t(Conj(Med1)) %*% Med1,
                    Euclidean = Med1)
    }
  }
  return(Med)
}


