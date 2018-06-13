#' Intrinsic 1D wavelet-based clustering of multivariate spectra.
#'
#' \code{pdSpecClust1D} performs clustering of multivariate spectral matrices via a two-step fuzzy
#' clustering algorithm in the intrinsic manifold wavelet domain of curves in the space of HPD matrices
#' equipped with a metric, e.g. the Riemannian metric, specified by the user.
#'
#' The input array \code{P} contains initial noisy HPD spectral estimates of the
#' (\eqn{d,d})-dimensional spectral matrices at \eqn{n} different frequencies for \eqn{S}
#' different subjects, where \eqn{n} is a dyadic number. The initial spectral estimates can
#' be e.g. the tapered HPD periodograms given as output by \code{\link{pdPgram}}. \cr
#' For each subject \eqn{s}, thresholded wavelet coefficients in the intrinsic manifold wavelet domain are
#' calculated by \code{\link{pdSpecEst1D}}.\cr
#' The maximum wavelet scale taken into account in the clustering procedure is determined by the arguments
#' \code{jmax} and \code{d.jmax}. The maximum scale is set to the minimum of \code{jmax} and the wavelet
#' scale \eqn{j} for which the proportion of nonzero thresholded wavelet coefficients (averaged
#' across subjects) is smaller than \code{d.jmax}.\cr
#' The \eqn{S} subjects are assigned to \eqn{K} different clusters in a probabilistic fashion according to a
#' two-step procedure:
#' \enumerate{
#' \item In the first step, an intrinsic fuzzy c-medoids algorithm, with fuzziness parameter \eqn{m} is applied to the
#' \eqn{S} coarsest midpoints at scale \code{j = 0} in the subject-specific midpoints pyramids. Note that the intrinsic
#' c-medoids algorithm crucially relies on the metric that the space of HPD matrices gets equipped with.
#' \item In the second step, a weighted fuzzy c-means algorithm based on the Euclidean
#' distance function, also with fuzziness parameter \eqn{m}, is applied to the nonzero thresholded wavelet
#' coefficients of the \eqn{S} different subjects. The tuning parameter \code{tau} controls the weight given
#' to the cluster assignments obtained in the first step of the clustering algorithm.
#' }
#' If \code{return.centers = T}, the function also returns the \code{K} HPD spectral curves corresponding to
#' the cluster centers based on the given metric by applying the intrinsic inverse 1D wavelet transform (
#' \code{\link{InvWavTransf1D}}) to the cluster centers in the wavelet domain.
#'
#' @param P a (\eqn{d,d,n,S})-dimensional array of \eqn{n}-dimensional sequences of HPD matrices for \code{S}
#' different subjects, with \eqn{n} a dyadic number.
#' @param K the number of clusters, should be a integer larger than 1.
#' @param jmax an upper bound on the maximum wavelet scale to be considered in the clustering procedure. If
#' \code{jmax} is not specified, it is set equal to the largest (i.e. finest) possible wavelet scale.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can also be one of: \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
#' \code{"Euclidean"}. The intrinsic distance measures in the clustering algorithm fundamentally rely on the chosen metric.
#' @param m the fuzziness parameter for both the fuzzy c-medoids and the weighted fuzzy c-means algorithm. \code{m}
#' should be larger or equal to \eqn{1}. If \eqn{m = 1} the cluster assignments are no longer fuzzy (i.e. the procedure
#' performs hard clustering).
#' @param d.jmax a proportion that is used to determine the maximum wavelet scale to be considered in the clustering
#' procedure. A larger value \code{d.jmax} leads to less wavelet coefficients being taken into account, and therefore
#' lower computational effort in the procedure. If \code{d.jmax} is not specified, by default \code{d.jmax = 0.1}.
#' @param eps an optional vector with two elements determining the stopping criterion. The fuzzy c-medoids algorithm
#' (i.e. first clustering step) terminates if the (integrated) intrinsic distance between cluster centers is smaller than
#' \code{eps[1]}. The weighted fuzzy c-means (i.e. second clustering step) terminates if the (integrated) distance between
#' cluster centers is smaller than \code{eps[2]}. If \code{eps} is not specified, by default \code{eps = c(1E-4, 1E-4)}.
#' @param tau an optional argument tuning the weight given to the cluster assignments obtained in the first step of
#' the clustering algorithm. If \code{tau} is not specified, by default \code{tau = 0.5}.
#' @param max.iter an optional argument tuning the maximum number of iterations in both the first and second step of the
#' clustering algorithm, defaults to \code{max.iter = 50}.
#' @param return.centers should the cluster centers transformed back the space of HPD matrices also be returned?
#' Defaults to \code{return.centers = F}.
#' @param ... additional arguments for internal use.
#'
#' @return The function returns a list with 6 components:
#' \describe{
#'   \item{cl.prob }{an (\eqn{S,K})-dimensional matrix, where the value at position (\eqn{s,k}) in the
#'   matrix corresponds to the probabilistic cluster membership assignment of subject \eqn{s} with respect
#'   to cluster \eqn{k}.}
#'   \item{cl.centers.D }{a list of \code{K} wavelet coefficient pyramids, where each pyramid of wavelet
#'   coefficients is associated to a cluster center.}
#'   \item{cl.centers.M0 }{a list \code{K} arrays of coarse-scale midpoints at scale \code{j = 0}, where each
#'   array is associated to a cluster center.}
#'   \item{cl.centers.f }{if \code{return.centers = T} returns a list of \code{K} \eqn{(d,d,n)}-dimensional arrays,
#'   where each array corresponds to a discretized curve of HPD matrices associated to a cluster center. If
#'   \code{return.centers = F}, \code{cl.centers.f} returns \code{NULL}.}
#'   \item{cl.jmax }{the maximum wavelet scale taken into account in the clustering procedure determined by
#'   the input arguments \code{jmax} and \code{d.jmax}.}
#'   \item{iters }{the number of iterations in respectively the first and second step of the clustering procedure.}
#' }
#'
#' @examples
#' ## ARMA(1,1) process: Example 11.4.1 in (Brockwell and Davis, 1991)
#'
#' Phi1 <- array(c(0.5, 0, 0, 0.6, rep(0, 4)), dim = c(2, 2, 2))
#' Phi2 <- array(c(0.7, 0, 0, 0.4, rep(0, 4)), dim = c(2, 2, 2))
#' Theta <- array(c(0.5, -0.7, 0.6, 0.8, rep(0, 4)), dim = c(2, 2, 2))
#' Sigma <- matrix(c(1, 0.71, 0.71, 2), nrow = 2)
#'
#' ## Generate periodogram data for 10 subjects
#' pgram <- function(Phi) pdPgram(rARMA(2^9, 2, Phi, Theta, Sigma)$X)$P
#' P <- array(c(replicate(5, pgram(Phi1)), replicate(5, pgram(Phi2))), dim=c(2,2,2^8,10))
#'
#' cl <- pdSpecClust1D(P, K = 2, metric = "logEuclidean")
#'
#' @seealso \code{\link{pdSpecEst1D}}, \code{\link{WavTransf1D}}, \code{\link{pdDist}}, \code{\link{pdPgram}}
#'
#' @references Chau, J. and von Sachs, R. (2017). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#' @references Brockwell, P.J. and Davis, R.A. (1991). \emph{Time series: Theory and Methods}. New York: Springer.
#'
#' @export
pdSpecClust1D <- function(P, K, jmax, metric = "Riemannian", m = 2, d.jmax = 0.1,
                          eps = c(1e-04, 1e-04), tau = 0.5, max_iter = 50, return.centers = F, ...) {

  ## Set variables
  d = dim(P)[1]
  S = dim(P)[4]
  metric = match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))
  jmax = (if(missing(jmax)) log2(dim(P)[3]) - 2 else jmax)

  ## Compute denoised wavelet coefficients
  D <- D.est <- list()
  M0 <- array(dim = c(d, d, S))
  D.nzero <- matrix(NA, S, jmax + 1)
  for(s in 1:S) {
    D.s <- pdSpecEst1D(P[, , , s], metric = metric, return = "D", periodic = F,
                       jmax = jmax, return.D = "D.white", ...)
    D[[s]] <- D.s$D.white
    M0[, , s] <- D.s$M0
    D.nzero[s, ] <- sapply(0:jmax, function(j) ifelse(j == 0, T, sum(D.s$tree.weights[[j]])/2^j))
    D.est[[s]] <- D.s$D
  }

  ## C++ c-medoids algorithm coarsest midpoints
  clust <- cMeans_C(M0, M0[, , 1:K], S, K, m, eps[1], max_iter,
                    ifelse(metric == "Riemannian", "Riemannian", "Euclidean"), matrix(1, S, K))

  ## Set up variables for weighted c-means in wavelet domain
  jmax <- min(jmax, sum(colMeans(D.nzero) > d.jmax) - 1)
  n_jmax <- 2^(jmax + 1) - 1
  DD <- array(aperm(array(sapply(D, function(D) unlist(D[1:(jmax + 1)])),
                          dim = c(d, d, n_jmax, S)), c(1, 2, 4, 3)), dim = c(d, d, n_jmax * S))
  cent <- array(dim = c(d, d, n_jmax * K))
  for(k in 1:K) {
    w <- clust[, k, 1]^m / sum(clust[, k, 1]^m)
    for(i in 1:n_jmax) {
      cent[, , (k - 1) * n_jmax + i] <- apply(sweep(DD[, , (i - 1) * S + 1:S], 3, w, "*"), c(1, 2), sum)
    }
  }
  dist_weights <- (1 - exp(-tau * clust[, , 2])) / (1 + exp(-tau * clust[, , 2]))

  ## C++ weighted c-means algorithm wavelet domain
  clust <- cMeans_C(DD, cent, S, K, m, eps[2], max_iter, "Euclidean", dist_weights)[, , 1]
  rownames(clust) <- paste0("Subject", 1:S)
  colnames(clust) <- paste0("Cluster", 1:K)

  ## Compute wavelet domain centroids
  cent.M0 <- cent.D <- list()
  for(k in 1:K) {
    w <- clust[, k]^m / sum(clust[, k]^m)
    cent.M0[[k]] <- array(pdMean(M0, w, ifelse(metric == "Riemannian", "Riemannian", "Euclidean"),
                                 grad_desc = T), dim = c(d, d, 1))
    cent.D[[k]] <- lapply(1:length(D.est[[1]]), function(j) {
      apply(sweep(sapply(D.est, "[[", j, simplify = "array"), 4, w, "*"), c(1, 2, 3), sum)
    })
  }
  names(cent.D) <- names(cent.M0) <- paste0("cluster.center.", 1:K)

  ## Inverse wavelet transforms and return output
  if(isTRUE(return.centers)){
    cent.f <- lapply(1:K, function(k) InvWavTransf1D(cent.D[[k]], cent.M0[[k]], periodic = F, metric = metric, ...))
    names(cent.f) <- paste0("cluster.center.", 1:K)
    res <- list(cl.prob = clust, cl.centers.D = cent.D, cl.centers.M0 = cent.M0, cl.centers.f = cent.f, cl.jmax = jmax)
  } else {
    res <- list(cl.prob = clust, cl.centers.D = cent.D, cl.centers.M0 = cent.M0, cl.jmax = jmax)
  }
  return(res)
}

#' Intrinsic 2D wavelet-based clustering of multivariate time-varying spectra.
#'
#' \code{pdSpecClust2D} performs clustering of multivariate time-varying spectral matrices via a two-step fuzzy
#' clustering algorithm in the intrinsic manifold wavelet domain of surface in the space of HPD matrices
#' equipped with a metric, e.g. the Riemannian metric, specified by the user. This function extends
#' \code{pdSpecClust2D} for clustering surfaces instead of curves of HPD matrices.
#'
#' The input array \code{P} corresponds to initial noisy HPD time-varying spectral estimates of the (\eqn{d, d})-
#' dimensional spectral matrices at \eqn{m_1 \times m_2} different time-frequency points for \eqn{S} different
#' subjects, with \eqn{m_1, m_2} dyadic numbers. The initial spectral estimates can be e.g. the tapered HPD
#' time-varying  periodograms given as output by \code{\link{pdPgram2D}}. \cr
#' For each subject \eqn{s}, thresholded wavelet coefficients in the intrinsic 2D manifold wavelet domain are
#' calculated by \code{\link{pdSpecEst2D}}.\cr
#' The maximum wavelet scale taken into account in the clustering procedure is determined by the arguments
#' \code{jmax} and \code{d.jmax}. The maximum scale is set to the minimum of \code{jmax} and the wavelet
#' scale \eqn{j} for which the proportion of nonzero thresholded wavelet coefficients (averaged
#' across subjects) is smaller than \code{d.jmax}.\cr
#' The \eqn{S} subjects are assigned to \eqn{K} different clusters in a probabilistic fashion according to a
#' two-step procedure:
#' \enumerate{
#' \item In the first step, an intrinsic fuzzy c-medoids algorithm, with fuzziness parameter \eqn{m} is applied to the
#' \eqn{S} coarsest midpoints at scale \code{j = 0} in the subject-specific 2D midpoints pyramids. Note that the intrinsic
#' c-medoids algorithm crucially relies on the metric that the space of HPD matrices gets equipped with.
#' \item In the second step, a weighted fuzzy c-means algorithm based on the Euclidean
#' distance function, also with fuzziness parameter \eqn{m}, is applied to the nonzero thresholded wavelet
#' coefficients of the \eqn{S} different subjects. The tuning parameter \code{tau} controls the weight given
#' to the cluster assignments obtained in the first step of the clustering algorithm.
#' }
#' If \code{return.centers = T}, the function also returns the \code{K} HPD time-varying spectral matrices corresponding
#' to the cluster centers based on the given metric by applying the intrinsic inverse 2D wavelet transform (
#' \code{\link{InvWavTransf2D}}) to the cluster centers in the wavelet domain.
#'
#' @param P a (\eqn{d,d,n_1,n_2,S})-dimensional array corresponding to discretized surfaces of HPD matrices for \code{S}
#' different subjects at \eqn{n_1 \times n_2} different time-frequency points (on a rectangular tensor grid), with \eqn{n_1}
#' and \eqn{n_2} dyadic numbers.
#' @param K the number of clusters, should be a integer larger than 1.
#' @param jmax an upper bound on the maximum wavelet scale to be considered in the clustering procedure. If
#' \code{jmax} is not specified, it is set equal to the largest (i.e. finest) possible wavelet scale.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can also be one of: \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
#' \code{"Euclidean"}. The intrinsic distance measures in the clustering algorithm fundamentally rely on the chosen metric.
#' @param m the fuzziness parameter for both the fuzzy c-medoids and the weighted fuzzy c-means algorithm. \code{m}
#' should be larger or equal to \eqn{1}. If \eqn{m = 1} the cluster assignments are no longer fuzzy (i.e. the procedure
#' performs hard clustering).
#' @param d.jmax a proportion that is used to determine the maximum wavelet scale to be considered in the clustering
#' procedure. A larger value \code{d.jmax} leads to less wavelet coefficients being taken into account, and therefore
#' lower computational effort in the procedure. If \code{d.jmax} is not specified, by default \code{d.jmax = 0.1}.
#' @param eps an optional vector with two elements determining the stopping criterion. The fuzzy c-medoids algorithm
#' (i.e. first clustering step) terminates if the (integrated) intrinsic distance between cluster centers is smaller than
#' \code{eps[1]}. The weighted fuzzy c-means (i.e. second clustering step) terminates if the (integrated) distance between
#' cluster centers is smaller than \code{eps[2]}. If \code{eps} is not specified, by default \code{eps = c(1E-4, 1E-4)}.
#' @param tau an optional argument tuning the weight given to the cluster assignments obtained in the first step of
#' the clustering algorithm. If \code{tau} is not specified, by default \code{tau = 0.5}.
#' @param max.iter an optional argument tuning the maximum number of iterations in both the first and second step of the
#' clustering algorithm, defaults to \code{max.iter = 50}.
#' @param return.centers should the cluster centers transformed back the space of HPD matrices also be returned?
#' Defaults to \code{return.centers = F}.
#' @param ... additional arguments for internal use.
#'
#' @return The function returns a list with 6 components:
#' \describe{
#'   \item{cl.prob }{an (\eqn{S,K})-dimensional matrix, where the value at position (\eqn{s,k}) in the
#'   matrix corresponds to the probabilistic cluster membership assignment of subject \eqn{s} with respect
#'   to cluster \eqn{k}.}
#'   \item{cl.centers.D }{a list of \code{K} wavelet coefficient pyramids, where each 2D pyramid of wavelet
#'   coefficients is associated to a cluster center.}
#'   \item{cl.centers.M0 }{an array \code{K} \eqn{(d,d)}-dimensional coarse-scale midpoints at scale \code{j = 0},
#'   where each midpoint is associated to a cluster center.}
#'   \item{cl.centers.f }{if \code{return.centers = T} returns a list of \code{K} \eqn{(d,d,n_1,n_2)}-dimensional arrays,
#'   where each array corresponds to a discretized surface of HPD matrices associated to a cluster center. If
#'   \code{return.centers = F}, \code{cl.centers.f} returns \code{NULL}.}
#'   \item{cl.jmax }{the maximum wavelet scale taken into account in the clustering procedure determined by
#'   the input arguments \code{jmax} and \code{d.jmax}.}
#'   \item{iters }{the number of iterations in respectively the first and second step of the clustering procedure.}
#' }
#'
#' @examples
#' ## Generate periodogram data for 4 subjects
#' pgram <- function(rescale) rescale * rExamples2D(c(2^5, 2^5), 2, example = "smiley")$P
#' P <- array(c(replicate(2, pgram(1)), replicate(2, pgram(2))), dim=c(2,2,2^5,2^5,4))
#'
#' cl <- pdSpecClust2D(P, K = 2, metric = "logEuclidean")
#'
#' @seealso \code{\link{pdSpecEst2D}}, \code{\link{WavTransf2D}}, \code{\link{pdDist}}, \code{\link{pdPgram2D}}
#'
#' @export
pdSpecClust2D <- function(P, K, jmax, metric = "Riemannian", m = 2, d.jmax = 0.1,
                          eps = c(1e-04, 1e-04), tau = 0.5, max_iter = 50, return.centers = F, ...) {

  ## Set variables
  d = dim(P)[1]
  S = dim(P)[5]
  jmax = (if(missing(jmax)) max(log2(dim(P)[3]), log2(dim(P)[4])) - 2 else jmax)
  metric = match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))

  ## Compute denoised wavelet coefficients
  D <- D.est <- list()
  M0 <- array(dim = c(d, d, S))
  D.nzero <- D.nzero1 <- matrix(NA, S, jmax + 1)
  for(s in 1:S) {
    D.s <- pdSpecEst2D(P[, , , , s], metric = metric, return = "D", jmax = jmax, return.D = "D.white", ...)
    D[[s]] <- D.s$D.white
    M0[, , s] <- D.s$M0[, , 1, 1]
    D.nzero[s, ] <- sapply(0:jmax, function(j) ifelse(j == 0, 1, mean(D.s$tree.weights[[j]])))
    D.est[[s]] <- D.s$D
  }

  ## C++ c-medoids algorithm coarsest midpoints
  clust <- cMeans_C(M0, M0[, , sample(1:S, K)], S, K, m, eps[1], max_iter,
                    ifelse(metric == "Riemannian", "Riemannian", "Euclidean"), matrix(1, S, K))

  ## Set up variables for weighted c-means in wavelet domain
  jmax <- min(jmax, sum(colMeans(D.nzero) > d.jmax) - 1)
  n_jmax <- sum(apply(sapply(D[[1]], dim)[c(3,4), 1:(jmax + 1)], 2, prod))
  DD <- array(aperm(array(sapply(D, function(D) unlist(D[1:(jmax + 1)])),
                          dim = c(d, d, n_jmax, S)), c(1, 2, 4, 3)), dim = c(d, d, n_jmax * S))
  cent <- array(dim = c(d, d, n_jmax * K))
  for(k in 1:K) {
    w <- clust[, k, 1]^m / sum(clust[, k, 1]^m)
    for(i in 1:n_jmax) {
      cent[, , (k - 1) * n_jmax + i] <- apply(sweep(DD[, , (i - 1) * S + 1:S], 3, w, "*"), c(1, 2), sum)
    }
  }
  dist_weights <- (1 - exp(-tau * clust[, , 2])) / (1 + exp(-tau * clust[, , 2]))

  ## C++ weighted c-means algorithm wavelet domain
  clust <- cMeans_C(DD, cent, S, K, m, eps[2], max_iter, "Euclidean", dist_weights)[, , 1]
  rownames(clust) <- paste0("Subject", 1:S)
  colnames(clust) <- paste0("Cluster", 1:K)

  ## Compute wavelet domain centroids
  cent.M0 <- cent.D <- list()
  for(k in 1:K) {
    w <- clust[, k]^m / sum(clust[, k]^m)
    cent.M0[[k]] <- array(pdMean(M0, w, ifelse(metric == "Riemannian", "Riemannian", "Euclidean"),
                                 grad_desc = T), dim = c(d, d, 1, 1))
    cent.D[[k]] <- lapply(1:length(D.est[[1]]), function(j) {
      apply(sweep(sapply(D.est, "[[", j, simplify = "array"), 5, w, "*"), 1:4, sum)
    })
  }
  names(cent.D) <- names(cent.M0) <- paste0("cluster.center.", 1:K)

  ## Inverse wavelet transforms and return output
  if(isTRUE(return.centers)){
    cent.f <- lapply(1:K, function(k) InvWavTransf2D(cent.D[[k]], cent.M0[[k]], metric = metric, ...))
    names(cent.f) <- paste0("cluster.center.", 1:K)
    res <- list(cl.prob = clust, cl.centers.D = cent.D, cl.centers.M0 = cent.M0, cl.centers.f = cent.f, cl.jmax = jmax)
  } else {
    res <- list(cl.prob = clust, cl.centers.D = cent.D, cl.centers.M0 = cent.M0, cl.jmax = jmax)
  }
  return(res)
}



