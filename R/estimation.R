#' Wavelet-thresholded multivariate spectral estimator
#'
#' \code{pdSpecEst} calculates a \eqn{(d \times d)}-dimensional Hermitian PD wavelet-denoised multivariate
#' spectral estimator by thresholding wavelet coefficients in the manifold wavelet domain. The
#' estimation procedure is described in detail (Chau and von Sachs, 2017a).
#'
#' The input array \code{P} corresponds to an initial noisy Hermitian PD spectral estimate of the
#' (\eqn{d \times d})-dimensional spectral matrix at \code{m} different frequencies, with \eqn{m = 2^J} for some
#' \eqn{J > 0}. This can be e.g. the curve of averaged periodograms given as output in \code{\link{pdPgram}}.\cr
#' \code{P} is transformed to the manifold wavelet domain by the function \code{\link{WavTransf}}, and the
#' wavelet coefficients are decomposed in terms of an orthonormal basis of the space of Hermitian matrices. \cr
#' The variances of the components of the wavelet coefficients are standardized across scales via weighted
#' log-linear regression, where the weights increase exponentially across scales at a rate depending on the
#' tuning parameter \code{alpha}. Without prior knowledge of the sparsity of the signal, a choice \code{alpha}
#' in \eqn{[0.5, 1)} is reasonable in most settings. If \code{alpha} is not specified, by default \code{alpha = 0.75} \cr
#' The components of the wavelet coefficients are thresholded based on the hard (keep-or-kill)
#' threshold \code{lam}. If \code{lam} is unspecified, the threshold is determined in a data-adaptive manner
#' by a twofold cross-validation procedure, which is described in detail in (Chau and von Sachs, 2017). \cr
#' If \code{return == 'f'} the thresholded wavelet coefficients are transformed back to the frequency domain by
#' the function \code{\link{InvWavTransf}} giving the wavelet-denoised Hermitian PD spectral estimate.
#'
#' @param P a (\eqn{d,d,m})-dimensional array of Hermitian PD matrices, with \eqn{m} a dyadic number.
#' @param lam an optional argument specifying the wavelet threshold, if \code{lam}
#'  is not specified the threshold is calculated by a twofold cross-validation procedure.
#' @param order an odd integer between 1 and 9 corresponding to the refinement order of the MI wavelet transform.
#' @param return an optional argument that specifies whether the denoised spectral estimator
#'  is returned or not.
#' @param alpha an optional argument tuning the weights used to normalize the variances
#' of the wavelet coefficients across scales. By default, \code{alpha} is set to \code{0.75}.
#'
#' @return The function returns a list with four components:
#' \item{f }{a (\eqn{d,d,m})-dimensional array corresponding to the wavelet-denoised Hermitian PD (\eqn{d \times d})-dimensional
#' spectral estimate at the \code{m} different frequencies. If \code{!(return == 'f')}, the inverse wavelet transform
#' of the thresholded wavelet coefficients is not computed and \code{f} is set equal to \code{NULL}.}
#' \item{D }{a list of arrays, each (\eqn{d, d, 2^j})-dimensional array contains the thresholded
#' (\eqn{d \times d})-dimensional wavelet coefficients at the \eqn{2^j} different locations in the given wavelet scale
#' \eqn{j}. The first list element contains the midpoints at the coarsest scale in the
#' midpoint pyramid \eqn{j=1}, see (Chau and von Sachs, 2017) for more details.}
#' \item{lam }{the hard threshold used to threshold the components of the wavelet coefficients.}
#' \item{components }{a list with two elements. The first element \code{thresholded} is a list of arrays; each (\eqn{d^2, 2^j})-dimensional
#' array contains the components of the thresholded (\eqn{d \times d})-dimensional wavelet coefficients in terms of an
#' orthonormal basis of the space of (\eqn{d \times d})-dimensional Hermitian matrices. The columns correspond to the
#' \eqn{d^2} basis components at each of the \eqn{2^j} different locations at wavelet scale \eqn{j}.  Analogous, for the second element
#' \code{non-thresholded}, but containing the components of the non-thresholded wavelet coefficients.}
#'
#' @examples
#' ## ARMA(1,1) process: Example 11.4.1 in (Brockwell and Davis, 1991)
#'
#' Phi <- array(c(0.7, 0, 0, 0.6, rep(0, 4)), dim = c(2, 2, 2))
#' Theta <- array(c(0.5, -0.7, 0.6, 0.8, rep(0, 4)), dim = c(2, 2, 2))
#' Sigma <- matrix(c(1, 0.71, 0.71, 2), nrow = 2)
#' ts.sim <- rARMA(2^10, 2, Phi, Theta, Sigma)
#' ts.plot(ts.sim$X) # plot generated time series traces.
#'
#' pgram <- pdPgram(ts.sim$X)
#' f <- pdSpecEst(pgram$P)
#'
#' @seealso \code{\link{pdPgram}}, \code{\link{WavTransf}}, \code{\link{InvWavTransf}}
#'
#' @references Chau, J. and von Sachs, R. (2017a). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#' @references Brockwell, P.J. and Davis, R.A. (1991). \emph{Time series: Theory and Methods}. New York: Springer.
#'
#' @importFrom stats mad
#' @export
pdSpecEst <- function(P, lam = NULL, order = 5, return = "f", alpha = 0.75) {
  .Deprecated("pdSpecEst1D")

  ## Define variables
  J <- log2(dim(P)[3])
  if (!isTRUE(all.equal(as.integer(J), J))) {
    warning(paste0("Input length is non-dyadic, please change length ", dim(P)[3],
                   " to dyadic number."))
  }
  stopifnot(isTRUE(all.equal(as.integer(J), J)))
  if (!(order %in% c(1, 3, 5, 7, 9))) {
    warning("Refinement order should be an odd integer between 1 and 9, by default set to 5")
    order <- 5
  }
  dim <- dim(P)[1]

  ## Find optimal threshold
  if (is.null(lam)) {

    ## Wavelet transforms
    P.half <- list(odd = P[, , c(T, F)], even = P[, , c(F, T)])
    D.half <- list(odd = WavTransf(P.half$odd, order)$D, even = WavTransf(P.half$even, order)$D)
    d <- list(odd = list(), even = list())
    for (j in 1:(J - 2)) {
      d$odd[[j]] <- sapply(1:2^j, function(k) E_coeff(D.half$odd[[j + 1]][, , k]))
      d$even[[j]] <- sapply(1:2^j, function(k) E_coeff(D.half$even[[j + 1]][, , k]))
    }

    ## Normalize variances wavelet coefficients
    W <- diag(exp(alpha * (1:(J - 2))))
    X <- matrix(c(rep(1, J - 2), 1:(J - 2)), nrow = J - 2)
    mads <- lapply(1:2, function(l) sapply(1:(J - 2), function(j) apply(d[[l]][[j]], 1, stats::mad)))
    Y <- colMeans(rbind(log(mads[[1]]), log(mads[[2]])))
    beta <- ((solve((t(X) %*% W) %*% X) %*% t(X)) %*% W) %*% Y
    sd.j <- function(j) exp(sum(beta * c(1, j)))

    d.vec <- list(odd = NULL, even = NULL)
    d.new <- list(odd = list(), even = list())
    for (j in 1:(J - 2)) {
      d.new$odd[[j]] <- d$odd[[j]]/sd.j(j)
      d.vec$odd <- cbind(d.vec$odd, d.new$odd[[j]])
      d.new$even[[j]] <- d$even[[j]]/sd.j(j)
      d.vec$even <- cbind(d.vec$even, d.new$even[[j]])
    }

    ## Cross-validation scores
    cv <- function(lam) {
      D.lam <- D.half
      d.lam <- d
      for (j in 3:(J - 2)) {
        zero.odd <- sapply(1:2^j, function(k) (abs(d.new$odd[[j]][, k]) <
                                                 lam) | (d.lam$odd[[j - 1]][, ceiling(k/2)] == 0))
        zero.even <- sapply(1:2^j, function(k) (abs(d.new$even[[j]][, k]) <
                                                  lam) | (d.lam$even[[j - 1]][, ceiling(k/2)] == 0))
        d.lam$odd[[j]][zero.odd] <- 0
        d.lam$even[[j]][zero.even] <- 0
      }
      for (j in 1:(J - 2)) {
        D.lam$odd[[j + 1]] <- sapply(1:2^j, function(k) E_coeff_inv(d.lam$odd[[j]][, k]),
                                     simplify = "array")
        D.lam$even[[j + 1]] <- sapply(1:2^j, function(k) E_coeff_inv(d.lam$even[[j]][, k]),
                                      simplify = "array")
      }
      f.hat <- list(odd = InvWavTransf(D.lam$odd, order), even = InvWavTransf(D.lam$even, order))

      ## Predicted points
      f.t.even <- sapply(1:(2^(J - 1) - 1), function(k) Mid(f.hat$odd[, , k], f.hat$odd[, , k + 1]),
                         simplify = "array")
      f.t.even <- array(c(f.t.even, f.hat$odd[, , 2^(J - 1)]), dim = c(dim, dim, 2^(J - 1)))
      f.t.odd <- sapply(1:(2^(J - 1) - 1), function(k) Mid(f.hat$even[, , k], f.hat$even[, , k + 1]),
                        simplify = "array")
      f.t.odd <- array(c(f.hat$even[, , 1], f.t.odd), dim = c(dim, dim, 2^(J - 1)))

      return(sum(sapply(1:2^(J - 1), function(k) RiemmDist(f.t.even[, , k], P.half$even[, , k])^2 +
                          RiemmDist(f.t.odd[, , k], P.half$odd[, , k])^2)))
    }

    ## Golden section search
    lam.range <- sort(abs(c(d.vec$even)), decreasing = T)[c(10, round(0.25 * length(c(d.vec$even))))]
    lam.cv <- gss(lam.range, cv)
  }

  ## Rescale threshold to twice number of observations
  lam.cv <- ifelse(is.null(lam), 1/sqrt((1 - log(2)/log(2^J * dim^2))) * lam.cv, lam)

  ## Transform original noisy data
  D <- WavTransf(P, order)$D
  d <- list()
  for (j in 1:(J - 1)) {
    d[[j]] <- sapply(1:2^j, function(k) E_coeff(D[[j + 1]][, , k]))
  }

  ## Normalize variances wavelet coefficients
  W <- diag(exp(alpha * (1:(J - 1))))
  X <- matrix(c(rep(1, J - 1), 1:(J - 1)), nrow = J - 1)
  Y <- colMeans(log(sapply(1:(J - 1), function(j) apply(d[[j]], 1, stats::mad))))
  beta <- ((solve((t(X) %*% W) %*% X) %*% t(X)) %*% W) %*% Y
  sd.j <- function(j) exp(sum(beta * c(1, j)))

  d.vec <- NULL
  d.new <- list()
  for (j in 1:(J - 1)) {
    d.new[[j]] <- d[[j]]/sd.j(j)
    d.vec <- cbind(d.vec, d.new[[j]])
  }
  d1 <- d.new

  ## Threshold coefficients
  for (j in 3:(J - 1)) {
    zero <- sapply(1:2^j, function(k) (abs(d.new[[j]][, k]) < lam.cv) |
                     (d[[j - 1]][, ceiling(k/2)] == 0))
    d.new[[j]][zero] <- 0
    d[[j]][zero] <- 0
  }

  ## Inverse transform denoised data
  for (j in 1:(J - 1)) {
    D[[j + 1]] <- sapply(1:2^j, function(k) E_coeff_inv(d[[j]][, k]), simplify = "array")
  }

  if (return == "f") {
    f <- InvWavTransf(D, order)
  } else {
    f <- NULL
  }
  return(list(f = f, D = D, lam = lam.cv, components = list(thresholded = d.new, not_thresholded = d1)))
}

#' Tree-structured wavelet-thresholded 2D spectral estimator
#'
#' @export
pdCART <- function(D, tree = T, lam = NA, alpha = 1, periodic = T) {

  J <- length(D)
  d <- dim(D[[1]])[1]
  N <- dim(D[[1]])[3]
  L <- (N - 1) / 2
  L_b <- (if(periodic) ceiling(L/2) else 0)
  is_2D <- ifelse(length(dim(D[[1]])) == 4, T, F)
  D_trace_full <- if(is_2D){
    lapply(2:J, function(j) apply(D[[j]], c(3, 4), function(A) Re(sum(diag(A)))))
  } else {
    lapply(2:J, function(j) apply(D[[j]], 3, function(A) Re(sum(diag(A)))))
  }
  D_trace <- lapply(1:(J-1), function(j) D_trace_full[[j]][L_b + 1:2^j])
  s_e <- mad(c(D_trace[[J-1]]))
  if(is.na(lam)){
    lam <- alpha * s_e * sqrt(2 * log(length(unlist(D_trace))))
    # lam <- alpha * s_e
  }

  if(tree){

    ## Dyadic CART
    w <- D_trace
    # w[[1]] <- (if(is_2D) array(T, dim = dim(D_trace[[1]])) else rep(T, length(D_trace[[1]])))
    for(j in (J-1):1){
      if(j == (J-1)){
        w[[j]] <- ifelse(abs(D_trace[[j]]) > lam, T, F)
        R <- pmin(D_trace[[j]]^2, lam^2)
        V <- D_trace[[j]]^2
      } else{
        if(is_2D){
          dims <- dim(D_trace[[j]])
          if(dims[1] > 1){
            l1 <- t(sapply(1:dims[1], function(i) c(1, 2) + 2 * (i - 1)))
          } else {
            l1 <- matrix(1, ncol = 1, nrow = 1)
          }
          if(dims[2] > 1){
            l2 <- t(sapply(1:dims[2], function(i) c(1, 2) + 2 * (i - 1)))
          } else {
            l2 <- matrix(1, ncol = 1, nrow = 1)
          }
          grid <- expand.grid(1:dims[1], 1:dims[2])
          V <- array(c(mapply(function(i1, i2) sum(V[l1[i1, ], l2[i2, ]]), grid$Var1, grid$Var2)),
                     dim = dims) + D_trace[[j]]^2
          R <- array(c(mapply(function(i1, i2) sum(R[l1[i1, ], l2[i2, ]]), grid$Var1, grid$Var2)),
                     dim = dims) + lam^2
        } else {
          dims <- length(D_trace[[j]])
          V <- sapply(1:dims, function(i) V[2 * i - 1] + V[2 * i]) + D_trace[[j]]^2
          R <- sapply(1:dims, function(i) R[2 * i - 1] + R[2 * i]) + lam^2
        }
        w[[j]] <- ifelse(R < V, T, F)
        R <- pmin(V, R)
      }
    }
  } else{
    w <- lapply(1:(J - 1), function(j) abs(D_trace[[j]]) > lam)
  }

  D_w <- D
  for(j in 2:J){
    if(is_2D){ ## NEEDS TO BE UPDATED!
      dims <- dim(w[[j]])
      roots <- matrix(rep(matrix(rep(t(w[[j - 1]]), each = ifelse(dims[2] == 1, 1, 2)), byrow = T,
                                 ncol = ncol(w[[j]])), each = ifelse(dims[1] == 1, 1, 2)), nrow = nrow(w[[j]]))
      w[[j]] <- w[[j]] & roots
      D0 <- array(D_w[[j]], dim = c(d, d, dim(D_w[[j]])[3] * dim(D_w[[j]])[4]))
      D0[, , !(w[[j]])] <- 0
      D_w[[j]] <- array(D0, dim = c(d, d, dim(D_w[[j]])[3], dim(D_w[[j]])[4]))
    } else {
      w[[j - 1]] <- (if(tree) w[[j - 1]] & rep((if(j == 2) T else w[[j - 2]]), each = 2) else w[[j - 1]])
      if(periodic & (L_b > 0)){
        zeros <- !(c(abs(D_trace_full[[j - 1]][1:L_b]) > lam, w[[j - 1]], abs(D_trace_full[[j - 1]][2^(j - 1) + L_b + 1:L_b]) > lam))
      } else{
        zeros <- !(w[[j - 1]])
      }
      D_w[[j]][, , zeros] <- 0
    }
  }
  return(list(w = w, D_w = D_w))
}

#' Automatic tree-structured wavelet regression for curves of HPD matrices
#'
#' @export
pdSpecEst1D <- function(P, order = 5, policy = c("universal", "cv"), metric = "Riemannian", periodic = T,
                        alpha = 1, return = "f", ...) {

  ## Set variables
  dots = list(...)
  tol = (if(is.null(dots$tol)) 0.01 else dots$tol)
  alpha.range = (if(is.null(dots$alpha.range)) c(0.5, 2) else dots$alpha.range)
  tree = (if(is.null(dots$tree)) T else dots$tree)

  policy = match.arg(policy, c("universal", "cv"))
  metric = match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))
  periodic = isTRUE(periodic)
  J = log2(dim(P)[3])
  d = dim(P)[1]
  B = (if(is.null(dots$B)) d else dots$B)
  progress = (if(is.null(dots$progress)) F else dots$progress)
  J.out = (if(is.null(dots$J.out)) J else dots$J.out)
  jmax = min((if(is.null(dots$jmax)) J - 3 else dots$jmax), J.out - 1)
  jmax.cv = min((if(is.null(dots$jmax.cv)) J - 3 else dots$jmax.cv), J.out - 1)

  # Manifold bias-correction
  P = (if(metric == "Riemannian" | metric == "logEuclidean") {
          B * exp(-1/d * sum(digamma(B - (d - 1:d)))) * P } else P)

  if(policy == "cv"){

    ## Two-fold cross-validation (Nason, 1996)
    Pmid <- (if(metric == "logEuclidean") {
      sapply(1:2^J, function(i) Logm(diag(d), P[, , i]), simplify = "array")
    } else if(metric == "Cholesky") {
      sapply(1:2^J, function(i) t(Conj(Chol(P[, , i]))), simplify = "array")
    } else if(metric == "rootEuclidean"){
      sapply(1:2^J, function(i) Sqrt(P[, , i]), simplify = "array")
    } else P)

    for (j in J:(jmax.cv + 1)) {
      Pmid <- sapply(1:(dim(Pmid)[3]/2), function(i) (if(metric == "Riemannian"){
                  Mid(Pmid[, , 2 * i - 1], Pmid[, , 2 * i]) } else { 0.5 * (Pmid[, , 2 * i - 1] +
                  Pmid[, , 2 * i]) }), simplify = "array")
    }

    P1 <- list(odd = Pmid[, , c(T, F)], even = Pmid[, , c(F, T)])

    if(metric == "Riemannian"){
    for(m in c(1,2)){
      P1[[m]] <- sapply(1:2^(jmax.cv - 1), function(i) Logm(diag(d), P1[[m]][, , i]), simplify = "array")
      }
    }

    coeff.odd <- WavTransf(P1$odd, order, periodic = periodic, metric = "Euclidean", progress = F)
    coeff.even <- WavTransf(P1$even, order, periodic = periodic, metric = "Euclidean", progress = F)

    cv <- function(alpha){

      D <- list(odd = pdCART(coeff.odd$D, tree = F, alpha = alpha, periodic = periodic)$D_w,
                even = pdCART(coeff.even$D, tree = F, alpha = alpha, periodic = periodic)$D_w)

      f.hat <- list(odd = InvWavTransf(D$odd, coeff.odd$M[[1]], order,
                                       periodic = periodic, metric = "Euclidean", return_val = "tangent"),
                    even = InvWavTransf(D$even, coeff.even$M[[1]], order,
                                        periodic = periodic, metric = "Euclidean", return_val = "tangent"))
      ## Predicted points
      f.pred <- list(even = array(c(sapply(1:(2^(jmax.cv - 1) - 1), function(k) 0.5 * (f.hat$odd[, , k] + f.hat$odd[, , k + 1]),
                                  simplify = "array"), f.hat$odd[, , 2^(jmax.cv - 1)]), dim = c(d, d, 2^(jmax.cv - 1))),
                     odd = array(c(f.hat$even[, , 1], sapply(1:(2^(jmax.cv - 1) - 1), function(k) 0.5 * (f.hat$even[, , k] +
                                 f.hat$even[, , k + 1]), simplify = "array")), dim = c(d, d, 2^(jmax.cv - 1))))

      return(mean(sapply(1:2^(jmax.cv - 1), function(k) NormF(f.pred$even[, , k] - P1$even[, , k])^2 +
                                                      NormF(f.pred$odd[, , k] - P1$odd[, , k])^2)))
    }

    ## Find minimum and rescale twice number of data points
    alpha.opt <- gss(alpha.range, cv, tol)

  } else {
    alpha.opt <- alpha
  }

  ## Threshold full data using 'alpha.opt'
  coeff <- (if(policy == "cv"){
              WavTransf(P, order, jmax = jmax.cv - 1, periodic = periodic, metric = metric, progress = progress)
            } else  WavTransf(P, order, jmax = jmax - 1, periodic = periodic, metric = metric, progress = progress))
  coeff.opt <- pdCART(coeff$D, alpha = alpha.opt, tree = tree, periodic = periodic)
  f <- (if(return == "f"){
    InvWavTransf(coeff.opt$D_w, coeff$M[[1]], order, jmax = J.out, periodic = periodic, metric = metric, progress = progress)
  } else NULL)

  return(list(f = f, D = coeff.opt$D_w, M0 = coeff$M[[1]], tree.weights = coeff.opt$w, alpha.opt = alpha.opt))
}
