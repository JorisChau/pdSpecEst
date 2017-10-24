#' Tapered HPD periodogram matrix
#'
#' \code{pdPgram} calculates a tapered Hermitian PD periodogram matrix based on averaging raw
#' Hermitian PSD periodogram matrices of tapered multivariate time series segments. The tapered Hermitian
#' PD periodogram matrix is rescaled by the \emph{manifold} bias-correction in (Chau and von Sachs, 2017).
#'
#' If \code{method == "multitaper"}, \code{pdPgram} calculates a \eqn{(d \times d)}-dimensional multitaper
#' periodogram matrix based on \eqn{B} Discrete Prolate Spheroidal (i.e. Slepian) orthogonal tapering functions
#' as in \code{\link[multitaper]{dpss}} applied to the \eqn{d}-dimensional time series \code{X}. If \code{method == "bartlett"}, \code{pdPgram} calculates
#' Bartlett's spectral estimator by averaging the periodogram matrices of \code{B} non-overlapping segments of the
#' \eqn{d}-dimensional time series \code{X}.Note that Bartlett's spectral estimator is a
#' specific (trivial) case of a multitaper spectral estimator with uniform orthogonal tapering windows.
#'
#' @param X  a multivariate time series, with the \code{d} columns corresponding to the components
#'  of the time series.
#' @param B  depending on \code{method}, either the number of orthogonal Slepian tapers, or the number of
#' non-overlapping segments to compute Bartlett's averaged periodogram. By default,
#' \code{B = d}, such that the averaged periodogram is  positive definite.
#' @param method the tapering method, either \code{"multitaper"} or \code{"bartlett"} explained in the Details
#' section below. Defaults to \code{"bartlett"}.
#' @param bias.corr should the manifold bias-correction be applied to the Hermitian PD periodogram matrix?
#' Defaults to \code{TRUE}.
#'
#' @return A list containing two components:
#'    \item{\code{freq} }{ vector of frequencies at which the periodogram is computed.}
#'    \item{\code{P} }{ a \code{(d, d, length(freq))}-dimensional array containing the
#'      (\eqn{d \times d})-dimensional tapered periodogram matrices at frequencies corresponding
#'      to \code{freq}.}
#'
#' @references Bartlett, M.S. (1950). \emph{Periodogram analysis and continuous spectra}.
#' Biometrika 37 (1-2): 1-16.
#' @references Chau, J. and von Sachs, R. (2017). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#' @references Brockwell, P.J. and Davis, R.A. (1991). \emph{Time series: Theory and Methods}. New York: Springer.
#'
#' @note The curve of HPD periodogram matrices obtained from \code{pdPgram(X)$P} can be
#' used as input in the functions \code{\link{WavTransf}} or \code{\link{pdSpecEst}}.
#'
#' @seealso \code{\link{rARMA}}, \code{\link[multitaper]{dpss}}
#'
#' @examples
#' ## ARMA(1,1) process: Example 11.4.1 in (Brockwell and Davis, 1991)
#'
#' Phi <- array(c(0.7, 0, 0, 0.6, rep(0, 4)), dim = c(2, 2, 2))
#' Theta <- array(c(0.5, -0.7, 0.6, 0.8, rep(0, 4)), dim = c(2, 2, 2))
#' Sigma <- matrix(c(1, 0.71, 0.71, 2), nrow = 2)
#' ts.sim <- rARMA(200, 2, Phi, Theta, Sigma)
#' ts.plot(ts.sim$X) # plot generated time series traces
#' pgram <- pdPgram(ts.sim$X)
#'
#' @importFrom multitaper dpss
#' @export
pdPgram <- function(X, B, method = c("multitaper", "bartlett"), bias.corr = T) {

  if(missing(method)){
    method <- "bartlett"
  }
  method <- match.arg(method, c("multitaper", "bartlett"))
  d <- ncol(X)
  n <- nrow(X)
  if (missing(B)) {
    B <- d
  }
  pgram <- function(Y){
    dft <- 1/sqrt(nrow(Y)) * mvfft(Y)
    return(sapply(1:floor(nrow(Y)/2), function(i) dft[i, ] %*% t(Conj(dft[i, ])), simplify = "array"))
  }

  if (method == "bartlett") {
    Per <- sapply(1:B, function(b) 1/(2 * pi) * pgram(X[floor(n/B) * (b - 1) + 1:floor(n/B), ]), simplify = "array")
  } else if (method == "multitaper") {
    h <- multitaper::dpss(n, B, 1, returnEigenvalues = F)$v * sqrt(n)
    Per <- sapply(1:B, function(k) 1/(2 * pi) * pgram(h[, k] * X), simplify = "array")
  }
  if(bias.corr){
    P <- B * exp(-1/d * sum(digamma(B - (d - 1:d)))) * apply(Per, c(1, 2, 3), mean)
  } else{
    P <- apply(Per, c(1, 2, 3), mean)
  }
  freq <- pi * (0:(dim(P)[3]-1))/dim(P)[3]

  return(list(freq = freq, P = P))
}

#' Tapered HPD time-varying periodogram matrix
#'
#' @export
pdPgram2D <- function(X, B, method = c("DPSS", "Hermite"), bias.corr = F, time, freq){

  d <- ncol(X)
  if(missing(method)) {
    method <- "DPSS"
  }
  method <- match.arg(method, c("DPSS", "Hermite"))
  if (missing(B)) {
    B <- d
  }
  n <- nrow(X)
  L <- length(time)
  x <- head(seq(0, 1, len = n + 1), -1)
  bias.corr <- ifelse(bias.corr, B * exp(-1/d * sum(digamma(B - (d - 1:d)))), 1)

  ## Periodogram matrix pointwise in time
  Per_t <- function(Y){
    dft <- 1/sqrt(nrow(Y)) * mvfft(Y)
    return(sapply(ceiling(freq * nrow(Y)), function(i) dft[i, ] %*% t(Conj(dft[i, ])), simplify = "array"))
  }

  if(method == "Hermite"){

    ## Hermite functions
    Hermite <- function(k, t){
      # if(!all.equal(k, as.integer(k)) | (k < 0)){
      #   stop(paste0("The Hermite order k = ", k, " is not a non-negative integer!"))
      # }
      Herm_poly <- function(ti){
        H <- numeric(k + 1)
        for(i in 0:k){
          if(identical(i, as.integer(0))){
            H[i + 1] <- 1
          } else if(identical(i, as.integer(1))){
            H[i + 1] <- 2 * ti
          } else{
            H[i + 1] <- 2 * ti * H[i] - 2 * (i - 1) * H[i - 1]
          }
        }
        H[k + 1]
      }
      sapply(t, function(ti) exp(-ti^2 / 2) * Herm_poly(ti) / sqrt(sqrt(pi) * 2^k * factorial(k)))
    }

    ## Time-varying periodogram matrix
    Per <- sapply(time, function(ti) bias.corr * apply(sapply(0:(B - 1), function(b) 1/(2 * pi) *
                          Per_t(Hermite(b, (n^0.2 * L * (x - ti))[(round(ti * n)  - floor(n / (2 * L) - 1)):
                          (round(ti * n) + floor(n / (2 * L)))]) * X[(round(ti * n)  - floor(n / (2 * L) - 1)):
                          (round(ti * n) + floor(n / (2 * L))),]), simplify = "array"),
                          c(1, 2, 3), mean), simplify = "array")
  } else if(method == "DPSS") {

    h <- multitaper::dpss(2 * floor(n / (2 * L)), B, 1, returnEigenvalues = F)$v * sqrt(2 * floor(n / (2 * L)))

    ## Sliding DPSS multitaper
    Per <- sapply(time, function(ti) bias.corr * apply(sapply(1:B, function(b) 1/(2 * pi) *
                  Per_t(h[, b] * X[(round(ti * n)  - floor(n / (2 * L) - 1)):(round(ti * n) + floor(n / (2 * L))),]),
                  simplify = "array"), c(1, 2, 3), mean), simplify = "array")
  }

  return(list(P = aperm(Per, c(1, 2, 4, 3)), time = time, freq = freq))

}
