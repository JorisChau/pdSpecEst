#' Simulate vARMA(2,2) process.
#'
#' \code{rARMA} generates \code{d}-dimensional time series observations from a vector ARMA(2,2)
#' (autoregressive moving average) process based on Gaussian white noise for testing and simulation
#' purposes.
#'
#' @param n  number of time series observations to be generated.
#' @param d  dimension of the multivariate time series.
#' @param Phi a (\eqn{d, d, 2})-dimensional array, with \code{Phi[, , 1]} and \code{Phi[, , 2]}
#'  the autoregressive parameter matrices.
#' @param Theta a (\eqn{d, d, 2})-dimensional array, with \code{Theta[, , 1]} and \code{Theta[, , 2]}
#'  the moving-average parameter matrices.
#' @param Sigma the covariance matrix of the Gaussian white noise component.
#' @param burn  a burn-in period when generating the time series observations, by default \code{burn = 100}.
#' @param freq  an optional vector of frequencies, if \code{!is.null(freq)} the function also returns the
#' underlying spectral matrix of the stationary generating process at the frequencies corresponding to \code{freq}.
#'
#' @return The function returns a list with two components:
#'    \item{\code{X} }{ generated time series observations, the \code{d} columns correspond to the components of
#'     the multivariate time series.}
#'    \item{\code{f} }{ if \code{!is.null(freq)}, \code{f} is a \code{(d, d, length(freq))}-dimensional array containing
#'     the underlying spectral matrix of the process at frequencies corresponding to \code{freq}. If
#'     \code{is.null(freq)}, \code{f} is set to \code{NULL}.}
#'
#' @references Brockwell, P.J. and Davis, R.A. (1991). \emph{Time series: Theory and Methods}. New York: Springer.
#'
#' @examples
#' ## ARMA(1,1) process: Example 11.4.1 in (Brockwell and Davis, 1991)
#'
#' freq <- seq(from = pi / 100, to = pi, length = 100)
#' Phi <- array(c(0.7, 0, 0, 0.6, rep(0, 4)), dim = c(2, 2, 2))
#' Theta <- array(c(0.5, -0.7, 0.6, 0.8, rep(0, 4)), dim = c(2, 2, 2))
#' Sigma <- matrix(c(1, 0.71, 0.71, 2), nrow = 2)
#' ts.sim <- rARMA(200, 2, Phi, Theta, Sigma, freq=freq)
#' ts.plot(ts.sim$X) # plot generated time series traces.
#'
#' @useDynLib pdSpecEst, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm
#'
#' @export
rARMA <- function(n, d, Phi, Theta, Sigma, burn = 100, freq = NULL) {

  ## Check arguments
  if (missing(Phi)) {
    warning("Phi is not specified. By default the AR-components are equal to the zero matrices")
    Phi <- array(0, c(d, d, 2))
  } else if (!isTRUE(all.equal(dim(Phi), c(d, d, 2)))) {
    warning("Phi is incorrectly specified. By default the AR-components are equal to the zero matrices")
    Phi <- array(0, c(d, d, 2))
  }
  if (missing(Theta)) {
    warning("Theta is not specified. By default the MA-components are equal to the zero matrices")
    Theta <- array(0, c(d, d, 2))
  } else if (!isTRUE(all.equal(dim(Theta), c(d, d, 2)))) {
    warning("Theta is incorrectly specified. By default the MA-components are equal to the zero matrices")
    Theta <- array(0, c(d, d, 2))
  }
  if (missing(Sigma)) {
    warning("Sigma is not specified. By default Sigma is equal to the diagonal matrix")
    Sigma <- diag(d)
  } else if(!isTRUE(all.equal(dim(Sigma), c(d, d)))) {
    warning("Sigma is incorrectly specified. By default Sigma is equal to the diagonal matrix")
    Sigma <- diag(d)
  }
  if (!isTRUE(all(eigen(Sigma, symmetric = T, only.values = T)$values >= 0))) {
    stop("Sigma is not a PSD matrix")
  }

  ## Generate time series
  Se <- Sqrt(Sigma)
  Z <- replicate(n + burn, Se %*% stats::rnorm(d), simplify = T)
  X <- t(ARMA(Phi, Theta, Z, n + burn)[, (burn + 1):(n + burn)])

  ## Compute spectrum
  f <- NULL
  if (!is.null(freq)) {
    f.nu <- function(nu) {
      PhiB <- diag(d) - Phi[, , 1] * exp(complex(imaginary = -nu)) - Phi[, , 2] *
        exp(complex(imaginary = -2 * nu))
      ThetaB <- diag(d) + Theta[, , 1] * exp(complex(imaginary = -nu)) + Theta[, , 2] *
        exp(complex(imaginary = -2 * nu))
      return((((solve(PhiB) %*% ThetaB) %*% Sigma) %*% t(Conj(ThetaB))) %*% t(solve(Conj(PhiB))))
    }
    f <- sapply(freq, function(freq) 1/(2 * pi) * f.nu(freq), simplify = "array")
  }
  return(list(X = X, f = f))
}

#' Test functions
#'
#' @export
rExamples <- function(n, example = c("heaviSine", "bumps", "two-cats", "sine", "gaussian-mix"), ...){

  ## Set variables
  example = match.arg(example, c("heaviSine", "bumps", "two-cats", "sine", "gaussian-mix"))
  if(n %% 2 != 0){
    warning("'n' is odd and cannot be divided by 2, instead 'n' is replaced by 'n + 1'")
    n = n + 1
  }
  m = n/2
  d = ifelse(example == "gaussian-mix", 2, 3)
  dots = list(...)
  B = (if(is.null(dots$B)) d else dots$B)
  method = (if(is.null(dots$method)) "multitaper" else dots$method)
  bias.corr = (if(is.null(dots$bias.corr)) F else dots$bias.corr)
  breaks = (if(is.null(dots$breaks)) 2 else dots$breaks)
  v0 = (if(is.null(dots$v0)) NULL else dots$v0)
  v1 = (if(is.null(dots$v1)) NULL else dots$v1)
  x = seq(0, 1, length.out = n / 2)

  if(example == "heaviSine"){
    P0 <- sapply(x, function(t) E_coeff_inv(sin(HeaviSine[2, ] * pi * t + HeaviSine[1, ])), simplify = "array")
    P1 <- sapply(x, function(t) E_coeff_inv(sin(HeaviSine[4, ] * pi * t + HeaviSine[3, ])), simplify = "array")

      if(breaks == 3){
        P2 <- sapply(x, function(t) E_coeff_inv(HeaviSine[2, ] * pi * t + HeaviSine[1, ]), simplify = "array")
        P <- sapply(1:m, function(i) (if(i <= m/2) 3 * Expm(diag(2, d), P0[, , i]) else if(i <= 3/4 * m){
          Expm(diag(2, 3), P1[, , i]) } else 2 * Expm(diag(2, d), P2[, , i])), simplify = "array")
      } else if(breaks == 2){
        P <- sapply(1:m, function(i) (if(i <= m/2) Expm(diag(2, d), P0[, , i]) else 3 * Expm(diag(2, d), P1[, , i])),
                    simplify = "array")
      }
    } else if(example == "sine") {
      if(!is.null(v0)){
        P0 <- sapply(x, function(t) E_coeff_inv(v0[1, ] * sin(v0[2, ] * 2 * pi * t + v0[3, ] * pi)), simplify = "array")
      } else {
        P0 <- sapply(x, function(t) E_coeff_inv(Sine1[1, ] * sin(Sine1[2, ] * 2 * pi * t + Sine1[3, ] * pi)), simplify = "array")
      }
      P <- sapply(1:m, function(i) 3 * Expm(diag(2, d), P0[, , i]), simplify = "array")
  } else if(example == "two-cats"){
    v1 <- (if(is.null(v1)) c(2.0000000, 0.4961828, -0.8595843, -0.8290600, 2.5000000, -1.3037704, -1.4214866,
                             1.7449149, 3.0000000) else v1)
    f <- cat_fun[round(seq(from = 1, to = 2^10, length = m))]
    P0 <- sapply(1:d^2, function(i) if(i %in% c(1,5,9)) v1[i] * f + 0.05 else v1[i] * f)
    P <- sapply(1:m, function(i) t(Conj(E_coeff_inv(P0[i,]))) %*% E_coeff_inv(P0[i,]), simplify = "array")
  } else if(example == "bumps"){
    P0 <- sapply(x, function(t)  (1 - 0.4 * t) * E_coeff_inv(sqrt(t * (1 - t) + 1) *
                                                               sin(pi/(0.4 * t + 0.1)) * (1 + Bumps)), simplify = "array")
    P <- sapply(1:m, function(i) Expm(diag(3), P0[,,i]), simplify = "array")
  } else if(example == "gaussian-mix"){
    v0 <- c(7.21935, 13.05307, 12.70734, 14.73554)
    # v1 <- c(8.645193, 13.572339, 10.166949, 13.223771)
    P0 <- sapply(x, function(t) E_coeff_inv(v0/sqrt(2*pi) * exp(-((t - 0.5)/0.3)^2/2))
                                            # +v1/sqrt(2*pi) * exp(-((t - 0.25)/0.15)^2/2))
                                            , simplify = "array")
    P <- sapply(1:m, function(i) Expm(diag(2, d), P0[,,i]), simplify = "array")
  }

  P.sqrt <- sapply(1:m, function(i) Sqrt(P[,,i]), simplify = "array")

  ## Generate time series via Cramer
  chi <- matrix(nrow = d, ncol = 2 * m)
  chi[, 1:(m - 1)] <- replicate(m - 1, complex(d, rnorm(d, sd = sqrt(1/2)), rnorm(d, sd = sqrt(1/2))))
  chi[, c(m, 2 * m)] <- replicate(2, rnorm(d))
  chi[, (m + 1):(2 * m - 1)] <- Conj(chi[, 1:(m - 1)])

  P.sqrt1 <- array(c(P.sqrt, Conj(P.sqrt[, , m:1])), dim = c(d, d, 2 * m))
  ts <- sqrt(2 * pi) / sqrt(2 * m) * mvfft(t(sapply(1:(2 * m), function(i) P.sqrt1[, , i] %*% chi[, i])), inverse = T)

  ## Compute pre-smoothed periodogram
  per <- pdPgram(ts, B = B, method = method, bias.corr = bias.corr)

  return(list(f = P, freq = per$freq, per = per$P, ts = ts))

}



