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
rExamples <- function(n, example = c("heaviSine", "bumps", "sine1", "sine2"), ...){

  example <- match.arg(example, c("heaviSine", "bumps", "sine1", "sine2"))
  if(n %% 2 != 0){
    warning("'n' is odd and cannot be divided by 2, instead 'n' is replaced by 'n + 1'")
    n <- n + 1
  }
  dots <- list(...)

  B <- dots$B
  if (is.null(B)){
    B <- 3
  }
  method <- dots$method
  if (is.null(method)) {
    method <- "multitaper"
  }
  bias.corr <- dots$bias.corr
  if (is.null(bias.corr)) {
    bias.corr <- F
  }

  x <- seq(0, 1, length.out = n / 2)

  if(example %in% c("heaviSine", "sine1", "sine2")){
    P0 <- sapply(x, function(t) E_coeff_inv(sin(HeaviSine[2, ] * pi * t + HeaviSine[1, ])), simplify = "array")
    P1 <- sapply(x, function(t) E_coeff_inv(sin(HeaviSine[4, ] * pi * t + HeaviSine[3, ])), simplify = "array")
    P <- (if(example == "heaviSine"){
      sapply(1:(n/2), function(i) (if(i <= n / 4) Expm(diag(2, 3), P0[, , i]) else 3 * Expm(diag(2, 3), P1[, , i])),
             simplify = "array")
    } else{
      sapply(1:(n/2), function(i) Expm(diag(2, 3), (if(example == "sine1") P0[, , i] else P1[, , i])),
             simplify = "array")
    })
  } else if(example == "bumps"){
    P0 <- sapply(x, function(t)  (1 - 0.4 * t) * E_coeff_inv(sqrt(t * (1 - t) + 1) *
                                                               sin(pi/(0.4 * t + 0.1)) * (1 + Bumps)), simplify = "array")
    P <- sapply(1:(n/2), function(i) Expm(diag(3), P0[,,i]), simplify = "array")
  }

  P.sqrt <- sapply(1:(n/2), function(i) Sqrt(P[,,i]), simplify = "array")

  ## Generate time series via Cramer
  chi <- matrix(nrow = 3, ncol = n)
  chi[, 1:(n/2 - 1)] <- replicate(n/2 - 1, complex(3, rnorm(3, sd = sqrt(1/2)), rnorm(3, sd = sqrt(1/2))))
  chi[, c(n/2, n)] <- replicate(2, rnorm(3))
  chi[, (n/2 + 1):(n - 1)] <- Conj(chi[, 1:(n/2 - 1)])

  P.sqrt1 <- array(c(P.sqrt, Conj(P.sqrt[, , (n/2):1])), dim = c(3, 3, n))
  ts <- sqrt(2 * pi) / sqrt(n) * mvfft(t(sapply(1:n, function(i) P.sqrt1[, , i] %*% chi[, i])), inverse = T)

  ## Compute pre-smoothed periodogram
  per <- pdPgram(ts, B = B, method = method, bias.corr = bias.corr)

  return(list(f = P, freq = per$freq, per = per$P, ts = ts))

}



