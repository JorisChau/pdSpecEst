#' Simulate vARMA(2,2) time series observations
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

#' Several example spectral matrices
#'
#' \code{rExamples()} generates stationary time series observations from several example HPD spectral matrices
#' for testing and simulation purposes. For more details, we refer to the simulation studies in (Chau and von Sachs, 2017).
#' The examples include: (i) a \eqn{(3 \times 3)} heaviSine HPD spectral matrix consisting of smooth sinosoids with a break,
#' (ii) a \eqn{(3 \times 3)} bumps HPD spectral matrix containing peaks and bumps of various smoothness degrees, (iii) a
#' \eqn{(3 \times 3)} two-cats HPD spectral matrix visualizing the contour of two side-by-side cats, with inhomogeneous
#' smoothness across frequency, and (iv) a \eqn{(3 \times 3)} Gaussian HPD spectral matrix consisting of smooth random Gaussian
#' functions. The time series observations are generated via the Cram\'er representation based on the transfer function of
#' the example spectral matrix and complex normal random variates.
#'
#' @param n number of time series observations to be generated.
#' @param example the example spectral matrix, one of \code{'heaviSine'}, \code{'bumps'}, \code{'two-cats'} or \code{'gaussian'}.
#' @param ... additional arguments passed on to \code{\link{pdPgram}}.
#'
#' @return Returns a list with four components:
#'    \item{\code{f} }{ example spectral matrix, \code{f} is a \code{(d, d, length(freq))}-dimensional array, corresponding
#'    to a \eqn{(d \times d)}-dimensional spectral matrix at the Fourier frequencies \code{freq}.}
#'    \item{\code{freq} }{ vector of Fourier frequencies from zero to \eqn{\pi}.}
#'    \item{\code{per} }{ multitaper HPD periodogram of the generated time series, by default pre-smoothed using \code{B = d}
#'    DPSS tapering functions, see \code{\link{pdPgram}} for details.}
#'    \item{\code{ts} }{ generated time series observations, the \code{d} columns correspond to the components of
#'     the multivariate time series.}
#'
#' @examples
#' example <- rExamples(100, example = "bumps")
#' ts.plot(example$ts) # plot generated time series.
#'
#' @seealso \code{\link{pdPgram}}
#'
#' @references Chau, J. and von Sachs, R. (2017). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#'
#' @export
rExamples <- function(n, example = c("heaviSine", "bumps", "two-cats", "gaussian"), ...){

  ## Set variables
  example = match.arg(example, c("heaviSine", "bumps", "two-cats", "gaussian"))
  if(n %% 2 != 0){
    warning("'n' is odd and cannot be divided by 2, instead 'n' is replaced by 'n + 1'")
    n = n + 1
  }
  m = n/2
  d = ifelse(example == "gaussian", 2, 3)
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
    P0 <- sapply(1:d^2, function(i) if(i %in% c(1, 5, 9)) v1[i] * f + 0.05 else v1[i] * f)
    P <- sapply(1:m, function(i) t(Conj(E_coeff_inv(P0[i,]))) %*% E_coeff_inv(P0[i,]), simplify = "array")
  } else if(example == "bumps"){
    P0 <- sapply(x, function(t)  (1 - 0.4 * t) * E_coeff_inv(sqrt(t * (1 - t) + 1) *
                                                               sin(pi/(0.4 * t + 0.1)) * (1 + Bumps)), simplify = "array")
    P <- sapply(1:m, function(i) Expm(diag(3), P0[,,i]), simplify = "array")
  } else if(example == "gaussian"){
    v0 <- c(7.21935, 13.05307, 12.70734, 14.73554)
    P0 <- sapply(x, function(t) E_coeff_inv(v0/sqrt(2*pi) * exp(-((t - 0.5)/0.3)^2/2))
                                            , simplify = "array")
    P <- sapply(1:m, function(i) Expm(diag(2, d), P0[, , i]), simplify = "array")
  }

  P.sqrt <- sapply(1:m, function(i) Sqrt(P[, , i]), simplify = "array")

  ## Generate time series via Cramer representation
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

#' 2D Examples
#'
#' @export
rExamples2D <- function(n, d, B, example = c("smiley", "arma")){

  d = (if(example == "smiley") 3 else d)
  B = (if(missing(B)) d else B)
  bias.corr = B * exp(-1/d * sum(digamma(B - (d - 1:d))))

  ## Create spectrum
  if(example == "smiley"){

    f <- array(dim = c(d, d, n[1], n[2]))
    center <- c(n[1]/2, n[2]/2)
    eye1 <- c(1/3 * n[1], 2/3 * n[2])
    eye2 <- c(2/3 * n[1], 2/3 * n[2])
    p_i <- replicate(4, Expm(2*diag(d), Re(H.coeff(rnorm(d^2), inverse = T))))

    for(i1 in 1:n[1]){
      for(i2 in 1:n[2]){
        if(((i1 - center[1])^2/(n[1] * 3/7)^2 + (i2 - center[2])^2/(n[2] * 2/5)^2) <= 1){
          f[,,i1,i2] <- p_i[,,1]
        } else {
          f[,,i1,i2] <- p_i[,,3]
        }
        if(((i1 - center[1])^2/(n[1] * 1/sqrt(10))^2 + (i2 - center[2])^2 /
            (n[2] * 1/sqrt(10))^2) <= 1 &
           ((i1 - center[1])^2/(n[1] * 1/sqrt(25))^2 + (i2 - center[2])^2 /
            (n[2] * 1/sqrt(25))^2) >= 1 &
           (i2 - center[2]) < -0.25 * abs(i1 - center[1])) {
          f[,,i1,i2] <- p_i[,,2]
        }
        if(((i1 - eye1[1])^2/(n[1] * 1/sqrt(100))^2 + (i2 - eye1[2])^2 /
            (n[2] * 1/sqrt(100))^2) <= 1 |
           ((i1 - eye2[1])^2/(n[1] * 1/sqrt(100))^2 + (i2 - eye2[2])^2 /
            (n[2] * 1/sqrt(100))^2) <= 1) {
          f[,,i1,i2] <- p_i[,,4]
        }
      }
    }
  }

    if(example == "arma"){

      Sigma <- Expm(diag(d), H.coeff(rnorm(d^2), inverse = T))
      x <- head(seq(0, 1, len = n[2] + 1), -1)
      v0 <- matrix(0.1, nrow=d, ncol=d) + diag(0.5, nrow = d)
      v1 <- 2*rnorm(d^2)
      v2 <- 2+rnorm(d^2)
      Phi.t <- sapply(x, function(t) matrix(c(v0 * sin(v2*pi*t + v1)), nrow = d), simplify = "array")

      f.nu_t <- function(nu, t) {
        PhiB <- solve(diag(d) - Phi.t[, , t] * exp(-1i * pi * nu / n[1]))
        return(1/(2 * pi) * ((PhiB %*% Sigma) %*% t(Conj(PhiB))))
      }

      grid <- expand.grid(1:n[1], 1:n[2])
      f <- array(c(mapply(function(nu, t) f.nu_t(nu, t), grid$Var1, grid$Var2)),
                 dim = c(d, d, n[1], n[2]))
    }

    ## Create iid Wishart periodogram
    grid.n <- expand.grid(1:n[1], 1:n[2])
    ast <- function(A, B) (A %*% B) %*% t(Conj(A))
    f.sqrt <- array(c(apply(f, c(3, 4), function(f) pdSpecEst:::Sqrt(f))), dim = dim(f))
    X0 <- array(c(replicate(n[1] * n[2] * d, complex(d, rnorm(d, sd = sqrt(1/2)), rnorm(d, sd = sqrt(1/2))))),
                  dim = c(d, n[1] * d, n[2]))
    W0 <- array(c(apply(X0, c(2,3), function(X) X %*% t(Conj(X)))), dim = c(d,d,n[1]*d,n[2]))
    W <- array(c(mapply(function(i1, i2) apply(W0[, , (i1 - 1) * d + 1:d, i2], c(1, 2), mean), grid.n$Var1,
                 grid.n$Var2)), dim = c(d, d, n[1], n[2]))
    X <- array(c(mapply(function(i1, i2) bias.corr * ast(f.sqrt[, , i1, i2], W[, , i1, i2]), grid.n$Var1,
                        grid.n$Var2)), dim = dim(W))

    return(list(f = f, per = X))
}

