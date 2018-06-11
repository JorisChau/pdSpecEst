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
#' \code{rExamples1D()} generates stationary time series observations from several example HPD spectral matrices
#' for testing and simulation purposes. For more details, we refer to the simulation studies in (Chau and von Sachs, 2017).
#'
#' The examples include: (i) a \eqn{(3 \times 3)} heaviSine HPD spectral matrix consisting of smooth sinosoids with a break,
#' (ii) a \eqn{(3 \times 3)} bumps HPD spectral matrix containing peaks and bumps of various smoothness degrees, (iii) a
#' \eqn{(3 \times 3)} two-cats HPD spectral matrix visualizing the contour of two side-by-side cats, with inhomogeneous
#' smoothness across frequency, and (iv) a \eqn{(2 \times 2)} Gaussian HPD spectral matrix consisting of smooth random Gaussian
#' functions. The time series observations are generated via the Cram\'er representation based on the transfer function of
#' the example spectral matrix and complex normal random variates.
#'
#' @param n number of sampled matrices to be generated.
#' @param d row- (resp. column-)dimension of the generated matrices.
#' @param example the example HPD spectral matrix curve, one of \code{'heaviSine'}, \code{'bumps'}, \code{'two-cats'},
#'     \code{'gaussian'} or \code{'mix-gaussian'}.
#' @param return.ts if \code{isTRUE(return.ts)} the function returns time series observations generated via the Cram\'er representation
#' based on the transfer function of the example spectral matrix and complex normal random variates.
#' @param noise noise distribution for the generated noisy spectral matrices, one of \code{'riem-gaussian'},
#' \code{'periodogram'}, \code{'log-gaussian'}, \code{'wishart'} or \code{'log-wishart'}. Additional details are given below.
#' @param noise.level parameter to tune the signal-to-noise ratio for the generated noisy spectral matrices.
#' @param df.wishart optional parameter to specify the degrees of freedom in the case of a Wishart noise distribution or the
#' number of DPSS taper functions in the case of generated periodogram matrices. Only used if \code{noise} is equal to \code{'periodogram'},
#' \code{'wishart'} or \code{'log-wishart'}. By default \code{df.wishart} is equal to the dimension \code{d} to
#' guarantee positive definiteness of the generated matrices.
#'
#' @return Depending on the input arguments returns a list with two or three components:
#'    \item{\code{f} }{ example spectral matrix, \code{f} is a \code{(d, d, length(freq))}-dimensional array, corresponding
#'    to a \eqn{(d \times d)}-dimensional spectral matrix at the Fourier frequencies \code{freq}.}
#'    \item{\code{P} }{ multitaper HPD periodogram of the generated time series, by default pre-smoothed using \code{df.wishart = d}
#'    DPSS tapering functions, see \code{\link{pdPgram}} for details.}
#'    \item{\code{ts}}{ generated \eqn{d}-dimensional time series observations, only available if \code{isTRUE(return.ts)}.}
#'
#' @examples
#' example <- rExamples1D(100, example = "bumps", return.ts = TRUE)
#' ts.plot(Re(example$ts)) # plot generated time series
#'
#' @seealso \code{\link{pdPgram}}
#'
#' @references Chau, J. and von Sachs, R. (2017). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#'
#' @export
rExamples1D <- function(n, d = 3, example = c("heaviSine", "bumps", "two-cats", "gaussian", "mix-gaussian"),
                        return.ts = FALSE, noise = "riem-gaussian", noise.level = 1, df.wishart = NULL){

  ## Set variables
  example <- match.arg(example, c("heaviSine", "bumps", "two-cats", "gaussian", "mix-gaussian"))
  noise <- match.arg(noise, c("riem-gaussian", "periodogram", "log-gaussian", "wishart", "log-wishart"))
  d <- switch(example,
              "heaviSine" = 3,
              "bumps" = 3,
              "two-cats" = 3,
              "gaussian" = 2,
              "mix-gaussian" = d)
  df.wishart <- (if(is.null(df.wishart)) d else df.wishart)
  if(df.wishart < d & noise %in% c("wishart", "log-wishart")) {
    warning("'df.wishart' is smaller than the dimension 'd';
             the generated wishart matrices are not positive definite!")
    if(noise == "log-wishart") {
      stop("Matrix logarithm of generated wishart matrices fails due to zero eigenvalues;
            increase value of 'df.wishart' to resolve problem.")
    }
  }
  x <- seq(0, 1, length.out = n)

  if(example == "heaviSine"){
    ## Discontinuous heavisine matrix curve (d = 3)
    f <- sapply(x, function(t){
      Expm(diag(2, 3), H.coeff(sin(HeaviSine[ifelse(t < 0.5, 2, 4),] * pi * t +
                                   HeaviSine[ifelse(t < 0.5, 1, 3),]), inverse = T))
    }, simplify = "array")
  } else if(example == "two-cats"){
    ## Locally smooth two-cats matrix curve (d = 3)
    v1 <- c(2.0000000, 0.4961828, -0.8595843, -0.8290600, 2.5000000, -1.3037704, -1.4214866, 1.7449149, 3.0000000)
    f0 <- sapply(1:n, function(i){
      H.coeff(v1 * cat_fun[ceiling(1024 / n * i)] + c(0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05), inverse = T)
    }, simplify = "array")
    f <- Ptransf2D_C(f0, T, F, "rootEuclidean")
  } else if(example == "bumps"){
    ## Locally smooth bumps matrix curve (d = 3)
    f <- sapply(x, function(t){
      Expm(diag(3), (1 - 0.4 * t) * H.coeff(sqrt(t * (1 - t) + 1) *
                                        sin(pi / (0.4 * t + 0.1)) * (1 + Bumps), inverse = T))
    }, simplify = "array")
  } else if(example == "gaussian"){
    ## Gaussian matrix curve (d = 2)
    v0 <- c(7.21935, 13.05307, 12.70734, 14.73554)
    f <- sapply(x, function(t){
      Expm(diag(2, 2), H.coeff(v0 / sqrt(2 * pi) * exp(-((t - 0.5) / 0.3)^2 / 2), inverse = T))
    }, simplify = "array")
  } else if(example == "mix-gaussian"){
    ## Random mixture of gaussian matrix curves
    pars <- list(mu = runif(4, min = -3, max = 3),
                 sd = abs(rnorm(4)),
                 w = abs(rnorm(4)),
                 w.comp = replicate(4, rnorm(d^2, mean = 1)))
    pars$w <- pars$w / sum(pars$w)
    f0 <- apply(sapply(1:4, function(i){
            pars$w[i] * sapply(2 * pi * x - pi, function(t) {
            H.coeff(pars$w.comp[, i] * pars$sd[i] / sqrt(2 * pi) *
                      exp(-(t - pars$mu[i])^2/(2 * pars$sd[i]^2)), inverse = T)
            }, simplify = "array")
      }, simplify = "array"), c(1,2,3), sum)
    f <- sapply(1:n, function(i) Expm(diag(0.25, d), f0[, , i]), simplify = "array")
  }

  ## Standardize int. Frob. norm of f
  f <- f / mean(sapply(1:n, function(i) NormF(f[,,i])))
  f1 <- Ptransf2D_C(f, F, F, ifelse(noise %in% c("riem-gaussian", "wishart", "periodogram"),
                                    "rootEuclidean", "logEuclidean"))

  if(return.ts | noise == "periodogram") {

    ## Generate time series via Cramer representation
    chi <- matrix(nrow = d, ncol = 2 * n)
    chi[, 1:(n - 1)] <- replicate(n - 1, complex(d, rnorm(d, sd = sqrt(1/2)), rnorm(d, sd = sqrt(1/2))))
    chi[, c(n, 2 * n)] <- replicate(2, rnorm(d))
    chi[, (n + 1):(2 * n - 1)] <- Conj(chi[, 1:(n - 1)])
    f2 <- array(c(f1, Conj(f1[, , n:1])), dim = c(d, d, 2 * n))
    ts <- sqrt(2 * pi) / sqrt(2 * n) * mvfft(t(sapply(1:(2 * n), function(i) f2[, , i] %*% chi[, i])), inverse = T)

    ## Generate multitaper periodogram from time series
    if(noise == "periodogram") {
      P <- pdPgram(ts, B = df.wishart, method = "multitaper")$P
    }
  }

  if(noise %in% c("riem-gaussian", "log-gaussian")) {
    ## Generate iid zero mean (intrinsic) gaussian noise
    W0 <- replicate(n, H.coeff(rnorm(d^2, sd = sqrt(noise.level)), inverse = T))
    if(noise == "riem-gaussian") {
      W <- Ptransf2D_C(W0, T, F, "logEuclidean")
      P <- sapply(1:n, function(i) (t(Conj(f1[, , i])) %*% W[, , i]) %*% f1[, , i], simplify = "array")
    } else if(noise == "log-gaussian") {
      P <- Ptransf2D_C(f1 + W0, T, F, "logEuclidean")
    }
  } else if(noise %in% c("wishart", "log-wishart")) {
    ## Generate iid zero mean (intrinsic) wishart noise
    X0 <- replicate(n * df.wishart, complex(d, rnorm(d, sd = sqrt(1/2)), rnorm(d, sd = sqrt(1/2))))
    W0 <- array(apply(X0, 2, function(x) x %*% t(Conj(x))), dim = c(d, d, n * df.wishart))
    W <- sapply(1:n, function(i) apply(W0[, , (i - 1) * df.wishart + 1:df.wishart], c(1, 2), mean), simplify = "array")
    if(noise.level != 1 & df.wishart >= d) {
      W <- Ptransf2D_C(sqrt(noise.level) * Ptransf2D_C(W, F, F, "logEuclidean"), T, F, "logEuclidean")
    }
    if(noise == "wishart") {
      P <- sapply(1:n, function(i) (t(Conj(f1[, , i])) %*% W[, , i]) %*% f1[, , i], simplify = "array")
    } else {
      P <- Ptransf2D_C(f1 + Ptransf2D_C(W, F, F, "logEuclidean"), T, F, "logEuclidean")
    }
  }

  ## return list with objects
  if(return.ts) {
    res <-  list(f = f, P = P, ts = ts)
  } else {
    res <- list(f = f, P = P)
  }

  return(res)
}

#' Several example time-varying spectral matrices
#'
#' \code{rExamples2D()} generates several example HPD time-varying spectral matrices for testing and simulation purposes.
#'
#' The examples include: (i) a \eqn{(d \times d)} smiley HPD spectral matrix consisting of constant surfaces of random HPD matrices in
#' the shape of a smiley face, (ii) a \eqn{(d \times d)} tvar HPD spectral matrix generated from a time-varying vector auto-
#' regressive process of order 1 with a random time-varying coefficient matrix (\eqn{\Phi}), (iii) a \eqn{(d \times d)} generally smooth
#' HPD spectral matrix containing a pronounced peak in the center of the discretized time-frequency grid and (iv) a \eqn{(d \times d)}
#' facets HPD spectral matrix consisting of several facets generated from random geodesic surfaces. In addition to the HPD spectral matrices,
#' the function also generates pseudo HPD periodogram observations as independent complex random HPD Wishart matrices centered around
#' the generating HPD spectral matrix with \eqn{B} degrees of freedom. Informally, such random matrix behavior corresponds to the asymptotic
#' distribution of actual HPD periodogram observations (as obtained via \code{\link{pdPgram2D}}) of a multivariate time series with the given
#' generating HPD spectral matrix.
#'
#' @param n numeric vector with two components corresponding to the size of the discretized time-frequency grid at which the spectral matrix
#' is generated.
#' @param d row- (resp. column-)dimension of the generated spectral matrix.
#' @param example the example spectral matrix, one of \code{'smiley'}, \code{'tvar'}, \code{'peak'} or \code{'facets'}.
#' @param noise noise distribution for the generated noisy time-varying spectral matrices, one of \code{'riem-gaussian'},
#'    \code{'log-gaussian'}, \code{'wishart'} or \code{'log-wishart'}. Additional details are given below.
#' @param noise.level parameter to tune the signal-to-noise ratio for the generated noisy spectral matrices.
#' @param df.wishart complex Wishart distribution degrees of freedom, defaults to \code{df.wishart = d}, such that the generated pseudo periodogram matrix
#' observations are positive definite. Note that \code{df.wishart} coincides with the number of tapering functions in the (time-varying) multitaper
#' periodogram in \code{\link{pdPgram2D}}.
#'
#' @return Returns a list with two components:
#'    \item{\code{f} }{ example spectral matrix, \code{f} is a \code{(d, d, n[1], n[2])}-dimensional array, corresponding
#'    to a \eqn{(d \times d)}-dimensional spectral matrix at a \eqn{(n_1,n_2)}-dimensional grid of time-frequency points.}
#'    \item{\code{P} }{ pseudo HPD periodogram observations generated as random complex Wishart matrices centered around \code{f}
#'    with \code{df.wishart} degrees of freedom at \eqn{(n_1,n_2)}-dimensional grid of time-frequency points.}
#'
#' @examples
#' example <- rExamples2D(n = c(32, 32), example = "smiley")
#'
#' @seealso \code{\link{pdPgram2D}}
#'
#' @references Chau, J. and von Sachs, R. (2017). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#'
#' @export
rExamples2D <- function(n, d = 2, example = c("smiley", "tvar", "peak", "facets"),
                        noise = "riem-gaussian", noise.level = 1, df.wishart = NULL){
  ## Set variables
  example <- match.arg(example, c("smiley", "tvar", "peak", "facets"))
  noise <- match.arg(noise, c("riem-gaussian", "log-gaussian", "wishart", "log-wishart"))
  df.wishart <- (if(is.null(df.wishart)) d else df.wishart)
  if(df.wishart < d & noise %in% c("wishart", "log-wishart")) {
    warning("'df.wishart' is smaller than the dimension 'd';
             the generated wishart matrices are not positive definite!")
    if(noise == "log-wishart") {
      stop("Matrix logarithm of generated wishart matrices fails due to zero eigenvalues;
            increase value of 'df.wishart' to resolve problem.")
    }
  }
  grid.n <- expand.grid(1:n[1], 1:n[2])
  ast <- function(A, B) (t(Conj(A)) %*% B) %*% A
  ## wrapper for Ptransf2D_C
  Ptransf2D <- function(X, inv, metric){
      array(Ptransf2D_C(array(X, dim = c(dim(X)[1], dim(X)[2], dim(X)[3] * dim(X)[4])), inv, F, metric),
            dim = dim(X))
  }

  ## Create spectrum
  if(example == "smiley"){

    ## Create 2D smiley spectrum
    f <- array(dim = c(d, d, n))
    center <- n / 2
    eye1 <- c(1/3, 2/3) * n
    eye2 <- c(2/3, 2/3) * n
    p_i <- replicate(4, Expm(diag(2, d), H.coeff(rnorm(d^2, sd = 2), inverse = T)))
    for(i1 in 1:n[1]){
      for(i2 in 1:n[2]){
        if(((i1 - center[1])^2 / (n[1] * 3/7)^2 + (i2 - center[2])^2 / (n[2] * 2/5)^2) <= 1){
          f[, , i1, i2] <- p_i[, , 1]
        } else {
          f[, , i1, i2] <- p_i[, , 3]
        }
        if(((i1 - center[1])^2 / (n[1] / sqrt(10))^2 + (i2 - center[2])^2 / (n[2] / sqrt(10))^2) <= 1 &
           ((i1 - center[1])^2/(n[1] / sqrt(25))^2 + (i2 - center[2])^2 / (n[2] / sqrt(25))^2) >= 1 &
           (i2 - center[2]) < -0.25 * abs(i1 - center[1])) {
          f[, , i1, i2] <- p_i[, , 2]
        }
        if(((i1 - eye1[1])^2/(n[1] / 10)^2 + (i2 - eye1[2])^2 / (n[2] / 10)^2) <= 1 |
           ((i1 - eye2[1])^2/(n[1] / 10)^2 + (i2 - eye2[2])^2 / (n[2] / 10)^2) <= 1) {
          f[, , i1, i2] <- p_i[, , 4]
        }
      }
    }

  } else if(example == "tvar"){

    ## Create 2D (smooth) arma spectrum
    Sigma <- Expm(diag(d), H.coeff(rnorm(d^2), inverse = T))
    x <- head(seq(0, 1, length.out = n[2] + 1), -1)
    v0 <- matrix(0.1, d, d) + diag(0.5, d)
    v1 <- rnorm(d^2, sd = 2)
    v2 <- 2 + rnorm(d^2)
    Phi.t <- sapply(x, function(t) matrix(c(v0 * sin(v2 * pi * t + v1)), d), simplify = "array")
    f.nu_t <- function(nu, t) {
      PhiB <- solve(diag(d) - Phi.t[, , t] * exp(-1i * pi * nu / n[1]))
      return(1/(2 * pi) * ((PhiB %*% Sigma) %*% t(Conj(PhiB))))
    }
    f <- array(c(mapply(function(nu, t) f.nu_t(nu, t), grid.n$Var1, grid.n$Var2)), dim = c(d, d, n))

  } else if(example == "peak"){

    ## Create 2D peak spectrum
    v1 <- abs(rnorm(d^2, 0.5, 0.5))
    grid.x <- expand.grid(seq(-1, 1, length.out = n[1]), seq(-1, 1, length.out = n[2]))
    mu <- Expm(diag(d), H.coeff(rnorm(d^2), inverse = T))
    f <- array(c(mapply(function(x, y) Expm(mu, H.coeff(exp(-(abs(x)^v1 + abs(y)^v1)),
                            inverse = T)), grid.x$Var1, grid.x$Var2)), dim = c(d, d, n))

  } else if(example == "facets"){

    ## Create 2D facets spectrum
    p <- replicate(2, Expm(diag(d), H.coeff(rnorm(d^2), inverse = T)))
    facets <- sapply(1:2, function(i) pdNeville(array(replicate(4, Expm(i * diag(d), H.coeff(rnorm(d^2), inverse = T))),
                                      dim=c(d, d, 2, 2)), X = list(x = c(0, 1), y = c(0, 1)),
                                      x = list(x = seq(0, 1, length.out = n[1]), y = seq(0, 1, length.out = n[2]))), simplify = "array")
    select <- function(i1, i2){
      if((n[2] / n[1] * i1) < i2 & i2 <= (n[2] - n[2] / n[1] * i1)){
        res <- facets[, , i1, i2, 1]
      } else if((n[2] / n[1] * i1) < i2 & i2 > (n[2] - n[2] / n[1] * i1)){
        res <- facets[, , i1, i2, 2]
      } else if((n[2] / n[1] * i1) >= i2 & i2 >= (n[2] - n[2] / n[1] * i1)){
        res <- ast(p[, , 1], facets[, , i1, i2, 1])
      } else if((n[2] / n[1] * i1) >= i2 & i2 < (n[2] - n[2] / n[1] * i1)){
        res <- ast(p[, , 2], facets[, , i1, i2, 2])
      }
      res
    }
    f <- array(c(mapply(select, grid.n$Var1, grid.n$Var2)), dim = c(d, d, n))

  }

  ## Standardize int. Frob. norm of f
  f <- f / mean(apply(f, c(3, 4), NormF))
  f1 <- Ptransf2D(f, F, ifelse(noise %in% c("riem-gaussian","wishart"), "rootEuclidean", "logEuclidean"))

  if(noise %in% c("riem-gaussian", "log-gaussian")) {
    ## Generate iid zero mean (intrinsic) gaussian noise
    W0 <- replicate(n[2], replicate(n[1], H.coeff(rnorm(d^2, sd = sqrt(noise.level)), inverse = T)))
    if(noise == "riem-gaussian") {
      W <- Ptransf2D(W0, T, "logEuclidean")
      P <- array(c(mapply(function(i, j) ast(f1[, , i, j], W[, , i, j]), grid.n$Var1, grid.n$Var2)), dim = dim(W))
    } else if(noise == "log-gaussian") {
      P <- Ptransf2D(f1 + W0, T, "logEuclidean")
    }
  } else if(noise %in% c("wishart", "log-wishart")) {
    ## Generate iid zero mean (intrinsic) wishart noise
    X0 <- array(c(replicate(n[1] * n[2] * df.wishart, complex(d, rnorm(d, sd = sqrt(1/2)),
                  rnorm(d, sd = sqrt(1/2))))), dim = c(d, n[1] * df.wishart, n[2]))
    W0 <- array(c(apply(X0, c(2,3), function(x) x %*% t(Conj(x)))), dim = c(d,d, n[1] * df.wishart, n[2]))
    W <- array(c(mapply(function(i1, i2) apply(W0[, , (i1 - 1) * df.wishart + 1:df.wishart, i2], c(1, 2), mean),
                        grid.n$Var1, grid.n$Var2)), dim = c(d, d, n))
    if(noise.level != 1 & df.wishart >= d) {
      W <- Ptransf2D(sqrt(noise.level) * Ptransf2D(W, F, "logEuclidean"), T, "logEuclidean")
    }
    if(noise == "wishart") {
      P <- array(c(mapply(function(i1, i2) ast(f1[, , i1, i2], W[, , i1, i2]), grid.n$Var1,
                          grid.n$Var2)), dim = dim(W))
    } else {
      P <- Ptransf2D(f1 + Ptransf2D(W, F, "logEuclidean"), T, "logEuclidean")
    }
  }
  return(list(f = f, P = P))
}

