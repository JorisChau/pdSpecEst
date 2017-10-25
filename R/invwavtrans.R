#' Inverse MI wavelet transform
#'
#' \code{InvWavTransf} computes the inverse \emph{midpoint-interpolation} (MI) wavelet
#' transform of an array of coarse-scale Hermitian PD midpoints combined with a pyramid of Hermitian
#' matrix-valued wavelet coefficients as described in (Chau and von Sachs, 2017).
#'
#' @param D a list of arrays containing coarse-scale midpoints and Hermitian matrix-valued wavelet
#'  coefficients in the same format as the \code{$D} component given as output by the function
#'  \code{\link{WavTransf}}.
#' @param order an odd integer between 1 and 9 corresponding to the order of the MI refinement scheme.
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
#' D <- WavTransf(pgram$P)$D
#' P <- InvWavTransf(D)
#' all.equal(pgram$P, P)
#'
#' @return Returns a (\eqn{d, d, m})-dimensional array corresponding to a curve of length \eqn{m} of
#' (\eqn{d \times d})-dimensional Hermitian PD matrices.
#'
#' @seealso \code{\link{WavTransf}}, \code{\link{pdSpecEst}}
#'
#' @references Chau, J. and von Sachs, R. (2017) \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#' @references Brockwell, P.J. and Davis, R.A. (1991). \emph{Time series: Theory and Methods}. New York: Springer.
#'
#' @export
InvWavTransf <- function(D, M0, order = 5, jmax, periodic = T, metric = "Riemannian", ...) {

  dots <- list(...)
  return_val <- (if(is.null(dots$return_val)) "manifold" else dots$return_val)
  progress <- (if(is.null(dots$progress)) F else dots$progress)

  if (!(order %% 2 == 1)) {
    warning("Refinement order should be an odd integer, by default set to 5")
    order <- 5
  }
  metric <- match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))
  if (isTRUE(order > 9) & !isTRUE(metric == "Riemannian")) {
    warning(paste0(metric ," metric: refinement order should be an odd integer <= 9, by default set to 5"))
    order <- 5
  }

  L <- (order - 1) / 2
  L_round <- 2 * ceiling(L / 2)
  d <- nrow(D[[1]][, , 1])
  if(missing(jmax)){
    J <- length(D)
  } else{
    J <- jmax
  }
  m1 <- M0

  if(progress){
    pb <- utils::txtProgressBar(1, 100, style = 3)
  }
  for (j in 0:(J - 1)) {
    tm1 <- Impute1D(m1, L, method = ifelse(order <= 9, "weights", "neville"), metric = metric)
    tm1 <- tm1[, , ifelse(j > 0, L_round, 2 * floor(L / 2)) + 1:(2^(j + 1) + 2 * L_round)]

    reconstr_even <- function(i) {
      if((j + 1) <= length(D)){
        if (any(c(D[[j + 1]][, , i]) != 0)) {
          if(metric == "Riemannian"){
            Sqrt_tm1 <- Sqrt(tm1[, , 2 * i])
            m1_i <- (Sqrt_tm1 %*% Expm(diag(d), 2^(j/2) * D[[j + 1]][, , i])) %*% Sqrt_tm1
          } else{
            m1_i <- 2^(j/2) * D[[j + 1]][, , i] + tm1[, , 2 * i]
          }
        } else {
          m1_i <- tm1[, , 2 * i]
        }
      } else {
        m1_i <- tm1[, , 2 * i]
      }
      return(m1_i)
    }

    grid_j <- (if((j + 1) <= length(D)) 1:dim(D[[j + 1]])[3] else 1:(dim(tm1)[3]/2))
    m2 <- array(dim = c(d, d, dim(tm1)[3]))
    m2[, , c(F, T)] <- sapply(grid_j, reconstr_even, simplify = "array")
    if(metric == "Riemannian"){
      m2[, , c(T, F)] <- sapply(grid_j, function(i) (m1[, , i + L_round / 2 +
          ifelse((j > 0) | (L %% 2  == 0), 0, -1)] %*% solve(m2[, , 2 * i])) %*%
          m1[, , i + L_round / 2 + ifelse((j > 0) | (L %% 2 == 0), 0, -1)],
          simplify = "array")
    } else{
      m2[, , c(T, F)] <- sapply(grid_j, function(i) 2 * m1[, , i + L_round / 2 +
          ifelse((j > 0) | (L %% 2  == 0), 0, -1)] - m2[, , 2 * i], simplify = "array")
    }
    m1 <- m2

    if(progress){
      utils::setTxtProgressBar(pb, round(100 * (j + 1) / J))
    }
  }
  if(progress){
    close(pb)
  }
  if(return_val == "manifold"){
    m1 <- (if(metric == "logEuclidean"){
      sapply(1:dim(m1)[3], function(i) Expm(diag(d), m1[, , i]), simplify = "array")
    } else if(metric == "Cholesky"){
      sapply(1:dim(m1)[3], function(i) Chol_inv(m1[, , i]), simplify = "array")
    } else if(metric == "rootEuclidean") {
      sapply(1:dim(m1)[3], function(i) t(Conj(m1[, , i])) %*% m1[, , i], simplify = "array")
    } else m1)
  }

  return((if(periodic) m1[, , L_round + 1:2^J] else m1))
}

#' 2D Inverse AI wavelet transform
#'
#' @export
InvWavTransf2D <- function(D, M0, order = c(3,3), jmax, cores = NULL) {

  ## Set variables
  d <- dim(D[[1]])[1]
  if(missing(jmax)){
    J <- length(D)
  } else{
    J <- jmax
  }
  m1 <- M0
  J0_2D <- sum(sapply(1:length(D), function(j) any(dim(D[[j]]) == 1)))

  # pb <- tcltk::tkProgressBar(max = 100)
  for (j in 0:(J - 1)) {

    if(dim(D[[j + 1]])[3] == 1){
      tm1 <- array(Impute1D(array(m1[, , 1, ], dim = c(d, d, 2^j)),
                            (order[2] - 1) / 2), dim = c(d, d, 1, 2^(j + 1)))
    } else if (dim(D[[j + 1]])[4] == 1){
      tm1 <- array(Impute1D(array(m1[, , , 1], dim = c(d, d, 2^j)),
                            (order[1] - 1) / 2), dim = c(d, d, 2^(j + 1), 1))
    } else{
      # if(is.null(cores)){
        tm1 <- Impute2D(m1, (order - 1) / 2)
      # } else{
        # tm1 <- Impute2D_multicore(m1, (order - 1) / 2, cores)
      # }
    }

    reconstr <- function(i1, i2) {
      if((j + 1) <= length(D)){
        if (any(c(D[[j + 1]][, , i1, i2]) != 0)) {
          Sqrt_tm1 <- Sqrt(tm1[, , i1, i2])
          m1_i <- (Sqrt_tm1 %*% Expm(diag(d), ifelse(any(dim(D[[j + 1]]) == 1),
                                                     (J0_2D - j) * 2^(j/2), 2^j) * D[[j + 1]][, , i1, i2])) %*% Sqrt_tm1
        } else {
          m1_i <- tm1[, , i1, i2]
        }
      } else {
        m1_i <- tm1[, , i1, i2]
      }
      return(m1_i)
    }

    if((j + 1) <= length(D)){
      grid_j <- expand.grid(1:dim(D[[j + 1]])[3], 1:dim(D[[j + 1]])[4])
    } else{
      grid_j <- expand.grid(1:(dim(D[[length(D)]])[3] * 2^((j + 1) - length(D))),
                            1:(dim(D[[length(D)]])[4] * 2^((j + 1) - length(D))))
    }
    m1 <- array(c(mapply(function(i1, i2) reconstr(i1, i2), grid_j$Var1, grid_j$Var2, SIMPLIFY = "array")),
                dim = c(d, d, attributes(grid_j)$out.attrs$dim[1], attributes(grid_j)$out.attrs$dim[2]))

    # tcltk::setTkProgressBar(pb, value = round(j/(J-1) * 100), label =
                              # paste0("Computed up to scale ", j + 1, ", (", round(sum(4^(0:j))/sum(4^(0:(J - 1))) * 100),"% done)"))
  }
  # close(pb)

  return(m1)

}

