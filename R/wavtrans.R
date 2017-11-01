#' Forward AI wavelet transform
#'
#' \code{WavTransf} computes the forward \emph{average-interpolation} (AI) wavelet transform of a
#' curve of length \eqn{m} of (\eqn{d \times d})-dimensional Hermitian PD matrices as described in
#' (Chau and von Sachs, 2017).
#'
#' @param P a (\eqn{d,d,m})-dimensional array of Hermitian PD matrices, with \eqn{m = 2^J} for some \eqn{J > 0}.
#' @param order an odd integer between 1 and 9 corresponding to the order of the MI refinement scheme.
#' @param jmax the maximum scale up to which the wavelet coefficients are computed. If \code{jmax} is not
#' specified it is set equal to the maximum possible scale \code{jmax = J-1}.
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
#' D <- WavTransf(pgram$P)
#'
#' @return The function returns a list with two components:
#' \item{D }{a list of arrays, where each (\eqn{d, d, 2^j})-dimensional array contains the (\eqn{d \times d})-
#' dimensional wavelet coefficients at the \eqn{2^j} different locations in the given wavelet scale. The first
#' list element is a (\eqn{d, d, 2})-dimensional array containing the (\eqn{d \times d})-dimensional midpoints
#' at the coarsest scale (\eqn{j = 1}) in the midpoint pyramid.}
#' \item{M }{a list of arrays, where each (\eqn{d, d, 2^j})-dimensional array contains the (\eqn{d \times d})-
#' dimensional midpoints at the \eqn{2^j} different locations in the given midpoint scale. The first list
#' element is equivalent to the first list element in the \code{$D} component.}
#'
#' @seealso \code{\link{InvWavTransf}}, \code{\link{pdSpecEst}}
#'
#' @references Chau, J. and von Sachs, R. (2017). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#' @references Brockwell, P.J. and Davis, R.A. (1991). \emph{Time series: Theory and Methods}. New York: Springer.
#'
#' @importFrom utils txtProgressBar
#'
#' @export
WavTransf <- function(P, order = 5, jmax, periodic = T, metric = "Riemannian", progress = F) {

  ## Set variables
  J <- log2(dim(P)[3])
  if (!isTRUE(all.equal(as.integer(J), J))) {
    stop(paste0("Input length is non-dyadic, please change length ", dim(P)[3],
                " to dyadic number."))
  }
  if (!isTRUE(order %% 2 == 1)) {
    warning("Refinement order should be an odd integer, by default set to 5")
    order <- 5
  }
  metric <- match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))
  if (isTRUE(order > 9) & !isTRUE(metric == "Riemannian")) {
    warning(paste0(metric ," metric: refinement order should be an odd integer <= 9, by default set to 5"))
    order <- 5
  }

  d <- dim(P)[1]
  L <- (order - 1) / 2
  L_round <- 2 * ceiling(L / 2)

  P <- (if(metric == "logEuclidean"){
    sapply(1:2^J, function(i) Logm(diag(d), P[, , i]), simplify = "array")
  } else if(metric == "Cholesky"){
    sapply(1:2^J, function(i) t(Conj(Chol(P[, , i]))), simplify = "array")
  } else if(metric == "rootEuclidean"){
    sapply(1:2^J, function(i) Sqrt(P[, , i]), simplify = "array")
  } else P)

  if(periodic & (order > 1)){
    P_per <- array(if(L %% 2 == 0) {
      c(rep(c(P, P[, , c(2^J, 2^J:2)]), times = L), P)
    } else {
      c(P[, , c(2^J, 2^J:2)], rep(c(P, P[, , c(2^J, 2^J:2)]), times = L))
    }, dim = c(d, d, (2 * L + 1) * 2^J))
  } else{
    P_per <- P
  }

  M <- list()
  for (j in J:0) {
    if (j == J) {
      Mper <- P_per
    } else {
      if(!(metric == "Riemannian")){
        Mper <- sapply(1:(dim(Mper)[3]/2), function(i) 0.5 * (Mper[, , 2 * i - 1] +  Mper[, , 2 * i]),
                       simplify = "array")
      } else {
        Mper <- sapply(1:(dim(Mper)[3]/2), function(i) Mid(Mper[, , 2 * i - 1], Mper[, , 2 * i]),
                       simplify = "array")
      }
    }
    M[[j + 1]] <- if(periodic & (j > 0)){
      Mper[, , L * (2^j) - L_round + 1:(2^j + 2 * L_round)]
    } else{
      Mper
    }
  }
  names(M) <- paste0("M.scale", 0:J)

  D <- list()
  tM <- list()
  if (missing(jmax)) {
    jmax <- J - 1
  }
  if (jmax > J - 1) {
    warning(paste0("jmax cannot exceed maximum scale j=", J - 1))
    jmax <- J - 1
  }

  ## Compute wavelet transform
  if(progress){
    pb <- utils::txtProgressBar(1, 100, style = 3)
  }
  for (j in 0:jmax) {
    tm1 <- pdSpecEst:::Impute1D(M[[j + 1]], L, ifelse(order <= 9, "weights", "neville"), inverse = F, metric = metric)
    tM[[j + 1]] <- (if(periodic) tm1[, , L_round / 2 + ifelse(j > 0 | L %% 2 == 0, 0, -1) + 1:(2^j + L_round), drop = F] else tm1)
    if(!(metric == "Riemannian")){
      D[[j + 1]] <- sapply(1:dim(tM[[j + 1]])[3], function(l) 2^(-j/2) * (M[[j + 2]][, , 2 * l] - tM[[j + 1]][, , l]),
                           simplify = "array")
    } else{
      iSqrt_tm1 <- sapply(1:dim(tM[[j + 1]])[3], function(l) pdSpecEst:::iSqrt(tM[[j + 1]][, , l]), simplify = "array")
      D[[j + 1]] <- sapply(1:dim(tM[[j + 1]])[3], function(l) 2^(-j/2) * Logm(diag(d), (iSqrt_tm1[, , l] %*%
                         M[[j + 2]][, , 2 * l]) %*% iSqrt_tm1[, , l]), simplify = "array")
    }
    names(D)[j + 1] <- paste0("D.scale", j + 1)
    if(progress){
      utils::setTxtProgressBar(pb, round(100 * (j + 1) / (jmax + 1)))
    }
  }
  if(progress){
    close(pb)
  }
  return(list(D = D, M = M[1:J], tM = tM))
}

#' 2D Forward AI wavelet transform
#'
#' @export
WavTransf2D <- function(P, order = c(3,3), jmax, cores = NULL) {

  # #' @importFrom foreach "%dopar%"

  ## Set variables
  J1 <- log2(dim(P)[3])
  J2 <- log2(dim(P)[4])
  J <- max(J1, J2)
  J0_2D <- abs(J1 - J2)
  if (!isTRUE(all.equal(as.integer(J1), J1) & all.equal(as.integer(J2), J2))) {
    stop(paste0("Input length is non-dyadic, please change length ", dim(P)[3], " or ",
                dim(P)[4], " to dyadic number."))
  }
  d <- dim(P)[1]
  M <- list()
  Mper_2D <- function(Mper, j1, j2) {
    l1 <- t(sapply(1:2^j1, function(i) c(1, 2) + 2 * (i - 1)))
    l2 <- t(sapply(1:2^j2, function(i) c(1, 2) + 2 * (i - 1)))
    grid <- expand.grid(1:2^j1, 1:2^j2)
    Mper_new <- array(c(mapply(function(i1, i2) KarchMean(array(c(Mper[,,l1[i1,],l2[i2,]]), dim = c(d,d,4))),
                               grid$Var1, grid$Var2, SIMPLIFY = "array")), dim = c(d, d, 2^j1, 2^j2))
    return(Mper_new)
  }

  grid_j <- cbind((J1:0)[1:(J + 1)], (J2:0)[1:(J + 1)])

  M <- list()
  for (j in J:0) {
    if(j == J){
      Mper <- unname(P)
    } else if(j >= J0_2D) {
      Mper <- Mper_2D(Mper, grid_j[J + 1 - j, 1], grid_j[J + 1 - j, 2])
    } else {
      if(is.na(grid_j[J + 1 - j, 1])) {
        j2 <- grid_j[J + 1 - j, 2]
        Mper <- array(c(sapply(1:2^j2, function(i) Mid(Mper[, , , 2 * i - 1],
                                                       Mper[, , , 2 * i]))), dim = c(d, d, 1, 2^j2))
      } else if(is.na(grid_j[J + 1 - j, 2])) {
        j1 <- grid_j[J + 1 - j, 1]
        Mper <- array(c(sapply(1:2^j1, function(i) Mid(Mper[, , 2 * i - 1, ],
                                                       Mper[, , 2 * i, ]))), dim = c(d, d, 2^j1, 1))
      }
    }
    M[[j + 1]] <- Mper
  }
  names(M) <- paste0("M.scale", 0:J)

  D <- list()
  tM <- list()
  if (missing(jmax)) {
    jmax <- J - 1
  }
  if (jmax > J - 1) {
    warning(paste0("jmax cannot exceed maximum scale j=", J - 1))
    jmax <- J - 1
  }

  ## Compute 2D wavelet transform
  # pb <- tcltk::tkProgressBar(max = 100)
  for (j in 0:jmax) {

    if(is.na(grid_j[J + 1 - j, 1])) {
      tm1 <- array(Impute1D(array(M[[j + 1]][, , 1, ], dim = c(d, d, 2^j)),
                            (order[2] - 1) / 2), dim = c(d, d, 1, 2^(j + 1)))
      iSqrt_tm1 <- sapply(1:2^(j + 1), function(i) iSqrt(tm1[, , 1, i]), simplify = "array")
      D[[j + 1]] <- array(c(sapply(1:2^(j + 1), function(i) ifelse(any(dim(M[[min(length(M), j + 2)]]) == 1),
                                                                   1/(J0_2D - j) * 2^(-j/2), 2^(-j)) * Logm(diag(d), (iSqrt_tm1[, , i] %*% M[[j + 2]][, , 1, i]) %*%
                                                                                                              iSqrt_tm1[, , i]))), dim = c(d, d, 1, 2^(j + 1)))
    } else if(is.na(grid_j[J + 1 - j, 2])) {
      tm1 <- array(Impute1D(array(M[[j + 1]][, , , 1], dim = c(d, d, 2^j)),
                            (order[1] - 1) / 2), dim = c(d, d, 2^(j + 1), 1))
      iSqrt_tm1 <- sapply(1:2^(j + 1), function(i) iSqrt(tm1[, , i, 1]), simplify = "array")
      D[[j + 1]] <- array(c(sapply(1:2^(j + 1), function(i) ifelse(any(dim(M[[min(length(M), j + 2)]]) == 1),
                                                                   1/(J0_2D - j) * 2^(-j/2), 2^(-j)) * Logm(diag(d), (iSqrt_tm1[, , i] %*% M[[j + 2]][, , i, 1]) %*%
                                                                                                              iSqrt_tm1[, , i]))), dim = c(d, d, 2^(j + 1), 1))
    } else{
      # if(is.null(cores)){
        tm1 <- Impute2D(M[[j + 1]], (order - 1) / 2)
      # } else{
      #   tm1 <- Impute2D_multicore(M[[j + 1]], (order - 1) / 2, cores)
      # }
      grid <- expand.grid(1:dim(tm1)[3], 1:dim(tm1)[4])
      iSqrt_tm1 <- array(c(mapply(function(i1, i2) iSqrt(tm1[, , i1, i2]), grid$Var1, grid$Var2, SIMPLIFY = "array")),
                         dim = c(d, d, dim(tm1)[3], dim(tm1)[4]))
      D[[j + 1]] <- array(c(mapply(function(i1, i2) 2^(-j) * Logm(diag(d), (iSqrt_tm1[, , i1, i2] %*%
                                                                              M[[j + 2]][, , i1, i2]) %*% iSqrt_tm1[, , i1, i2]),
                                   grid$Var1, grid$Var2, SIMPLIFY = "array")), dim = c(d, d, dim(tm1)[3], dim(tm1)[4]))
    }
    tM[[j + 1]] <- tm1
    names(D)[j + 1] <- paste0("D.scale", j + 1)

    # tcltk::setTkProgressBar(pb, value=round(j/jmax*100), label=paste0("Computed up to scale ",
                                                                      # j + 1, ", (", round(sum(4^(0:j))/sum(4^(0:jmax))*100),"% done)"))
  }
  # close(pb)

  return(list(M = M, D = D, tM = tM))
}
