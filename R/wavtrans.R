#' MI forward wavelet transform.
#'
#' \code{WavTransf} computes the \emph{midpoint-interpolation} (MI) forward wavelet
#' transform of a curve of length \eqn{m} of (\eqn{d x d})-dimensional HPD matrices
#' (e.g. a multivariate spectrum) as described in (Chau and von Sachs, 2017).
#'
#' @param P a (\eqn{d,d,m})-dimensional array, with \eqn{m} a dyadic number.
#' @param order an odd integer taking values between 1 and 9 corresponding to the
#'  order of the MI refinement scheme.
#' @param jmax the maximum scale at which the wavelet coefficients are computed, such
#' that \eqn{jmax < log2(m)}.
#'
#' @examples
#' m <- function(){
#'  X <- matrix(complex(real=rnorm(9), imaginary=rnorm(9)), nrow=3)
#'  t(Conj(X)) %*% X
#' }
#' M <- replicate(100, m())
#' z <- rnorm(100)
#' w <- abs(z)/sum(abs(z))
#' KarchMean(M, w)
#'
#' @return Returns an integer sequence with cluster assignments.
#'
#' @seealso \code{\link{InvWavTransf}}
#'
#' @references Chau, J. and von Sachs, R. \emph{Positive-definite multivariate spectral
#' estimation: a geometric wavelet approach}. (Submitted)
#'
#' @export
WavTransf <- function(P, order = 5, jmax)
{
  J <- log(dim(P)[3], 2)
  if(!isTRUE(all.equal(as.integer(J), J))){
    print(paste0("Input length is non-dyadic, please change length ", dim(P)[3],
      " to dyadic number."))
    result <- NA
  } else {
    if(!(order %in% c(1,3,5,7,9))){
      print("Refinement order should be an odd integer between 1 and 9, by default set to 5")
      order <- 5
    }
    d <- dim(P)[1]
    M <- list()
    for (j in J:1){
      if (j == J){
        Mper <- unname(P)
      } else{
        Mper <- sapply(1:(dim(Mper)[3]/2), function(i) Mid(Mper[, , 2 * i - 1], Mper[, , 2 * i]), simplify = "array")
      }
      M[[j]] <- array(Mper[, , 1:2^j], dim = c(d, d, 2^j))
    }
    names(M) <- paste0("M.scale", 1:J)
    D <- list(M.scale1 = M[[1]])
    Nw <- list(N1 = 1, N3 = c(-1, 8, 1)/8, N5 = c(3, -22, 128, 22, -3)/128,
               N7 = c(-5, 44, -201, 1024, 201, -44, 5)/1024,
               N9 = c(35, -370, 1898, -6922, 32768, 6922, -1898, 370, -35)/32768)
    if(missing(jmax)) jmax <- J-1
    if(jmax > J-1){
      print(paste0("jmax cannot exceed maximum scale j=", J-1))
      jmax <- J-1
    }
    for (j in 1:jmax){
      tm1 <- Impute_man(M[[j]], (order-1)/2, Nw)
      iSqrt_tm1 <- sapply(1:2^j, function(l) iSqrt(tm1[, , l]), simplify = "array")
      D[[j + 1]] <- sapply(1:2^j, function(l) Logm(diag(d), (iSqrt_tm1[ , , l] %*%
                                M[[j + 1]][ , , 2 * l]) %*% iSqrt_tm1[ , , l]), simplify = "array")
      names(D)[j + 1] <- paste0("D.scale", j)
    }
    result <- list(D=D, M=M)
  }
  return(result)
}
