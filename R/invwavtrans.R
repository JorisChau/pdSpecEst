#' Computes MI inverse wavelet transform.
#'
#' @param wt list, pyramid of wavelet coefficients and midpoints on scale \eqn{j_0}.
#' @param D integer, number of neighbors included in MI scheme.
#' @param start numeric, left-boundary of frequency range, defaults to \code{0}.
#' @param end numeric, right-boundary of frequency range, defaults to \code{pi}.
#' @return Returns a reconstruction of the curve of Hermitian PD matrices in the form
#'  of a list with values:
#'  \itemize{
#'    \item \code{t} the frequency points of the sampled curve.
#'    \item \code{p} an array of Hermitian PD matrices.
#'  }
#'  If \code{coh == TRUE} returns the squared coherence of the reconstructed
#'  curve of Hermitian PD matrices.
#' @export
InvWavTrans <- function(wt, D, start = 0, end = pi)
{
  Dnew <- 0
  d <- nrow(wt$b1$m[, , 1])
  J <- length(wt)
  Nw <- list(N1 = 1, N3 = c(-1, 8, 1)/8, N5 = c(3, -22, 128, 22, -3)/128,
             N7 = c(-5, 44, -201, 1024, 201, -44, 5)/1024,
             N9 = c(35, -370, 1898, -6922, 32768, 6922, -1898, 370, -35)/32768)
  m1 <- wt[[1]]

  for (j in 1:(J-1))
  {
    tm1 <- Impute_man(m1, D, j, Nw)
    m2 <- list(k = sort(c(tm1$k - 1, tm1$k)), m = array(NA, c(d, d, 2^(j+1))))

    reconstr <- function(i){
      if(any(c(wt[[j+1]]$coeff[,,i]) != 0)){
        Sqrt_tm1 <- Sqrt(tm1$m[,,i])
        m2_even <- (Sqrt_tm1 %*% Expm(diag(d), wt[[j+1]]$coeff[,,i])) %*% Sqrt_tm1
      } else{
        m2_even <- tm1$m[,,i]
      }
      return(m2_even)
    }

    m2$m[,,c(F,T)] <- sapply(1:2^j, reconstr, simplify="array")
    m2$m[,,c(T,F)] <- sapply(1:2^j, function(i) solveMid(m2$m[,,2*i], m1$m[,,ceiling(Dnew/2)+i]), simplify="array")
    m1 <- m2
  }

  P <- list(t = m1$k[2 * ceiling(Dnew/2) + 1:2^J] * 2^(-J) * (end - start),
            p = m1$m[, , 2 * ceiling(Dnew/2) + 1:2^J])

  return(P)
}
