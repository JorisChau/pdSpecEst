#' MI inverse wavelet transform.
#'
#' \code{InvWavTransf} computes the \emph{midpoint-interpolation} (MI) inverse wavelet
#' transform of a pyramid of Hermitian matrix-valued wavelet coefficients combined with an
#' array of coarse-scale Hermitian positive-definite midpoints as described in (Chau and von Sachs, 2017).
#'
#' @param D a list of wavelet coefficients and coarse-scale midpoints as obtained from
#'  the function \code{WavTransf}.
#' @param order an odd integer taking values between 1 and 9 corresponding to the
#'  order of the MI refinement scheme.
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
#' @seealso \code{\link{WavTransf}}
#'
#' @references Chau, J. and von Sachs, R. \emph{Positive-definite multivariate spectral
#' estimation: a geometric wavelet approach}. (Submitted)
#'
#' @export
InvWavTransf <- function(D, order=5)
{
  if(!(order %in% c(1,3,5,7,9))){
    print("Refinement order should be an odd integer between 1 and 9, by default set to 5")
    order <- 5
  }
  d <- nrow(D[[1]][, , 1])
  J <- length(D)
  Nw <- list(N1 = 1, N3 = c(-1, 8, 1)/8, N5 = c(3, -22, 128, 22, -3)/128,
             N7 = c(-5, 44, -201, 1024, 201, -44, 5)/1024,
             N9 = c(35, -370, 1898, -6922, 32768, 6922, -1898, 370, -35)/32768)
  m1 <- D[[1]]
  for (j in 1:(J-1)) {
    tm1 <- Impute_man(m1, (order-1)/2, Nw)
    m2 <- array(dim=c(d, d, 2^(j+1)))
    reconstr <- function(i){
      if(any(c(D[[j+1]][,,i]) != 0)){
        Sqrt_tm1 <- Sqrt(tm1[,,i])
        m2_even <- (Sqrt_tm1 %*% Expm(diag(d), D[[j+1]][,,i])) %*% Sqrt_tm1
      } else{
        m2_even <- tm1[,,i]
      }
      return(m2_even)
    }
    m2[,,c(F,T)] <- sapply(1:2^j, reconstr, simplify="array")
    m2[,,c(T,F)] <- sapply(1:2^j, function(i) solveMid(m2[,,2*i], m1[,,i]), simplify="array")
    m1 <- m2
  }
  return(m1)
}
