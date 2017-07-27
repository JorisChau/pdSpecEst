#' Generate intrinsic HPD polynomials
#'
#' @export
pdPolynomial <- function(p0, v0, delta.t=0.01, steps = 100) {

  d <- dim(p0)[1]
  k <- dim(v0)[3]
  p <- array(dim=c(d,d,steps))
  p[,,1] <- p0
  t <- 0
  vi <- v0

    for(ti in 1:(steps-1)){
      w <- vi[,,1]
      vi[,,1:(k-1)] <- sapply(1:(k-1), function(i) pdSpecEst:::ParTrans(p[,,ti], delta.t * w, vi[,,i] +
                                                              delta.t * vi[,,i+1]), simplify = "array")
      vi[,,k] <- pdSpecEst:::ParTrans(p[,,ti], delta.t * w, vi[,,k])
      p[,,ti+1] <- Expm(p[,,ti], delta.t * w)
    }

  return(p)
}
