#' Computes MI forward wavelet transform of a curve of Hermitian PD matrices.
#'
#' @param P list, sampled curve of Hermitian PD matrices, defaults to \code{NULL}.
#' @param D integer, number of neighbors included in MI scheme.
#' @param M list, midpoint pyramid, defaults to \code{NULL}.
#' @return Returns a list with values:
#'  \itemize{
#'    \item \code{wt} pyramid of wavelet coefficients (Hermitian matrices).
#'    \item \code{M} pyramid of midpoints.
#'  }
#' If \code{M == NULL} is \code{FALSE}, the function computes the pyramid of wavelet
#' coefficients directly from the given midpoint pyramid \code{M}.
#' @export
WavTrans <- function(P = NULL, D, M = NULL, uptoscale = NULL)
{
  if (!is.null(P))
  {
    J <- log(dim(P$p)[3], 2)
    tper <- P$t
    pper <- P$p
    d <- dim(P$p)[1]
    M <- lapply(1:J, function(j) list(k = (-Dnew):(2^j - 1 + Dnew),
                                      m = array(NA, dim = c(d, d, 2^j + 2 * Dnew))))
    for (j in J:1)
    {
      if (j == J)
      {
        Mper <- unname(pper)
      } else
      {
        Mper <- sapply(1:(dim(Mper)[3]/2), function(i) Mid(Mper[, , 2 * i - 1], Mper[, , 2 * i]), simplify = "array")
      }
      M[[j]]$m <- array(Mper[, , Dnew * (2^j - 1) + 1:(2^j + 2 * Dnew)], dim = c(d, d, 2^j + 2 * Dnew))
    }
    names(M) <- paste0("m", 1:J)
  } else
  {
    J <- length(M)
    d <- dim(M$m1$m)[1]
  }
  wt <- list(b1 = M$m1)
  Nw <- list(N1 = 1, N3 = c(-1, 8, 1)/8, N5 = c(3, -22, 128, 22, -3)/128,
             N7 = c(-5, 44, -201, 1024, 201, -44, 5)/1024,
             N9 = c(35, -370, 1898, -6922, 32768, 6922, -1898, 370, -35)/32768)

  uptoJ <- ifelse(is.null(uptoscale), J-1, uptoscale)

  for (j in 1:uptoJ)
  {
    tm1 <- Impute_man(M[[j]], D, j, Nw)
    iSqrt_tm1 <- sapply(1:2^j, function(l) iSqrt(tm1$m[, , l]), simplify = "array")
    wt[[j + 1]] <- list(k = (tm1$k - 1)/2, coeff = sapply(1:2^j, function(l) Logm(diag(d),
                            (iSqrt_tm1[ , , l] %*% M[[j + 1]]$m[ , , 2 * l]) %*% iSqrt_tm1[ , , l]), simplify = "array"))
    names(wt)[j + 1] <- paste0("a", j)
  }
  return(list(wt = wt, M = M))
}
