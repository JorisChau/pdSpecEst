#' Impute midpoints at a finer scale (internal use only).
#'
#' @param M_j0 list of midpoints.
#' @param D integer.
#' @param j0 integer.
#' @param E list of basis elements.
#' @param Nw list of weights.
#' @return Returns a list of imputed midpoints at the next finer scale \code{j0 + 1}.
#' @export
Impute_man <- function(M_j0, D, j0, Nw)
{
  d <- dim(M_j0$m)[1]
  K_j0 <- M_j0$k
  K_j1 <- 2*K_j0 + 1
  tM_j1 <- list(k = K_j1, m = array(NA, dim = c(d, d, 2^j0)))
  P0 <- M_j0$m[, , 2^(j0-1)]
  M_j0.log <- sapply(1:2^j0, function(i) Logm(P0, M_j0$m[,,i]), simplify="array")

  if (D == 0)
  {
    for(i in 1:2^(j0-1)){
      P0 <- M_j0$m[, , 2 * i - 1]
      tM_j1$m[, , 2 * i - 1] <- Expm(P0, (1 / 4) * Logm(P0, M_j0$m[, , 2 * i]))
      tM_j1$m[, , 2 * i] <- Expm(P0, (5 / 4) * Logm(P0, M_j0$m[, , 2 * i]))
    }
  } else
  {
    for (i in 1:2^j0)
    {
      if(i == 1){
        tM_j1$m[,,i] <- Expm(M_j0$m[,,i], (1/4)*Logm(M_j0$m[,,i], M_j0$m[,,i+1]))
      } else if(i == 2^j0){
        tM_j1$m[,,i] <-   Expm(M_j0$m[,,i], (-1/4)*Logm(M_j0$m[,,i], M_j0$m[,,i-1]))
      } else{
        nbrs0 <- abs(K_j0 - K_j0[i]) <= D
        Di <- min(sum(nbrs0[1:(i-1)]), sum(nbrs0[(i+1):2^j0]))
        nbrs <- abs(K_j0 - K_j0[i]) <= Di
        # P0 <- M_j0$m[,,i]
        tM_j1$m[,,i] <- Expm(P0, apply(array(rep(Nw[[Di+1]], each=d^2), dim=c(d,d,2*Di+1)) *
                                         M_j0.log[,,(i-Di):(i+Di)], c(1,2), sum))
        # tM_j1$m[,,i] <- Expm(P0, apply(sapply((1:(2*Di+1))[-(Di+1)], function(l) Nw[[Di+1]][l] *
                                    # Logm(P0, M_j0$m[,,i-(Di+1)+l]), simplify="array"), c(1,2), sum))
      }
    }
  }
  return(tM_j1)
}
