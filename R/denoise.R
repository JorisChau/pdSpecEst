#' Applies nonlinear tree-structured wavelet thresholding based on cross-validation.
#'
#' @param P list, signal.
#' @return Returns a list containing the wavelet denoised signal..
#' @export
Denoise <- function(P = NULL, D, alpha = 0.75, lam = NULL, return=c("P", "wt"))
{
  J <- log(dim(P$p)[3], 2) - 1
  d <- dim(P$p)[1]
  E <- E_basis(d)

  if(is.null(lam))
  {
    P0 <- list(odd = list(t = P$t[c(T,F)], p = P$p[,,c(T,F)]), even = list(t = P$t[c(F,T)], p = P$p[,,c(F,T)]))
    wt <- list(odd = WavTrans(P0$odd, D)$wt, even = WavTrans(P0$even, D)$wt)

    coeff_E <- list(odd = list(), even = list())
    for (j in 1:(J-1))
    {
      coeff_E$odd[[j]] <- sapply(1:2^j, function(i) E_transf(wt$odd[[j + 1]]$coeff[, , i], E))
      coeff_E$even[[j]] <- sapply(1:2^j, function(i) E_transf(wt$even[[j + 1]]$coeff[, , i], E))
    }

    W <- diag(exp(alpha*(1:(J-1))))
    X <- matrix(c(rep(1, J-1), 1:(J-1)), nrow=J-1)
    mads <- lapply(1:2, function(l) sapply(1:(J-1), function(j) apply(coeff_E[[l]][[j]], 1, stats::mad)))
    Y <- colMeans(rbind(log(mads[[1]]), log(mads[[2]])))
    beta <- ((solve((t(X) %*% W) %*% X) %*% t(X)) %*% W) %*% Y
    sd.j <- function(j) exp(sum(beta * c(1, j)))

    coeff_tot <- list(odd = NULL, even = NULL)
    stand_coeff <- list(odd = list(), even = list())
    for (j in 1:(J-1))
    {
      stand_coeff$odd[[j]] <- coeff_E$odd[[j]]/sd.j(j)
      coeff_tot$odd <- cbind(coeff_tot$odd, stand_coeff$odd[[j]])
      stand_coeff$even[[j]] <- coeff_E$even[[j]]/sd.j(j)
      coeff_tot$even <- cbind(coeff_tot$even, stand_coeff$even[[j]])
    }

    CV_score <- function(lam){
      wt_lam <- list(odd = wt$odd, even = wt$even)
      coeff_E_lam <- list(odd = coeff_E$odd, even = coeff_E$even)
      for(j in 3:(J-1)){
        zeros_odd <- sapply(1:2^j, function(k) (abs(stand_coeff$odd[[j]][,k]) < lam) |
                              (coeff_E_lam$odd[[j-1]][,ceiling(k/2)] == 0))
        zeros_even <- sapply(1:2^j, function(k) (abs(stand_coeff$even[[j]][,k]) < lam) |
                               (coeff_E_lam$even[[j-1]][,ceiling(k/2)] == 0))
        coeff_E_lam$odd[[j]][zeros_odd] <- 0
        coeff_E_lam$even[[j]][zeros_even] <- 0
      }
      for (j in 1:(J-1)){
        wt_lam$odd[[j+1]]$coeff <- sapply(1:2^j, function(i) E_transf_inv(coeff_E_lam$odd[[j]][,i],
                                                                          E), simplify = "array")
        wt_lam$even[[j+1]]$coeff <- sapply(1:2^j, function(i) E_transf_inv(coeff_E_lam$even[[j]][,i],
                                                                           E), simplify = "array")
      }
      P_den <- list(odd = InvWavTrans(wt_lam$odd, D)$p, even = InvWavTrans(wt_lam$even, D)$p)
      P_pred_even <- sapply(1:(2^J-1), function(i) Mid(P_den$odd[,,i], P_den$odd[,,i+1]), simplify="array")
      P_pred_even <- unname(abind::abind(P_pred_even, P_den$odd[,,2^J], along = 3))
      P_pred_odd <- sapply(1:(2^J-1), function(i) Mid(P_den$even[,,i], P_den$even[,,i+1]), simplify="array")
      P_pred_odd <- unname(abind::abind(P_den$even[,,1], P_pred_odd, along = 3))

      return(sum(sapply(1:2^J, function(i) Dist(P_pred_even[,,i], P0$even$p[,,i])^2 + Dist(P_pred_odd[,,i], P0$odd$p[,,i])^2)))
    }

    gss <- function(x1, x4, tol = 0.01){
      x2 <- x4 + (x1-x4)/((1+sqrt(5))/2)
      x3 <- x1 + (x4-x1)/((1+sqrt(5))/2)
      fx2 <- CV_score(x2)
      fx3 <- CV_score(x3)
      i <- 0
      while(!isTRUE(all.equal(fx2, fx3)) & (abs(x2-x3) > tol) & i <= 20){
        i <- i+1
        # print(i)
        if(fx2 < fx3){
          x4 <- x3
          x3 <- x2
          fx3 <- fx2
          x2 <- x4 + (x1-x4)/((1+sqrt(5))/2)
          fx2 <- CV_score(x2)
        } else {
          x1 <- x2
          x2 <- x3
          fx2 <- fx3
          x3 <- x1 + (x4-x1)/((1+sqrt(5))/2)
          fx3 <- CV_score(x3)
        }
      }
      return(mean(c(x2,x3)))
    }

    lam1 <- sort(abs(c(coeff_tot$even)), decreasing = T)[5]
    lam0 <- sort(abs(c(coeff_tot$even)), decreasing = T)[round(0.2 * length(c(coeff_tot$even)))]
    lam.cv <- gss(lam0, lam1)
  }

  J <- J+1
  lam.cv <- ifelse(is.null(lam), 1/sqrt((1-log(2)/log(2^J*d^2))) * lam.cv, lam)

  wt <- WavTrans(P, D)$wt
  coeff_E <- list()
  for (j in 1:(J-1))
  {
    coeff_E[[j]] <- sapply(1:length(wt[[j + 1]]$k), function(i) E_transf(wt[[j + 1]]$coeff[, , i], E))
  }
  W <- diag(exp(alpha*(1:(J-1))))
  X <- matrix(c(rep(1, J-1), 1:(J-1)), nrow=J-1)
  Y <-  colMeans(log(sapply(1:(J-1), function(j) apply(coeff_E[[j]], 1, stats::mad))))
  beta <- ((solve((t(X) %*% W) %*% X) %*% t(X)) %*% W) %*% Y
  sd.j <- function(j) exp(sum(beta * c(1, j)))
  coeff_tot <- NULL
  stand_coeff <- list()
  for (j in 1:(J-1))
  {
    stand_coeff[[j]] <- coeff_E[[j]]/sd.j(j)
    coeff_tot <- cbind(coeff_tot, stand_coeff[[j]])
  }
  for(j in 3:(J-1)){
    zeros<- sapply(1:2^j, function(k) (abs(stand_coeff[[j]][,k]) < lam.cv) |
                     (coeff_E[[j-1]][,ceiling(k/2)] == 0))
    coeff_E[[j]][zeros] <- 0
  }
  for (j in 1:(J-1)){
    wt[[j+1]]$coeff <- sapply(1:2^j, function(i) E_transf_inv(coeff_E[[j]][,i], E), simplify = "array")
  }
  if(return == "P"){
    P_den <- InvWavTrans(wt, D, start = P$t[1], end = utils::tail(P$t, 1) + diff(c(P$t[1], P$t[2])))
  } else{
    P_den <- NULL
  }

  # den_full <- function(lam){
  #   wt_lam <- wt
  #   coeff_E_lam <- coeff_E
  #   for(j in 3:(J-1)){
  #     zeros<- sapply(1:2^j, function(k) (abs(stand_coeff[[j]][,k]) < lam) |
  #                           (coeff_E_lam[[j-1]][,ceiling(k/2)] == 0))
  #     coeff_E_lam[[j]][zeros] <- 0
  #   }
  #   for (j in 1:(J-1)){
  #     wt_lam[[j+1]]$coeff <- sapply(1:2^j, function(i) E_transf_inv(coeff_E_lam[[j]][,i], E), simplify = "array")
  #   }
  #   P_den <- InvWavTrans(wt_lam, D)$p
  #
  #   return(sum(sapply(1:2^J, function(i) Dist(P_den[,,i], P2$p[,,i])^2)))
  # }

  # t <- seq(from = lam0, to = lam1, length = 100)
  # test1 <- rep(NA, 100)
  # pb <- winProgressBar(title="Progress bar", label="0% done", min=0, max=100, initial=0)
  #
  # for(i in 1:100){
  #   test[i] <- score(t[i])
  #   info <- sprintf("%d%% done", round(i))
  #   setWinProgressBar(pb, round(i), label=info)
  # }
  # close(pb)
  # plot(t, test, type="l")

  return(list(P = P_den, wt = wt, lam = lam.cv, coeff = coeff_E))
}
