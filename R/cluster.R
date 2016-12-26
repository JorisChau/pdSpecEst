#' Wavelet-based clustering of multivariate spectra.
#'
#' \code{SpecClust} performs clustering of multivariate spectral matrices via a two-step fuzzy
#' clustering algorithm in the manifold wavelet domain as detailed in (Chau and von Sachs, 2017).
#'
#' This is a two-step clustering algorithm.
#'
#' @param P list, curves of PD matrices to be clustered.
#' @param wt list, wavelet representations of curves of PD matrices to be clustered.
#' @param D integer, number of neighbors used in wavelet transform (only necessary when \code{is.null(wt)}).
#' @param K integer, number of clusters.
#' @param uptoscale integer, number of wavelet scales to be considered in iterative clustering algorithm.
#' @param anytime logical, output cluster assigments (yes/no) after each step in iterative clustering algorithm.
#' @param alpha numeric, hyperparameter for weighting function (cluster preference) in iterative clustering algorithm.
#' @param tau numeric, hyperparameter for weighting function (cluster preference) in iterative clustering algorithm.
#' @param m numeric, fuzziness parameter (\code{m > 1}).
#' @param d.uptoscale numeric, if fraction of components of wavelet coefficients at scale \code{j} is smaller than
#' \code{d.uptoscale} the iterative clustering algorithm terminates at scale \code{j}.
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
#' @references Chau, J. and von Sachs, R. \emph{Positive-definite multivariate spectral
#' estimation: a geometric wavelet approach}. (Submitted)
#'
#' @export
SpecClust <- function(P = NULL, wt = NULL, D = 2, K = 2, jmax = NULL, denoise = F,
                      m = 2, d.jmax=0.1){

  if(is.null(wt)){
    d <- dim(P$p[[1]])[1]
    S <- length(P$p)
    d.nonzero <- matrix(1, ncol = uptoscale, nrow = S)
    if(!denoise){
      wt <- lapply(1:S, function(s) WavTrans(list(t=P$P0$t, p=P$p[[s]]), D, uptoscale=uptoscale)$wt)
    } else{
      wt <- lapply(1:S, function(s) list())
      for(s in 1:S){
        wt_den <- Denoise(P = list(t=P$P0$t, p=P$p[[s]]), D, return = "wt")
        wt[[s]] <- wt_den$wt
        d.nonzero[s,] <- sapply(1:uptoscale, function(j) sum(as.logical(wt_den$coeff[[j]])) / (d^2 * 2^j))
      }
      uptoscale <- min(uptoscale, sum(colMeans(d.nonzero) > d.uptoscale))
      # print(paste0("Finest scale set to (j_max = ", uptoscale, "), as most wavelet coefficients at finer scales are zero..."))
    }
  } else {
    d <- dim(wt[[1]]$wt$b1$m)[1]
    S <- length(wt)
    wt_den <- wt
    wt <- lapply(1:S, function(s) list())
    d.nonzero <- matrix(ncol = uptoscale, nrow = S)
    for(s in 1:S){
      wt[[s]] <- wt_den[[s]]$wt
      d.nonzero[s,] <- sapply(1:uptoscale, function(j) sum(as.logical(wt_den[[s]]$coeff[[j]])) / (d^2 * 2^j))
    }
    uptoscale <- min(uptoscale, sum(colMeans(d.nonzero) > d.uptoscale))
    # print(paste0("Finest scale set to (j_max = ", uptoscale, "), as most wavelet coefficients at finer scales are zero..."))
  }

  M <- sapply(1:S, function(s) wt[[s]]$b1$m, simplify="array")
  centers <- M[,,,sample(1:S, K)]
  stopit <- F
  i <- 0

  distM <- function(M1, M2){
    Dist(M1[,,1], M2[,,1])^2 + Dist(M1[,,2], M2[,,2])^2
  }
  distj <- function(D1, D2, j){
    sum(sapply(1:2^j, function(i) 2^(uptoscale-j) * NormF(D1[,,i] - D2[,,i])^2))
  }

  # while((!stopit) & (i < 20)){
  while((i < 10)){
    distances <- t(sapply(1:S, function(s) t(sapply(1:K, function(k) distM(M[,,,s], centers[,,,k])))))
    mu_s <- function(s){
      if(!any(distances[s,] < .Machine$double.eps)){
        sapply(1:K, function(k) 1 / sum((distances[s,k] / distances[s,])^(1 / (m - 1))))
      } else{
        as.numeric(distances[s,] < .Machine$double.eps) / sum(distances[s,] < .Machine$double.eps)
      }
    }
    clusters <- t(sapply(1:S, mu_s))
    centers.new <- array(dim=c(d,d,2,K))
    for(k in 1:K){
      weights <- clusters[,k]^m / sum(clusters[,k]^m)
      centers.new[,,,k] <- sapply(1:2, function(i) kMean(M = do.call(rbind, lapply(1:S, function(s) M[,,i,s])),
                                                         mu = weights[2:S]/ifelse(cumsum(weights)[2:S] == 0, 1,
                                                                                  cumsum(weights)[2:S])), simplify="array")
    }
    # stopit <- (sum(sapply(1:K, function(k) distM(centers[,,,k], centers.new[,,,k]))) < 1E-4)
    # stopit <- ifelse(isTRUE(!stopit), F, T)
    # print(sum(sapply(1:K, function(k) distM(centers[,,,k], centers.new[,,,k]))))
    centers <- centers.new
    i <- i+1
  }

  dist.prior <- distances
  clust.prior <- clusters

  D <- lapply(1:uptoscale, function(j) sapply(1:S, function(s) wt[[s]][[j+1]]$coeff, simplify="array"))
  centers <- lapply(1:uptoscale, function(j) array(dim=c(d,d,2^j,K)))
  for(k in 1:K){
    weights <- clusters[,k]^m / sum(clusters[,k]^m)
    for(j in 1:uptoscale){
      centers[[j]][,,,k] <- sapply(1:2^j, function(i) apply(array(rep(weights, each=d^2), dim=c(d,d,S)) *
                                                              D[[j]][,,i,], c(1,2), sum), simplify="array")
    }
  }

  stopit <- F
  i <- 0
  # while((!stopit) & (i < 50)){
  while((i < 10)){
    distances <- t(sapply(1:S, function(s) t(sapply(1:K, function(k) sum(sapply(1:uptoscale, function(j)
      distj(D[[j]][,,,s], centers[[j]][,,,k], j)))))))
    mu_s <- function(s){
      if(!any(distances[s,] < .Machine$double.eps)){
        sapply(1:K, function(k) 1 / sum((distances[s,k] / distances[s,])^(1 / (m - 1))))
      } else{
        as.numeric(distances[s,] < .Machine$double.eps) / sum(distances[s,] < .Machine$double.eps)
      }
    }
    clusters <- t(sapply(1:S, mu_s))
    centers.new <- lapply(1:uptoscale, function(j) array(dim=c(d,d,2^j,K)))
    for(k in 1:K){
      weights <- clusters[,k]^m / sum(clusters[,k]^m)
      for(j in 1:uptoscale){
        centers.new[[j]][,,,k] <- sapply(1:2^j, function(i) apply(array(rep(weights, each=d^2), dim=c(d,d,S)) *
                                                                    D[[j]][,,i,], c(1,2), sum), simplify="array")
      }
    }
    # stopit <- sum(sapply(1:K, function(k) sum(sapply(1:uptoscale, function(j) distj(centers[[j]][,,,k],
                                                                                    # centers.new[[j]][,,,k], j)))))
    # print(stopit)
    # stopit <- ifelse(isTRUE(stopit > 1E-3), F, T)
    centers <- centers.new
    i <- i+1
  }

  return(clusters)
}




