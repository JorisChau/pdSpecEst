#' HPD depth-based confidence intervals for wavelet-based spectral matrix estimators
#'
#' @export
pdConfInt <- function(f, alpha = 0.05, conf_indices, N_samples = 250,
                      depth = c("gdd", "zonoid", "spatial"), ...) {

  ## Set variables
  J = log2(dim(f))[3]
  n = 2 * dim(f)[3]
  d = dim(f)[1]
  depth_vec = depth

  if(missing(conf_indices)) {
    conf_indices <- 1:(n/2)
    warning("'conf_indices' is missing, by default it is set to the entire domain")
  }
  if(!isTRUE(all(alpha > 0 & alpha < 1))) {
   stop(paste0("alpha = ", alpha, " should be a numeric vector (of quantiles) between 0 and 1"))
  }
  if (!isTRUE(all.equal(as.integer(J), J))) {
    stop(paste0("Input length is non-dyadic, please change length ", dim(f)[3],
                " to dyadic number."))
  }

  dots = list(...)
  B = (if(is.null(dots$B)) d else dots$B)
  method = (if(is.null(dots$method)) "multitaper" else dots$method)
  bias.corr = (if(is.null(dots$bias.corr)) F else dots$bias.corr)
  policy = (if(is.null(dots$policy)) "universal" else dots$policy)
  order = (if(is.null(dots$order)) 5 else dots$order)
  lam = (if(is.null(dots$lam)) 1 else dots$lam)
  metric = (if(is.null(dots$metric)) "Riemannian" else dots$metric)
  jmax.cv = (if(is.null(dots$jmax.cv)) J - 3 else dots$jmax.cv)
  jmax = (if(is.null(dots$jmax)) J - 3 else dots$jmax)
  n_ts = (if(is.null(dots$n)) n else dots$n_ts)
  # Jmax_out = (if(is.null(dots$Jmax_out)) J else dots$Jmax_out)
  return_f = (if(is.null(dots$return_f)) F else dots$return_f)

  ## Generate bootstrap samples
  if(!(n_ts == n)){
    f <- f[, , round(seq(from = 1, to = n/2, length = n_ts/2))]
  }
  f.sqrt <- sapply(1:(n_ts/2), function(i) pdSpecEst:::Sqrt(f[, , i]), simplify = "array")
  f_m <- array(dim = c(d, d, (n_ts/2), N_samples))
  cat("1. Generating bootstrap samples...")
  pb <- utils::txtProgressBar(1, 100, style = 3)

  for(m in 1:N_samples){
    ## Generate time series via Cramer
    chi <- matrix(nrow = 3, ncol = n_ts)
    chi[, 1:(n_ts/2 - 1)] <- replicate(n_ts/2 - 1, complex(3, rnorm(3, sd = sqrt(1/2)), rnorm(3, sd = sqrt(1/2))))
    chi[, c(n_ts/2, n_ts)] <- replicate(2, rnorm(3))
    chi[, (n_ts/2 + 1):(n_ts - 1)] <- Conj(chi[, 1:(n_ts/2 - 1)])

    f.sqrt1 <- array(c(f.sqrt, Conj(f.sqrt[, , (n_ts/2):1])), dim = c(3, 3, n_ts))
    ts <- sqrt(2 * pi) / sqrt(n_ts) * mvfft(t(sapply(1:n_ts, function(i) f.sqrt1[, , i] %*% chi[, i])), inverse = T)

    ## Compute m-th pre-smoothed periodogram
    per_m <- pdPgram(ts, B = B, method = method, bias.corr = bias.corr)$P

    ## Compute m-th spectral estimator
    coeffs <- pdSpecEst1D(per_m, order = order, policy = policy, metric = metric, alpha = lam,
                                jmax = jmax, jmax.cv = jmax.cv, return = "coeff")
    f_m[, , , m] <- InvWavTransf(coeffs$D, coeffs$M0, order = order, jmax = J, metric = metric)

    utils::setTxtProgressBar(pb, round(100 * m / N_samples))
  }
  close(pb)

  depth_CI <- rad_CI <- f_CI <- lapply(1:length(depth_vec), function(j) NA)
  names(depth_CI) <- names(rad_CI) <- names(f_CI) <- depth_vec

  for(depth in depth_vec){

    ## Generate integrated depth values
    if(depth == 'gdd'){

      dist <- matrix(0, nrow = N_samples, ncol = N_samples)
      grid <- expand.grid(1:N_samples, 1:N_samples)
      grid <- grid[grid$Var1 > grid$Var2, ]

      cat("2. Computing integrated depth values...")
      pb <- utils::txtProgressBar(1, 100, style = 3)

      for(l in 1:nrow(grid)){
        dist[grid[l, 1], grid[l, 2]] <- mean(sapply(conf_indices, function(t) pdDist(f_m[, , t, grid[l, 1]],
                                                                                f_m[, , t, grid[l, 2]])))
        utils::setTxtProgressBar(pb, round(100 * l / nrow(grid)))
      }
      close(pb)

      dist[upper.tri(dist)] <- t(dist)[upper.tri(dist)]
      dd <- exp(-colMeans(dist))

    } else {

      if(depth == "zonoid"){

        intDepth <- function(y, X) {
          depth.t <- numeric(dim(y)[3])
          for (t in 1:dim(y)[3]) {
            # E_y <- T_basis(E, y[, , t])
            depth.t[t] <- ddalpha::depth.zonoid(t(as.matrix(rep(0, d^2))), t(sapply(1:N_samples,
                                 function(m) pdSpecEst:::E_coeff(Logm(y[, , t], X[, , t, m])))))
          }
          return(mean(depth.t))
        }

      } else if(depth == "spatial"){

      intDepth <- function(y, X) {
        depth.t <- numeric(dim(y)[3])
        for (t in 1:dim(y)[3]) {
          y.isqrt <- iSqrt(y[, , t])
          log.yx <- sapply(1:N_samples, function(m) Logm(diag(d), (y.isqrt %*% X[, , t, m]) %*% y.isqrt), simplify = "array")
          depth.t[t] <- 1 - NormF(apply(sapply(1:N_samples,
                                  function(m) log.yx[, , m]/NormF(log.yx[, , m]), simplify = "array"), c(1, 2), mean))
        }
        return(mean(depth.t))
      }

    }

    cat("2. Computing integrated depth values...")
    pb <- utils::txtProgressBar(1, 100, style = 3)
    dd <- numeric(N_samples)
    for(m in 1:N_samples){
      dd[m] <- intDepth(f_m[, , conf_indices, m], f_m[, , conf_indices, , drop = F])
      utils::setTxtProgressBar(pb, round(100 * m / N_samples))
    }
    close(pb)

    }

    ## Construct confidence regions
    depth_CI[[which(depth_vec == depth)]] <- t(sapply(1:length(alpha), function(i) c(1, unname(quantile(dd, alpha[i])))))
    colnames(depth_CI[[which(depth_vec == depth)]]) <- c("max-depth", "min-depth", "radius")
    rownames(depth_CI[[which(depth_vec == depth)]]) <- sapply(1:length(alpha), function(i) paste0(100 * (1 - alpha[i]), "%-CI"))
    if(return_f){
      f_CI[[which(depth_vec == depth)]] <- f_m[, , , which(dd > depth_CI[[which(depth_vec == depth)]][1, 2])]
    }
    depth_CI[[which(depth_vec == depth)]][, 3] <- sapply(1:length(alpha), function(j) mean(sapply(conf_indices, function(i)
                            pdDist(f[, , i], f_m[, , i, which.min(abs(dd - depth_CI[[which(depth_vec == depth)]][j, 2]))]))))
  }
  return(list(depth_CI = depth_CI, f_CI = (if(return_f) f_CI else NULL)))
}


