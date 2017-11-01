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
  J.out = (if(is.null(dots$J.out)) J else dots$J.out)
  return_f = (if(is.null(dots$return_f)) F else dots$return_f)
  jmax = min((if(is.null(dots$jmax)) J - 3 else dots$jmax), J.out - 1)
  jmax.cv = min(if(is.null(dots$jmax.cv)) J - 3 else dots$jmax.cv, J.out - 1)
  f.0 = dots$f.0

  if(missing(conf_indices)) {
    conf_indices <- list(1:2^J.out)
    warning("'conf_indices' is missing, by default it is set to the entire domain")
  } else {
    conf_indices <- lapply(1:length(conf_indices), function(L)
      max(round(conf_indices[[L]][1] * 2^J.out), 1):min(round(conf_indices[[L]][2] * 2^J.out), 2^J.out))
  }

  if(J.out < J){
    f.hat.J <- f
    if(!is.null(f.0)) f.0.J <- f.0
    for (j in J:(J.out + 1)) {
      f.hat.J <- sapply(1:(dim(f.hat.J)[3]/2), function(i) Mid(f.hat.J[, , 2 * i - 1], f.hat.J[, , 2 * i]),
                        simplify = "array")
      if(!is.null(f.0)){
        f.0.J <- sapply(1:(dim(f.0.J)[3]/2), function(i) Mid(f.0.J[, , 2 * i - 1], f.0.J[, , 2 * i]),
                        simplify = "array")
      }
    }
  }

  ## Generate bootstrap samples
  f.sqrt <- sapply(1:(n/2), function(i) pdSpecEst:::Sqrt(f[, , i]), simplify = "array")
  f_m <- array(dim = c(d, d, 2^J.out, N_samples))
  cat("1. Generating bootstrap samples...")
  pb <- utils::txtProgressBar(1, 100, style = 3)

  for(m in 1:N_samples){
    ## Generate time series via Cramer
    chi <- matrix(nrow = d, ncol = n)
    chi[, 1:(n/2 - 1)] <- replicate(n/2 - 1, complex(d, rnorm(d, sd = sqrt(1/2)), rnorm(d, sd = sqrt(1/2))))
    chi[, c(n/2, n)] <- replicate(2, rnorm(d))
    chi[, (n/2 + 1):(n - 1)] <- Conj(chi[, 1:(n/2 - 1)])

    f.sqrt1 <- array(c(f.sqrt, Conj(f.sqrt[, , (n/2):1])), dim = c(d, d, n))
    ts <- sqrt(2 * pi) / sqrt(n) * mvfft(t(sapply(1:n, function(i) f.sqrt1[, , i] %*% chi[, i])), inverse = T)

    ## Compute m-th pre-smoothed periodogram
    per_m <- pdPgram(ts, B = B, method = method, bias.corr = bias.corr)$P

    ## Compute m-th spectral estimator
    coeffs <- pdSpecEst1D(per_m, order = order, policy = policy, alpha = lam,
                          jmax = jmax, jmax.cv = jmax.cv, return = "coeff")
    f_m[, , , m] <- InvWavTransf(coeffs$D, coeffs$M0, order = order, jmax = J.out)

    utils::setTxtProgressBar(pb, round(100 * m / N_samples))
  }
  close(pb)

  depth_CI <- lapply(1:length(depth_vec), function(j) array(dim = c(length(alpha), 3, length(conf_indices))))
  f_CI <- lapply(1:length(depth_vec), function(j) lapply(1:length(conf_indices), function(I) NA))
  cover_f.0 <- lapply(1:length(depth_vec), function(j) matrix(nrow = length(alpha), ncol = length(conf_indices)))
  depth_f.0 <- lapply(1:length(depth_vec), function(j) numeric(length(conf_indices)))
  names(depth_CI) <- names(f_CI) <- names(cover_f.0) <- names(depth_f.0) <- depth_vec

  for(depth in depth_vec){
    dimnames(depth_CI[[which(depth_vec == depth)]]) <- list(sapply(1:length(alpha), function(i) paste0(100 * (1 - alpha[i]), "%-CI")),
                                                            c('max-depth', 'min-depth', 'radius'),
                                                            sapply(1:length(conf_indices), function(i) paste0("indices_", i)))
    dimnames(cover_f.0[[which(depth_vec == depth)]]) <- dimnames(depth_CI[[which(depth_vec == depth)]])[c(1, 3)]
    names(depth_f.0[[which(depth_vec == depth)]]) <- dimnames(depth_CI[[which(depth_vec == depth)]])[[3]]
    names(f_CI[[which(depth_vec == depth)]]) <- dimnames(depth_CI[[which(depth_vec == depth)]])[[3]]
    dd <- matrix(nrow = length(conf_indices), ncol = N_samples)

    for(I in 1:length(conf_indices)){

      ## Generate integrated depth values
      if(depth == 'gdd'){

        dist <- matrix(0, nrow = N_samples, ncol = N_samples)
        grid <- expand.grid(1:N_samples, 1:N_samples)
        grid <- grid[grid$Var1 > grid$Var2, ]

        cat("2. Computing integrated depth values...")
        pb <- utils::txtProgressBar(1, 100, style = 3)

        for(l in 1:nrow(grid)){
          dist[grid[l, 1], grid[l, 2]] <- mean(sapply(conf_indices[[I]], function(t) pdDist(f_m[, , t, grid[l, 1]],
                                                                                            f_m[, , t, grid[l, 2]])))
          utils::setTxtProgressBar(pb, round(100 * l / nrow(grid)))
        }
        close(pb)

        dist[upper.tri(dist)] <- t(dist)[upper.tri(dist)]
        dd[I, ] <- exp(-colMeans(dist))

      } else {

        if(depth == "zonoid"){

          intDepth <- function(y, X) {
            depth.t <- numeric(dim(y)[3])
            for (t in 1:dim(y)[3]) {
              # E_y <- T_basis(E, y[, , t])
              depth.t[t] <- ddalpha::depth.zonoid(t(as.matrix(rep(0, d^2))), t(sapply(1:N_samples,
                                                                                      function(m) E_coeff(Logm(y[, , t], X[, , t, m])))))
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
        for(m in 1:N_samples){
          dd[I, m] <- intDepth(f_m[, , conf_indices[[I]], m], f_m[, , conf_indices[[I]], , drop = F])
          utils::setTxtProgressBar(pb, round(100 * m / N_samples))
        }
        close(pb)

      }

      ## Construct confidence regions
      depth_CI[[which(depth_vec == depth)]][, 1:2, I] <- t(sapply(1:length(alpha), function(i) c(1, unname(quantile(dd[I, ], alpha[i])))))
      if(return_f){
        f_CI[[which(depth_vec == depth)]][[I]] <- f_m[, , , which(dd[I, ] > depth_CI[[which(depth_vec == depth)]][1, 2, I])]
      }
      depth_CI[[which(depth_vec == depth)]][, 3, I] <- sapply(1:length(alpha), function(j) mean(sapply(conf_indices[[I]], function(i)
        pdDist(f.hat.J[, , i], f_m[, , i, which.min(abs(dd[I, ] - depth_CI[[which(depth_vec == depth)]][j, 2, I]))]))))

      if(!is.null(f.0)){
        depth_f.0[[which(depth_vec == depth)]][I] <- if(depth == 'gdd'){
          pdDepth(f.0.J[, , conf_indices[[I]]], f_m[, , conf_indices[[I]], , drop = F], method = 'gdd') } else {
            intDepth(f.0.J[, , conf_indices[[I]]], f_m[, , conf_indices[[I]], , drop = F]) }
        cover_f.0[[which(depth_vec == depth)]][, I] <- (depth_f.0[[which(depth_vec == depth)]][I] > depth_CI[[which(depth_vec == depth)]][, 2, I])
      }
    }
  }
  return(list(depth_CI = depth_CI, f_CI = (if(return_f) f_CI else NULL), f_J.out = f.hat.J,
              cover_f = (if(!is.null(f.0)) cover_f.0 else NULL), depth_f = (if(!is.null(f.0)) depth_f.0 else NULL)))
}


