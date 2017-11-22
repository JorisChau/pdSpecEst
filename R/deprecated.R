## Wavelet-based spectral matrix estimation
pdSpecEst <- function(P, lam = NULL, order = 5, return = "f", alpha = 0.75) {
  .Deprecated("pdSpecEst1D")

  ## Define variables
  J <- log2(dim(P)[3])
  if (!isTRUE(all.equal(as.integer(J), J))) {
    warning(paste0("Input length is non-dyadic, please change length ", dim(P)[3],
                   " to dyadic number."))
  }
  stopifnot(isTRUE(all.equal(as.integer(J), J)))
  if (!(order %in% c(1, 3, 5, 7, 9))) {
    warning("Refinement order should be an odd integer between 1 and 9, by default set to 5")
    order <- 5
  }
  dim <- dim(P)[1]

  ## Find optimal threshold
  if (is.null(lam)) {

    ## Wavelet transforms
    P.half <- list(odd = P[, , c(T, F)], even = P[, , c(F, T)])
    D.half <- list(odd = WavTransf(P.half$odd, order)$D, even = WavTransf(P.half$even, order)$D)
    d <- list(odd = list(), even = list())
    for (j in 1:(J - 2)) {
      d$odd[[j]] <- sapply(1:2^j, function(k) E_coeff(D.half$odd[[j + 1]][, , k]))
      d$even[[j]] <- sapply(1:2^j, function(k) E_coeff(D.half$even[[j + 1]][, , k]))
    }

    ## Normalize variances wavelet coefficients
    W <- diag(exp(alpha * (1:(J - 2))))
    X <- matrix(c(rep(1, J - 2), 1:(J - 2)), nrow = J - 2)
    mads <- lapply(1:2, function(l) sapply(1:(J - 2), function(j) apply(d[[l]][[j]], 1, stats::mad)))
    Y <- colMeans(rbind(log(mads[[1]]), log(mads[[2]])))
    beta <- ((solve((t(X) %*% W) %*% X) %*% t(X)) %*% W) %*% Y
    sd.j <- function(j) exp(sum(beta * c(1, j)))

    d.vec <- list(odd = NULL, even = NULL)
    d.new <- list(odd = list(), even = list())
    for (j in 1:(J - 2)) {
      d.new$odd[[j]] <- d$odd[[j]]/sd.j(j)
      d.vec$odd <- cbind(d.vec$odd, d.new$odd[[j]])
      d.new$even[[j]] <- d$even[[j]]/sd.j(j)
      d.vec$even <- cbind(d.vec$even, d.new$even[[j]])
    }

    ## Cross-validation scores
    cv <- function(lam) {
      D.lam <- D.half
      d.lam <- d
      for (j in 3:(J - 2)) {
        zero.odd <- sapply(1:2^j, function(k) (abs(d.new$odd[[j]][, k]) <
                                                 lam) | (d.lam$odd[[j - 1]][, ceiling(k/2)] == 0))
        zero.even <- sapply(1:2^j, function(k) (abs(d.new$even[[j]][, k]) <
                                                  lam) | (d.lam$even[[j - 1]][, ceiling(k/2)] == 0))
        d.lam$odd[[j]][zero.odd] <- 0
        d.lam$even[[j]][zero.even] <- 0
      }
      for (j in 1:(J - 2)) {
        D.lam$odd[[j + 1]] <- sapply(1:2^j, function(k) E_coeff_inv(d.lam$odd[[j]][, k]),
                                     simplify = "array")
        D.lam$even[[j + 1]] <- sapply(1:2^j, function(k) E_coeff_inv(d.lam$even[[j]][, k]),
                                      simplify = "array")
      }
      f.hat <- list(odd = InvWavTransf(D.lam$odd, order), even = InvWavTransf(D.lam$even, order))

      ## Predicted points
      f.t.even <- sapply(1:(2^(J - 1) - 1), function(k) Mid(f.hat$odd[, , k], f.hat$odd[, , k + 1]),
                         simplify = "array")
      f.t.even <- array(c(f.t.even, f.hat$odd[, , 2^(J - 1)]), dim = c(dim, dim, 2^(J - 1)))
      f.t.odd <- sapply(1:(2^(J - 1) - 1), function(k) Mid(f.hat$even[, , k], f.hat$even[, , k + 1]),
                        simplify = "array")
      f.t.odd <- array(c(f.hat$even[, , 1], f.t.odd), dim = c(dim, dim, 2^(J - 1)))

      return(sum(sapply(1:2^(J - 1), function(k) RiemmDist(f.t.even[, , k], P.half$even[, , k])^2 +
                          RiemmDist(f.t.odd[, , k], P.half$odd[, , k])^2)))
    }

    ## Golden section search
    lam.range <- sort(abs(c(d.vec$even)), decreasing = T)[c(10, round(0.25 * length(c(d.vec$even))))]
    lam.cv <- gss(lam.range, cv)
  }

  ## Rescale threshold to twice number of observations
  lam.cv <- ifelse(is.null(lam), 1/sqrt((1 - log(2)/log(2^J * dim^2))) * lam.cv, lam)

  ## Transform original noisy data
  D <- WavTransf(P, order)$D
  d <- list()
  for (j in 1:(J - 1)) {
    d[[j]] <- sapply(1:2^j, function(k) E_coeff(D[[j + 1]][, , k]))
  }

  ## Normalize variances wavelet coefficients
  W <- diag(exp(alpha * (1:(J - 1))))
  X <- matrix(c(rep(1, J - 1), 1:(J - 1)), nrow = J - 1)
  Y <- colMeans(log(sapply(1:(J - 1), function(j) apply(d[[j]], 1, stats::mad))))
  beta <- ((solve((t(X) %*% W) %*% X) %*% t(X)) %*% W) %*% Y
  sd.j <- function(j) exp(sum(beta * c(1, j)))

  d.vec <- NULL
  d.new <- list()
  for (j in 1:(J - 1)) {
    d.new[[j]] <- d[[j]]/sd.j(j)
    d.vec <- cbind(d.vec, d.new[[j]])
  }
  d1 <- d.new

  ## Threshold coefficients
  for (j in 3:(J - 1)) {
    zero <- sapply(1:2^j, function(k) (abs(d.new[[j]][, k]) < lam.cv) |
                     (d[[j - 1]][, ceiling(k/2)] == 0))
    d.new[[j]][zero] <- 0
    d[[j]][zero] <- 0
  }

  ## Inverse transform denoised data
  for (j in 1:(J - 1)) {
    D[[j + 1]] <- sapply(1:2^j, function(k) E_coeff_inv(d[[j]][, k]), simplify = "array")
  }

  if (return == "f") {
    f <- InvWavTransf(D, order)
  } else {
    f <- NULL
  }
  return(list(f = f, D = D, lam = lam.cv, components = list(thresholded = d.new, not_thresholded = d1)))
}

## Basis Hermitian matrices
E_basis <- function(d) {
  .Deprecated("pdSpecEst:::E_coeff or pdSpecEst:::E_coeff_inv")
  E <- function(i, j) {
    E_ij <- matrix(0, nrow = d, ncol = d)
    if (i == j) {
      E_ij[i, j] <- 1
    } else if (i < j) {
      E_ij[c((j - 1) * d + i, (i - 1) * d + j)] <- 1/sqrt(2)
    } else {
      E_ij[c((j - 1) * d + i, (i - 1) * d + j)] <- complex(imaginary = c(-1, 1))/sqrt(2)
    }
    E_ij
  }

  indices <- expand.grid(1:d, 1:d)
  return(mapply(E, indices$Var1, indices$Var2, SIMPLIFY = "array"))
}

## Basis Tangent space
T_basis <- function(E, y) {
  .Deprecated("pdSpecEst:::T_coeff or pdSpecEst:::T_coeff_inv")
  d <- nrow(E)
  y.sqrt <- Sqrt(y)
  return(array(c(apply(E, 3, function(E) (y.sqrt %*% E) %*% y.sqrt)), dim = c(d, d, d^2)))
}

## Wavelet refinement scheme
Impute_man <- function(M_j0, D, Nw) {
  .Deprecated("pdSpecEst:::Impute1D")
  d <- dim(M_j0)[1]
  j0 <- log(dim(M_j0)[3], 2)
  tM_j1 <- array(dim = c(d, d, 2^(j0+1)))

  if (D == 0) {
    for (i in 1:2^(j0 - 1)) {
      P0 <- M_j0[, , 2 * i - 1]
      tM_j1[, , 2 * i - 1] <- Expm(P0, (1/4) * Logm(P0, M_j0[, , 2 * i]))
      tM_j1[, , 2 * i] <- Expm(P0, (5/4) * Logm(P0, M_j0[, , 2 * i]))
    }
  } else {
    for (i in 1:2^j0) {
      if (i == 1) {
        tM_j1[, , i] <- Expm(M_j0[, , i], (1/4) * Logm(M_j0[, , i], M_j0[, , i + 1]))
      } else if (i == 2^j0) {
        tM_j1[, , i] <- Expm(M_j0[, , i], (-1/4) * Logm(M_j0[, , i], M_j0[, , i - 1]))
      } else {
        nbrs0 <- abs((0:2^j0) - (0:2^j0)[i]) <= D
        Di <- min(sum(nbrs0[1:(i - 1)]), sum(nbrs0[(i + 1):2^j0]))
        # nbrs <- abs((0:2^j0) - (0:2^j0)[i]) <= Di
        tM_j1[, , i] <- KarchMean(M = M_j0[, , (i - Di):(i + Di)], w = Nw[[Di + 1]])
      }
    }
  }
  return(tM_j1)
}

