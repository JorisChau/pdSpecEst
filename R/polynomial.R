#' Generate intrinsic HPD polynomials
#'
#' @export
pdPolynomial <- function(p0, v0, delta.t=0.01, steps = 100) {

  d <- dim(p0)[1]
  k <- dim(v0)[3]
  p <- array(dim = c(d, d, steps))
  p[ , , 1] <- p0
  vi <- v0

  for(ti in 1:(steps - 1)){
    w <- vi[ , , 1]
    vi[ , , 1:(k - 1)] <- sapply(1:(k - 1), function(i) ParTrans(p[ , , ti], delta.t * w, vi[ , , i] +
                                                            delta.t * vi[ , , i+1]), simplify = "array")
    vi[ , , k] <- ParTrans(p[ , , ti], delta.t * w, vi[ , , k])
    p[ , , ti+1] <- Expm(p[ , , ti], delta.t * w)
  }

  return(p)
}

#' HPD cubic smoothing spline regression
#'
#' @export
pdSplineReg <- function(P, f.hat0, lam, initial_step = 1, max.iter = 100, eps = 1E-3, ...) {

  dots <- list(...)
  tau <- dots$tau
  if (is.null(tau)) {
    tau <- 0.5
  }
  sigma <- dots$sigma
  if (is.null(sigma)) {
    sigma <- 0.5
  }

  d <- dim(P)[1]
  n <- dim(P)[3]
  if(missing(lam)){
    lam <- 1
  }

  ast <- function(A, B) (A %*% B) %*% t(Conj(A))

  DLog <- function(A, H){
    e <- eigen(A, symmetric = T)
    e_vec <- e$vectors
    e_val <- e$values
    H1 <- (t(Conj(e_vec)) %*% H) %*% e_vec
    grid <- expand.grid(1:d, 1:d)
    Z1 <- array(mapply(function(i1, i2) if(i1 == i2) 1/e_val[i1] else (log(e_val[i1]) - log(e_val[i2])) /
                   (e_val[i1] - e_val[i2]), grid$Var1, grid$Var2), dim = c(d, d))
    return(e_vec %*% (H1 * Z1) %*% t(Conj(e_vec)))
  }

  grad_fA <- function(A, B) -2 * Logm(A, B)

  grad_gA <- function(A, B, C){
    B1 <- pdSpecEst:::iSqrt(B)
    C1 <- pdSpecEst:::iSqrt(C)
    return(ast(A, ast(B1, DLog(ast(B1, A), Logm(diag(d),
            ast(pdSpecEst:::Sqrt(solve(ast(B1, C))), ast(B1, A)))))) +
           ast(A, ast(C1, DLog(ast(C1, A), Logm(diag(d),
            ast(pdSpecEst:::Sqrt(ast(C1, A)), solve(ast(C1, B))))))))
  }

  grad_gB <- function(A, B, C){
    A1 <- pdSpecEst:::iSqrt(A)
    return(ast(B, ast(A1, DLog(ast(A1, B), Logm(diag(d), ast(A1, C))))))
  }

  grad_E <- function(gamma){
    return(sapply(1:n, function(j) 0.5 / n * grad_fA(gamma[, , j], P[, , j]) +
              n^3 / (n - 2) * lam / 2 * (if(j %in% 2:(n-1)){
              2 * grad_fA(gamma[, , j], gamma[, , j + 1]) +
              2 * grad_fA(gamma[, , j], gamma[, , j - 1])  +
              2 * grad_gA(gamma[, , j], gamma[, , j + 1], gamma[, , j - 1]) +
              (if(j > 2) 2 * grad_gB(gamma[, , j - 1], gamma[, , j], gamma[, , j - 2]) else 0) +
              (if(j < n-1) 2 * grad_gB(gamma[, , j + 1], gamma[, , j], gamma[, , j + 2]) else 0)
              } else 0),
           simplify = "array"))
  }

  E <- function(gamma){
    gamma.isqrt <- sapply(1:n, function(i) pdSpecEst:::iSqrt(gamma[, , i]), simplify = "array")
    return((0.5 * mean(sapply(1:n, function(i) pdDist(P[,,i], gamma[,,i])^2)) +
              n^3 * lam / 2 * mean(sapply(2:(n-1), function(i) pdSpecEst:::NormF(ast(gamma.isqrt[, , i],
                Logm(gamma[, , i], gamma[, , i + 1]) + Logm(gamma[, , i], gamma[, , i - 1])))^2))))
  }

  backtrack <- function(p, gamma, E0){
    alpha <- initial_step
    t <- mean(sapply(1:n, function(i) sigma * pdSpecEst:::NormF(ast(pdSpecEst:::iSqrt(gamma[,,i]), p[,,i]))^2))
    E1 <- NULL
    while(is.null(E1) | isTRUE((E0 - E1) < (alpha * t))){
      E1 <- tryCatch({ E(sapply(1:dim(p)[3], function(i) Expm(gamma[, , i], alpha * p[, , i]),
                 simplify = "array")) }, error = function(e) return(NULL))
      alpha <- tau * alpha
    }
    return(alpha)
  }

  p <- -grad_E(f.hat0)
  cost <- E(f.hat0)
  gamma_0 <- f.hat0
  k <- 0
  cost_diff <- -1

  while(isTRUE(abs(cost_diff) > eps) & isTRUE(cost_diff < 0) & isTRUE(k < max.iter)){
    alpha <- backtrack(p, gamma_0, tail(cost, 1))
    gamma_1 <- sapply(1:n, function(i) Expm(gamma_0[, , i], alpha * p[, , i]),
                      simplify = "array")
    grad_1 <- grad_E(gamma_1)
    beta_1 <- sapply(1:n, function(i) pdSpecEst:::NormF(grad_1[, , i])^2/pdSpecEst:::NormF(p[, , i])^2)
    p <- -grad_1 + array(rep(beta_1, each = d^2), dim = c(d, d, n)) * p
    cost <- c(cost, E(gamma_1))
    cost_diff <- tail(cost, 1) - tail(cost, 2)[1]
    gamma_0 <- gamma_1
    k <- k + 1
    if(k == max.iter){
      message("Reached maximum number of iterations")
    }
  }

  return(list(f.hat = gamma_0, cost = cost, total.iter = k))

}

