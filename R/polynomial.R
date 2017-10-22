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
pdSplineReg <- function(P, f.hat0, lam, K, initial_step = 1, max.iter = 100, eps = 1E-3, ...) {

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
  if(missing(K)){
    K <- n/2
  }

  K.seq <- round(seq(from = 1, to = n, length = K))
  K.ind <- sapply(1:n, function(j) which.min(abs(K.seq - j)))
  delta.t <- mean(diff(K.seq))

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
    return(sapply(1:K, function(k) 0.5 * apply(sapply(which(K.ind == k),
            function(ki) grad_fA(gamma[, , k], P[, , ki]), simplify = "array"),
              c(1, 2), sum) + delta.t^3  * lam / 2 * (if(k %in% 2:(K-1)){
              2 * grad_fA(gamma[, , k], gamma[, , k + 1]) +
              2 * grad_fA(gamma[, , k], gamma[, , k - 1])  +
              2 * grad_gA(gamma[, , k], gamma[, , k + 1], gamma[, , k - 1]) +
              (if(k > 2) 2 * grad_gB(gamma[, , k - 1], gamma[, , k], gamma[, , k - 2]) else 0) +
              (if(k < K-1) 2 * grad_gB(gamma[, , k + 1], gamma[, , k], gamma[, , k + 2]) else 0)
              } else 0),
           simplify = "array"))
  }

  E <- function(gamma){
    gamma.isqrt <- sapply(1:K, function(k) iSqrt(gamma[, , k]), simplify = "array")
    return((0.5 * sum(sapply(1:n, function(i) pdDist(P[,,i], gamma[,,which.min(abs(K.seq - i))])^2)) +
              delta.t^3 * lam / 2 * sum(sapply(2:(K-1), function(k) NormF(ast(gamma.isqrt[, , k],
                Logm(gamma[, , k], gamma[, , k + 1]) + Logm(gamma[, , k], gamma[, , k - 1])))^2))))
  }

  backtrack <- function(p, gamma, E0){
    alpha <- initial_step
    t <- mean(sapply(1:K, function(k) sigma * pdSpecEst:::NormF(ast(pdSpecEst:::iSqrt(gamma[,,k]), p[,,k]))^2))
    E1 <- NULL
    while(is.null(E1) | isTRUE((E0 - E1) < (alpha * t))){
      E1 <- tryCatch({ E(sapply(1:dim(p)[3], function(i) Expm(gamma[, , i], alpha * p[, , i]),
                 simplify = "array")) }, error = function(e) return(NULL))
      alpha <- tau * alpha
    }
    return(alpha)
  }

  p <- -grad_E(f.hat0[, , K.seq])
  cost <- E(f.hat0[, , K.seq])
  gamma_0 <- f.hat0[, , K.seq]
  iter <- 0
  cost_diff <- -1

  while(isTRUE(abs(cost_diff) > eps) & isTRUE(cost_diff < 0) & isTRUE(iter < max.iter)){
    alpha <- backtrack(p, gamma_0, tail(cost, 1))
    gamma_1 <- sapply(1:K, function(k) Expm(gamma_0[, , k], alpha * p[, , k]),
                      simplify = "array")
    grad_1 <- grad_E(gamma_1)
    beta_1 <- sapply(1:K, function(k) pdSpecEst:::NormF(grad_1[, , k])^2/pdSpecEst:::NormF(p[, , k])^2)
    p <- -grad_1 + array(rep(beta_1, each = d^2), dim = c(d, d, K)) * p
    cost <- c(cost, E(gamma_1))
    cost_diff <- tail(cost, 1) - tail(cost, 2)[1]
    gamma_0 <- gamma_1
    iter <- iter + 1
    if(iter == max.iter){
      message("Reached maximum number of iterations")
    }
  }

  return(list(f.hat = gamma_0, cost = cost, total.iter = iter))

}

