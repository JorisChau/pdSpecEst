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

## Iterative Cholesky with bias-correction (Dai & Guo, 2004)

Chol_inv <- function(R, bias.corr = T){
  d <- dim(R)[1]
  b <- (if(bias.corr) gamma(d - 1:d + 3/2) / (sqrt(d) * gamma(d - 1:d + 1)) else rep(1, d))
  S <- matrix(nrow = d, ncol = d)
  S[1, 1] <- R[1, 1]^2 / b[1]^2
  for(k in 1:(d - 1)){
    S[k + 1, 1:k] <- (t(R[k + 1, 1:k]) %*% solve(R[1:k, 1:k])) %*% S[1:k, 1:k]
    S[1:k, k + 1] <- Conj(S[k + 1, 1:k])
    S[k + 1, k + 1] <- R[k + 1, k + 1]^2 / b[k + 1]^2 + (t(S[k + 1, 1:k]) %*%
                                        solve(S[1:k, 1:k])) %*% Conj(S[k + 1, 1:k])
  }
  return(S)
}

## Basis Cholesky tangent space

E_chol <- function(R, inverse = F){
  d <- (if(!inverse) dim(R)[1] else round(sqrt(length(R))))

  if(!inverse){
    P <- c(Re(R)[lower.tri(R, diag = T)], Im(R)[lower.tri(R)])
  } else {
    P1 <- array(0, dim = c(d, d))
    P2 <- array(0i, dim = c(d, d))
    P1[lower.tri(P1, diag = T)] <- R[1:(d * (d + 1)/2)]
    P2[lower.tri(P2)] <- complex(imaginary = tail(R, -d * (d + 1)/2))
    P <- P1 + P2
  }
  return(P)
}



