## Basis Hermitian matrices

E_basis <- function(d) {
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
  V <- E
  d <- nrow(V)
  y_inv <- solve(y)
  Inner_y <- function(e_i, e_j) sum(diag((y_inv %*% e_i) %*% (y_inv %*% e_j)))
  Proj <- function(v,u) Inner_y(u, v) / Inner_y(u, u) * u

  U <- array(dim=c(d,d,d^2))
  U[,,1] <- V[,,1] / sqrt(Inner_y(V[,,1], V[,,1]))

  for(k in 2:d^2){
      U_k <- V[,,k] - apply(sapply(1:(k-1), function(j) Proj(V[,,k], U[,,j]),
                                   simplify = "array"), c(1,2), sum)
      U[,,k] <- U_k / sqrt(Inner_y(U_k, U_k))
  }
  return(U)
}






