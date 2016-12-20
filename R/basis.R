#' Orthonormal basis for space of Hermitian matrices.
#'
#' \code{E_basis} constructs an orthonormal basis (with respect to the Frobenius norm) of matrix elements
#' for the real-valued vector space of Hermitian matrices.
#'
#' @param d the row-dimension of the matrix elements.
#' @return \code{E_basis(d)} returns a \eqn{(d, d, d^2)}-dimensional array, containing the \eqn{d^2} orthonormal basis elements
#'  (i.e. \eqn{d x d}-dimensional matrices).
#' @examples
#'  E_basis(3)
#'
#' @useDynLib pdSpecEst
#' @importFrom Rcpp sourceCpp
#' @export
E_basis <- function(d) {
  E <- function(i, j) {
    E_ij <- matrix(0, nrow = d, ncol = d)
    if (i == j) {
      E_ij[i, j] <- 1
    } else if (i < j) {
      E_ij[c((j - 1) * d + i, (i - 1) * d + j)] <- 1 / sqrt(2)
    } else{
      E_ij[c((j - 1) * d + i, (i - 1) * d + j)] <- complex(imaginary = c(-1, 1)) / sqrt(2)
    }
    E_ij
  }
  indices <- expand.grid(1:d, 1:d)
  return(mapply(E, indices$Var1, indices$Var2, SIMPLIFY = "array"))
}

