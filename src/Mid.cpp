#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

//' Geodesic midpoint between HPD-matrices
//'
//' \code{Mid} calculates the geodesic midpoint between two Hermitian PD matrices as in
//' (Bhatia, 2009, Chapter 6).
//'
//' @param A,B square Hermitian positive-definite matrices (of equal dimension).
//'
//' @examples
//'  a <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  A <- t(Conj(a)) %*% a
//'  b <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  B <- t(Conj(b)) %*% b
//'  Mid(A, B)
//' @references Bhatia, R. (2009). \emph{Positive Definite Matrices}. New Jersey: Princeton University Press.
//'
//' @seealso \code{\link{KarchMean}}
//'
//' @export
// [[Rcpp::export()]]

Eigen::MatrixXcd Mid(Eigen::MatrixXcd A, Eigen::MatrixXcd B) {

  Eigen::MatrixXcd A1 = A.sqrt();

  Eigen::MatrixXcd A2 = A1.inverse();

  Eigen::MatrixXcd C = A2 * B * A2;

  return A1 * C.sqrt() * A1;

}
