#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

//' Exponential map
//'
//' \code{Expm(P, H)} computes the projection of a Hermitian matrix \code{H} from the tangent space at a Hermitian
//' PD matrix \code{P} to the Riemannian manifold of Hermitian PD matrices via the
//' exponential map as in (Pennec, 2006). This is the unique inverse of the logarithmic map \code{\link{Logm}}.
//'
//' @param P a square Hermitian positive-definite matrix.
//' @param H a square Hermitian matrix (of equal dimension as \code{P}).
//'
//' @examples
//'  H <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  diag(H) <- rnorm(3)
//'  H[lower.tri(H)] <- t(Conj(H))[lower.tri(H)]
//'  p <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  P <- t(Conj(p)) %*% p
//'  Expm(P, H)
//'
//' @references
//' Pennec, X. (2006). Intrinsic statistics on Riemannian manifolds: Basic tools for geometric
//' measurements. \emph{Journal of Mathematical Imaging and Vision} 25(1), 127-154.
//'
//' @seealso \code{\link{Logm}}
//'
//' @export
// [[Rcpp::export()]]

Eigen::MatrixXcd Expm(Eigen::MatrixXcd P, Eigen::MatrixXcd H) {

  Eigen::MatrixXcd P1 = P.sqrt();

  Eigen::MatrixXcd P2 = P1.inverse();

  Eigen::MatrixXcd P3 = (P2 * H * P2).exp();

  return P1 * P3 * P1;

}


