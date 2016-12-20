#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

//' @export
// [[Rcpp::export()]]

Eigen::MatrixXcd Mid(Eigen::MatrixXcd A, Eigen::MatrixXcd B) {

  Eigen::MatrixXcd A1 = A.sqrt();

  Eigen::MatrixXcd A2 = A1.inverse();

  Eigen::MatrixXcd C = A2 * B * A2;

  return A1 * C.sqrt() * A1;

}
