#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

//' @export
// [[Rcpp::export()]]

Eigen::MatrixXcd solveMid(Eigen::MatrixXcd B, Eigen::MatrixXcd C) {

  return C * B.inverse() * C;

}
