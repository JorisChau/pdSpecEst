#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


//' @export
// [[Rcpp::export()]]

Eigen::MatrixXcd iSqrt(Eigen::MatrixXcd M) {

  Eigen::MatrixXcd M1 = M.sqrt();

  return M1.inverse();

}
