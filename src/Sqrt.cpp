#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export()]]

Eigen::MatrixXcd Sqrt(Eigen::MatrixXcd M) {

  return M.sqrt();

}
