#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


//' @export
// [[Rcpp::export()]]

Eigen::MatrixXcd Sqrt(Eigen::MatrixXcd M) {

  return M.sqrt();

}
