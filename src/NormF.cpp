# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' @export
// [[Rcpp::export()]]

double NormF(arma::cx_mat M) {

  return norm(M, "fro");

}


