# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Orthonormal basis recomposition.
//'
//' \code{E_coeff_inv} computes the \eqn{(d x d)}-dimensional Hermitian matrix with (real-valued) coefficients
//' \code{coeff} in the orthonormal basis \code{E_basis(d)}.
//'
//' @param coeff  a numeric vector.
//' @param E  an array of orthonormal basis elements obtained by \code{E_basis(d)}.
//'
//' @examples
//'  E <- E_basis(3)
//'  coeff <- rnorm(9)
//'  E_coeff_inv(coeff, E)
//'
//' @seealso \code{\link{E_coeff}}, \code{\link{E_basis}}
//'
//' @export
// [[Rcpp::export()]]

arma::cx_mat E_coeff_inv(arma::vec coeff, arma::cx_cube E) {

 int d2 = coeff.size();

 int d = E.n_rows;

 arma::cx_mat M = arma::zeros<arma::cx_mat>(d, d);

 for(int i=0; i < d2; i++){

   arma::cx_mat Ei = E.slice(i);

   M += coeff[i] * Ei;

 }

 return M;

}


