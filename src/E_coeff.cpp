# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


//' Orthonormal basis decomposition.
//'
//' \code{E_coeff} computes the \eqn{d^2} real-valued basis coefficients with respect to an orthonormal
//' (in terms of the Frobenius metric) basis of the real vector space of \eqn{(d x d)}-dimensional
//' Hermitian matrices.
//'
//' @param H a \eqn{(d x d)}-dimensional Hermitian matrix.
//' @param E an array of orthonormal \eqn{(d x d)}-dimensional basis elements obtained by \code{E_basis(d)}.
//'
//' @examples
//'  E <- E_basis(3)
//'  h <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  H <- t(Conj(h)) %*% h
//'  E_coeff(H, E)
//'
//' @seealso \code{\link{E_coeff_inv}}, \code{\link{E_basis}}
//'
//' @export
// [[Rcpp::export()]]

arma::vec E_coeff(arma::cx_mat H, arma::cx_cube E) {

 int d = H.n_rows;

 arma::vec coeff(d * d);

 for(int i=0; i < (d * d); i++){

   arma::cx_mat Ei = E.slice(i);

   arma::cx_double coeff0 = accu(H % conj(Ei));

   coeff[i] = real(coeff0);

 }

 return coeff;

}



