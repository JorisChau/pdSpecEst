#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Geodesic midpoint between HPD-matrices
//'
//' \code{Mid} calculates the geodesic midpoint between two Hermitian PD matrices as in
//' (Bhatia, 2009, Chapter 6).
//'
//' @param A,B Hermitian positive definite matrices (of equal dimension).
//'
//' @examples
//'  a <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  A <- t(Conj(a)) %*% a
//'  b <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  B <- t(Conj(b)) %*% b
//'  Mid(A, B)
//' @references Bhatia, R. (2009). \emph{Positive Definite Matrices}. New Jersey: Princeton University Press.
//'
//' @seealso \code{\link{pdMean}}
//'
//' @export
// [[Rcpp::export()]]

arma::cx_mat Mid(arma::cx_mat A, arma::cx_mat B) {
  arma::cx_mat A1 = arma::sqrtmat_sympd(A);
  arma::cx_mat A2 = arma::inv_sympd(A1);
  arma::cx_mat C = A2 * B * A2;
  return A1 * arma::sqrtmat_sympd(C) * A1;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat Sqrt(arma::cx_mat M) {
  return arma::sqrtmat_sympd(M);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat Chol_C(arma::cx_mat M) {
  return arma::chol(M);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat iSqrt(arma::cx_mat M) {
  arma::cx_mat M1 = arma::sqrtmat_sympd(M);
  return arma::inv_sympd(M1);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

double NormF(arma::cx_mat M) {
  return arma::norm(M, "fro");
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_cube wavPyr_C(arma::cx_cube P, int L, int J, arma::ivec Nj, std::string metric) {

  // Set parameters
  int d = P.n_rows;
  int n = P.n_slices;
  int N = (2 * L + 1) * n;
  arma::cx_cube P1 = P;
  arma::cx_cube P_rev(d, d, n);
  arma::cx_cube P_per(d, d, N);
  arma::cx_cube M_per(d, d, arma::sum(Nj));

  // Transform according to metric
  for(int i = 0; i < n; ++i) {
    if(metric == "logEuclidean") {
      P1.slice(i) = arma::logmat_sympd(P.slice(i));
    }
    else if(metric == "Cholesky") {
      P1.slice(i) = arma::chol(P.slice(i));
    }
    else if(metric == "rootEuclidean") {
      P1.slice(i) = arma::sqrtmat_sympd(P.slice(i));
    }
  }
  // Construct reverse cube
  P_rev.slice(0) = P1.slice(n - 1);
  for(int i = 1; i < n; ++i) {
    P_rev.slice(i) = P1.slice(n - i);
  }
  // Construct periodized cube
  for(int l = 0; l < (2 * L + 1); ++l) {
    if((L % 2 == 0) && (l % 2 == 0)) {
      P_per.slices(l * n, (l + 1) * n - 1) = P1;
    }
    else if((L % 2 == 0) && (l % 2 != 0)) {
      P_per.slices(l * n, (l + 1) * n - 1) = P_rev;
    }
    else if((L % 2 != 0) && (l % 2 == 0)) {
      P_per.slices(l * n, (l + 1) * n - 1) = P_rev;
    }
    else if((L % 2 != 0) && (l % 2 != 0)) {
      P_per.slices(l * n, (l + 1) * n - 1) = P1;
    }
  }

  M_per.slices(0, N - 1) = P_per;
  for(int j = 1; j <= J; ++j) {
      int len = arma::sum(Nj.head(j));
      for(int k = 0; k < Nj(j); ++k) {
        if(metric == "Riemannian") {
          M_per.slice(len + k) = Mid(M_per.slice(len - Nj(j - 1) + 2 * k),
                                   M_per.slice(len - Nj(j - 1) + 2 * k + 1));
        }
        else {
          M_per.slice(len + k) = (M_per.slice(len - Nj(j - 1) + 2 * k) +
                                  M_per.slice(len - Nj(j - 1) + 2 * k + 1)) / 2;
        }
      }
    }
  return M_per;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_cube wavCoeff_C(arma::cx_cube M1, arma::cx_cube M0, double j, std::string metric) {

  // Set parameters
  int d = M1.n_cols;
  int n = M1.n_slices;
  arma::cx_cube W(d, d, 2 * n);
  arma::cx_mat Msqrt(d, d);
  arma::cx_mat Misqrt(d, d);

  for(int k = 0; k < n; ++k) {
    // Compute (non-)whitened wavelet coeff's
    if(metric == "Riemannian") {
      // Riemannian coefficients
      Msqrt = arma::sqrtmat_sympd(M1.slice(k));
      Misqrt = arma::inv_sympd(Msqrt);
      W.slice(k) = std::pow((double)2, (double)(-j / 2)) * arma::logmat_sympd(Misqrt * M0.slice(k) * Misqrt);
      W.slice(k + n) = Msqrt * W.slice(k) * Msqrt;
    }
    else {
      // Euclidean coefficients
      W.slice(k) = std::pow((double)2, (double)(-j / 2)) * (M0.slice(k) - M1.slice(k));
      W.slice(k + n) = W.slice(k);
    }
  }
  return W;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_cube reconstr_C(arma::cx_cube M1, arma::cx_cube M0, arma::cx_cube Dj, double j,
                         int nj, bool in_sample, int L1, std::string metric) {
  // Initialize variables
  int d = M1.n_cols;
  arma::cx_cube M2(d, d, 2 * nj);
  arma::cx_mat Msqrt(d, d);
  arma::cx_mat Misqrt(d, d);

  for(int k = 0; k < nj; ++k) {
    // Reconstruct odd midpoints from non-zero wav. coeffs
    if((arma::norm(Dj.slice(k), "inf") > 1E-10) && in_sample) {
      if(metric == "Riemannian") {
        Msqrt = arma::sqrtmat_sympd(M1.slice(2 * k + 1));
        Misqrt = arma::inv_sympd(Msqrt);
        M2.slice(2 * k + 1) = Msqrt * arma::expmat_sym(std::pow((double)2, (double)(j/2)) *
                              Misqrt * Dj.slice(k) * Misqrt) * Msqrt;
      }
      else {
        M2.slice(2 * k + 1) = std::pow((double)2, (double)(j/2)) * Dj.slice(k) + M1.slice(2 * k + 1);
      }
    }
    else {
      // Reconstruct odd midpoints from zero wav.coeffs
      M2.slice(2 * k + 1) = M1.slice(2 * k + 1);
    }
    // Reconstruct even midpoints from odd midpoints + coarse midpoints
    if(metric == "Riemannian") {
        M2.slice(2 * k) = M0.slice(k + L1) * arma::inv_sympd(M2.slice(2 * k + 1)) * M0.slice(k + L1);
    }
    else {
        M2.slice(2 * k) = 2 * M0.slice(k + L1) - M2.slice(2 * k + 1);
    }
  }
  return M2;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_cube reconstr2D_C(arma::cx_cube M1, arma::cx_cube D, double j,
                           arma::ivec n, bool in_sample, std::string metric) {
  // Initialize variables
  int d = M1.n_cols;
  arma::cx_cube M2(d, d, n(0) * n(1));
  arma::cx_mat Msqrt(d, d);
  arma::cx_mat Misqrt(d, d);

  for(int k1 = 0; k1 < n(0); ++k1) {
    for(int k2 = 0; k2 < n(1); ++k2) {
      // Reconstruct midpoints from non-zero wav. coeffs
      if((arma::norm(D.slice(k2 * n(0) + k1), "inf") > 1E-10) && in_sample) {
        if(metric == "Riemannian") {
        // Riemannian metric
        Msqrt = arma::sqrtmat_sympd(M1.slice(k2 * n(0) + k1));
        Misqrt = arma::inv_sympd(Msqrt);
        M2.slice(k2 * n(0) + k1) = Msqrt * arma::expmat_sym(std::pow((double)2, (double)(j/2)) *
          Misqrt * D.slice(k2 * n(0) + k1) * Misqrt) * Msqrt;
        }
        else {
          // Euclidean metric
          M2.slice(k2 * n(0) + k1) = std::pow((double)2, (double)(j/2)) * D.slice(k2 * n(0) + k1) +
            M1.slice(k2 * n(0) + k1);
        }
      }
      else {
        // Reconstruct midpoints from zero wav. coeffs
        M2.slice(k2 * n(0) + k1) = M1.slice(k2 * n(0) + k1);
      }
    }
  }
  return M2;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

double pdDist_C(arma::cx_mat A, arma::cx_mat B, std::string method) {
  try{
    arma::cx_mat A2;
    // compute distance
    if(method == "Riemannian") {
      arma::cx_mat A1 = arma::inv_sympd(arma::sqrtmat_sympd(A));
      A2 = arma::logmat_sympd(A1 * B * A1);
    }
    else if (method == "logEuclidean") {
      A2 = arma::logmat_sympd(A) - arma::logmat_sympd(B);
    }
    else if (method == "Cholesky") {
      A2 = arma::chol(A) - arma::chol(B);
    }
    else if (method == "rootEuclidean") {
      A2 = arma::sqrtmat_sympd(A) - arma::sqrtmat_sympd(B);
    }
    else {
      A2 = A - B;
    }
    if(A2.has_nan()) {
      Rcpp::stop("c++ function logmat_sympd() failed, matrix possibly not positive definite");
    }
    return arma::norm(A2, "fro");
    // catch exceptions
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::exception("c++ exception (unknown reason)");
  }
}
