# define ARMA_DONT_PRINT_ERRORS
# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat Mid_w(arma::cx_mat A, arma::cx_mat B, double w) {
  arma::cx_mat B1 = arma::sqrtmat_sympd(B);
  arma::cx_mat B2 = arma::inv_sympd(B1);
  arma::vec eigval;
  arma::cx_mat eigvec;
  arma::eig_sym(eigval, eigvec, B2 * A * B2);
  return B1 * (eigvec * arma::diagmat(arma::pow(eigval, w)) * eigvec.t()) * B1;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat pdMean_C_approx(arma::cx_cube M, arma::vec mu) {
  try{
    // initialize params
    int d = M.n_cols;
    int n = mu.size();
    arma::cx_mat M1 = M.slice(0);
    arma::cx_mat Mi(d, d);
    arma::cx_mat M1sq(d, d);
    arma::cx_mat M1isq(d, d);
    arma::vec eigval;
    arma::cx_mat eigvec;
    double mu1;
    // compute approx mean
    for(int i = 1; i < n; ++i) {
      if (i % 100 == 0) {
        Rcpp::checkUserInterrupt();
      }
      Mi = M.slice(i);
      M1sq = arma::sqrtmat_sympd(M1);
      M1isq = arma::inv_sympd(M1sq);
      mu1 = arma::sum(mu.head(i + 1));
      arma::eig_sym(eigval, eigvec, M1isq * Mi * M1isq);
      arma::cx_mat M11 = M1sq * eigvec;
      if(mu1 == 0){
        M1 = M11 * arma::diagmat(arma::pow(eigval, mu[i])) * M11.t();
      } else {
        M1 = M11 * arma::diagmat(arma::pow(eigval, mu[i] / mu1)) * M11.t();
      }
    }
    return M1;
    // catch exceptions
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::exception("c++ exception (unknown reason)");
  }
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat pdMean_C(arma::cx_mat M0, arma::cx_cube M, arma::vec mu,
                      int maxit, double reltol) {
  try{
    // initialize params
    int d = M0.n_cols;
    int n = mu.size();
    arma::cx_mat M1 = M0;
    arma::cx_mat M0sqrt;
    arma::cx_mat M0isqrt;
    arma::cx_mat M1log;
    double reltol_m = reltol + 1;
    int i = 0;
    // run gradient descent
    while(reltol_m > reltol && i < maxit) {
      if (i % 100 == 0) {
        Rcpp::checkUserInterrupt();
      }
      M0 = M1;
      M0sqrt = arma::sqrtmat_sympd(M0);
      M0isqrt = arma::inv_sympd(M0sqrt);
      M1log = arma::zeros<arma::cx_mat>(d, d);
      for(int j = 0; j < n; ++j) {
        M1log += mu(j) * arma::logmat_sympd(M0isqrt * M.slice(j) * M0isqrt);
        if(M1log.has_nan()) {
          Rcpp::stop("c++ matrix functions failed, matrix possibly not positive definite");
        }
      }
      M1 = M0sqrt * arma::expmat_sym(M1log) * M0sqrt;
      reltol_m = arma::norm(arma::logmat_sympd(M0isqrt * M1 * M0isqrt), "fro");
      ++ i;
    }
    return M1;
    // catch errors
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::exception("c++ exception (unknown reason)");
  }
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat pdMedian_C(arma::cx_mat M0, arma::cx_cube M, arma::vec mu,
                        int maxit, double reltol) {
  try{
    // initialize params
    int d = M0.n_cols;
    int n = mu.size();
    arma::cx_mat M1 = M0;
    arma::cx_mat M0sqrt;
    arma::cx_mat M0isqrt;
    arma::cx_mat M1log;
    arma::cx_mat Mlog;
    double dist;
    double W;
    double reltol_m = reltol + 1;
    int i = 0;
    // run Weiszfeld algorithm (gradient descent)
    while(reltol_m > reltol && i < maxit) {
      if (i % 100 == 0) {
        Rcpp::checkUserInterrupt();
      }
      M0 = M1;
      M0sqrt = arma::sqrtmat_sympd(M0);
      M0isqrt = arma::inv_sympd(M0sqrt);
      M1log = arma::zeros<arma::cx_mat>(d, d);
      W = 0;
      for(int j = 0; j < n; ++j) {
        Mlog = arma::logmat_sympd(M0isqrt * M.slice(j) * M0isqrt);
        if(Mlog.has_nan()) {
          Rcpp::stop("c++ matrix functions failed, matrix possibly not positive definite");
        }
        dist = arma::norm(Mlog, "fro");
        if(dist > 1E-10) {
          M1log += mu(j) / dist * Mlog;
          W += mu(j) / dist;
        }
      }
      M1 = M0sqrt * arma::expmat_sym(M1log / W) * M0sqrt;
      reltol_m = arma::norm(arma::logmat_sympd(M0isqrt * M1 * M0isqrt), "fro");
      ++ i;
    }
    return M1;
    // catch errors
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::exception("c++ exception (unknown reason)");
  }
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat Euclid_Median_C(arma::cx_mat M0, arma::cx_cube M, arma::vec mu,
                             int maxit, double reltol) {
  try{
    // initialize params
    int d = M0.n_cols;
    int n = mu.size();
    arma::cx_mat M1 = M0;
    double dist;
    double W;
    double reltol_m = reltol + 1;
    int i = 0;
    // run Weiszfeld algorithm (gradient descent)
    while(reltol_m > reltol && i < maxit) {
      if (i % 100 == 0) {
        Rcpp::checkUserInterrupt();
      }
      M0 = M1;
      M1 = arma::zeros<arma::cx_mat>(d, d);
      W = 0;
      for(int j = 0; j < n; ++j) {
        dist = arma::norm(M.slice(j) - M0, "fro");
        if(dist > 1E-10) {
          M1 += mu(j) / dist * M.slice(j);
          W += mu(j) / dist;
        }
      }
      M1 = M1 / W;
      reltol_m = arma::norm(M1 - M0, "fro");
      ++ i;
    }
    return M1;
    // catch errors
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::exception("c++ exception (unknown reason)");
  }
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_cube pdNeville_C(arma::cx_cube P, arma::vec X, arma::vec x, std::string method) {

  int n_P = X.size() - 1;
  int n_x = x.size();
  int d = P.n_cols;
  arma::cx_cube res(d, d, n_x);
  for(int j = 0; j < n_x; ++j) {
    arma::cx_cube p = P;
    for(int k = 0; k < n_P; ++k) {
      for(int i = 0; i < (n_P - k); ++i) {
        if(method == "Riemannian"){
          p.slice(i) = Mid_w(p.slice(i + 1), p.slice(i), (x(j) - X(i)) / (X(i + k + 1) - X(i)));
        }
        else {
          p.slice(i) = (X(i + k + 1) - x(j)) / (X(i + k + 1) - X(i)) * p.slice(i) +
            (x(j) - X(i)) / (X(i + k + 1) - X(i)) * p.slice(i + 1);
        }
      }
    }
    res.slice(j) = p.slice(0);
  }
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_cube impute_C(arma::cx_cube M0, arma::mat W, int L, bool inverse,
                       std::string metric, std::string method) {
  // Set parameters
  int n = M0.n_slices;
  int d = M0.n_rows;
  int N = 2 * L + 1;
  arma::vec w = arma::ones<arma::vec>(N);
  arma::cx_cube M(d, d, N);
  arma::cx_cube M1 = arma::zeros<arma::cx_cube>(d, d, 2 * n);

  if(L == 0 || n == 1) {
    // Trivial refinement
    for(int k = 0; k < n; ++k) {
      M1.slice(2 * k + 1) = M0.slice(k);
      if(inverse) {
        M1.slice(2 * k) = M0.slice(k);
      }
    }
  }
  else if (L <= 4){
    // Nontrivial refinement with filter weights
    for(int k = 0; k < n; ++k) {
      // Coarse-scale midpoints and filter weights
      if((k - L) < 0) {
        M = M0.slices(0, N - 1);
        w = arma::conv_to<arma::vec>::from(W.row(2 * k + 1));
      }
      else if((k + L) > (n - 1)) {
        M = M0.slices(n - N, n - 1);
        w = arma::conv_to<arma::vec>::from(W.row(2 * (N - (n - k)) + 1));
      }
      else {
        M = M0.slices(k - L, k + L);
        w = arma::conv_to<arma::vec>::from(W.row(N));
      }
      // Compute predicted midpoints
      if(metric == "Riemannian") {
        // Riemannian weighted mean (odd locations)
        M1.slice(2 * k + 1) = pdMean_C_approx(M, w);
        if(method != "fast") {
          M1.slice(2 * k + 1) = pdMean_C(M1.slice(2 * k + 1), M, w, 1000, 1E-10);
        }
        // Prediction even locations
        if(inverse) {
          M1.slice(2 * k) = M0.slice(k) * arma::inv_sympd(M1.slice(2 * k + 1)) * M0.slice(k);
        }
      } else {
        // Euclidean weighted mean (odd locations)
        arma::cx_mat M_w = arma::zeros<arma::cx_mat>(d, d);
        for(int i = 0; i < N; ++i){
          M_w += w(i) * M.slice(i);
        }
        M1.slice(2 * k + 1) = M_w;
        // Prediction even locations
        if(inverse) {
          M1.slice(2 * k) = 2 * M0.slice(k) - M1.slice(2 * k + 1);
        }
      }
    }
  }
  else if(L > 4) {
    // Nontrivial refinement with Neville's algorithm
    arma::cx_cube M_bar(d, d, N);
    arma::cx_cube M_neville(d, d, 1);
    arma::vec N_seq = arma::linspace(1, N, N);

    for(int k = 0; k < n; ++k) {
      // Coarse-scale midpoints
      if((k - L) < 0) {
        M = M0.slices(0, N - 1);
      }
      else if((k + L) > (n - 1)) {
        M = M0.slices(n - N, n - 1);
      }
      else {
        M = M0.slices(k - L, k + L);
      }
      // Compute predicted midpoints (Neville's algorithm)
      arma::cx_mat M_w = arma::zeros<arma::cx_mat>(d, d);
      for(int i = 0; i < N; ++i) {
        // Cumulative intrinsic means
        if(metric == "Riemannian") {
          M_bar.slice(i) = pdMean_C_approx(M.head_slices(i + 1), w.head(i + 1));
          // if(method != "fast") {
          //   M_bar.slice(i) = pdMean_C(M_bar.slice(i), M.head_slices(i + 1), w.head(i + 1), 1000, 1E-10);
          //   }
        }
        else {
          M_bar.slice(i) = mean(M.head_slices(i + 1), 2);
        }
      }
      if((k - L) < 0) {
        // Prediction at left boundary, uneven locations
        M = M0.slices(0, N - 1);
        arma::vec k_seq(1);
        k_seq(0) = k + 0.5;
        M_neville = pdNeville_C(M_bar, N_seq, k_seq, metric);
        if(metric == "Riemannian") {
          M1.slice(2 * k + 1) = Mid_w(M_neville.slice(0), M_bar.slice(k), -(2 * k + 1));
        }
        else {
          M1.slice(2 * k + 1) = (2 * k + 2) * M_bar.slice(k) - (2 * k + 1) * M_neville.slice(0);
        }
      }
      else if((k + L) > (n - 1)) {
        // Prediction at right boundary, uneven locations
        M = M0.slices(n - N, n - 1);
        int k1 = N - (n - k);
        arma::vec k1_seq(1);
        k1_seq(0) = k1 + 0.5;
        M_neville = pdNeville_C(M_bar, N_seq, k1_seq, metric);
        if(metric == "Riemannian") {
          M1.slice(2 * k + 1) = Mid_w(M_neville.slice(0), M_bar.slice(k1), -(2 * k1 + 1));
        }
        else {
          M1.slice(2 * k + 1) = (2 * k1 + 2) * M_bar.slice(k1) - (2 * k1 + 1) * M_neville.slice(0);
        }
      }
      else {
        // Prediction away from boundary, uneven locations
        M = M0.slices(k - L, k + L);
        arma::vec D_seq(1);
        D_seq(0) = L + 0.5;
        M_neville = pdNeville_C(M_bar, N_seq, D_seq, metric);
        if(metric == "Riemannian") {
          M1.slice(2 * k + 1) = Mid_w(M_neville.slice(0), M_bar.slice(L), -N);
        }
        else {
          M1.slice(2 * k + 1) = (1 + N) * M_bar.slice(L) - N * M_neville.slice(0);
        }
      }
      // Prediction at even locations
      if(inverse) {
        if(metric == "Riemannian") {
          M1.slice(2 * k) = M0.slice(k) * arma::inv_sympd(M1.slice(2 * k + 1)) * M0.slice(k);
        }
        else {
          M1.slice(2 * k) = 2 * M0.slice(k) - M1.slice(2 * k + 1);
        }
      }
    }
  }
  return M1;
}
