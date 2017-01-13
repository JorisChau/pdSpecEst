# include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export()]]

Eigen::MatrixXcd kMean(Eigen::MatrixXcd M, Eigen::VectorXd mu) {

  int d = M.cols();

  int n = M.rows()/d;

  Eigen::MatrixXcd M1 = M.block(0, 0, d, d);

  Eigen::MatrixXcd Mi(d, d);

  Eigen::MatrixXcd M1sq(d, d);

  Eigen::MatrixXcd M1isq(d, d);

  double mu1;

  for(int i = 1; i < n; ++i) {

    Mi = M.block(i * d, 0, d, d);

    M1sq = M1.sqrt();

    M1isq = M1sq.inverse();

    mu1 = mu.head(i + 1).sum();

    if(mu1 == 0){

      M1 = M1sq * (M1isq * Mi * M1isq).pow(mu[i]) * M1sq;

    } else {

      M1 = M1sq * (M1isq * Mi * M1isq).pow(mu[i] / mu1) * M1sq;

    }

  }

  return M1;
}
