// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ARMA
arma::mat ARMA(arma::cube Phi, arma::cube Theta, arma::mat Z, int len);
RcppExport SEXP pdSpecEst_ARMA(SEXP PhiSEXP, SEXP ThetaSEXP, SEXP ZSEXP, SEXP lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    rcpp_result_gen = Rcpp::wrap(ARMA(Phi, Theta, Z, len));
    return rcpp_result_gen;
END_RCPP
}
// E_coeff
arma::vec E_coeff(arma::cx_mat H, arma::cx_cube E);
RcppExport SEXP pdSpecEst_E_coeff(SEXP HSEXP, SEXP ESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cx_mat >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::cx_cube >::type E(ESEXP);
    rcpp_result_gen = Rcpp::wrap(E_coeff(H, E));
    return rcpp_result_gen;
END_RCPP
}
// E_coeff_inv
arma::cx_mat E_coeff_inv(arma::vec coeff, arma::cx_cube E);
RcppExport SEXP pdSpecEst_E_coeff_inv(SEXP coeffSEXP, SEXP ESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type coeff(coeffSEXP);
    Rcpp::traits::input_parameter< arma::cx_cube >::type E(ESEXP);
    rcpp_result_gen = Rcpp::wrap(E_coeff_inv(coeff, E));
    return rcpp_result_gen;
END_RCPP
}
// Expm
arma::cx_mat Expm(arma::cx_mat P, arma::cx_mat H);
RcppExport SEXP pdSpecEst_Expm(SEXP PSEXP, SEXP HSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cx_mat >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::cx_mat >::type H(HSEXP);
    rcpp_result_gen = Rcpp::wrap(Expm(P, H));
    return rcpp_result_gen;
END_RCPP
}
// iSqrt
arma::cx_mat iSqrt(arma::cx_mat M);
RcppExport SEXP pdSpecEst_iSqrt(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cx_mat >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(iSqrt(M));
    return rcpp_result_gen;
END_RCPP
}
// kMean
arma::cx_mat kMean(arma::cx_mat M, arma::vec mu);
RcppExport SEXP pdSpecEst_kMean(SEXP MSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cx_mat >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(kMean(M, mu));
    return rcpp_result_gen;
END_RCPP
}
// Logm
arma::cx_mat Logm(arma::cx_mat P, arma::cx_mat Q);
RcppExport SEXP pdSpecEst_Logm(SEXP PSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cx_mat >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::cx_mat >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(Logm(P, Q));
    return rcpp_result_gen;
END_RCPP
}
// Mid
arma::cx_mat Mid(arma::cx_mat A, arma::cx_mat B);
RcppExport SEXP pdSpecEst_Mid(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cx_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::cx_mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(Mid(A, B));
    return rcpp_result_gen;
END_RCPP
}
// NormF
double NormF(arma::cx_mat M);
RcppExport SEXP pdSpecEst_NormF(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cx_mat >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(NormF(M));
    return rcpp_result_gen;
END_RCPP
}
// RiemmDist
double RiemmDist(arma::cx_mat A, arma::cx_mat B);
RcppExport SEXP pdSpecEst_RiemmDist(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cx_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::cx_mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(RiemmDist(A, B));
    return rcpp_result_gen;
END_RCPP
}
// solveMid
arma::cx_mat solveMid(arma::cx_mat B, arma::cx_mat C);
RcppExport SEXP pdSpecEst_solveMid(SEXP BSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cx_mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::cx_mat >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(solveMid(B, C));
    return rcpp_result_gen;
END_RCPP
}
// Sqrt
arma::cx_mat Sqrt(arma::cx_mat M);
RcppExport SEXP pdSpecEst_Sqrt(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cx_mat >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(Sqrt(M));
    return rcpp_result_gen;
END_RCPP
}
