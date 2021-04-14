// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// get_lambda
arma::cube get_lambda(double lam0, arma::cube distSq, double kappa, arma::icube farOut, arma::imat notOper, int nthreads);
RcppExport SEXP _scrmove_get_lambda(SEXP lam0SEXP, SEXP distSqSEXP, SEXP kappaSEXP, SEXP farOutSEXP, SEXP notOperSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lam0(lam0SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type distSq(distSqSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::icube >::type farOut(farOutSEXP);
    Rcpp::traits::input_parameter< arma::imat >::type notOper(notOperSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_lambda(lam0, distSq, kappa, farOut, notOper, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// get_lambda_nou
arma::cube get_lambda_nou(double lam0, arma::mat distSq, double sigma, arma::imat farOut, arma::imat notOper, int nthreads);
RcppExport SEXP _scrmove_get_lambda_nou(SEXP lam0SEXP, SEXP distSqSEXP, SEXP sigmaSEXP, SEXP farOutSEXP, SEXP notOperSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lam0(lam0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type distSq(distSqSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::imat >::type farOut(farOutSEXP);
    Rcpp::traits::input_parameter< arma::imat >::type notOper(notOperSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_lambda_nou(lam0, distSq, sigma, farOut, notOper, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// get_lambda_noK
arma::mat get_lambda_noK(double lam0, arma::mat distSq, double kappa, arma::mat farOut, arma::ivec trapOper, int nthreads);
RcppExport SEXP _scrmove_get_lambda_noK(SEXP lam0SEXP, SEXP distSqSEXP, SEXP kappaSEXP, SEXP farOutSEXP, SEXP trapOperSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lam0(lam0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type distSq(distSqSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type farOut(farOutSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type trapOper(trapOperSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_lambda_noK(lam0, distSq, kappa, farOut, trapOper, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// get_ld_y_n
arma::vec get_ld_y_n(arma::icube y, arma::cube lambda, int nthreads);
RcppExport SEXP _scrmove_get_ld_y_n(SEXP ySEXP, SEXP lambdaSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::icube >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::cube >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_ld_y_n(y, lambda, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// get_ld_y_nK
arma::mat get_ld_y_nK(arma::icube y, arma::cube lambda, int nthreads);
RcppExport SEXP _scrmove_get_ld_y_nK(SEXP ySEXP, SEXP lambdaSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::icube >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::cube >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_ld_y_nK(y, lambda, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// get_ld_y
arma::vec get_ld_y(arma::icube y, arma::cube lambda, arma::ivec z, int M);
RcppExport SEXP _scrmove_get_ld_y(SEXP ySEXP, SEXP lambdaSEXP, SEXP zSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::icube >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::cube >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(get_ld_y(y, lambda, z, M));
    return rcpp_result_gen;
END_RCPP
}
// get_ld_y_noK
arma::vec get_ld_y_noK(arma::imat y, arma::mat lambda, int nthreads);
RcppExport SEXP _scrmove_get_ld_y_noK(SEXP ySEXP, SEXP lambdaSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::imat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_ld_y_noK(y, lambda, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// get_distSq_ux
arma::cube get_distSq_ux(arma::cube u, arma::mat x, int nthreads);
RcppExport SEXP _scrmove_get_distSq_ux(SEXP uSEXP, SEXP xSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_distSq_ux(u, x, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// get_qstar_cam
double get_qstar_cam(int nsim, arma::mat x, int K, double sigma, double rho, double lam0, double kappa, double buffer, double maxDistSq, arma::imat notOper, int nthreads);
RcppExport SEXP _scrmove_get_qstar_cam(SEXP nsimSEXP, SEXP xSEXP, SEXP KSEXP, SEXP sigmaSEXP, SEXP rhoSEXP, SEXP lam0SEXP, SEXP kappaSEXP, SEXP bufferSEXP, SEXP maxDistSqSEXP, SEXP notOperSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type lam0(lam0SEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type buffer(bufferSEXP);
    Rcpp::traits::input_parameter< double >::type maxDistSq(maxDistSqSEXP);
    Rcpp::traits::input_parameter< arma::imat >::type notOper(notOperSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_qstar_cam(nsim, x, K, sigma, rho, lam0, kappa, buffer, maxDistSq, notOper, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// get_ld_ydet_aug
double get_ld_ydet_aug(int nsim, double s1, double s2, arma::mat x, int K, double sigma, double rho, double lam0, double kappa, double maxDistSq, arma::imat notOper, int nthreads);
RcppExport SEXP _scrmove_get_ld_ydet_aug(SEXP nsimSEXP, SEXP s1SEXP, SEXP s2SEXP, SEXP xSEXP, SEXP KSEXP, SEXP sigmaSEXP, SEXP rhoSEXP, SEXP lam0SEXP, SEXP kappaSEXP, SEXP maxDistSqSEXP, SEXP notOperSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< double >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type lam0(lam0SEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type maxDistSq(maxDistSqSEXP);
    Rcpp::traits::input_parameter< arma::imat >::type notOper(notOperSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_ld_ydet_aug(nsim, s1, s2, x, K, sigma, rho, lam0, kappa, maxDistSq, notOper, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// get_qstar_cam_nomove
double get_qstar_cam_nomove(int nsim, arma::mat x, int K, double lam0, double kappa, double buffer, double maxDistSq, arma::ivec trapOper, int nthreads);
RcppExport SEXP _scrmove_get_qstar_cam_nomove(SEXP nsimSEXP, SEXP xSEXP, SEXP KSEXP, SEXP lam0SEXP, SEXP kappaSEXP, SEXP bufferSEXP, SEXP maxDistSqSEXP, SEXP trapOperSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type lam0(lam0SEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type buffer(bufferSEXP);
    Rcpp::traits::input_parameter< double >::type maxDistSq(maxDistSqSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type trapOper(trapOperSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_qstar_cam_nomove(nsim, x, K, lam0, kappa, buffer, maxDistSq, trapOper, nthreads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scrmove_get_lambda", (DL_FUNC) &_scrmove_get_lambda, 6},
    {"_scrmove_get_lambda_nou", (DL_FUNC) &_scrmove_get_lambda_nou, 6},
    {"_scrmove_get_lambda_noK", (DL_FUNC) &_scrmove_get_lambda_noK, 6},
    {"_scrmove_get_ld_y_n", (DL_FUNC) &_scrmove_get_ld_y_n, 3},
    {"_scrmove_get_ld_y_nK", (DL_FUNC) &_scrmove_get_ld_y_nK, 3},
    {"_scrmove_get_ld_y", (DL_FUNC) &_scrmove_get_ld_y, 4},
    {"_scrmove_get_ld_y_noK", (DL_FUNC) &_scrmove_get_ld_y_noK, 3},
    {"_scrmove_get_distSq_ux", (DL_FUNC) &_scrmove_get_distSq_ux, 3},
    {"_scrmove_get_qstar_cam", (DL_FUNC) &_scrmove_get_qstar_cam, 11},
    {"_scrmove_get_ld_ydet_aug", (DL_FUNC) &_scrmove_get_ld_ydet_aug, 12},
    {"_scrmove_get_qstar_cam_nomove", (DL_FUNC) &_scrmove_get_qstar_cam_nomove, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_scrmove(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
