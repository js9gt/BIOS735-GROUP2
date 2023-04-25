// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// c_pi
arma::vec c_pi(const arma::mat& X, const arma::vec& gammak);
RcppExport SEXP _ZINB_c_pi(SEXP XSEXP, SEXP gammakSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type gammak(gammakSEXP);
    rcpp_result_gen = Rcpp::wrap(c_pi(X, gammak));
    return rcpp_result_gen;
END_RCPP
}
// c_mu
arma::vec c_mu(const arma::mat& X, const arma::vec& betak);
RcppExport SEXP _ZINB_c_mu(SEXP XSEXP, SEXP betakSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type betak(betakSEXP);
    rcpp_result_gen = Rcpp::wrap(c_mu(X, betak));
    return rcpp_result_gen;
END_RCPP
}
// Qik
arma::vec Qik(const arma::vec& pi_k, const arma::vec& mu_k, const arma::vec& y_k, double thetak);
RcppExport SEXP _ZINB_Qik(SEXP pi_kSEXP, SEXP mu_kSEXP, SEXP y_kSEXP, SEXP thetakSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type pi_k(pi_kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_k(mu_kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y_k(y_kSEXP);
    Rcpp::traits::input_parameter< double >::type thetak(thetakSEXP);
    rcpp_result_gen = Rcpp::wrap(Qik(pi_k, mu_k, y_k, thetak));
    return rcpp_result_gen;
END_RCPP
}
// score
arma::vec score(const arma::mat& X, const arma::vec& y_k, const arma::vec& z_k, const arma::vec& mu_k, const arma::vec& pi_k, double thetak);
RcppExport SEXP _ZINB_score(SEXP XSEXP, SEXP y_kSEXP, SEXP z_kSEXP, SEXP mu_kSEXP, SEXP pi_kSEXP, SEXP thetakSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y_k(y_kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type z_k(z_kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_k(mu_kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pi_k(pi_kSEXP);
    Rcpp::traits::input_parameter< double >::type thetak(thetakSEXP);
    rcpp_result_gen = Rcpp::wrap(score(X, y_k, z_k, mu_k, pi_k, thetak));
    return rcpp_result_gen;
END_RCPP
}
// In
arma::mat In(const arma::mat& X, const arma::vec& y_k, const arma::vec& z_k, const arma::vec& mu_k, const arma::vec& pi_k, const arma::vec& gammak, double thetak);
RcppExport SEXP _ZINB_In(SEXP XSEXP, SEXP y_kSEXP, SEXP z_kSEXP, SEXP mu_kSEXP, SEXP pi_kSEXP, SEXP gammakSEXP, SEXP thetakSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y_k(y_kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type z_k(z_kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_k(mu_kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pi_k(pi_kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type gammak(gammakSEXP);
    Rcpp::traits::input_parameter< double >::type thetak(thetakSEXP);
    rcpp_result_gen = Rcpp::wrap(In(X, y_k, z_k, mu_k, pi_k, gammak, thetak));
    return rcpp_result_gen;
END_RCPP
}
// M_estimator
double M_estimator(double thetak, const arma::vec& yik, const arma::vec& muik, const arma::vec& zik, int p);
RcppExport SEXP _ZINB_M_estimator(SEXP thetakSEXP, SEXP yikSEXP, SEXP muikSEXP, SEXP zikSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type thetak(thetakSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type yik(yikSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type muik(muikSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type zik(zikSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(M_estimator(thetak, yik, muik, zik, p));
    return rcpp_result_gen;
END_RCPP
}
// M_deriv
double M_deriv(double thetak, const arma::vec& yik, const arma::vec& muik, const arma::vec& zik);
RcppExport SEXP _ZINB_M_deriv(SEXP thetakSEXP, SEXP yikSEXP, SEXP muikSEXP, SEXP zikSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type thetak(thetakSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type yik(yikSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type muik(muikSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type zik(zikSEXP);
    rcpp_result_gen = Rcpp::wrap(M_deriv(thetak, yik, muik, zik));
    return rcpp_result_gen;
END_RCPP
}
// theta_update
double theta_update(double thetak, int numIter, const arma::vec& yik, const arma::vec& muik, const arma::vec& zik, int p);
RcppExport SEXP _ZINB_theta_update(SEXP thetakSEXP, SEXP numIterSEXP, SEXP yikSEXP, SEXP muikSEXP, SEXP zikSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type thetak(thetakSEXP);
    Rcpp::traits::input_parameter< int >::type numIter(numIterSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type yik(yikSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type muik(muikSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type zik(zikSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(theta_update(thetak, numIter, yik, muik, zik, p));
    return rcpp_result_gen;
END_RCPP
}
// dznbinom
arma::vec dznbinom(const arma::vec& x, double size, const arma::vec& prob, arma::vec& infl);
RcppExport SEXP _ZINB_dznbinom(SEXP xSEXP, SEXP sizeSEXP, SEXP probSEXP, SEXP inflSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type infl(inflSEXP);
    rcpp_result_gen = Rcpp::wrap(dznbinom(x, size, prob, infl));
    return rcpp_result_gen;
END_RCPP
}
// par_EM
List par_EM(const arma::mat& X, const arma::vec& y_k, double tol, int maxIter, int inflation, Nullable<List> initial);
RcppExport SEXP _ZINB_par_EM(SEXP XSEXP, SEXP y_kSEXP, SEXP tolSEXP, SEXP maxIterSEXP, SEXP inflationSEXP, SEXP initialSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y_k(y_kSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< int >::type inflation(inflationSEXP);
    Rcpp::traits::input_parameter< Nullable<List> >::type initial(initialSEXP);
    rcpp_result_gen = Rcpp::wrap(par_EM(X, y_k, tol, maxIter, inflation, initial));
    return rcpp_result_gen;
END_RCPP
}
// LRT1G
List LRT1G(const arma::mat& X, const arma::vec& y_k, double tol, int maxIter, int inflation, Nullable<List> initial);
RcppExport SEXP _ZINB_LRT1G(SEXP XSEXP, SEXP y_kSEXP, SEXP tolSEXP, SEXP maxIterSEXP, SEXP inflationSEXP, SEXP initialSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y_k(y_kSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< int >::type inflation(inflationSEXP);
    Rcpp::traits::input_parameter< Nullable<List> >::type initial(initialSEXP);
    rcpp_result_gen = Rcpp::wrap(LRT1G(X, y_k, tol, maxIter, inflation, initial));
    return rcpp_result_gen;
END_RCPP
}
// LRTnG
List LRTnG(const arma::mat& X, const arma::mat& Y, const arma::vec& inflation, double tol, int maxIter, Nullable<List> initial);
RcppExport SEXP _ZINB_LRTnG(SEXP XSEXP, SEXP YSEXP, SEXP inflationSEXP, SEXP tolSEXP, SEXP maxIterSEXP, SEXP initialSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type inflation(inflationSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< Nullable<List> >::type initial(initialSEXP);
    rcpp_result_gen = Rcpp::wrap(LRTnG(X, Y, inflation, tol, maxIter, initial));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ZINB_c_pi", (DL_FUNC) &_ZINB_c_pi, 2},
    {"_ZINB_c_mu", (DL_FUNC) &_ZINB_c_mu, 2},
    {"_ZINB_Qik", (DL_FUNC) &_ZINB_Qik, 4},
    {"_ZINB_score", (DL_FUNC) &_ZINB_score, 6},
    {"_ZINB_In", (DL_FUNC) &_ZINB_In, 7},
    {"_ZINB_M_estimator", (DL_FUNC) &_ZINB_M_estimator, 5},
    {"_ZINB_M_deriv", (DL_FUNC) &_ZINB_M_deriv, 4},
    {"_ZINB_theta_update", (DL_FUNC) &_ZINB_theta_update, 6},
    {"_ZINB_dznbinom", (DL_FUNC) &_ZINB_dznbinom, 4},
    {"_ZINB_par_EM", (DL_FUNC) &_ZINB_par_EM, 6},
    {"_ZINB_LRT1G", (DL_FUNC) &_ZINB_LRT1G, 6},
    {"_ZINB_LRTnG", (DL_FUNC) &_ZINB_LRTnG, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_ZINB(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
