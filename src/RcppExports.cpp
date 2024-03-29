// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// CTCRWNLL
Rcpp::List CTCRWNLL(const arma::mat& y, const arma::mat& Hmat, const arma::vec& beta, const arma::vec& sig2, const arma::vec& delta, const arma::vec& noObs, const arma::vec& active, const arma::colvec& a, const arma::mat& P);
RcppExport SEXP _crawl_CTCRWNLL(SEXP ySEXP, SEXP HmatSEXP, SEXP betaSEXP, SEXP sig2SEXP, SEXP deltaSEXP, SEXP noObsSEXP, SEXP activeSEXP, SEXP aSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Hmat(HmatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type noObs(noObsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type active(activeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(CTCRWNLL(y, Hmat, beta, sig2, delta, noObs, active, a, P));
    return rcpp_result_gen;
END_RCPP
}
// CTCRWNLL_DRIFT
Rcpp::List CTCRWNLL_DRIFT(const arma::mat& y, const arma::mat& Hmat, const arma::vec& beta, const arma::vec& beta_drift, const arma::vec& sig2, const arma::vec& sig2_drift, const arma::vec& delta, const arma::vec& noObs, const arma::vec& active, const arma::colvec& a, const arma::mat& P);
RcppExport SEXP _crawl_CTCRWNLL_DRIFT(SEXP ySEXP, SEXP HmatSEXP, SEXP betaSEXP, SEXP beta_driftSEXP, SEXP sig2SEXP, SEXP sig2_driftSEXP, SEXP deltaSEXP, SEXP noObsSEXP, SEXP activeSEXP, SEXP aSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Hmat(HmatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta_drift(beta_driftSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sig2_drift(sig2_driftSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type noObs(noObsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type active(activeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(CTCRWNLL_DRIFT(y, Hmat, beta, beta_drift, sig2, sig2_drift, delta, noObs, active, a, P));
    return rcpp_result_gen;
END_RCPP
}
// CTCRWPREDICT
Rcpp::List CTCRWPREDICT(const arma::mat& y, const arma::mat& Hmat, const arma::vec& beta, const arma::vec& sig2, const arma::vec& delta, const arma::vec& noObs, const arma::vec& active, const arma::colvec& a, const arma::mat& P);
RcppExport SEXP _crawl_CTCRWPREDICT(SEXP ySEXP, SEXP HmatSEXP, SEXP betaSEXP, SEXP sig2SEXP, SEXP deltaSEXP, SEXP noObsSEXP, SEXP activeSEXP, SEXP aSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Hmat(HmatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type noObs(noObsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type active(activeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(CTCRWPREDICT(y, Hmat, beta, sig2, delta, noObs, active, a, P));
    return rcpp_result_gen;
END_RCPP
}
// CTCRWPREDICT_DRIFT
Rcpp::List CTCRWPREDICT_DRIFT(const arma::mat& y, const arma::mat& Hmat, const arma::vec& beta, const arma::vec& beta_drift, const arma::vec& sig2, const arma::vec& sig2_drift, const arma::vec& delta, const arma::vec& noObs, const arma::vec& active, const arma::colvec& a, const arma::mat& P);
RcppExport SEXP _crawl_CTCRWPREDICT_DRIFT(SEXP ySEXP, SEXP HmatSEXP, SEXP betaSEXP, SEXP beta_driftSEXP, SEXP sig2SEXP, SEXP sig2_driftSEXP, SEXP deltaSEXP, SEXP noObsSEXP, SEXP activeSEXP, SEXP aSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Hmat(HmatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta_drift(beta_driftSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sig2_drift(sig2_driftSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type noObs(noObsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type active(activeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(CTCRWPREDICT_DRIFT(y, Hmat, beta, beta_drift, sig2, sig2_drift, delta, noObs, active, a, P));
    return rcpp_result_gen;
END_RCPP
}
// CTCRWSAMPLE
Rcpp::List CTCRWSAMPLE(const arma::mat& y, const arma::mat& Hmat, const arma::vec& beta, const arma::vec& sig2, const arma::vec& delta, const arma::vec& noObs, const arma::vec& active, const arma::colvec& a, const arma::mat& P);
RcppExport SEXP _crawl_CTCRWSAMPLE(SEXP ySEXP, SEXP HmatSEXP, SEXP betaSEXP, SEXP sig2SEXP, SEXP deltaSEXP, SEXP noObsSEXP, SEXP activeSEXP, SEXP aSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Hmat(HmatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type noObs(noObsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type active(activeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(CTCRWSAMPLE(y, Hmat, beta, sig2, delta, noObs, active, a, P));
    return rcpp_result_gen;
END_RCPP
}
// CTCRWSAMPLE_DRIFT
Rcpp::List CTCRWSAMPLE_DRIFT(const arma::mat& y, const arma::mat& Hmat, const arma::vec& beta, const arma::vec& beta_drift, const arma::vec& sig2, const arma::vec& sig2_drift, const arma::vec& delta, const arma::vec& noObs, const arma::vec& active, const arma::colvec& a, const arma::mat& P);
RcppExport SEXP _crawl_CTCRWSAMPLE_DRIFT(SEXP ySEXP, SEXP HmatSEXP, SEXP betaSEXP, SEXP beta_driftSEXP, SEXP sig2SEXP, SEXP sig2_driftSEXP, SEXP deltaSEXP, SEXP noObsSEXP, SEXP activeSEXP, SEXP aSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Hmat(HmatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta_drift(beta_driftSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sig2_drift(sig2_driftSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type noObs(noObsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type active(activeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(CTCRWSAMPLE_DRIFT(y, Hmat, beta, beta_drift, sig2, sig2_drift, delta, noObs, active, a, P));
    return rcpp_result_gen;
END_RCPP
}
// makeT
arma::mat makeT(const double& b, const double& delta, const double& active);
RcppExport SEXP _crawl_makeT(SEXP bSEXP, SEXP deltaSEXP, SEXP activeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const double& >::type active(activeSEXP);
    rcpp_result_gen = Rcpp::wrap(makeT(b, delta, active));
    return rcpp_result_gen;
END_RCPP
}
// makeQ
arma::mat makeQ(const double& b, const double& sig2, const double& delta, const double& active);
RcppExport SEXP _crawl_makeQ(SEXP bSEXP, SEXP sig2SEXP, SEXP deltaSEXP, SEXP activeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< const double& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const double& >::type active(activeSEXP);
    rcpp_result_gen = Rcpp::wrap(makeQ(b, sig2, delta, active));
    return rcpp_result_gen;
END_RCPP
}
// makeT_drift
arma::mat makeT_drift(const double& b, const double& b_drift, const double& delta, const double& active);
RcppExport SEXP _crawl_makeT_drift(SEXP bSEXP, SEXP b_driftSEXP, SEXP deltaSEXP, SEXP activeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type b_drift(b_driftSEXP);
    Rcpp::traits::input_parameter< const double& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const double& >::type active(activeSEXP);
    rcpp_result_gen = Rcpp::wrap(makeT_drift(b, b_drift, delta, active));
    return rcpp_result_gen;
END_RCPP
}
// makeQ_drift
arma::mat makeQ_drift(const double& b, const double& b_drift, const double& sig2, const double& sig2_drift, const double& delta, const double& active);
RcppExport SEXP _crawl_makeQ_drift(SEXP bSEXP, SEXP b_driftSEXP, SEXP sig2SEXP, SEXP sig2_driftSEXP, SEXP deltaSEXP, SEXP activeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type b_drift(b_driftSEXP);
    Rcpp::traits::input_parameter< const double& >::type sig2(sig2SEXP);
    Rcpp::traits::input_parameter< const double& >::type sig2_drift(sig2_driftSEXP);
    Rcpp::traits::input_parameter< const double& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const double& >::type active(activeSEXP);
    rcpp_result_gen = Rcpp::wrap(makeQ_drift(b, b_drift, sig2, sig2_drift, delta, active));
    return rcpp_result_gen;
END_RCPP
}
