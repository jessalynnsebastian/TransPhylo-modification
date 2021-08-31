// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// probSubtree
double probSubtree(NumericMatrix tab, double rate);
RcppExport SEXP _TransPhylo_probSubtree(SEXP tabSEXP, SEXP rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type tab(tabSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    rcpp_result_gen = Rcpp::wrap(probSubtree(tab, rate));
    return rcpp_result_gen;
END_RCPP
}
// probPTreeGivenTTree
double probPTreeGivenTTree(NumericMatrix ctree, double neg, IntegerVector w);
RcppExport SEXP _TransPhylo_probPTreeGivenTTree(SEXP ctreeSEXP, SEXP negSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ctree(ctreeSEXP);
    Rcpp::traits::input_parameter< double >::type neg(negSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(probPTreeGivenTTree(ctree, neg, w));
    return rcpp_result_gen;
END_RCPP
}
// coalescent
double coalescent(NumericVector leaves, NumericVector nodes, double alpha);
RcppExport SEXP _TransPhylo_coalescent(SEXP leavesSEXP, SEXP nodesSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type leaves(leavesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(coalescent(leaves, nodes, alpha));
    return rcpp_result_gen;
END_RCPP
}
// log_sum_exp
double log_sum_exp(double u, double v);
RcppExport SEXP _TransPhylo_log_sum_exp(SEXP uSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(log_sum_exp(u, v));
    return rcpp_result_gen;
END_RCPP
}
// log_subtract_exp
double log_subtract_exp(double u, double v);
RcppExport SEXP _TransPhylo_log_subtract_exp(SEXP uSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(log_subtract_exp(u, v));
    return rcpp_result_gen;
END_RCPP
}
// log_sum_exp_vec
double log_sum_exp_vec(NumericVector w);
RcppExport SEXP _TransPhylo_log_sum_exp_vec(SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(log_sum_exp_vec(w));
    return rcpp_result_gen;
END_RCPP
}
// wbar
NumericVector wbar(double tinf, double dateT, double rOff, double pOff, double pi, double shGen, double scGen, double shSam, double scSam, double delta_t);
RcppExport SEXP _TransPhylo_wbar(SEXP tinfSEXP, SEXP dateTSEXP, SEXP rOffSEXP, SEXP pOffSEXP, SEXP piSEXP, SEXP shGenSEXP, SEXP scGenSEXP, SEXP shSamSEXP, SEXP scSamSEXP, SEXP delta_tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type tinf(tinfSEXP);
    Rcpp::traits::input_parameter< double >::type dateT(dateTSEXP);
    Rcpp::traits::input_parameter< double >::type rOff(rOffSEXP);
    Rcpp::traits::input_parameter< double >::type pOff(pOffSEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    Rcpp::traits::input_parameter< double >::type shGen(shGenSEXP);
    Rcpp::traits::input_parameter< double >::type scGen(scGenSEXP);
    Rcpp::traits::input_parameter< double >::type shSam(shSamSEXP);
    Rcpp::traits::input_parameter< double >::type scSam(scSamSEXP);
    Rcpp::traits::input_parameter< double >::type delta_t(delta_tSEXP);
    rcpp_result_gen = Rcpp::wrap(wbar(tinf, dateT, rOff, pOff, pi, shGen, scGen, shSam, scSam, delta_t));
    return rcpp_result_gen;
END_RCPP
}
// probTTree
double probTTree(NumericMatrix ttree, double rOff, double pOff, double pi, double shGen, double scGen, double shSam, double scSam, double dateT, double delta_t);
RcppExport SEXP _TransPhylo_probTTree(SEXP ttreeSEXP, SEXP rOffSEXP, SEXP pOffSEXP, SEXP piSEXP, SEXP shGenSEXP, SEXP scGenSEXP, SEXP shSamSEXP, SEXP scSamSEXP, SEXP dateTSEXP, SEXP delta_tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ttree(ttreeSEXP);
    Rcpp::traits::input_parameter< double >::type rOff(rOffSEXP);
    Rcpp::traits::input_parameter< double >::type pOff(pOffSEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    Rcpp::traits::input_parameter< double >::type shGen(shGenSEXP);
    Rcpp::traits::input_parameter< double >::type scGen(scGenSEXP);
    Rcpp::traits::input_parameter< double >::type shSam(shSamSEXP);
    Rcpp::traits::input_parameter< double >::type scSam(scSamSEXP);
    Rcpp::traits::input_parameter< double >::type dateT(dateTSEXP);
    Rcpp::traits::input_parameter< double >::type delta_t(delta_tSEXP);
    rcpp_result_gen = Rcpp::wrap(probTTree(ttree, rOff, pOff, pi, shGen, scGen, shSam, scSam, dateT, delta_t));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TransPhylo_probSubtree", (DL_FUNC) &_TransPhylo_probSubtree, 2},
    {"_TransPhylo_probPTreeGivenTTree", (DL_FUNC) &_TransPhylo_probPTreeGivenTTree, 3},
    {"_TransPhylo_coalescent", (DL_FUNC) &_TransPhylo_coalescent, 3},
    {"_TransPhylo_log_sum_exp", (DL_FUNC) &_TransPhylo_log_sum_exp, 2},
    {"_TransPhylo_log_subtract_exp", (DL_FUNC) &_TransPhylo_log_subtract_exp, 2},
    {"_TransPhylo_log_sum_exp_vec", (DL_FUNC) &_TransPhylo_log_sum_exp_vec, 1},
    {"_TransPhylo_wbar", (DL_FUNC) &_TransPhylo_wbar, 10},
    {"_TransPhylo_probTTree", (DL_FUNC) &_TransPhylo_probTTree, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_TransPhylo(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
