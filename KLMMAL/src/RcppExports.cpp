// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// LMEMIXEDQUICK
Rcpp::List LMEMIXEDQUICK(Rcpp::List MyData, bool kinIn, bool covVarRead, double epsilon, int remliter);
RcppExport SEXP _KLMMAL_LMEMIXEDQUICK(SEXP MyDataSEXP, SEXP kinInSEXP, SEXP covVarReadSEXP, SEXP epsilonSEXP, SEXP remliterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type MyData(MyDataSEXP);
    Rcpp::traits::input_parameter< bool >::type kinIn(kinInSEXP);
    Rcpp::traits::input_parameter< bool >::type covVarRead(covVarReadSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type remliter(remliterSEXP);
    rcpp_result_gen = Rcpp::wrap(LMEMIXEDQUICK(MyData, kinIn, covVarRead, epsilon, remliter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_KLMMAL_LMEMIXEDQUICK", (DL_FUNC) &_KLMMAL_LMEMIXEDQUICK, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_KLMMAL(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
