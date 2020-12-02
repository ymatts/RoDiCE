// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rcpp_pobs
NumericMatrix rcpp_pobs(NumericMatrix x);
RcppExport SEXP _RoDiCE_rcpp_pobs(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_pobs(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_cvm
double rcpp_cvm(NumericMatrix mat1, NumericMatrix mat2);
RcppExport SEXP _RoDiCE_rcpp_cvm(SEXP mat1SEXP, SEXP mat2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat1(mat1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mat2(mat2SEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_cvm(mat1, mat2));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_coptest
List rcpp_coptest(NumericMatrix mat1, NumericMatrix mat2, int nperm);
RcppExport SEXP _RoDiCE_rcpp_coptest(SEXP mat1SEXP, SEXP mat2SEXP, SEXP npermSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat1(mat1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mat2(mat2SEXP);
    Rcpp::traits::input_parameter< int >::type nperm(npermSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_coptest(mat1, mat2, nperm));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RoDiCE_rcpp_pobs", (DL_FUNC) &_RoDiCE_rcpp_pobs, 1},
    {"_RoDiCE_rcpp_cvm", (DL_FUNC) &_RoDiCE_rcpp_cvm, 2},
    {"_RoDiCE_rcpp_coptest", (DL_FUNC) &_RoDiCE_rcpp_coptest, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_RoDiCE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}