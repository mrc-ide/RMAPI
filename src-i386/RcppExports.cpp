// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// calc_hex_values_cpp
Rcpp::List calc_hex_values_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress);
RcppExport SEXP _RMAPI_calc_hex_values_cpp(SEXP argsSEXP, SEXP args_functionsSEXP, SEXP args_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type args_functions(args_functionsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type args_progress(args_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_hex_values_cpp(args, args_functions, args_progress));
    return rcpp_result_gen;
END_RCPP
}
// calc_intersections_cpp
Rcpp::List calc_intersections_cpp(Rcpp::List args);
RcppExport SEXP _RMAPI_calc_intersections_cpp(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_intersections_cpp(args));
    return rcpp_result_gen;
END_RCPP
}
// hexbarrier01
Rcpp::List hexbarrier01(Rcpp::List args_h);
RcppExport SEXP _RMAPI_hexbarrier01(SEXP args_hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args_h(args_hSEXP);
    rcpp_result_gen = Rcpp::wrap(hexbarrier01(args_h));
    return rcpp_result_gen;
END_RCPP
}
// rmapi_analysis_cpp
Rcpp::List rmapi_analysis_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress);
RcppExport SEXP _RMAPI_rmapi_analysis_cpp(SEXP argsSEXP, SEXP args_functionsSEXP, SEXP args_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type args_functions(args_functionsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type args_progress(args_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(rmapi_analysis_cpp(args, args_functions, args_progress));
    return rcpp_result_gen;
END_RCPP
}
// sim_falciparum_cpp
Rcpp::List sim_falciparum_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress);
RcppExport SEXP _RMAPI_sim_falciparum_cpp(SEXP argsSEXP, SEXP args_functionsSEXP, SEXP args_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type args_functions(args_functionsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type args_progress(args_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_falciparum_cpp(args, args_functions, args_progress));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RMAPI_calc_hex_values_cpp", (DL_FUNC) &_RMAPI_calc_hex_values_cpp, 3},
    {"_RMAPI_calc_intersections_cpp", (DL_FUNC) &_RMAPI_calc_intersections_cpp, 1},
    {"_RMAPI_hexbarrier01", (DL_FUNC) &_RMAPI_hexbarrier01, 1},
    {"_RMAPI_rmapi_analysis_cpp", (DL_FUNC) &_RMAPI_rmapi_analysis_cpp, 3},
    {"_RMAPI_sim_falciparum_cpp", (DL_FUNC) &_RMAPI_sim_falciparum_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_RMAPI(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}