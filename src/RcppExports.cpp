// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pf_barnes2014
NumericMatrix pf_barnes2014(NumericMatrix dem);
RcppExport SEXP _flowdem_pf_barnes2014(SEXP demSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dem(demSEXP);
    rcpp_result_gen = Rcpp::wrap(pf_barnes2014(dem));
    return rcpp_result_gen;
END_RCPP
}
// pf_eps_barnes2014
NumericMatrix pf_eps_barnes2014(NumericMatrix dem);
RcppExport SEXP _flowdem_pf_eps_barnes2014(SEXP demSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dem(demSEXP);
    rcpp_result_gen = Rcpp::wrap(pf_eps_barnes2014(dem));
    return rcpp_result_gen;
END_RCPP
}
// pf_basins_barnes2014
List pf_basins_barnes2014(NumericMatrix dem);
RcppExport SEXP _flowdem_pf_basins_barnes2014(SEXP demSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dem(demSEXP);
    rcpp_result_gen = Rcpp::wrap(pf_basins_barnes2014(dem));
    return rcpp_result_gen;
END_RCPP
}
// comp_breach_lindsay2016
NumericMatrix comp_breach_lindsay2016(NumericMatrix dem);
RcppExport SEXP _flowdem_comp_breach_lindsay2016(SEXP demSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dem(demSEXP);
    rcpp_result_gen = Rcpp::wrap(comp_breach_lindsay2016(dem));
    return rcpp_result_gen;
END_RCPP
}
// d8_flow_directions
IntegerMatrix d8_flow_directions(NumericMatrix dem);
RcppExport SEXP _flowdem_d8_flow_directions(SEXP demSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dem(demSEXP);
    rcpp_result_gen = Rcpp::wrap(d8_flow_directions(dem));
    return rcpp_result_gen;
END_RCPP
}
// d8_flow_accum
NumericMatrix d8_flow_accum(IntegerMatrix flowdirs);
RcppExport SEXP _flowdem_d8_flow_accum(SEXP flowdirsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type flowdirs(flowdirsSEXP);
    rcpp_result_gen = Rcpp::wrap(d8_flow_accum(flowdirs));
    return rcpp_result_gen;
END_RCPP
}
// d8_watershed_nested
IntegerMatrix d8_watershed_nested(IntegerMatrix flowdirs, NumericMatrix target_rc, bool nested);
RcppExport SEXP _flowdem_d8_watershed_nested(SEXP flowdirsSEXP, SEXP target_rcSEXP, SEXP nestedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type flowdirs(flowdirsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type target_rc(target_rcSEXP);
    Rcpp::traits::input_parameter< bool >::type nested(nestedSEXP);
    rcpp_result_gen = Rcpp::wrap(d8_watershed_nested(flowdirs, target_rc, nested));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_flowdem_pf_barnes2014", (DL_FUNC) &_flowdem_pf_barnes2014, 1},
    {"_flowdem_pf_eps_barnes2014", (DL_FUNC) &_flowdem_pf_eps_barnes2014, 1},
    {"_flowdem_pf_basins_barnes2014", (DL_FUNC) &_flowdem_pf_basins_barnes2014, 1},
    {"_flowdem_comp_breach_lindsay2016", (DL_FUNC) &_flowdem_comp_breach_lindsay2016, 1},
    {"_flowdem_d8_flow_directions", (DL_FUNC) &_flowdem_d8_flow_directions, 1},
    {"_flowdem_d8_flow_accum", (DL_FUNC) &_flowdem_d8_flow_accum, 1},
    {"_flowdem_d8_watershed_nested", (DL_FUNC) &_flowdem_d8_watershed_nested, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_flowdem(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
