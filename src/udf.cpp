#define USE_FC_LEN_T
#define STRICT_R_HEADER
#include <Rcpp.h>
using namespace Rcpp;

Function loadNamespace("loadNamespace", R_BaseNamespace);
//Function requireNamespace("requireNamespace", R_BaseNamespace);

extern "C" SEXP rxode2parse_getUdf(const char *fun) {
BEGIN_RCPP
  Environment rxode2parseNS = loadNamespace("rxode2parse");
  Function rxode2parse_getUdf_ = as<Function>(rxode2parseNS[".getUdfInfo"]);
  return rxode2parse_getUdf_(fun);
END_RCPP
}

extern "C" double _rxode2parse_evalUdf(const char *fun, int n, const double *args) {
BEGIN_RCPP
  Environment rxode2parseNS = loadNamespace("rxode2parse");
  Function rxode2parse_evalUdf = as<Function>(rxode2parseNS[".udfCall"]);
  List retL(n);
  CharacterVector funC(1);
  funC = fun;
  for (unsigned int i = 0; i < n; ++i) {
    NumericVector nv(1);
    nv[0] = args[i];
    retL[i] = nv;
  }
  NumericVector ret = rxode2parse_evalUdf(funC, retL);
  return ret[0];
VOID_END_RCPP
  return NA_REAL;
}

extern "C" void _rxode2parse_resetUdf() {
BEGIN_RCPP
  Environment rxode2parseNS = loadNamespace("rxode2parse");
  Function resetUdf = as<Function>(rxode2parseNS[".udfReset"]);
  resetUdf();
VOID_END_RCPP
}

extern "C" SEXP _rxode2parse_getUdf() {
BEGIN_RCPP
  Environment rxode2parseNS = loadNamespace("rxode2parse");
  Function getUdf = as<Function>(rxode2parseNS[".udfInfo"]);
  return getUdf();
END_RCPP
}
