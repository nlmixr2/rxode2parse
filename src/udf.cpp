#define USE_FC_LEN_T
#define STRICT_R_HEADER
#include <Rcpp.h>
using namespace Rcpp;

Function loadNamespace("loadNamespace", R_BaseNamespace);
//Function requireNamespace("requireNamespace", R_BaseNamespace);
Environment rxode2parseNS = loadNamespace("rxode2parse");
Function rxode2parse_getUdf_ = as<Function>(rxode2parseNS[".getUdfInfo"]);
Function rxode2parse_evalUdf = as<Function>(rxode2parseNS[".udfCall"]);

extern "C" SEXP rxode2parse_getUdf(const char *fun) {
  return rxode2parse_getUdf_(fun);
}

extern "C" double _rxode2parse_evalUdf(const char *fun, int n, const double *args) {
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
}
