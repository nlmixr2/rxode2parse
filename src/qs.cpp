#define USE_FC_LEN_T
#define STRICT_R_HEADERS
#include <Rcpp.h>

using namespace Rcpp;

Rcpp::Function loadNamespaceQsParse("loadNamespace", R_BaseNamespace);
Rcpp::Environment qsParseNs;
bool loadQsParseC = false;

static void loadQsParse() {
  if (!loadQsParseC) {
    qsParseNs = loadNamespaceQsParse("qs");
    loadQsParseC = true;
  }
}

//[[Rcpp::export]]
Rcpp::CharacterVector rxQsParse(SEXP const x) {
  loadQsParse();
  Rcpp::Function base91_encode = Rcpp::as<Rcpp::Function>(qsParseNs["base91_encode"]);
  Rcpp::Function qsParseerialize = Rcpp::as<Rcpp::Function>(qsParseNs["qserialize"]);
  return base91_encode(qsParseerialize(x, Rcpp::CharacterVector::create("high"), Rcpp::CharacterVector::create("zstd"),
				    Rcpp::IntegerVector::create(22),
				    Rcpp::IntegerVector::create(15), Rcpp::LogicalVector::create(true)));
}

//[[Rcpp::export]]
SEXP rxQrParse(const std::string& encoded_string) {
  loadQsParse();
  Rcpp::Function base91_decode = Rcpp::as<Rcpp::Function>(qsParseNs["base91_decode"]);
  Rcpp::Function qdeserialize = Rcpp::as<Rcpp::Function>(qsParseNs["qdeserialize"]);
  return qdeserialize(base91_decode(Rcpp::wrap(encoded_string)), false, false);
}
