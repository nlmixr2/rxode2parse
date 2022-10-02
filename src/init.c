#define USE_FC_LEN_T
#define STRICT_R_HEADERS
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "../inst/include/rxode2parse.h"
#define __DOINIT__
#include "tran.h"

SEXP _rxode2parse_codeLoaded();

SEXP _rxode2parse_codegen(SEXP c_file, SEXP prefix, SEXP libname, SEXP pMd5, SEXP timeId, SEXP lastMv);
SEXP _rxode2parse_parseModel(SEXP type);
SEXP _rxode2parse_isLinCmt();

void transIniNull();

void R_init_rxode2parse(DllInfo *info){
  R_CallMethodDef callMethods[]  = {
    {"_rxode2parse_trans", (DL_FUNC) &_rxode2parse_trans, 8},
    {"_rxode2parse_codegen", (DL_FUNC) &_rxode2parse_codegen, 6},
    {"_rxode2parse_codeLoaded", (DL_FUNC) &_rxode2parse_codeLoaded, 0},
    {"_rxode2parse_parseModel", (DL_FUNC) &_rxode2parse_parseModel, 1},
    {"_rxode2parse_isLinCmt", (DL_FUNC) &_rxode2parse_isLinCmt, 0},
    {NULL, NULL, 0} 
  };
  // C callable to assign environments.
  R_RegisterCCallable("rxode2parse","_rxode2parse_trans", (DL_FUNC) &_rxode2parse_trans);
  R_RegisterCCallable("rxode2parse","_rxode2parse_codegen", (DL_FUNC) &_rxode2parse_codegen);
  R_RegisterCCallable("rxode2parse","_rxode2parse_codeLoaded", (DL_FUNC) &_rxode2parse_codeLoaded);
  R_RegisterCCallable("rxode2parse","_rxode2parse_parseModel", (DL_FUNC) &_rxode2parse_parseModel);
  R_RegisterCCallable("rxode2parse","_rxode2parse_isLinCmt", (DL_FUNC) &_rxode2parse_isLinCmt); 
  // log likelihoods used in calculations
  static const R_CMethodDef cMethods[] = {
    {NULL, NULL, 0, NULL}
  };

  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  transIniNull();
}

void R_unload_rxode2parse(DllInfo *info){
  parseFree(1);
}
