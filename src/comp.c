#define USE_FC_LEN_T
#define STRICT_R_HEADERS
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>   /* dj: import intptr_t */
//#include "ode.h"
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("rxode2parse", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

#define op_global _rxode2parse_op_global
#define rx_global _rxode2parse_rx_global
#define AMT _rxode2parse_AMT
#define LAG _rxode2parse_LAG
#define RATE _rxode2parse_RATE
#define DUR _rxode2parse_DUR
#define calc_mtime _rxode2parse_calc_mtime
#define getTime_ _rxode2parse_getTime_
#define getTime _rxode2parse_getTime
#define _locateTimeIndex _rxode2parse_locateTimeIndex
#define handle_evidL _rxode2parse_handle_evidL

#include "../inst/include/rxode2parse.h"
#define _calcDerived _rxode2parse_calcDerived

extern rx_solving_options _rxode2parse_op_global;
extern rx_solve _rxode2parse_rx_global;
extern t_handle_evidL _rxode2parse_handle_evidL;
extern t_getDur _rxode2parse_getDur;
#define _getDur _rxode2parse_getDur

#include "../inst/include/rxode2parseHandleEvid.h"
#include "../inst/include/rxode2parseGetTime.h"
#include "../inst/include/rxode2parseHandleSs.h"
#ifndef max2
#define max2( a , b )  ( (a) > (b) ? (a) : (b) )
#endif
#include "solComp.h"
#include "compSSc.h"
#include "comp.h"
void solveWith1Pt_lin(double *yp,
                      double xout, double xp,
                      int *i,
                      int *istate,
                      rx_solving_options *op,
                      rx_solving_options_ind *ind,
                      t_update_inis u_inis,
                      void *ctx) {
  rx_solve *rx=(&rx_global);
  lin_context_c_t *lin =  (lin_context_c_t*)(ctx);
  int linCmt = ind->linCmt;
  /* int ret = 1; */
  /* int linNcmt = linCmtI[RxMvFlag_ncmt]; */
  /* int linKa = linCmtI[RxMvFlag_ka]; */
  /* rx->linKa = linKa; */
  /* rx->linNcmt = linNcmt; */
  int ret = 1;
  switch(rx->linNcmt) {
  case 3:
    ret = comp1solve3(yp + linCmt, &xout, &xp, // last time
                      lin->rate, // rate in central compartment
                      &(lin->ka),  &(lin->k10),
                      &(lin->k12), &(lin->k21),
                      &(lin->k13), &(lin->k31));
    break;
  case 2:
    ret = comp1solve2(yp + linCmt, &xout, &xp,
                      lin->rate, // rate in central compartment
                      &(lin->ka), &(lin->k10), &(lin->k12), &(lin->k21));
    break;
  case 1:
    ret = comp1solve1(yp + linCmt, &xout, &xp, // last time
                      lin->rate, &(lin->ka), &(lin->k10));
    break;
  }
  if (ret == 0) {
    for (int j = op->neq*(ind->n_all_times); j--;){ \
      ind->solve[j] = NA_REAL;                                          \
    }                                                                   \
    op->badSolve = 1;                                                   \
    *i = ind->n_all_times-1; // Get out of here!
  }
}

void handleSSbolus_lin(double *yp,
                       double *xout, double xp,
                       int *i,
                       int *istate,
                       rx_solving_options *op,
                       rx_solving_options_ind *ind,
                       t_update_inis u_inis,
                       void *ctx,
                       double *xout2,
                       double *xp2,
                       double *curIi,
                       int *canBreak,
                       solveWith1Pt_fn solveWith1Pt) {
  lin_context_c_t *lin =  (lin_context_c_t*)(ctx);
  double ii = getIi(ind, ind->ix[*i]);
  double dose = getDose(ind, ind->ix[*i]);
  int oral0 = lin->oral0;
  int linCmt = ind->linCmt;
  int cmtOff = ind->cmt- ind->linCmt;
  int ret = 1;
  rx_solve *rx=(&rx_global);
  switch(rx->linNcmt) {
  case 3:
    if (cmtOff == 0 && oral0 == 1) {
      comp3ssBolusDepot(yp + linCmt, &ii, &dose, &(lin->ka), &(lin->k10), &(lin->k12), &(lin->k21),
                        &(lin->k13), &(lin->k31));
    } else {
      comp3ssBolusCentral(yp + linCmt, &ii, &dose, &(lin->ka), &(lin->k10), &(lin->k12), &(lin->k21),
                          &(lin->k13), &(lin->k31));
    }
    break;
  case 2:
    if (cmtOff == 0 && oral0 == 1) {
      comp2ssBolusDepot(yp + linCmt, &ii, &dose, &(lin->ka), &(lin->k10), &(lin->k12), &(lin->k21));
    } else {
      comp2ssBolusCentral(yp + linCmt, &ii, &dose, &(lin->ka), &(lin->k10), &(lin->k12), &(lin->k21));
    }
    break;
  case 1:
    if (cmtOff == 0 && oral0 == 1) {
      comp1ssBolusDepot(yp + linCmt, &ii, &dose, &(lin->ka), &(lin->k10));
    } else {
      comp1ssBolusCentral(yp + linCmt, &ii, &dose, &(lin->ka), &(lin->k10));
    }
    break;
  }
  if (ret == 0) {
    for (int j = op->neq*(ind->n_all_times); j--;){ \
      ind->solve[j] = NA_REAL;                                          \
    }                                                                   \
    op->badSolve = 1;                                                   \
    *i = ind->n_all_times-1; // Get out of here!
  }
}
void solveSSinf_lin(double *yp,
                    double *xout, double xp,
                    int *i,
                    int *istate,
                    rx_solving_options *op,
                    rx_solving_options_ind *ind,
                    t_update_inis u_inis,
                    void *ctx,
                    double *xout2,
                    double *xp2,
                    int *infBixds,
                    int *bi,
                    int *infEixds,
                    int *ei,
                    double *curIi,
                    double *dur,
                    double *dur2,
                    int *canBreak,
                    solveWith1Pt_fn solveWith1Pt) {
  lin_context_c_t *lin =  (lin_context_c_t*)(ctx);
  int linCmt = ind->linCmt;
  rx_solve *rx=(&rx_global);
  switch(rx->linNcmt) {
  case 3:
    comp3ssInf(yp + linCmt, dur, curIi, lin->rate,
               &(lin->ka), &(lin->k10), &(lin->k12), &(lin->k21),
               &(lin->k13), &(lin->k31));
    break;
  case 2:
    comp2ssInf(yp + linCmt, dur, curIi, lin->rate,
               &(lin->ka), &(lin->k10), &(lin->k12), &(lin->k21));
    break;
  case 1:
    comp1ssInf(yp + linCmt, dur, curIi, lin->rate,
               &(lin->ka), &(lin->k10));
    break;
  }
}

void rxode2parse_sortInd0(rx_solving_options_ind *ind);

SEXP _rxode2parse_compC(SEXP in) {
  rx_solve *rx=(&rx_global);
  rx_solving_options *op = rx->op;
  rx_solving_options_ind *oldInd = rx->subjects;
  iniSolvingRx(rx);
  iniSolvingOptions(op);
  int pro = 0;
  SEXP dat = PROTECT(VECTOR_ELT(in, 0)); pro++;
  SEXP par = PROTECT(VECTOR_ELT(in, 1)); pro++;
  double rate[2];
  rate[0] = rate[1] = 0.0;
  int cnt = 0;
  rx_solving_options_ind ind;
  ind.bT = 0.0;
  ind.slvr_counter = &cnt;
  ind.dadt_counter = &cnt;
  ind.jac_counter = &cnt;
  ind.InfusionRate = &rate;
  //ind.BadDose
  // ind.nBadDose
  ind.HMAX = 0.0; // Determined by diff
  ind.curDose = NA_REAL;
  ind.dosenum = 0;
  double tlast[2];
  tlast[0] = tlast[1] = NA_REAL;
  ind.tlastS = tlast;
  double curDoseS[2];
  curDoseS[0] = curDoseS[1] = NA_REAL;
  ind.curDoseS = curDoseS;
  double tfirstS[2];
  tfirstS[0] = tfirstS[1] = NA_REAL;
  ind.tfirstS = tfirstS;
  ind.podo = 0.0;
  double podoS[2];
  podoS[0] = podoS[1] = 0.0;
  ind.podoS = podoS;
  // double *par_ptr; // both time changing and time invariant
  ind.par_ptr= NULL;
  // double *solve;
  ind.solve = NULL;
  // double *mtime;
  ind.mtime = NULL;
  // double *lhs;
  ind.lhs = NULL;
  // double *cov_ptr;
  ind.cov_ptr = NULL;
  // int *cov_sample;
  ind.cov_sample=NULL;
  // double *dv;
  ind.dv = NULL;
  //double *limit;
  ind.limit = NULL;
  //int *cens;
  ind.cens = NULL;


  // ..$ TIME: num [1:135] 0 0 1 2 3 4 5 6 7 8 ...
  ind.all_times  = REAL(VECTOR_ELT(dat, 1));
  ind.n_all_times = Rf_length(VECTOR_ELT(dat, 1));
  // ..$ EVID: int [1:135] 101 0 0 0 0 0 0 0 0 0 ...
  ind.evid = INTEGER(VECTOR_ELT(dat, 1)); // $EVID;
  // ..$ AMT : num [1:135] 100 NA NA NA NA NA NA NA NA NA ...
  ind.dose = REAL(VECTOR_ELT(dat, 2)); // $AMT
  // ..$ II  : num [1:135] 0 0 0 0 0 0 0 0 0 0 ...
  ind.ii = REAL(VECTOR_ELT(dat, 3)); // $II

  double solveSave[2];
  solveSave[0]   = solveSave[1] = 0.0;
  ind.solveSave  = solveSave;
  double solveLast[2];
  solveLast[0]   = solveLast[1] =0.0;
  ind.solveLast  = solveLast;
  double solveLast2[2];
  solveLast2[0]  = solveLast2[1] = 0.0;
  ind.solveLast2 = solveLast2;
  /* int ixds; */
  ind.ixds=0;
  /* int ndoses; */
  ind.ndoses=0;
  /* int nevid2; */
  ind.nevid2 = 0;
  /* int id; */
  /* int solveid; */
  /* int idReal; */
  /* int sim; */
  /* int idx; */
  /* int solveid; */
  /* int idReal; */
  /* int sim; */
  /* int idx; */
  /* int yj; */
  /* int wh; */
  /* int wh100; */
  /* int cmt; */
  /* int whI; */
  /* int wh0; */
  /* int doSS; */
  /* int allCovWarn; */
  /* int wrongSSDur; */
  /* int _newind; */
  /* int _rxFlag; */
  /* int err; */
  /* int solved; */
  /* int _update_par_ptr_in; */
  /* int cacheME; */
  /* int inLhs; */

  ind.id = ind.solveid = ind.idReal = ind.sim = ind.idx = ind.yj =
    ind.wh = ind.wh100 = ind.cmt = ind.whI = ind.wh0 = ind.allCovWarn =
    ind.wrongSSDur = ind._newind = ind._rxFlag = ind.err =
    ind._update_par_ptr_in = ind.cacheME = ind.inLhs =
    ind.isIni =  ind.linCmt = 0;

  ind.solved = -1;

  ind.curShift = 0.0;
  /* bool lastIsSs2; */
  ind.lastIsSs2 = false;
  int on[2];
  on[0] = on[1] = 1;
  ind.on = on;
  
  /* double solveTime; */
  /* double curShift; */

  int *idose = malloc(2*ind.n_all_times*sizeof(int));
  if (idose == NULL) {
    Rf_errorcall(R_NilValue, _("ran out of memory"));
  }
  /* int *ix; */
  ind.ix = idose + ind.n_all_times;
  /* int *idose; */
  ind.idose = idose;

  for (int i = 0; i < ind.n_all_times; ++i) {
    ind.ix[i] = i;
    if (ind.evid[i] == 2) {
      ind.nevid2++;
    } else if (!isObs(ind.evid[i])) {
      ind.idose[ind.ndoses] = i;
      ind.ndoses++;
    }
  }
  int *tmpI = (int*)malloc(EVID_EXTRA_SIZE* sizeof(int));
  if (tmpI == NULL) {
    free(idose);
    Rf_errorcall(R_NilValue, _("ran out of memory"));
  }
  ind.ignoredDoses = tmpI;
  tmpI = (int*)malloc(6*sizeof(int));
  if (tmpI == NULL) {
    free(idose);
    free(ind.ignoredDoses);
    Rf_errorcall(R_NilValue, _("ran out of memory"));
  }
  ind.ignoredDosesN = tmpI;
  ind.ignoredDosesN[0] = 0;
  ind.ignoredDosesAllocN = ind.ignoredDosesN + 1;
  ind.ignoredDosesAllocN[0] = EVID_EXTRA_SIZE;
  ind.pendingDosesN = ind.ignoredDosesAllocN + 1;
  ind.pendingDosesN[0] = 0;
  ind.pendingDosesAllocN = ind.pendingDosesN + 1;
  ind.pendingDosesAllocN[0] =EVID_EXTRA_SIZE;
  ind.extraDoseN = ind.pendingDosesAllocN + 1;
  ind.extraDoseN[0] = 0;
  ind.extraDoseAllocN = ind.extraDoseN + 1;
  ind.extraDoseAllocN[0] = EVID_EXTRA_SIZE;
  tmpI = (int*)malloc(EVID_EXTRA_SIZE* sizeof(int));
  if (tmpI == NULL) {
    free(idose);
    free(ind.ignoredDoses);
    free(ind.ignoredDosesN);
    Rf_errorcall(R_NilValue, _("ran out of memory"));
  }
  ind.pendingDoses = tmpI;
  tmpI = (int*)malloc(EVID_EXTRA_SIZE* sizeof(int));
  if (tmpI == NULL) {
    free(idose);
    free(ind.ignoredDoses);
    free(ind.ignoredDosesN);
    free(ind.pendingDoses);
    Rf_errorcall(R_NilValue, _("ran out of memory"));
  }
  ind.extraDoseTimeIdx = tmpI;
  tmpI = (int*)malloc(EVID_EXTRA_SIZE* sizeof(int));
  if (tmpI == NULL) {
    free(idose);
    free(ind.ignoredDoses);
    free(ind.ignoredDosesN);
    free(ind.pendingDoses);
    free(ind.extraDoseTimeIdx);
    Rf_errorcall(R_NilValue, _("ran out of memory"));
  }
  ind.extraDoseEvid = tmpI;
  double *tmpD = (double*)malloc(EVID_EXTRA_SIZE* sizeof(double));
  if (tmpD == NULL) {
    free(idose);
    free(ind.ignoredDoses);
    free(ind.ignoredDosesN);
    free(ind.pendingDoses);
    free(ind.extraDoseTimeIdx);
    Rf_errorcall(R_NilValue, _("ran out of memory"));
  }
  ind.extraDoseTime = tmpD;
  if (tmpD == NULL) {
    free(idose);
    free(ind.ignoredDoses);
    free(ind.ignoredDosesN);
    free(ind.pendingDoses);
    free(ind.extraDoseTimeIdx);
    free(ind.extraDoseTime);
    Rf_errorcall(R_NilValue, _("ran out of memory"));
  }
  ind.extraDoseDose = tmpD;
  // extra doses
  ind.extraDoseNewXout = NA_REAL;
  ind.idxExtra = 0;
  ind.extraSorted = 0;

  rxode2parse_sortInd0(&ind);
  ind.ixds = ind.idx=0;

  rx->subjects =  &ind;

  /* if (op->badSolve) return 0; */
  /* if (ncmt) ind->pendingDosesN[0] = 0; */
  /* return 1; */

  /* int *rc; */
  /* // a b */
  /* // 1 4 */
  /* // 2 5 */
  /* // 3 6 */
  /* // Cache alag */
  /* double *alag; */
  /* // Cache F */
  /* double *cF; */
  /* // Cache rate; */
  /* double *cRate; */
  /* // Cache duration */
  /* double *cDur; */
  /* double *simIni; */
  /* int badIni; */
  /* double *llikSave; */
  /* // Add pointers for drifting atol/rtol values during optimization */
  /* double *rtol2; */
  /* double *atol2; */
  /* double *ssRtol; */
  /* double *ssAtol; */
  /* // ignored and pending doses */
  /* //double *extraDoseIi; // ii doses unsupported */
  /* double *timeThread; */
  free(idose);
  free(ind.ignoredDoses);
  free(ind.ignoredDosesN);
  free(ind.pendingDoses);
  free(ind.extraDoseTimeIdx);
  free(ind.extraDoseTime);
  free(ind.extraDoseDose);
  rx->subjects = oldInd;
  UNPROTECT(pro);
  return R_NilValue;
}
