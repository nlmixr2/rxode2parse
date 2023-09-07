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

#ifndef max2
#define max2( a , b )  ( (a) > (b) ? (a) : (b) )
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
int rxode2parse_handleExtraDose0(double *yp, double xout, double xp, int *i, rx_solving_options *op, rx_solving_options_ind *ind);

SEXP _rxode2parse_compC(SEXP in, SEXP mv) {
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
  rx_solving_options_ind indR;
  indR.bT = 0.0;
  indR.slvr_counter = &cnt;
  indR.dadt_counter = &cnt;
  indR.jac_counter = &cnt;
  indR.InfusionRate = rate;
  //indR.BadDose
  // indR.nBadDose
  indR.HMAX = 0.0; // Determined by diff
  indR.curDose = NA_REAL;
  indR.dosenum = 0;
  double tlast[2];
  tlast[0] = tlast[1] = NA_REAL;
  indR.tlastS = tlast;
  double curDoseS[2];
  curDoseS[0] = curDoseS[1] = NA_REAL;
  indR.curDoseS = curDoseS;
  double tfirstS[2];
  tfirstS[0] = tfirstS[1] = NA_REAL;
  indR.tfirstS = tfirstS;
  indR.podo = 0.0;
  double podoS[2];
  podoS[0] = podoS[1] = 0.0;
  indR.podoS = podoS;
  // double *par_ptr; // both time changing and time invariant
  indR.par_ptr= NULL;
  // double *solve;
  indR.solve = NULL;
  // double *mtime;
  indR.mtime = NULL;
  // double *lhs;
  indR.lhs = NULL;
  // double *cov_ptr;
  indR.cov_ptr = NULL;
  // int *cov_sample;
  indR.cov_sample=NULL;
  // double *dv;
  indR.dv = NULL;
  //double *limit;
  indR.limit = NULL;
  //int *cens;
  indR.cens = NULL;


  // ..$ TIME: num [1:135] 0 0 1 2 3 4 5 6 7 8 ...
  indR.all_times  = REAL(VECTOR_ELT(dat, 1));
  indR.n_all_times = Rf_length(VECTOR_ELT(dat, 1));
  // ..$ EVID: int [1:135] 101 0 0 0 0 0 0 0 0 0 ...
  indR.evid = INTEGER(VECTOR_ELT(dat, 1)); // $EVID;
  rx->nall = Rf_length(VECTOR_ELT(dat, 1));
  rx->nobs = 0;
  rx->nobs2 = 0;
  rx->nevid9 = 0;
  for (int j = 0; j < rx->nall; ++j) {
    if (isObs(indR.evid[j])) rx->nobs++;
    if (indR.evid[j] == 0) rx->nobs2++;
    if (indR.evid[j] == 9) rx->nevid9++;
  }
  rx->nKeepF     = 0;
  rx->nCov0      = 0;
  rx->neps       = 0;
  rx->neps       = 0;
  rx->neta       = 0;
  rx->hasFactors = 0;
  rx->nsub       = 1;
  rx->nsim       = 1;

  op->ncoresRV = 0;
  op->do_par_cov=0;
  op->cores = 0;
  op->neq = 0; // no states by default
  op->nLlik = 0; // no likelihoods by default
  op->badSolve = 0;
  op->naTime = 0;
  op->abort = 0;
  op->ATOL = NA_REAL;
  op->RTOL = NA_REAL;
  op->indLinPhiTol = NA_REAL;
  op->indLinMatExpType = NA_REAL;
  op->indLinPhiM = NA_REAL;
  op->indLinMatExpOrder=0;
  op->doIndLin=0;
  op->H0 = NA_REAL;
  op->HMIN = NA_REAL;
  op->mxstep = -1;
  op->MXORDN = 0;
  op->MXORDS = 0;
  op->nlhs = 0;
  op->is_locf = 1;
  op->f2 = 0.0;
  op->f1 = 1.0;
  op->kind = 0;
  op->extraCmt= INTEGER(VECTOR_ELT(mv, RxMv_extraCmt))[0];
  op->nDisplayProgress = 0;
  SEXP flagsS  = VECTOR_ELT(mv, RxMv_flags);
  int *flags   = INTEGER(flagsS);
  rx->linKa    = flags[RxMvFlag_ka];
  // RxMvFlag_linB
  op->linBflag = flags[RxMvFlag_linCmtFlg];
  rx->linNcmt  = flags[RxMvFlag_ncmt];

  // ..$ AMT : num [1:135] 100 NA NA NA NA NA NA NA NA NA ...
  indR.dose = REAL(VECTOR_ELT(dat, 2)); // $AMT
  // ..$ II  : num [1:135] 0 0 0 0 0 0 0 0 0 0 ...
  indR.ii = REAL(VECTOR_ELT(dat, 3)); // $II

  double solveSave[2];
  solveSave[0]   = solveSave[1] = 0.0;
  indR.solveSave  = solveSave;
  double solveLast[2];
  solveLast[0]   = solveLast[1] =0.0;
  indR.solveLast  = solveLast;
  double solveLast2[2];
  solveLast2[0]  = solveLast2[1] = 0.0;
  indR.solveLast2 = solveLast2;
  /* int ixds; */
  indR.ixds=0;
  /* int ndoses; */
  indR.ndoses=0;
  /* int nevid2; */
  indR.nevid2 = 0;
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

  indR.id = indR.solveid = indR.idReal = indR.sim = indR.idx = indR.yj =
    indR.wh = indR.wh100 = indR.cmt = indR.whI = indR.wh0 = indR.allCovWarn =
    indR.wrongSSDur = indR._newind = indR._rxFlag = indR.err =
    indR._update_par_ptr_in = indR.cacheME = indR.inLhs =
    indR.isIni =  indR.linCmt = 0;

  indR.solved = -1;

  indR.curShift = 0.0;
  /* bool lastIsSs2; */
  indR.lastIsSs2 = false;
  int on[2];
  on[0] = on[1] = 1;
  indR.on = on;

  /* double solveTime; */
  /* double curShift; */

  int *idose = malloc(2*indR.n_all_times*sizeof(int));
  if (idose == NULL) {
    Rf_errorcall(R_NilValue, _("ran out of memory"));
  }
  /* int *ix; */
  indR.ix = idose + indR.n_all_times;
  /* int *idose; */
  indR.idose = idose;
  for (int i = 0; i < indR.n_all_times; ++i) {
    indR.ix[i] = i;
    if (indR.evid[i] == 2) {
      indR.nevid2++;
    } else if (!isObs(indR.evid[i])) {
      indR.idose[indR.ndoses] = i;
      indR.ndoses++;
    }
  }
  int *tmpI = (int*)malloc(EVID_EXTRA_SIZE* sizeof(int));
  if (tmpI == NULL) {
    free(idose);
    Rf_errorcall(R_NilValue, _("ran out of memory"));
  }
  indR.ignoredDoses = tmpI;
  tmpI = (int*)malloc(6*sizeof(int));
  if (tmpI == NULL) {
    free(idose);
    free(indR.ignoredDoses);
    Rf_errorcall(R_NilValue, _("ran out of memory"));
  }
  indR.ignoredDosesN = tmpI;
  indR.ignoredDosesN[0] = 0;
  indR.ignoredDosesAllocN = indR.ignoredDosesN + 1;
  indR.ignoredDosesAllocN[0] = EVID_EXTRA_SIZE;
  indR.pendingDosesN = indR.ignoredDosesAllocN + 1;
  indR.pendingDosesN[0] = 0;
  indR.pendingDosesAllocN = indR.pendingDosesN + 1;
  indR.pendingDosesAllocN[0] =EVID_EXTRA_SIZE;
  indR.extraDoseN = indR.pendingDosesAllocN + 1;
  indR.extraDoseN[0] = 0;
  indR.extraDoseAllocN = indR.extraDoseN + 1;
  indR.extraDoseAllocN[0] = EVID_EXTRA_SIZE;
  tmpI = (int*)malloc(EVID_EXTRA_SIZE* sizeof(int));
  if (tmpI == NULL) {
    free(idose);
    free(indR.ignoredDoses);
    free(indR.ignoredDosesN);
    Rf_errorcall(R_NilValue, _("ran out of memory"));
  }
  indR.pendingDoses = tmpI;
  tmpI = (int*)malloc(EVID_EXTRA_SIZE* sizeof(int));
  if (tmpI == NULL) {
    free(idose);
    free(indR.ignoredDoses);
    free(indR.ignoredDosesN);
    free(indR.pendingDoses);
    Rf_errorcall(R_NilValue, _("ran out of memory"));
  }
  indR.extraDoseTimeIdx = tmpI;
  tmpI = (int*)malloc(EVID_EXTRA_SIZE* sizeof(int));
  if (tmpI == NULL) {
    free(idose);
    free(indR.ignoredDoses);
    free(indR.ignoredDosesN);
    free(indR.pendingDoses);
    free(indR.extraDoseTimeIdx);
    Rf_errorcall(R_NilValue, _("ran out of memory"));
  }
  indR.extraDoseEvid = tmpI;
  double *tmpD = (double*)malloc(EVID_EXTRA_SIZE* sizeof(double));
  if (tmpD == NULL) {
    free(idose);
    free(indR.ignoredDoses);
    free(indR.ignoredDosesN);
    free(indR.pendingDoses);
    free(indR.extraDoseTimeIdx);
    Rf_errorcall(R_NilValue, _("ran out of memory"));
  }
  indR.extraDoseTime = tmpD;
  if (tmpD == NULL) {
    free(idose);
    free(indR.ignoredDoses);
    free(indR.ignoredDosesN);
    free(indR.pendingDoses);
    free(indR.extraDoseTimeIdx);
    free(indR.extraDoseTime);
    Rf_errorcall(R_NilValue, _("ran out of memory"));
  }
  indR.extraDoseDose = tmpD;
  // extra doses
  indR.extraDoseNewXout = NA_REAL;
  indR.idxExtra = 0;
  indR.extraSorted = 0;
  tmpD = (double*)calloc(indR.n_all_times*(rx->linNcmt+rx->linKa), sizeof(double));
  if (tmpD == NULL) {
    free(idose);
    free(indR.ignoredDoses);
    free(indR.ignoredDosesN);
    free(indR.pendingDoses);
    free(indR.extraDoseTimeIdx);
    free(indR.extraDoseTime);
    free(indR.extraDoseDose);
    Rf_errorcall(R_NilValue, _("ran out of memory"));
  }
  indR.solve = tmpD;

  indR.ixds = indR.idx=0;
  rx_solving_options_ind* ind = &indR;
  rx->subjects =  ind;
  rxode2parse_sortInd0(ind);
  double *yp;
  SEXP CcSxp = PROTECT(Rf_allocVector(REALSXP, indR.n_all_times)); pro++;

  double *Cc = REAL(CcSxp);
  double *p1 = REAL(VECTOR_ELT(par, 0));
  double *v1 = REAL(VECTOR_ELT(par, 1));
  double *p2 = REAL(VECTOR_ELT(par, 2));
  double *p3 = REAL(VECTOR_ELT(par, 3));
  double *p4 = REAL(VECTOR_ELT(par, 4));
  double *p5 = REAL(VECTOR_ELT(par, 5));
  double *lagdepot = REAL(VECTOR_ELT(par, 6));
  double *fdepot   = REAL(VECTOR_ELT(par, 7));
  double *ratedepot = REAL(VECTOR_ELT(par, 8));
  double *durdepot = REAL(VECTOR_ELT(par, 9));
  double *ka = REAL(VECTOR_ELT(par, 10));
  double *lagcentral = REAL(VECTOR_ELT(par, 11));
  double *fcentral = REAL(VECTOR_ELT(par, 12));
  double *ratecentral = REAL(VECTOR_ELT(par, 13));
  double *durcentral = REAL(VECTOR_ELT(par, 14));
  double xout = 0.0, xp= 0.0;
  int istate = 1;
  void *ctx = NULL;
  for(int i=0; i < indR.n_all_times; ++i) {
    ind->idx=i;
    yp = getSolve(i);
    xout = getTime_(ind->ix[i], ind);
    if (getEvid(ind, ind->ix[i]) != 3) {
      if (ind->err){
        //*rc = -1000;
        // Bad Solve => NA
        //badSolveExit(i);
      } else {
        if (rxode2parse_handleExtraDose0(yp, xout, xp, &i, op, ind)) {
          if (!isSameTime(ind->extraDoseNewXout, xp)) {
            /* F77_CALL(dlsoda)(dydt_lsoda, neq, yp, &xp, &ind->extraDoseNewXout, &gitol, &(op->RTOL), &(op->ATOL), &gitask, */
            /*                  &istate, &giopt, rwork, &lrw, iwork, &liw, jdum, &jt); */
            /* postSolve(&istate, ind->rc, &i, yp, err_msg_ls, 7, true, ind, op, rx); */
          }
          int idx = ind->idx;
          int ixds = ind->ixds;
          int trueIdx = ind->extraDoseTimeIdx[ind->idxExtra];
          ind->idx = -1-trueIdx;
          handle_evid(ind->extraDoseEvid[trueIdx], yp, xout, ind);
          istate = 1;
          ind->ixds = ixds; // This is a fake dose, real dose stays in place
          ind->idx = idx;
          ind->idxExtra++;
          if (!isSameTime(xout, ind->extraDoseNewXout)) {
            /* F77_CALL(dlsoda)(dydt_lsoda, neq, yp, &ind->extraDoseNewXout, &xout, &gitol, &(op->RTOL), &(op->ATOL), &gitask, */
            /*                  &istate, &giopt, rwork, &lrw, iwork, &liw, jdum, &jt); */
            /* postSolve(&istate, ind->rc, &i, yp, err_msg_ls, 7, true, ind, op, rx); */
          }
          xp =  ind->extraDoseNewXout;
        }
        if (!isSameTime(xout, xp)) {
          /* F77_CALL(dlsoda)(dydt_lsoda, neq, yp, &xp, &xout, &gitol, &(op->RTOL), &(op->ATOL), &gitask, */
          /*                  &istate, &giopt, rwork, &lrw, iwork, &liw, jdum, &jt); */
          /* postSolve(&istate, ind->rc, &i, yp, err_msg_ls, 7, true, ind, op, rx); */
        }
        xp = xout;
      }
    }
  }

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
  free(indR.ignoredDoses);
  free(indR.ignoredDosesN);
  free(indR.pendingDoses);
  free(indR.extraDoseTimeIdx);
  free(indR.extraDoseTime);
  free(indR.extraDoseDose);
  free(indR.solve);
  rx->subjects = oldInd;
  UNPROTECT(pro);
  return R_NilValue;
}
