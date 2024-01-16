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
#include "parTrans.h"

#include "../inst/include/rxode2parseIniSubject.h"

// does nothing since linCmt test is single threaded
void rxode2parse_setIndPointersByThread(rx_solving_options_ind *ind) {}

// both these functions require sorting from C++, these are interfaces
void rxode2parse_sortInd0(rx_solving_options_ind *ind);
void rxode2parse_sortRest0(rx_solving_options_ind *ind, int i0);
int rxode2parse_handleExtraDose0(double *yp, double xout, double xp, int *i, rx_solving_options *op, rx_solving_options_ind *ind);


// dummy function for ini reset
void lincmt_ini_dummy0(int cSub, double *x){
  return;
}

t_update_inis u_inis_lincmt = lincmt_ini_dummy0;

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

SEXP _rxode2parse_solve1ptLin(SEXP inp, SEXP tS, SEXP kaS, SEXP k10S,
                              SEXP k12S, SEXP k21S,
                              SEXP k13S, SEXP k31S,
                              SEXP vS, SEXP rateS) {
  lin_context_c_t lin;
  lin.ka = REAL(kaS)[0];
  rx_solve* rx = &rx_global;
  if (lin.ka > 0.0) {
    lin.oral0 = 1;
    rx->linKa = 1;
  } else {
    lin.oral0 = 0;
    rx->linKa = 0;
  }
  lin.k10 = REAL(k10S)[0];
  lin.k12 = REAL(k12S)[0];
  lin.k21 = REAL(k21S)[0];
  lin.k13 = REAL(k13S)[0];
  lin.k31 = REAL(k31S)[0];
  lin.v = REAL(vS)[0];
  lin.rate= REAL(rateS);
  rx_solving_options *op = rx->op;
  rx->linNcmt = Rf_length(inp) - lin.oral0;

  rx_solving_options_ind indR;
  indR.n_all_times = 1;
  double solve[1];
  solve[0] = 0.0;
  indR.solve = solve;
  indR.linCmt = 0;//rx->linNcmt; // linNcmt
  SEXP ret = PROTECT(Rf_allocVector(REALSXP, Rf_length(inp)));
  double* inpd = REAL(inp);
  double *yp = REAL(ret);
  for (int i = Rf_length(inp);i--;) {
    yp[i] = inpd[i];
  }
  int i =  0;
  solveWith1Pt_lin(yp,
                   REAL(tS)[0], 0.0,
                   &i,
                   &i, op,
                   &indR,
                   u_inis_lincmt,
                   (void *)(&lin));
  UNPROTECT(1);
  return ret;
}

void solveWith1Pt_lin(double *yp,
                      double xout, double xp,
                      int *i,
                      int *istate,
                      rx_solving_options *op,
                      rx_solving_options_ind *ind,
                      t_update_inis u_inis,
                      void *ctx);

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
  rx_solve *rx=(&rx_global);
  double ii = getIi(ind, ind->ix[*i]);
  double dose = getDose(ind, ind->ix[*i]);
  int linCmt = ind->linCmt;
  int cmtOff = ind->cmt- ind->linCmt;
  int ret = 1;
  switch(rx->linNcmt) {
  case 3:
    comp3ssBolus(&cmtOff, yp + linCmt, &ii, &dose, &(lin->ka), &(lin->k10), &(lin->k12), &(lin->k21),
                        &(lin->k13), &(lin->k31));
    break;
  case 2:
    comp2ssBolus(&cmtOff, yp + linCmt, &ii, &dose, &(lin->ka), &(lin->k10), &(lin->k12), &(lin->k21));
    break;
  case 1:
    comp1ssBolus(&cmtOff, yp + linCmt, &ii, &dose, &(lin->ka), &(lin->k10));
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

void solveSSinf8_lin(double *yp,
                     double *xout, double xp,
                     int *i,
                     int *istate,
                     rx_solving_options *op,
                     rx_solving_options_ind *ind,
                     t_update_inis u_inis,
                     void *ctx,
                     int *infBixds,
                     int *bi,
                     double *rateOn,
                     double *xout2,
                     double *xp2,
                     int *canBreak,
                     solveWith1Pt_fn solveWith1Pt) {
  lin_context_c_t *lin =  (lin_context_c_t*)(ctx);
  int linCmt = ind->linCmt;
  rx_solve *rx=(&rx_global);
  switch(rx->linNcmt) {
  case 3:
    comp3ssInf8(yp + linCmt, lin->rate,
                &(lin->ka), &(lin->k10), &(lin->k12), &(lin->k21),
                &(lin->k13), &(lin->k31));
    break;
  case 2:
    comp2ssInf8(yp + linCmt,lin->rate, &(lin->ka), &(lin->k10), &(lin->k12), &(lin->k21));
    break;
  case 1:
    comp1ssInf8(yp + linCmt, lin->rate, &(lin->ka), &(lin->k10));
    break;
  }
}

#define handleSS(neq, BadDose, InfusionRate, dose, yp, xout, xp, id, i, nx, istate, op, ind, u_inis, ctx) handleSSGen(neq, BadDose, InfusionRate, dose, yp, xout, xp, id, i, nx, istate, op, ind, u_inis, ctx, solveWith1Pt_lin, handleSSbolus_lin, solveSSinf_lin, solveSSinf8_lin)


#define badSolveExit(i) for (int j = (op->neq + op->extraCmt)*(ind->n_all_times); j--;){ \
    ind->solve[j] = NA_REAL;                                            \
  }                                                                     \
  op->badSolve = 1;                                                     \
  i = ind->n_all_times-1; // Get out of here!

double linCmtCompA(rx_solve *rx, unsigned int id, double _t, int linCmt,
                   int i_cmt, int trans,
                   double p1, double v1,
                   double p2, double p3,
                   double p4, double p5,
                   double d_tlag, double d_F, double d_rate1, double d_dur1,
                   // Oral parameters
                   double d_ka, double d_tlag2, double d_F2,  double d_rate2, double d_dur2) {
  lin_context_c_t lin;
  rx_solving_options_ind *ind = &(rx->subjects[id]);
  rx_solving_options *op = rx->op;
  double t = _t - ind->curShift;
  unsigned int ncmt=0;
  int neq[2];
  neq[0] = op->neq;
  neq[1] = ind->solveid;
  int istate=0; // dummy variable for consistent interface with lsoda
  if (!parTrans(&trans, &p1, &v1, &p2, &p3, &p4, &p5,
                &ncmt, &(lin.k10), &(lin.v), &(lin.k12),
                &(lin.k21), &(lin.k13), &(lin.k31))){
    return NA_REAL;
  }
  lin.ka = d_ka;
  lin.rate = ind->InfusionRate + op->neq;
  double xp, xout;
  double *ypLast, *yp;
  double Alast0[4] = {0, 0, 0, 0};
  int oral0 = rx->linKa;
  if (oral0) {
    /* double d_tlag, double d_F, double d_rate1, double d_dur1, */
    /*   // Oral parameters */
    /*   double d_ka, double d_tlag2, double d_F2,  double d_rate2, double d_dur2 */
    ind->linCmtLag[0] = d_tlag; // depot
    ind->linCmtLag[1] = d_tlag2; // central
    ind->linCmtF[0] = d_F;
    ind->linCmtF[1] = d_F2;
    ind->linCmtDur[0] = d_dur1;
    ind->linCmtDur[1] = d_dur2;
    ind->linCmtRate[0] = d_rate1;
    ind->linCmtRate[1] = d_rate2;
  } else {
    ind->linCmtLag[0] = d_tlag;
    ind->linCmtF[0]   = d_F;
    ind->linCmtDur[0] = d_dur1;
    ind->linCmtRate[0] = d_rate1;
  }
  rxode2parse_sortRest0(ind, ind->idx);
  void *ctx = &(lin);
  if (ind->idx == 0) {
    // initialization
    xp = xout = getTime__(ind->ix[ind->idx], ind, 0);
    yp = ypLast = Alast0;
  } else {
    xp = (ind->idx == 0 ? 0.0 : getTime_(ind->ix[ind->idx-1], ind));
    xout = getTime__(ind->ix[ind->idx], ind, 0);
    ypLast=getAdvan(ind->idx-1);
  }
  yp = getAdvan(ind->idx);
  if (ind->idx <= ind->solved) {
    // Pull from last solved value (cached)
    if (yp[oral0] == 0.0) {
      // it is zero, perhaps it wasn't solved, double check
      ind->solved = max2(ind->idx-1, 0);
    } else {
      if (trans == 10) {
        return(yp[oral0]*(v1+p3+p5));
      } else {
        return(yp[oral0]/v1);
      }
    }
  } else if (ind->idx != 0) {
    for (int j=0; j < rx->linNcmt + rx->linKa; ++j) {
      yp[j] = ypLast[j];
    }
  }
  int i = ind->idx;
  if (getEvid(ind, ind->ix[i]) != 3) {
    if (ind->err){
      printErr(ind->err, ind->id);
      // Bad Solve => NA
      badSolveExit(i);
      return NA_REAL;
    } else {
      if (rxode2parse_handleExtraDose0(yp, xout, xp, &i, op, ind)) {
        if (!isSameTimeOp(ind->extraDoseNewXout, xp)) {
          solveWith1Pt_lin(yp, xout, xp, &(ind->idx),
                           &istate, op, ind, u_inis_lincmt, (void *)(&lin));
          /* postSolve(&idid, rc, &i, yp, err_msg, 4, true, ind, op, rx); */
        }
        int idx = ind->idx;
        int ixds = ind->ixds;
        int trueIdx = ind->extraDoseTimeIdx[ind->idxExtra];
        ind->idx = -1-trueIdx;
        handle_evid(ind->extraDoseEvid[trueIdx], yp, xout, ind);
        ind->idx = idx;
        ind->ixds = ixds;
        ind->idxExtra++;
        if (!isSameTimeOp(xout, ind->extraDoseNewXout)) {
          solveWith1Pt_lin(yp, ind->extraDoseNewXout, xout, &(ind->idx),
                           &istate, op, ind, u_inis_lincmt, (void *)(&lin));
          /* postSolve(&idid, rc, &i, yp, err_msg, 4, true, ind, op, rx); */
        }
        xp = ind->extraDoseNewXout;
      }
      if (!isSameTimeOp(xout, xp)) {
        solveWith1Pt_lin(yp, xout, xp, &(ind->idx),
                         &istate, op, ind, u_inis_lincmt, (void *)(&lin));
        /* postSolve(&idid, rc, &i, yp, err_msg, 4, true, ind, op, rx); */
        xp = xout;
      }
      //dadt_counter = 0;
    }
  }
  if (!op->badSolve) {
    ind->idx = i;
    if (getEvid(ind, ind->ix[i]) == 3){
      ind->curShift -= rx->maxShift;
      for (int j = op->neq + op->extraCmt; j--;) {
        ind->InfusionRate[j] = 0;
        ind->on[j] = 1;
        ind->cacheME=0;
      }
      cancelInfusionsThatHaveStarted(ind, ind->solveid, xout);
      cancelPendingDoses(ind, neq[1]);
      /* memcpy(yp, op->inits, neq[0]*sizeof(double)); */
      u_inis_lincmt(neq[1], yp); // Update initial conditions @ current time */
      ind->ixds++;
      xp=xout;
    } else if (handleEvid1(&i, rx, neq, yp, &xout)){
      handleSS(neq, ind->BadDose, ind->InfusionRate, ind->dose, yp, xout,
               xp, ind->id, &i, ind->n_all_times, &istate, op, ind, u_inis_lincmt, ctx);
      if (ind->wh0 == EVID0_OFF) {
        yp[ind->cmt] = op->inits[ind->cmt];
      }
      xp = xout;
    }
    if (i+1 != ind->n_all_times) {
      //ypLast=getAdvan(ind->idx-1);
      double *ypNext = getAdvan(ind->idx+1);
      for (int j=0; j < rx->linNcmt + rx->linKa; ++j) {
        ypNext[j] = yp[j];
      }
    }
    //if (i+1 != nx) memcpy(getSolve(i+1), getSolve(i), neq[0]*sizeof(double));
    //calc_lhs(neq[1], xout, getSolve(i), ind->lhs);
    //updateExtraDoseGlobals(ind);
  }
  return(yp[oral0]/lin.v);
}


SEXP _rxode2parse_compC(SEXP in, SEXP mv) {
  rx_solve *rx=(&rx_global);
  rx_solving_options *op = &(op_global);
  rx->op = op;
  rx_solving_options_ind *oldInd = rx->subjects;
  iniSolvingRx(rx);
  iniSolvingOptions(op);
  int pro = 0;
  SEXP dat = PROTECT(VECTOR_ELT(in, 0)); pro++; // event table
  SEXP par = PROTECT(VECTOR_ELT(in, 1)); pro++; // parameter table
  int trans = INTEGER(VECTOR_ELT(in, 2))[0];
  double sm = REAL(VECTOR_ELT(in, 3))[0];
  double rate[2];
  rate[0] = rate[1] = 0.0;
  int cnt = 0;
  rx_solving_options_ind indR;
  indR.bT = 0.0;
  indR.slvr_counter = &cnt;
  indR.dadt_counter = &cnt;
  indR.jac_counter = &cnt;
  indR.InfusionRate = rate;
  int BadDose[2];
  BadDose[0] = BadDose[1] = 0;
  indR.BadDose = BadDose;
  double linCmtLag[2];
  linCmtLag[0] = linCmtLag[1] = 0.0;
  indR.linCmtLag = linCmtLag;
  double linCmtF[2];
  linCmtF[0] = linCmtF[1] = 1.0;
  indR.linCmtF = linCmtF;
  double linCmtDur[2];
  linCmtDur[0] = linCmtDur[1] = 0.0;
  indR.linCmtDur = linCmtDur;

  double linCmtRate[2];
  linCmtRate[0] = linCmtRate[1] = 0.0;
  indR.linCmtRate = linCmtRate;
  int rc[1];
  rc[0] = 0;
  indR.rc=rc;
  indR.nBadDose = 0;
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
  indR.all_times  = REAL(VECTOR_ELT(dat, 0));
  indR.n_all_times = Rf_length(VECTOR_ELT(dat, 0));
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
  double ssAtol[2];
  ssAtol[0] = ssAtol[1] = 1.0e-8;
  double ssRtol[2];
  ssRtol[0] = ssRtol[1] = 1.0e-6;
  op->ssAtol=ssAtol;
  op->ssRtol=ssRtol;
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

  indR.linCmt = flags[RxMvFlag_linCmt];

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
  for (int i = indR.ndoses; i< indR.n_all_times; ++i) {
    indR.idose[i] = 0;
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
  tmpD = (double*)malloc(EVID_EXTRA_SIZE* sizeof(double));
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
  tmpD = (double*)calloc(indR.n_all_times*(rx->linNcmt+rx->linKa+1), sizeof(double));
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
  indR.timeThread = tmpD + indR.n_all_times*(rx->linNcmt+rx->linKa);
  indR.ixds = indR.idx=0;
  rx_solving_options_ind* ind = &indR;
  rx->subjects =  ind;
  op->nlin = rx->linNcmt+rx->linKa;
  iniSubject(0, 0, ind, op, rx, u_inis_lincmt);
  double *yp;
  SEXP CcSxp = PROTECT(Rf_allocVector(REALSXP, indR.n_all_times)); pro++;
  SEXP TimeSxp = PROTECT(Rf_allocVector(REALSXP, indR.n_all_times)); pro++;
  SEXP EvidSxp = PROTECT(Rf_allocVector(INTSXP, indR.n_all_times)); pro++;
  SEXP DoseSxp = PROTECT(Rf_allocVector(REALSXP, indR.n_all_times)); pro++;
  SEXP IiSxp   = PROTECT(Rf_allocVector(REALSXP, indR.n_all_times)); pro++;
  double *Cc = REAL(CcSxp);
  double *time = REAL(TimeSxp);
  int *evidOut = INTEGER(EvidSxp);
  double *doseOut = REAL(DoseSxp);
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
  bool isOral  = rx->linKa;
  int linCmt = 0; // states before linCmt model, in this case 0
  for(int i=0; i < indR.n_all_times; ++i) {
    ind->idx=i;
    xout = getTime__(ind->ix[i], ind, 0);
    if (isOral) {
      Cc[i] = linCmtCompA(rx, 0, xout, linCmt, rx->linNcmt, trans,
                          p1[i], v1[i], p2[i], p3[i], p4[i], p5[i],
                          lagdepot[i], fdepot[i], ratedepot[i], durdepot[i],
                          ka[i], lagcentral[i], fcentral[i],  ratecentral[i], durcentral[i]);
    } else {
      Cc[i] = linCmtCompA(rx, 0, xout, linCmt, rx->linNcmt, trans,
                          p1[i], v1[i], p2[i], p3[i], p4[i], p5[i],
                         lagcentral[i], fcentral[i], ratecentral[i], durcentral[i],
                          0.0, 0.0, 1.0,  0.0, 0.0);
    }
    Cc[i]= Cc[i]*sm;
    time[i] = getTime__(ind->ix[i], ind, 0);
    evidOut[i] = getEvid(ind, ind->ix[i]);
    doseOut[i] = getDose(ind, ind->ix[i]);
  }
  SEXP dfNames = PROTECT(Rf_allocVector(STRSXP, 19)); pro++;
  SEXP dfVals = PROTECT(Rf_allocVector(VECSXP, 19)); pro++;
  SEXP rnVals = PROTECT(Rf_allocVector(INTSXP, 2)); pro++;
  int *rnI = INTEGER(rnVals);

  rnI[0] = NA_INTEGER;
  rnI[1] = -indR.n_all_times;

  //"TIME"
  SET_STRING_ELT(dfNames, 0, Rf_mkChar("TIME"));
  SET_VECTOR_ELT(dfVals, 0, TimeSxp);

  // "EVID"
  SET_STRING_ELT(dfNames, 1, Rf_mkChar("EVID"));
  SET_VECTOR_ELT(dfVals,  1, EvidSxp);

  // "AMT"
  SET_STRING_ELT(dfNames, 2, Rf_mkChar("AMT"));
  SET_VECTOR_ELT(dfVals,  2, DoseSxp);

  // "II"
  SET_STRING_ELT(dfNames, 3, Rf_mkChar("II"));
  SET_VECTOR_ELT(dfVals,  3, IiSxp);

  // "Cc"
  SET_STRING_ELT(dfNames, 3, Rf_mkChar("Cc"));
  SET_VECTOR_ELT(dfVals,  3, CcSxp);

  SET_STRING_ELT(dfNames, 4, Rf_mkChar("p1"));
  SET_VECTOR_ELT(dfVals,  4, VECTOR_ELT(par, 0));

  SET_STRING_ELT(dfNames, 5, Rf_mkChar("v1"));
  SET_VECTOR_ELT(dfVals,  5, VECTOR_ELT(par, 1));

  SET_STRING_ELT(dfNames, 6, Rf_mkChar("p2"));
  SET_VECTOR_ELT(dfVals,  6, VECTOR_ELT(par, 2));

  SET_STRING_ELT(dfNames, 7, Rf_mkChar("p3"));
  SET_VECTOR_ELT(dfVals,  7, VECTOR_ELT(par, 3));

  SET_STRING_ELT(dfNames, 8, Rf_mkChar("p4"));
  SET_VECTOR_ELT(dfVals,  8, VECTOR_ELT(par, 4));

  SET_STRING_ELT(dfNames, 9, Rf_mkChar("p5"));
  SET_VECTOR_ELT(dfVals,  9, VECTOR_ELT(par, 5));

  SET_STRING_ELT(dfNames, 10, Rf_mkChar("lagdepot"));
  SET_VECTOR_ELT(dfVals,  10, VECTOR_ELT(par, 6));

  SET_STRING_ELT(dfNames, 11, Rf_mkChar("fdepot"));
  SET_VECTOR_ELT(dfVals,  11, VECTOR_ELT(par, 7));

  SET_STRING_ELT(dfNames, 12, Rf_mkChar("ratedepot"));
  SET_VECTOR_ELT(dfVals,  12, VECTOR_ELT(par, 8));

  SET_STRING_ELT(dfNames, 13, Rf_mkChar("durdepot"));
  SET_VECTOR_ELT(dfVals,  13, VECTOR_ELT(par, 9));

  SET_STRING_ELT(dfNames, 14, Rf_mkChar("ka"));
  SET_VECTOR_ELT(dfVals,  14, VECTOR_ELT(par, 10));

  SET_STRING_ELT(dfNames, 15, Rf_mkChar("lagcentral"));
  SET_VECTOR_ELT(dfVals,  15, VECTOR_ELT(par, 11));

  SET_STRING_ELT(dfNames, 16, Rf_mkChar("fcentral"));
  SET_VECTOR_ELT(dfVals,  16, VECTOR_ELT(par, 12));

  SET_STRING_ELT(dfNames, 17, Rf_mkChar("ratecentral"));
  SET_VECTOR_ELT(dfVals,  17, VECTOR_ELT(par, 13));

  SET_STRING_ELT(dfNames, 18, Rf_mkChar("durcentral"));
  SET_VECTOR_ELT(dfVals,  18, VECTOR_ELT(par, 14));

  // R_NameSymbol
  SEXP cls = PROTECT(Rf_allocVector(STRSXP, 1)); pro++;
  SET_STRING_ELT(cls, 0, Rf_mkChar("data.frame"));
  Rf_setAttrib(dfVals, R_NamesSymbol, dfNames);
  Rf_setAttrib(dfVals, R_RowNamesSymbol, rnVals);
  Rf_setAttrib(dfVals, R_ClassSymbol, cls);
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
  return dfVals;
}
