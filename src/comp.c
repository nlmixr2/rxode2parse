#define USE_FC_LEN_T
#define STRICT_R_HEADERS
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>   /* dj: import intptr_t */
//#include "ode.h"
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
