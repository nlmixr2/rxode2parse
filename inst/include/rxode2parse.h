#pragma once
#define STRICT_R_HEADERS
#ifndef __rxode2parse_H__
#define __rxode2parse_H__
#define rxLlikSaveSize 9


#include <R.h>
#include <stdbool.h>

#include <float.h>
#include <stdio.h>
#include <stdarg.h>

#include "rxode2parse_control.h"
#include <stdint.h>    // for uint64_t rather than unsigned long long

#ifdef _isrxode2parse_
#define _linCmtParse _rxode2parse_linCmtParse
#define _rxode2_linCmtGen _rxode2parse_linCmtGen
#define rc_buf_read _rxode2parse_rc_buf_read
#define sIniTo _rxode2parse_sIniTo
#define sFree _rxode2parse_sFree
#define sFreeIni _rxode2parse_sFreeIni
#define sAppendN _rxode2parse_sAppendN
#define sAppend _rxode2parse_sAppend
#define sPrint _rxode2parse_sPrint
#define lineIni _rxode2parse_lineIni
#define lineFree _rxode2parse_lineFree
#define addLine _rxode2parse_addLine
#define curLineProp _rxode2parse_curLineProp
#define curLineType _rxode2parse_curLineType
#define doDot _rxode2parse_doDot
#define doDot2 _rxode2parse_doDot2
#define _setSilentErr _rxode2parse__setSilentErr
#define _isRstudio2 _rxode2parse_isRstudio2
#define setSilentErr _rxode2parse_setSilentErr
#define setRstudioPrint _rxode2parse_setRstudioPrint
#define getSilentErr _rxode2parse_getSilentErr
#define getRstudioPrint _rxode2parse_getRstudioPrint
#define RSprintf _rxode2parse_RSprintf
#define parseFree _rxode2parse_parseFree
#define parseFreeLast _rxode2parse_parseFreeLast
#endif

typedef struct sbuf {
  char *s;        /* curr print buffer */
  int sN;
  int o;                        /* offset of print buffer */
} sbuf;
  
typedef struct vLines {
  char *s;
  int sN;
  int o;
  int n;
  int nL;
  char **line;
  int *lProp;
  int *lType;
  int *os;
} vLines;


typedef struct {
  // These options should not change based on an individual solve
  int badSolve;
  int naTime;
  double ATOL; //absolute error
  double RTOL; //relative error
  double H0;
  double HMIN;
  int mxstep;
  int MXORDN;
  int MXORDS;
  //
  int nlhs;
  int neq;
  int stiff;
  int ncov;
  char modNamePtr[1000];
  int *par_cov;
  double *inits;
  double *scale;
  bool do_par_cov;
  // approx fun options
  double f1;
  double f2;
  int kind;
  int is_locf;
  int cores;
  int doesRandom;
  int extraCmt;
  double hmax2; // Determined by diff
  double *rtol2;
  double *atol2;
  double *ssRtol;
  double *ssAtol;
  int *indLin;
  int indLinN;
  double indLinPhiTol;
  int indLinPhiM;
  int indLinMatExpType;
  int indLinMatExpOrder;
  int nDisplayProgress;
  int ncoresRV;
  int isChol;
  int nsvar;
  int abort;
  int minSS;
  int maxSS;
  int doIndLin;
  int strictSS;
  double infSSstep;
  int mxhnil;
  double hmxi;
  int nlin;
  int nlin2;
  int nlinR;
  int linBflag;
  bool cTlag;
  double hTlag;
  bool cF;
  double hF;
  bool cRate;
  double hRate;
  bool cDur;
  double hDur;
  bool cTlag2;
  double hTlag2;
  bool cF2;
  double hF2;
  bool cRate2;
  double hRate2;
  bool cDur2;
  double hDur2;
  int nLlik;
} rx_solving_options;


typedef struct {
  double bT;
  int *slvr_counter;
  int *dadt_counter;
  int *jac_counter;
  double *InfusionRate;
  int *BadDose;
  int nBadDose;
  double HMAX; // Determined by diff
  double tlast;
  double curDose;
  int dosenum;
  double tfirst;
  double *tlastS;
  double *curDoseS;
  double *tfirstS;
  double podo;
  double *podoS;
  double *par_ptr; // both time changing and time invariant
  double *dose;
  double *ii;
  double *solve;
  double *mtime;
  double *solveSave;
  double *solveLast;
  double *solveLast2;
  double *lhs;
  int  *evid;
  int *rc;
  double *cov_ptr;
  int *cov_sample;
  // a b
  // 1 4
  // 2 5
  // 3 6
  int n_all_times;
  int nevid2;
  int ixds;
  int ndoses;
  double *all_times;
  int *ix;
  double *dv;
  double *limit;
  int *cens;
  int *idose;
  int *on;
  int idosen;
  int id;
  int idReal;
  int sim;
  int idx;
  double ylow;
  double yhigh;
  double logitHi;
  double logitLow;
  double lambda;
  double yj;
  // Saved info
  int wh;
  int wh100;
  int cmt;
  int whI;
  int wh0;
  int doSS;
  int allCovWarn;
  int wrongSSDur;
  int _newind;
  int _rxFlag;
  int err;
  int solved;
  double *linCmtAdvan;
  double *linCmtRate;
  int linCmt;
  int linCmtAdvanSetup;
  int cacheME;
  int inLhs;
  // Cache alag
  double *alag;
  // Cache F
  double *cF;
  // Cache rate;
  double *cRate;
  // Cache duration
  double *cDur;
  double solveTime;
  double curShift;
  double *simIni;
  int isIni;
  int _update_par_ptr_in;
  int badIni;
  double *llikSave;
} rx_solving_options_ind;

typedef struct {
  rx_solving_options_ind *subjects;
  rx_solving_options *op;
  int nsub;
  int nsim;
  int neta;
  int neps;
  int nIndSim;
  int simflg;
  int nall;
  int nevid9;
  int nobs;
  int nobs2;
  int nr;
  int add_cov;
  int matrix;
  int needSort;
  int nMtime;
  double stateTrimU;
  double stateTrimL;
  int *stateIgnore;
  int nCov0;
  int *cov0;
  int nKeepF;
  int istateReset;
  int cens;
  int limit;
  int safeZero;
  int sumType;
  int prodType;
  int sensType;
  vLines factors;
  vLines factorNames;
  int factorNs[500];
  int hasFactors;
  // For forder
  uint64_t minD;
  uint64_t maxD;
  int maxAllTimes;
  uint8_t ***keys;// = NULL; keys per thread
  int *TMP;
  int *ordId;
  uint8_t *UGRP;
  int *nradix;
  double *ypNA;
  bool sample;
  int *par_sample;
  double maxShift;
  int linKa;
  int linNcmt;
  int maxwhile;
  int whileexit;
  int *svar;
  int *ovar;
  int hasEvid2;
  int useStdPow;
} rx_solve;


static inline void sNull(sbuf *sbb) {
  sbb->s = NULL;
  sbb->sN=0;
  sbb->o=0;
}

static inline void lineNull(vLines *sbb) {
  sbb->s = NULL;
  sbb->lProp = NULL;
  sbb->lType = NULL;
  sbb->line = NULL;
  sbb->os = NULL;
  sbb->sN = 0;
  sbb->nL = 0;
  sbb->n  = 0;
  sbb->o  = 0;
}



#endif
