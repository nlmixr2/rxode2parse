
#ifndef __RXODE2PARSESTRUCT_H__
#define __RXODE2PARSESTRUCT_H__

#if defined(__cplusplus)
extern "C" {
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

#define rxode2naTimeInputIgnore 1
#define rxode2naTimeInputWarn   2
#define rxode2naTimeInputError  3

  typedef struct {
    // These options should not change based on an individual solve
    int badSolve;
    int naTime;
    int naTimeInput;
    int naTimeInputWarn;
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

  static inline void iniSolvingOptions(rx_solving_options *op) {
    op->badSolve = 0;
    op->naTime = 0;
    op->ATOL = 1e-8; //absolute error
    op->RTOL = 1e-8; //relative error
    op->H0  = 0;
    op->HMIN = 0;
    op->mxstep = 70000;
    op->MXORDN =12;
    op->MXORDS = 5;
    //
    op->nlhs = 0;
    op->neq = 0;
    op->stiff = 0;
    op->ncov = 0;
    op->par_cov = NULL;
    op->inits = NULL;
    op->scale = NULL;
    op->do_par_cov=false;
    // approx fun options
    op->f1   = 0.0;
    op->f2   = 1.0;
    op->kind = 0;
    op->is_locf = 1;
    op->cores = 0;
    op->extraCmt = 0;
    op->hmax2 = 0; // Determined by diff
    op->rtol2 = NULL;
    op->atol2 = NULL;
    op->ssRtol = NULL;
    op->ssAtol = NULL;
    op->indLin = NULL;
    op->indLinN = 0;
    op->indLinPhiTol = 1e-7;
    op->indLinPhiM = 0;
    op->indLinMatExpType = 2;
    op->indLinMatExpOrder = 6;
    op->nDisplayProgress = 10000;
    op->isChol = 0;
    op->nsvar = 0;
    op->abort = 0;
    op->minSS = 10;
    op->maxSS = 1000;
    op->doIndLin  = 0;
    op->strictSS = 1;
    op->infSSstep = 12;
    op->ncoresRV = 1;
    op->mxhnil = 0;
    op->hmxi = 0.0;
    op->nlin = 0;
    op->nlin2 = 0;
    op->nlinR = 0;
    op->linBflag = 0;
    op->cTlag = false;
    op->hTlag = 0;
    op->cF = false;
    op->hF = 0;
    op->cRate = false;
    op->hRate = 0;
    op->cDur = false;
    op->hDur = 0;
    op->cTlag2 = false;
    op->hTlag2 = 0;
    op->cF2 = false;
    op->hF2 = 0;
    op->cRate2 = false;
    op->hRate2 = 0;
    op->cDur2 = false;
    op->hDur2 = 0;
    op->nLlik = 0;
  }

  typedef struct {
    int idx;
    double p1;
    double v1;
    double p2;
    double p3;
    double p4;
    double p5;
    double d_tlag;
    double d_F;
    double d_rate1;
    double d_dur1;
    // Oral parameters
    double d_ka;
    double d_tlag2;
    double d_F2;
    double d_rate2;
    double d_dur2;
  } rx_solving_linCmt_dose_single;

  typedef struct {
    rx_solving_linCmt_dose_single *pars;
    int n;
    int nAlloc;
  } rx_solving_linCmt_doses;

  #define iniLinCmtAlloc 2

  static inline void iniLinCmtDoses(rx_solving_linCmt_doses *doses) {
    doses->n = 0;
    doses->pars = (rx_solving_linCmt_dose_single*)(malloc(iniLinCmtAlloc*sizeof(rx_solving_linCmt_dose_single)));
    if (doses->pars == NULL) {
      Rf_error("ran out of memory");
    }
    doses->nAlloc = iniLinCmtAlloc;
  }

  static inline rx_solving_linCmt_dose_single* getLinCmtDoses(rx_solving_linCmt_doses *doses, int idx) {
    rx_solving_linCmt_dose_single *single;
    for (int i = 0; i < doses->n; ++i) {
      single = &(doses->pars[i]);
      if (single->idx == idx) {
        return single;
      }
    }
    return NULL;
  }

  static inline rx_solving_linCmt_dose_single* pushLinCmtDose(rx_solving_linCmt_doses *doses, int idx,
                                                              double p1, double v1,
                                                              double p2, double p3,
                                                              double p4, double p5,
                                                              double d_tlag, double d_F, double d_rate1,
                                                              double d_dur1,
                                                              // Oral parameters
                                                              double d_ka, double d_tlag2,
                                                              double d_F2,  double d_rate2, double d_dur2) {
    rx_solving_linCmt_dose_single* single = getLinCmtDoses(doses, idx);
    if (single == NULL) {
      if (doses->nAlloc <= doses->n + 1) {
        int size = doses->n + iniLinCmtAlloc + 1;
        rx_solving_linCmt_dose_single* n = (rx_solving_linCmt_dose_single*)(realloc(doses->pars, size*sizeof(rx_solving_linCmt_dose_single)));
        if (n == NULL) {
          Rf_error("ran out of memory");
        }
        doses->pars = n;
        doses->nAlloc = size;
      }
      single = &(doses->pars[doses->n]);
      doses->n++;
    }
    single->idx = idx;
    single->p1 = p1;
    single->v1 = v1;
    single->p2 = p2;
    single->p3 = p3;
    single->p4 = p4;
    single->p4 = p5;
    single->d_tlag = d_tlag;
    single->d_F = d_F;
    single->d_rate1 = d_rate1;
    single->d_dur1 = d_dur1;
    single->d_ka = d_ka;
    single->d_tlag2 = d_tlag2;
    single->d_F2 = d_tlag2;
    single->d_rate2 = d_rate2;
    single->d_dur2 = d_dur2;
    return single;
  }

  static inline void resetLinCmtDose(rx_solving_linCmt_doses *doses) {
    doses->n = 0;
  }

  static inline void freeLinCmtDoses(rx_solving_linCmt_doses *doses) {
    doses->n = 0;
    free(doses->pars);
    doses->nAlloc = 0;
  }

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
    int id;
    int solveid;
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
    int linCmt;
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
    // Add pointers for drifting atol/rtol values during optimization
    double *rtol2;
    double *atol2;
    double *ssRtol;
    double *ssAtol;
    // ignored and pending doses
    int *ignoredDoses;
    int *ignoredDosesN;
    int *ignoredDosesAllocN;
    int *pendingDoses;
    int *pendingDosesN;
    int *pendingDosesAllocN;
    // extra doses
    int *extraDoseTimeIdx;
    int *extraDoseN;
    int *extraDoseAllocN;
    double *extraDoseTime;
    int *extraDoseEvid;
    double *extraDoseDose;
    double extraDoseNewXout;
    int idxExtra; // extra idx
    int extraSorted; // extra sorted?
    //double *extraDoseIi; // ii doses unsupported
    bool lastIsSs2;
    double *timeThread;
    // linear compartmental solve
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
    int maxAllTimes;
    int *ordId;
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
    bool ss2cancelAllPending;
  } rx_solve;

  static inline void iniSolvingRx(rx_solve *rx) {
    rx->hasEvid2 = 0;
    rx->nsub = 0;
    rx->nsim = 0;
    rx->neta=0;
    rx->neps=0;
    rx->nIndSim=0;
    rx->simflg = 0;
    rx->nall = 0;
    rx->nevid9 = 0;
    rx->nobs = 0;
    rx->nobs2 = 0;
    rx->nr = 0;
    rx->add_cov = 0;
    rx->matrix = 0;
    rx->needSort = 0;
    rx->nMtime = 0;
    rx->stateTrimU = R_PosInf;
    rx->stateTrimL = R_NegInf;
    rx->stateIgnore = NULL;
    rx->nCov0 = 0;
    rx->cov0 = NULL;
    rx->nKeepF = 0;
    rx->istateReset=1;
    rx->cens = 0;
    rx->limit = 0;
    rx->safeZero = 1;
    rx->useStdPow = 0;
    rx->ss2cancelAllPending = false;
    rx->sumType = 1; // pairwise
    rx->prodType = 1; // long double
    rx->sensType = 4; // advan
    rx->hasFactors = 0;
    rx->ordId = NULL;
    rx->ypNA = NULL;
    rx->sample = false;
    rx->par_sample = NULL;
    rx->maxShift = 0.0;
    rx->linKa  = 0;
    rx->linNcmt = 0;
    rx->maxwhile = 100000;
    rx->whileexit= 0;
  }

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

  typedef double (*t_F)(int _cSub,  int _cmt, double _amt, double t, double *y);
  typedef double (*t_LAG)(int _cSub,  int _cmt, double t);
  typedef double (*t_RATE)(int _cSub,  int _cmt, double _amt, double t);
  typedef double (*t_DUR)(int _cSub,  int _cmt, double _amt, double t);

  typedef void (*t_calc_mtime)(int cSub, double *mtime);

  typedef void (*t_ME)(int _cSub, double _t, double t, double *_mat, const double *__zzStateVar__);
  typedef void (*t_IndF)(int _cSub, double _t, double t, double *_mat);

  typedef double (*t_getTime)(int idx, rx_solving_options_ind *ind);
  typedef int (*t_locateTimeIndex)(double obs_time,  rx_solving_options_ind *ind);
  typedef int (*t_handle_evidL)(int evid, double *yp, double xout, rx_solving_options_ind *ind) ;

#ifdef _isrxode2parse_

  extern rx_solving_options _rxode2parse_op_global;
  extern rx_solve _rxode2parse_rx_global;
  extern t_handle_evidL _rxode2parse_handle_evidL;
  extern t_getTime _rxode2parse_getTime;
#define op_global _rxode2parse_op_global
#define rx_global _rxode2parse_rx_global
#define AMT _rxode2parse_AMT
#define LAG _rxode2parse_LAG
#define RATE _rxode2parse_RATE
#define DUR _rxode2parse_DUR
#define calc_mtime _rxode2parse_calc_mtime
#define getTime _rxode2parse_getTime
#define _locateTimeIndex _rxode2parse_locateTimeIndex

#endif

#if defined(__cplusplus)
}
#endif
#endif
