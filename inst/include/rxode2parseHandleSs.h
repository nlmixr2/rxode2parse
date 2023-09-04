// -*- mode: c++; c-basic-offset: 2; tab-width: 2; indent-tabs-mode: nil; -*-
#ifndef __RXODE2PARSEHANDLSS_H___
#define __RXODE2PARSEHANDLSS_H___
#include "rxode2parse.h"
#include "rxode2parseHandleEvid.h"
#include "rxode2parseGetTime.h"
#define isSameTimeOp(xout, xp) (op->stiff == 0 ? isSameTimeDop(xout, xp) : isSameTime(xout, xp))
#define max2( a , b )  ( (a) > (b) ? (a) : (b) )
#if defined(__cplusplus)
#include "timsort.h"
#define SORT gfx::timsort
static inline int handleExtraDose(double *yp,
                                  double xout, double xp,
                                  int *i,
                                  int *istate,
                                  rx_solving_options *op,
                                  rx_solving_options_ind *ind,
                                  t_update_inis u_inis,
                                  void *ctx) {
    if (ind->extraDoseN[0] > ind->idxExtra) {
      if (ind->extraSorted == 0) {
        // do sort
        SORT(ind->extraDoseTimeIdx + ind->idxExtra, ind->extraDoseTimeIdx + ind->extraDoseN[0],
             [ind](int a, int b){
               double timea = ind->extraDoseTime[a],
                 timeb = ind->extraDoseTime[b];
               if (timea == timeb) {
                 int evida = ind->extraDoseEvid[a],
                   evidb = ind->extraDoseEvid[b];
                 if (evida == evidb){
                   return a < b;
                 }
                 return evida < evidb;
               }
               return timea < timeb;
             });
        ind->extraSorted=1;
        ind->idxExtra=0;
      }
      // Use "real" xout for handle_evid functions.
      int idx = ind->idx;
      int ixds = ind->ixds;
      int trueIdx = ind->extraDoseTimeIdx[ind->idxExtra];
      ind->idx = -1-trueIdx;
      double time = getAllTimes(ind, ind->idx);
      while (!isSameTimeOp(time, xp) && time < xp && ind->idxExtra < ind->extraDoseN[0]) {
        ind->idxExtra++;
        trueIdx = ind->extraDoseTimeIdx[ind->idxExtra];
        ind->idx = -1-trueIdx;
        time = getAllTimes(ind, ind->idx);
      }
      if ((isSameTimeOp(time, xp) || time > xp) && (isSameTimeOp(time, xout) || time <= xout)) {
        bool ignore = true;
        while (ignore && time <= xout) {
          ignore=false;
          for (int i = 0; i < ind->ignoredDosesN[0]; ++i) {
            int curIdx = ind->ignoredDoses[i];
            if (curIdx < 0 && -1-curIdx == trueIdx) {
              ignore = true;
              break;
            }
          }
          if (ignore) {
            ind->idxExtra++;
            if (ind->idxExtra < ind->extraDoseN[0]) {
              trueIdx = ind->extraDoseTimeIdx[ind->idxExtra];
              ind->idx = -1-trueIdx;
              time = getAllTimes(ind, ind->idx);
            } else {
              ind->idxExtra--;
              break;
            }
          } else {
            break;
          }
        }
        if (ignore) {
          ind->idx = idx;
          ind->ixds = ixds;
          return 0;
        } else {
          ind->extraDoseNewXout = time;
          ind->idx = idx;
          ind->ixds = ixds;
          // REprintf("time: %f; xp: %f; xout: %f; handleExtra\n", time, xp, xout);
          return 1;
        }
      }
      ind->idx = idx;
      ind->ixds = ixds;
      return 0;
    }
    return 0;
  }
#undef SORT

extern "C" {
#endif

#define badSolveExit(i) for (int j = op->neq*(ind->n_all_times); j--;){ \
    ind->solve[j] = NA_REAL;                                            \
  }                                                                     \
  op->badSolve = 1;                                                     \
  i = ind->n_all_times-1; // Get out of here!

  typedef void (*solveWith1Pt_fn)(double *yp,
                                  double xout, double xp,
                                  int *i,
                                  int *istate,
                                  rx_solving_options *op,
                                  rx_solving_options_ind *ind,
                                  t_update_inis u_inis,
                                  void *ctx);

  typedef void (*handleSSbolus_fn)(double *yp,
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
                                   solveWith1Pt_fn solveWith1Pt);

  typedef void (*solveSSinf_fn)(double *yp,
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
                                solveWith1Pt_fn solveWith1Pt);

  typedef void (*solveSSinfLargeDur_fn)(double *yp,
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
                                        int *numDoseInf,
                                        double *offTime,
                                        double *addTime,
                                        int *canBreak,
                                        solveWith1Pt_fn solveWith1Pt);

  typedef void (*handleSSinf8_fn)(double *yp,
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
                                  solveWith1Pt_fn solveWith1Pt);

  typedef void (*updateExtraDoseGlobals_fn)(rx_solving_options_ind* ind);

  static inline const char *getId(int id) {
    rx_solve *rx = &rx_global;
    int curLen=  rx->factorNs[0];
    const char *unknownId = "Unknown";
    if (id < 0) {
      return unknownId; // Bad value
    }
    if (id < curLen){
      if (id >= rx->factors.n) {
        return unknownId;
      }
      return rx->factors.line[id];
    } else {
      return unknownId;
    }
  }

  static inline void printErr(int err, int id) {
    REprintf("Recovered solving errors for internal ID %s (%d):\n", getId(id), err);
    if (err & rxErrCorruptETSort){
      REprintf("  Corrupted event table during sort (1)\n");
    }
    if (err & rxErrRate0){
      REprintf("  Rate is zero/negative\n");
    }
    if (err & rxErrModelRateAbsent){
      REprintf("  Modeled rate requested in event table, but not in model; use 'rate(cmt) ='\n");
    }
    if (err & rxErrCorruptETSort2){
      REprintf("  Corrupted event table during sort (2)\n");
    }
    if (err & rxErrDurNeg0){
      REprintf("  Duration is zero/negative\n");
    }
    if (err & rxErrModelDurAbsent){
      REprintf("  Modeled duration requested in event table, but not in model; use 'dur(cmt) ='\n");
    }
    if (err & rxErrModelData686){
      REprintf("  Data error 686\n");
    }
    if (err & rxErrModelDataNeg6){
      REprintf("  Data Error -6\n");
    }
    if (err & rxErrModelDataErr8){
      REprintf("  Data Error 8\n");
    }
    if (err & rxErrModelDataErr886){
      REprintf("  Data error 886\n");
    }
    if (err & rxErrModelDataErr797){
      REprintf("  Data error 797\n");
    }
    if (err & rxErrModelDataNeg7){
      REprintf("  Data Error -7\n");
    }
    if (err & rxErrModelDataErr9){
      REprintf("  Data Error 9\n");
    }
    if (err & rxErrModelDataErr997){
      REprintf("  Data error 997\n");
    }
    if (err & rxErrCorruptETSort3){
      REprintf("  Corrupted event table during sort (3)\n");
    }
    if (err & rxErrCorruptET) {
      REprintf("  Corrupted event table\n");
    }
    if (err & rxErrCorruptET2){
      REprintf("  Corrupted events\n");
    }
    if (err & rxErrNegCmt){
      REprintf("  Supplied an invalid EVID\n");
    }
    if (err & rxErrSync){
      REprintf("  Corrupted event table (during sync)\n");
    }
    if (err & rxErrSync2){
      REprintf("  Corrupted event table (end of sync)\n");
    }
    if (err & rxErrModeledFss2){
      REprintf("  SS=2 & Modeled F does not work\n");
    }

    if (err & rxErrModeledFss2n2){
      REprintf("  SS=2 & Modeled F does not work\n");
    }
    if (err & rxErrModeledFss2n3){
      REprintf("  SS=2 & Modeled F does not work\n");
    }
    if (err & rxErrRate02){
      REprintf(" Rate is zero/negative\n");
    }
  }

  static inline void handleSSGen(int *neq,
                                 int *BadDose,
                                 double *InfusionRate,
                                 double *dose,
                                 double *yp,
                                 double xout, double xp, int id,
                                 int *i, int nx,
                                 int *istate,
                                 rx_solving_options *op,
                                 rx_solving_options_ind *ind,
                                 t_update_inis u_inis,
                                 void *ctx,
                                 solveWith1Pt_fn solveWith1Pt,
                                 handleSSbolus_fn handleSSbolus,
                                 solveSSinf_fn solveSSinf,
                                 solveSSinfLargeDur_fn solveSSinfLargeDur,
                                 handleSSinf8_fn handleSSinf8,
                                 updateExtraDoseGlobals_fn updateExtraDoseGlobalsI) {
    rx_solve *rx = &rx_global;
    int j;
    int doSS2=0;
    int doSSinf=0;
    int maxSS = op->maxSS;
    int minSS = op->minSS;
    int isSsLag = ind->wh0 == EVID0_SS20 || ind->wh0 == EVID0_SS0;
    bool skipDosingEvent = false, isRateDose = false;
    bool isModeled = ind->whI == EVIDF_MODEL_DUR_ON ||
      ind->whI == EVIDF_MODEL_RATE_ON;
    double curIi = ind->ixds == 0 ? 0.0 : getIiNumber(ind, ind->ixds-1);
    if (((ind->wh0 == EVID0_SS2  || isSsLag ||
          ind->wh0 == EVID0_SS) &&
         curIi > 0) || ind->wh0 == EVID0_SSINF) {
      int oldIdx = ind->idx, oldIxds = ind->ixds;
      int ignoreDoses[4];
      ignoreDoses[0] = ignoreDoses[1] = ignoreDoses[2] = ignoreDoses[3] = -1;
      int nIgnoredDoses = 0;
      if (isSsLag) {
        maxSS--; minSS--;
      }
      ind->doSS=1;
      ind->ixds--; // This dose stays in place; Reverse dose
      if (ind->wh0 == EVID0_SS2 || ind->wh0 == EVID0_SS20){
        doSS2=1;
      } else if (ind->wh0 == EVID0_SSINF) {
        doSSinf=1;
      }
      double dur = 0, dur2=0, rateOn=0.0, rateOff = 0.0;
      int infBixds =0, infBixds2=0, infEixds = 0, infFixds = 0,
        infSixds = 0,
        ei=0, wh, cmt, wh100, whI, wh0, oldI,
        bi = *i, fi = *i;
      if (doSSinf){
      } else if (ind->whI == EVIDF_INF_RATE || ind->whI == EVIDF_INF_DUR) {
        if (getDose(ind, ind->idose[ind->ixds]) < 0) return;
        oldI = ind->whI;
        isRateDose=true;
        infBixds = infBixds2 = infFixds = ind->ixds;
        // Find the next fixed length infusion that is turned off.
        if (isSsLag) {
          ignoreDoses[nIgnoredDoses++] = infBixds;
          int bEvid = getEvid(ind, ind->idose[infBixds2]);
          double bIi = getIiNumber(ind, infBixds2);
          double cIi = bIi;
          double bDose = getDoseNumber(ind,ind->ixds);
          double cDose = bDose;
          getWh(bEvid, &wh, &cmt, &wh100, &whI, &wh0);
          // This structure is
          // TIME  EVID AMT II
          //    0 10209  10 24
          //    0 10208 -10 24
          //    0 10201  10  0
          // ...
          //   10 10201 -10  0
          //
          // or with ss=2
          //
          // TIME  EVID AMT II
          //    0 10219  10 24 (note the 19 flag)
          //    0 10208 -10 24
          //    0 10201  10  0
          // ...
          //   10 10201 -10  0
          int evid8 = bEvid - ind->wh0 + EVID0_INFRM;
          bool foundEvid8 = false;
          int evid1 = bEvid - ind->wh0 + EVID0_REGULAR;
          bool foundEvid1 = false;
          while (!foundEvid1 && infBixds2 < ind->ndoses) {
            infBixds2++;
            if (infBixds2 == ind->ndoses) {
              infBixds = infBixds2 = infEixds = -1;
              break;
            } else {
              bEvid = getEvid(ind, ind->idose[infBixds2]);
              cDose = getDoseNumber(ind, infBixds2);
              cIi = getIiNumber(ind, infBixds2);
              if (!foundEvid8 && bEvid == evid8 && cDose == -bDose && cIi == bIi) {
                foundEvid8 = true;
                ignoreDoses[nIgnoredDoses++] = infBixds2;
              } else if (!foundEvid1 && bEvid == evid1 && cDose == bDose && cIi == 0.0) {
                foundEvid1 = true;
                ignoreDoses[nIgnoredDoses++] = infBixds2;
                break;
              }
            }
          }
          if (infEixds != -1) {
            handleInfusionGetEndOfInfusionIndex(infBixds2,
                                                &infEixds, rx, op, ind);
          }
          if (infEixds == -1 && infBixds2 != -1) {
            ind->wrongSSDur=1;
            // // Bad Solve => NA
            badSolveExit(*i);
            return;
          } else {
            ignoreDoses[nIgnoredDoses++] = infEixds;
            double f = 1.0;
            if (ind->whI == EVIDF_INF_RATE) {
              f = getAmt(ind, ind->id, ind->cmt, 1.0, getAllTimes(ind, ind->idose[infBixds2]), yp);
              rateOn = getDose(ind, ind->idose[infBixds2]);
              rateOff = getDose(ind, ind->idose[infEixds]);
            } else {
              f = 1.0;
              rateOn = getAmt(ind, ind->id, ind->cmt, getDose(ind, ind->idose[infBixds2]), getAllTimes(ind, ind->idose[infBixds2]), yp);
              rateOff = getAmt(ind, ind->id, ind->cmt, getDose(ind, ind->idose[infEixds]), getAllTimes(ind, ind->idose[infEixds]), yp);
            }
            // REprintf("get dur isSsLag: %d\n", isSsLag);
            dur = getAllTimes(ind, ind->idose[infEixds]);// -
            // REprintf("\ttime infEixds: %f %d\n", dur, infEixds);
            dur -= getAllTimes(ind, ind->idose[infBixds2]);
            // REprintf("\ttime infBixds: %f %d\n", getAllTimes(ind, ind->idose[infBixds2]), infBixds2);
            dur *= f;
            // REprintf("\tdur: %f\n", dur);
            dur2 = curIi - dur;
            // REprintf("\tdur2: %f\n", dur2);
            infSixds = infBixds;
          }
        } else {
          // This is the infusion structure:
          // TIME  EVID AMT II
          //    0 10201  10  0
          // ...
          //   10 10201 -10  0
          infBixds=infBixds2=infEixds=infFixds=ind->ixds;
          ignoreDoses[nIgnoredDoses++] = infBixds;
          handleInfusionGetEndOfInfusionIndex(infBixds, &infEixds, rx, op, ind);
          if (infEixds == -1) {
            ind->wrongSSDur=1;
            // // Bad Solve => NA
            badSolveExit(*i);
            return;
          } else {
            ignoreDoses[nIgnoredDoses++] = infEixds;
            double f = 1.0;
            if (ind->whI == EVIDF_INF_RATE) {
              f = getAmt(ind, ind->id, ind->cmt, 1.0, getAllTimes(ind, ind->idose[infBixds2]), yp);
              rateOn = getDose(ind, ind->idose[infBixds2]);
              rateOff = getDose(ind, ind->idose[infEixds]);
            } else {
              f = 1.0;
              rateOn = getAmt(ind, ind->id, ind->cmt, getDose(ind, ind->idose[infBixds2]), getAllTimes(ind, ind->idose[infBixds2]), yp);
              rateOff = getAmt(ind, ind->id, ind->cmt, getDose(ind, ind->idose[infEixds]), getAllTimes(ind, ind->idose[infEixds]), yp);
            }
            // REprintf("get dur2\n");
            dur = getAllTimes(ind, ind->idose[infEixds]);// -
            // REprintf("\ttime infEixds: %f\n", dur);
            dur -= getAllTimes(ind, ind->idose[ind->ixds]);
            dur *= f;
            // REprintf("\tdur: %f\n", dur);
            dur2 =  curIi - dur;
            infSixds = infBixds;
            // REprintf("\tdur2: %f\n", dur2);
          }
        }
      } else if (isModeled) {
        isRateDose=true;
        // These are typically right next to another.
        infBixds=infBixds2=infEixds=infFixds=ind->ixds;
        if (isSsLag) {
          // The structure of this modeled rate with a lag item is:
          //
          // SS=1
          // TIME  EVID AMT II
          //    0 90209 100 24
          //    0 90201 100  0
          //    0 70201 100  0
          //
          // OR
          //
          // SS=2
          // TIME  EVID AMT II
          //    0 90219 100 24
          //    0 90201 100  0
          //    0 70201 100  0

          // The structure of this modeled duration item is:
          //
          // SS=1
          // TIME  EVID AMT II
          //    0 80209 100 24
          //    0 80201 100  0
          //    0 60201 100  0
          //
          // OR
          //
          // SS=2
          // TIME  EVID AMT II
          //    0 80219 100 24
          //    0 80201 100  0
          //    0 60201 100  0
          ignoreDoses[nIgnoredDoses++] = infBixds;
          int bEvid = getEvid(ind, ind->idose[infBixds]);
          double cIi = 0.0;
          double bDose = getDoseNumber(ind,ind->ixds);
          double cDose = bDose;
          int evid1 = bEvid - ind->wh0 + EVID0_REGULAR;
          bool foundEvid1 = false;
          int evidOff = bEvid - ind->wh0 + EVID0_REGULAR - 2*10000;
          bool foundEvidOff = false;
          int curRec = infBixds;
          while (!foundEvidOff && curRec < ind->ndoses) {
            curRec++;
            if (curRec == ind->ndoses) {
              infBixds = infBixds2 = infEixds = -1;
              break;
            } else {
              bEvid = getEvid(ind, ind->idose[curRec]);
              cDose = getDoseNumber(ind, curRec);
              cIi = getIiNumber(ind, curRec);
              if (!foundEvid1 && bEvid == evid1 && cDose == bDose && cIi == 0.0) {
                foundEvid1 = true;
                infBixds = infBixds2 = curRec;
                ignoreDoses[nIgnoredDoses++] = infBixds;
              } else if (!foundEvidOff && bEvid == evidOff  && cDose < 0.0 &&cIi == 0.0) {
                // note that this record stores the calculated infusion rate (in amt) and duration (in time)
                foundEvidOff = true;
                infEixds = curRec;
                ignoreDoses[nIgnoredDoses++] = infEixds;
                break;
              }
            }
          }
          if (infEixds == -1) {
            ind->wrongSSDur=1;
            // // Bad Solve => NA
            badSolveExit(*i);
            return;
          } else {
            // These use the getTime_() to grab calculated duration
            // REprintf("getDur3\n");
            dur = getTime_(ind->idose[infEixds], ind);// -
            dur -= getTime_(ind->idose[infBixds],ind);
            // REprintf("\tdur: %f\n", dur);
            dur2 = curIi - dur;
            // REprintf("\tdur2: %f\n", dur2);
            while (ind->ix[bi] != ind->idose[infBixds] && bi < ind->n_all_times) {
              bi++;
            }
            infSixds = infBixds;
          }
        } else {
          // SS=1
          // TIME  EVID AMT II
          //    0 90210 100 24
          //    0 70201 100  0
          //
          // OR
          //
          // SS=2
          // TIME  EVID AMT II
          //    0 90201 100 24
          //    0 70201 100  0

          // The structure of this modeled duration item is:
          //
          // SS=1
          // TIME  EVID AMT II
          //    0 80201 100 24
          //    0 60201 100  0
          //
          // OR
          //
          // SS=2
          // TIME  EVID AMT II
          //    0 80201 100 24
          //    0 60201 100  0
          infFixds = infBixds = infBixds2 = ind->ixds;
          infEixds = infBixds+1;
          ignoreDoses[nIgnoredDoses++] = infBixds;
          ignoreDoses[nIgnoredDoses++] = infEixds;
          dur = getAllTimes(ind, ind->idose[infBixds]);
          dur2 = getAllTimes(ind, ind->idose[infBixds+1]);
          dur = dur2-dur;
          dur2 = curIi - dur;
          infSixds = infBixds;
        }
        if (dur > curIi) {
          infSixds = infFixds;
        }
        rateOn = -getDose(ind, ind->idose[infBixds2+1]);
        rateOff = -rateOn;
      }
      if (ind->wh0 == EVID0_SSINF) {
        infEixds=infBixds=infBixds2 = ind->ixds;
        if (ind->whI == EVIDF_INF_RATE) {
          rateOn = getDose(ind, ind->idose[infBixds]);
          rateOff = -getDose(ind, ind->idose[infEixds]);
        } else if (isModeled) {
          rateOn =getRate(ind, ind->id, ind->cmt, 0.0,
                          getAllTimes(ind, ind->idose[ind->ixds]));
          rateOff = -rateOn;
        } else {
          // shouldn't ever get here modeled duration would have to be
          // infinite for this to occur... but I don't think that is
          // supported yet
          rateOn = getAmt(ind, ind->id, ind->cmt, getDose(ind, ind->idose[infBixds2]), getAllTimes(ind, ind->idose[infBixds2]), yp);
          rateOff = getAmt(ind, ind->id, ind->cmt, getDose(ind, ind->idose[infEixds]), getAllTimes(ind, ind->idose[infEixds]), yp);
        }
      } else if (ind->whI == EVIDF_INF_RATE ||
                 ind->whI == EVIDF_INF_DUR ||
                 isModeled) {
        ei = *i;
        while(ind->ix[ei] != ind->idose[infEixds] && ei < ind->n_all_times) {
          ei++;
        }
        if (ind->ix[ei] != ind->idose[infEixds]){
          /* Rf_errorcall(R_NilValue, "Cannot figure out infusion end time."); */
          if (!(ind->err & 8388608)){
            ind->err += 8388608;
            /* Rf_errorcall(R_NilValue, "Rate is zero/negative"); */
          }
          return;
        }
      }
      // REprintf("rateOn: %f rateOff: %f; dur: %f dur2: %f\n", rateOn, rateOff, dur, dur2);
      double startTimeD = 0.0;
      double curLagExtra = 0.0;
      int overIi = 0, regEvid = 0, extraEvid = 0;
      if (isRateDose) {
        startTimeD = getTime_(ind->idose[infSixds], ind);
        regEvid = getEvidClassic(ind->cmt+1, rateOn, rateOn, 0.0, 0.0, 1, 0);
        extraEvid = regEvid - EVID0_REGULAR + EVID0_RATEADJ;
      } else {
        startTimeD = getTime_(ind->idose[ind->ixds], ind);
      }
      if (isSsLag) {
        int wh0 = ind->wh0; ind->wh0=1;
        curLagExtra = getLag(ind, neq[1], ind->cmt, startTimeD) -
          startTimeD;
        ind->wh0 = wh0;
        overIi = floor(curLagExtra/curIi);
        curLagExtra = curLagExtra - overIi*curIi;
      }
      // First Reset
      for (j = neq[0]; j--;) {
        ind->InfusionRate[j] = 0;
        // ind->on[j] = 1; // nonmem doesn't reset on according to doc
      }
      // REprintf("reset & cancel pending doses\n");
      cancelInfusionsThatHaveStarted(ind, neq[1], startTimeD);
      if (!rx->ss2cancelAllPending && doSS2) {
      } else {
        cancelPendingDoses(ind, neq[1]);
      }
      ind->cacheME=0;
      // Reset LHS to NA
      ind->inLhs = 0;
      for (j = op->nlhs; j--;) ind->lhs[j] = NA_REAL;
      memcpy(yp,op->inits, neq[0]*sizeof(double));
      u_inis(neq[1], yp); // Update initial conditions @ current time
      if (rx->istateReset) *istate = 1;
      int k;
      double xp2, xout2;
      int canBreak=0;
      xp2 = xp;
      if (doSSinf || isSameTimeOp(curIi, dur)) {
        handleSSinf8(yp,
                     &xout, xp,
                     i,
                     istate,
                     op,
                     ind,
                     u_inis,
                     ctx,
                     &infBixds,
                     &bi,
                     &rateOn,
                     &xout2,
                     &xp2,
                     &canBreak, solveWith1Pt);
        if (isSameTimeOp(curIi, dur) && !isSameTimeOp(dur, 0.0)) {
          ind->InfusionRate[ind->cmt] = rateOn;
          pushDosingEvent(startTimeD+curIi, rateOff, extraEvid, ind);
          for (int ii = 0; ii < nIgnoredDoses; ++ii) {
            pushIgnoredDose(ignoreDoses[ii], ind);
          }
        } else {
          ind->InfusionRate[ind->cmt] = 0.0;
        }
        // REprintf("at ss: %f (inf: %f; rate: %f)\n", yp[ind->cmt],
        //          ind->InfusionRate[ind->cmt], rate);
        if (doSS2){
          // Add at the end
          for (j = neq[0];j--;) yp[j]+=ind->solveSave[j];
        }
        ind->doSS=0;
        updateExtraDoseGlobalsI(ind);
        return;
      } else if (dur == 0) {
        // Bolus
        handleSSbolus(yp,
                      &xout, xp,
                      i,
                      istate,
                      op,
                      ind,
                      u_inis,
                      ctx,
                      &xout2,
                      &xp2,
                      &curIi,
                      &canBreak, solveWith1Pt);
        if (isSsLag) {
          //advance the lag time
          ind->idx=*i;
          xout2 = xp2 + curIi - curLagExtra;
          // Use "real" xout for handle_evid functions.
          *istate=1;
          handle_evid(getEvid(ind, ind->ix[bi]),yp, xout, ind);
          // yp is last solve or y0
          solveWith1Pt(yp, xout2, xp2, i, istate, op, ind, u_inis, ctx);
          for (int cur = 0; cur < overIi; ++cur) {
            pushDosingEvent(startTimeD+curLagExtra+cur*curIi,
                            rateOn, regEvid, ind);
          }
          for (k = neq[0]; k--;){
            ind->solveLast[k] = yp[k];
          }
          ind->ixds--; // This dose stays in place
          // REprintf("ixds-- #3\n");
          xp2 = xout2;
        }
      } else {
        if (dur > curIi) {
          // in this case, the duration is greater than the inter-dose interval
          // number of doses before infusions turn off:
          int numDoseInf;
          double offTime, addTime;
          solveSSinfLargeDur(yp,
                             &xout, xp,
                             i,
                             istate,
                             op,
                             ind,
                             u_inis,
                             ctx,
                             &xout2,
                             &xp2,
                             &infBixds,
                             &bi,
                             &infEixds,
                             &ei,
                             &curIi,
                             &dur,
                             &numDoseInf,
                             &offTime,
                             &addTime,
                             &canBreak, solveWith1Pt);
          skipDosingEvent = true;
          // REprintf("Assign ind->ixds to %d (idx: %d) #1\n", indf->ixds, ind->idx);
          for (int cur = 0; cur < overIi; ++cur) {
            pushDosingEvent(startTimeD + offTime + cur*curIi + curLagExtra,
                            rateOff, extraEvid, ind);
            pushDosingEvent(startTimeD + (cur+1)*curIi + curLagExtra,
                            rateOn, extraEvid, ind);
          }
          for (int cur = 0; cur < numDoseInf+1; ++cur) {
            pushDosingEvent(startTimeD + offTime + (overIi+cur)*curIi + curLagExtra,
                            rateOff, extraEvid, ind);
          }
          if (curLagExtra > 0) {
            double solveTo=curIi - curLagExtra;
            if (solveTo > offTime) {
              // infusion where the lag time does not cause the infusion
              // to occur during the inter-dose interval
              xp2 = startTimeD;
              xout2 = xp2+offTime;
              ind->idx=bi;
              ind->ixds = infBixds;
              // REprintf("Assign ind->ixds to %d (idx: %d) #2\n", ind->ixds, ind->idx);
              handle_evid(getEvid(ind, ind->idose[infBixds]), yp, xout, ind);
              // yp is last solve or y0
              *istate=1;
              // yp is last solve or y0
              solveWith1Pt(yp,xout2, xp2, i, istate, op, ind, u_inis, ctx);
              xp2 = xout2;
              // Turn off Infusion, and solve to infusion stop time
              xout2 = xp2 + solveTo - offTime;
              ind->ixds = infEixds;
              ind->idx=ei;
              // REprintf("Assign ind->ixds to %d (idx: %d) #3\n", ind->ixds, ind->idx);
              handle_evid(getEvid(ind, ind->idose[infEixds]), yp, xout+dur, ind);
              *istate=1;
              solveWith1Pt(yp, xout2, xp2, i, istate, op, ind, u_inis, ctx);
              pushDosingEvent(startTimeD+curLagExtra,
                              rateOn, extraEvid, ind);
            } else {
              // infusion where the lag time occurs during the inter-dose interval split
              xp2 = startTimeD;
              xout2 = xp2 + solveTo;
              ind->idx=bi;
              ind->ixds = infBixds;
              // REprintf("Assign ind->ixds to %d (idx: %d) #4\n", ind->ixds, ind->idx);

              handle_evid(getEvid(ind, ind->idose[infBixds]), yp, xout, ind);
              // yp is last solve or y0
              *istate=1;
              // yp is last solve or y0
              solveWith1Pt(yp, xout2, xp2, i, istate, op, ind, u_inis, ctx);
              pushDosingEvent(startTimeD+offTime-solveTo,
                              rateOff, extraEvid, ind);
              pushDosingEvent(startTimeD+curLagExtra,
                              rateOn, extraEvid, ind);
            }
          } else {
            // infusion without a lag time.
            *istate=1;
            ind->idx = bi;
            ind->ixds = infBixds;
            // REprintf("Assign ind->ixds to %d (idx: %d) #5\n", ind->ixds, ind->idx);
            handle_evid(getEvid(ind, ind->idose[infBixds]), yp, xout, ind);
          }
          // REprintf("Assign ind->ixds to %d (idx: %d) #5a\n", ind->ixds, ind->idx);
          // yp is last solve or y0
          *istate=1;
          for (int ii = 0; ii < nIgnoredDoses; ++ii) {
            pushIgnoredDose(ignoreDoses[ii], ind);
          }
        } else if (ind->err) {
          printErr(ind->err, ind->id);
          badSolveExit(*i);
        } else {
          // Infusion
          dur2 = curIi-dur;
          if (isModeled && isSsLag) {
            // adjust start time for modeled w/ssLag
            startTimeD = getTime(ind->idose[infFixds],ind);
          }
          solveSSinf(yp,
                     &xout, xp,
                     i,
                     istate,
                     op,
                     ind,
                     u_inis,
                     ctx,
                     &xout2,
                     &xp2,
                     &infBixds,
                     &bi,
                     &infEixds,
                     &ei,
                     &curIi,
                     &dur,
                     &dur2,
                     &canBreak, solveWith1Pt);
          *istate=1;
          // REprintf("Assign ind->ixds to %d (idx: %d) #6\n", ind->ixds, ind->idx);
          for (k = neq[0]; k--;){
            ind->solveLast[k] = yp[k];
          }
          xp2 = xout2;
          bool doNoLag = false;
          if (isSsLag) {
            double totTime = xp2 + dur + dur2 - curLagExtra;
            if (curLagExtra > 0) {
              if (isSameTimeOp(curLagExtra+dur, curIi)) {
                // REprintf("isSameTimeOp(curLagExtra:%f+dur:%f, curIi:%f)\n",
                //          curLagExtra, dur, curIi);
                ind->idx  = bi;
                ind->ixds = infBixds;
                handle_evid(getEvid(ind, ind->idose[infBixds]), yp, xout, ind);
                xp2   = startTimeD;
                xout2 = startTimeD + dur;
                // yp is last solve or y0
                *istate=1;
                // yp is last solve or y0
                solveWith1Pt(yp, xout2, xp2, i, istate, op, ind, u_inis, ctx);
                ind->idx=ei;
                ind->ixds = infEixds;
                handle_evid(getEvid(ind, ind->idose[infEixds]), yp, xout, ind);
                ind->idx=ei;
                ind->ixds = infEixds;
                for (int cur = 0; cur < (overIi+1); ++cur) {
                  pushDosingEvent(startTimeD+curLagExtra+cur*curIi,
                                  rateOn, regEvid, ind);
                  pushDosingEvent(startTimeD+curLagExtra+dur+cur*curIi,
                                  rateOff, regEvid, ind);
                }
                ind->idx=fi;
                ind->ixds = infFixds;
              } else if (curLagExtra > dur2) {
                // REprintf("(curLagExtra: %f > dur2: %f)\n",
                //          curLagExtra, dur2);
                // dosing time occurs during the infusion
                double solveExtra=dur+dur2-curLagExtra;
                ind->idx=bi;
                ind->ixds = infBixds;
                handle_evid(getEvid(ind, ind->idose[infBixds]), yp, xout, ind);
                xp2   = startTimeD;
                xout2 = startTimeD + solveExtra;
                // yp is last solve or y0
                *istate=1;
                // yp is last solve or y0
                solveWith1Pt(yp, xout2, xp2, i, istate, op, ind, u_inis, ctx);

                for (int cur = 0; cur < (overIi+1); ++cur) {
                  pushDosingEvent(startTimeD+dur-solveExtra+cur*curIi,
                                  rateOff, (cur == 0) ? extraEvid : regEvid, ind);
                  pushDosingEvent(startTimeD+dur+dur2-solveExtra+cur*curIi,
                                  rateOn, regEvid, ind);
                }
                pushDosingEvent(startTimeD+dur-solveExtra+(overIi+1)*curIi,
                                rateOff, overIi == 0 ? extraEvid : regEvid, ind);
                ind->idx=fi;
                ind->ixds = infFixds;
              } else {
                if (xp2 + dur < totTime) {
                  xout2 = xp2 + dur;
                } else {
                  xout2 = totTime;
                }
                ind->idx=bi;
                ind->ixds = infBixds;
                handle_evid(getEvid(ind, ind->idose[infBixds]), yp, xout, ind);
                // yp is last solve or y0
                *istate=1;
                // yp is last solve or y0
                solveWith1Pt(yp, xout2, xp2, i, istate, op, ind, u_inis, ctx);
                if (!isSameTimeOp(totTime, xout2)) {
                  // don't give the infusion off dose
                  xp2 = xout2;
                  // Turn off Infusion, solve (dur-ii)
                  xout2 = totTime;
                  ind->ixds = infEixds;
                  ind->idx = ei;
                  handle_evid(getEvid(ind, ind->idose[infEixds]), yp, xout+dur, ind);
                  // yp is last solve or y0
                  *istate=1;
                  solveWith1Pt(yp, xout2, xp2, i, istate, op, ind, u_inis, ctx);
                }
                for (int cur = 0; cur < (overIi+1); ++cur) {
                  pushDosingEvent(startTimeD+curLagExtra+cur*curIi,
                                  rateOn, regEvid, ind);
                  pushDosingEvent(startTimeD+curLagExtra+dur+cur*curIi,
                                  rateOff, regEvid, ind);
                }
              }
            } else {
              doNoLag = true;
            }
          } else {
            doNoLag = true;
          }
          if (doNoLag) {
            ind->idx=bi;
            ind->ixds=infBixds;
            handle_evid(getEvid(ind, ind->idose[infBixds]), yp, xout, ind);
            pushDosingEvent(startTimeD+dur,
                            rateOff, extraEvid, ind);
          }
          *istate=1;
          for (k = neq[0]; k--;){
            ind->solveLast[k] = yp[k];
          }
          xp2 = xout2;
          ind->idx=fi;
          ind->ixds = infFixds;
          for (int ii = 0; ii < nIgnoredDoses; ++ii) {
            pushIgnoredDose(ignoreDoses[ii], ind);
          }
        }
      }
      if (doSS2){
        // Add at the end
        for (j = neq[0];j--;) yp[j]+=ind->solveSave[j];
      }
      if (!doSSinf && !isSsLag && !skipDosingEvent){
        // REprintf("handleEvid %d %d %d\n", doSSinf, isSsLag, skipDosingEvent);
        handle_evid(getEvid(ind, ind->ix[*i]), yp, xout, ind);
      }
      ind->doSS=0;
      ind->ixds=oldIxds; ind->idx=oldIdx;
    }
    updateExtraDoseGlobalsI(ind);
  }

#undef max2
#undef badSolveExit

#if defined(__cplusplus)
}
#endif

#endif
