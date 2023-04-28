// -*- mode: c++; c-basic-offset: 2; tab-width: 2; indent-tabs-mode: t; -*-
#ifndef __RXODE2PARSEGETTTIME_H__
#define __RXODE2PARSEGETTTIME_H__
#if defined(__cplusplus)

#include "rxode2parse.h"

extern t_F AMT;
extern t_LAG LAG;
extern t_RATE RATE;
extern t_DUR DUR;
extern t_calc_mtime calc_mtime;

#ifndef __DOINIT__

static inline double getLag(rx_solving_options_ind *ind, int id, int cmt, double time){
	if (ind->wh0 == EVID0_SS0 || ind->wh0 == EVID0_SS20) {
		return time;
	}
  double ret = LAG(id, cmt, time);
  if (ISNA(ret)) {
    rx_solving_options *op = &op_global;
    op->badSolve=1;
    op->naTime = 1;
  }
  return ret;
}

static inline double getRate(rx_solving_options_ind *ind, int id, int cmt, double dose, double t){
  double ret = RATE(id, cmt, dose, t);
  if (ISNA(ret)){
    rx_solving_options *op = &op_global;
    op->badSolve=1;
    op->naTime = 1;
  }
  return ret;
}

static inline double getDur(rx_solving_options_ind *ind, int id, int cmt, double dose, double t){
  double ret = DUR(id, cmt, dose, t);
  if (ISNA(ret)){
    rx_solving_options *op = &op_global;
    op->badSolve=1;
    op->naTime = 1;
  }
  return ret;
}


static inline int isEvidType(int evid, int type) {
  int wh, cmt, wh100, whI, wh0;
  getWh(evid, &wh, &cmt, &wh100, &whI, &wh0);
  return (whI == type);
}

#define isEvidModeledDurationStart(evid) isEvidType(evid, EVIDF_MODEL_DUR_ON)
#define isEvidModeledDurationStop(evid) isEvidType(evid, EVIDF_MODEL_DUR_OFF)
#define isEvidModeledRateStart(evid) isEvidType(evid, EVIDF_MODEL_RATE_ON)
#define isEvidModeledRateStop(evid) isEvidType(evid, EVIDF_MODEL_RATE_OFF)

static inline void updateDur(int idx, rx_solving_options_ind *ind, double *yp){
  double t = getAllTimes(ind, idx);
  double dur, amt;
  // The duration and f cannot depend on state values
  int oldIdx = ind->idx;
  ind->idx = idx;
  amt  = getAmt(ind, ind->id, ind->cmt, getDose(ind, idx), t, yp);
  dur  = getDur(ind, ind->id, ind->cmt, amt, t);
  ind->idx = oldIdx;
  if (dur > 0) {
		setDoseP1(ind, idx, -amt/dur);
		setAllTimesP1(ind, idx, t+dur);
  } else {
    rx_solve *rx = &rx_global;
    rx_solving_options *op = &op_global;
    if (ind->cmt < op->neq){
      if (rx->needSort & 4){
        if (!(ind->err & 16)){
          ind->err += 16;
        }
        return;
      } else {
        if (!(ind->err & 32)){
          ind->err += 32;
        }
        return;
      }
    }
  }
}

static inline void updateRate(int idx, rx_solving_options_ind *ind, double *yp) {
  double t = getAllTimes(ind, idx);
  int oldIdx = ind->idx;
  ind->idx=idx;
  double dur, rate, amt;
  amt  = getAmt(ind, ind->id, ind->cmt, getDose(ind,idx), t, yp);
  rate  = getRate(ind, ind->id, ind->cmt, amt, t);
  if (rate > 0){
    dur = amt/rate; // mg/hr
		setDoseP1(ind, idx, -rate);
		setAllTimesP1(ind, idx, t+dur);
    ind->idx=oldIdx;
  } else {
    rx_solve *rx;
    rx = &rx_global;
    rx_solving_options *op = &op_global;
    if (ind->cmt < op->neq){
      if (rx->needSort & 8){
        if (!(ind->err & 2)){
          ind->err += 2;
        }
        return;
      } else {
        // FIXME don't error out with linear compartmental model
        if (!(ind->err & 4)){
          ind->err += 4;
        }
        return;
      }
    }
  }
  ind->idx=oldIdx;
}

static inline void handleTurnOffModeledDuration(int idx, rx_solve *rx, rx_solving_options *op, rx_solving_options_ind *ind) {
  if (idx > 0){
    if (!isEvidModeledDurationStart(getEvidM1(ind, idx))) {
      if (!(ind->err & 64)){
        ind->err += 64;
      }
      return;
    }
  } else {
    if (!(ind->err & 128)){
      ind->err += 128;
    }
    return;
  }
}

static inline void handleTurnOnModeledDuration(int idx, rx_solve *rx, rx_solving_options *op, rx_solving_options_ind *ind) {
  // This calculates the rate and the duration and then assigns it to the next record
  if (idx >= ind->n_all_times){
    // error: Last record, can't be used.
    if (!(ind->err & 256)){
      ind->err += 256;
    }
    return;
  } else {
    if (!isEvidModeledDurationStop(getEvidP1(ind, idx))) {
      if (!(ind->err & 512)){
        ind->err += 512;
      }
      return;
    }
    updateDur(idx, ind, rx->ypNA);
  }
}

static inline void handleTurnOffModeledRate(int idx, rx_solve *rx, rx_solving_options *op, rx_solving_options_ind *ind) {
  if (idx > 0){
    if (!isEvidModeledRateStart(getEvidM1(ind, idx))) {
      if (!(ind->err & 1024)){
        ind->err += 1024;
      }
      return;
    }
  } else {
    if (!(ind->err & 2048)){
      ind->err += 2048;
    }
    return;
  }
}

static inline void handleTurnOnModeledRate(int idx, rx_solve *rx, rx_solving_options *op, rx_solving_options_ind *ind) {
  // This calculates the rate and the duration and then assigns it to the next record
  if (idx >= ind->n_all_times){
    // error: Last record, can't be used.
    if (!(ind->err & 4096)){
      ind->err += 4096;
    }
    /* Rf_errorcall(R_NilValue, "Data Error 9\n"); */
    return;
  } else {
    if (!isEvidModeledRateStop(getEvidP1(ind, idx))) {
      if (!(ind->err & 8192)){
        ind->err += 8192;
      }
      return;
    }
		setAllTimesP1(ind, idx, getAllTimes(ind, idx));
    updateRate(idx, ind, rx->ypNA);
  }
}

static inline double handleInfusionItem(int idx, rx_solve *rx, rx_solving_options *op, rx_solving_options_ind *ind) {
  double amt = getDose(ind, idx);
  if (amt > 0) {
		return getLag(ind, ind->id, ind->cmt, getAllTimes(ind, idx));
  } else if (amt < 0){
    int j = getDoseNumberFromIndex(ind, idx);
    if (j == -1){
      if (!(ind->err & 16384)){
        ind->err += 16384;
      }
      return 0.0;
      /* Rf_errorcall(R_NilValue, "Corrupted event table during sort (1)."); */
    }
		int k;
		if (ind->wh0 == EVID0_INFRM) {
			k = j+1;
			for (j = k; j < ind->ndoses; ++j) {
				if (getEvid(ind, ind->idose[j]) == getEvid(ind, ind->idose[k])) break;
				if (j == ind->ndoses-1) {
					if (!(ind->err & 32768)){
						ind->err += 32768;
					}
					return 0.0;
				}
			}
		}  else {
			for (k = j; k--;){
				if (getEvid(ind, ind->idose[j]) == getEvid(ind, ind->idose[k])) break;
				if (k == 0) {
					if (!(ind->err & 32768)){
						ind->err += 32768;
					}
					return 0.0;
				}
			}			
		}
    rx_solve *rx = &rx_global;
    double f = getAmt(ind, ind->id, ind->cmt, 1.0, getAllTimes(ind, ind->idose[j-1]), rx->ypNA);
    if (ISNA(f)){
      rx_solving_options *op = &op_global;
      op->badSolve=1;
      op->naTime = 1;
    }
    double durOld = (getAllTimes(ind, ind->idose[j]) -
                     getAllTimes(ind, ind->idose[k]));
    double dur = f*durOld;
    double t = getAllTimes(ind, ind->idose[k]) + dur;
    return getLag(ind, ind->id, ind->cmt, t);
  } else {
    /* Rf_errorcall(R_NilValue, "Corrupted events."); */
    if (!(ind->err & 131072)){
      ind->err += 131072;
    }
    return 0.0;
  }
}

static inline double getTimeCalculateInfusionTimes(int idx, rx_solve *rx, rx_solving_options *op, rx_solving_options_ind *ind) {
  switch(ind->whI){
  case EVIDF_MODEL_DUR_OFF:
    handleTurnOffModeledDuration(idx, rx, op, ind);
    break;
  case EVIDF_MODEL_DUR_ON:
    handleTurnOnModeledDuration(idx, rx, op, ind);
    break;
  case EVIDF_MODEL_RATE_OFF:
    handleTurnOffModeledRate(idx, rx, op, ind);
    break;
  case EVIDF_MODEL_RATE_ON:
    handleTurnOnModeledRate(idx, rx, op, ind);
    break;
  case EVIDF_INF_RATE:
    return handleInfusionItem(idx, rx, op, ind);
    break;
  }
  return getLag(ind, ind->id, ind->cmt, getAllTimes(ind,idx));
}

static inline double getTime__(int idx, rx_solving_options_ind *ind, int update) {
  rx_solving_options *op = &op_global;
  rx_solve *rx = &rx_global;
  int evid = getEvid(ind, idx);
  if (evid == 9) return 0.0;
  if (evid >= 10 && evid <= 99) return ind->mtime[evid-10];
  if (isObs(evid)) return getAllTimes(ind, idx);
  getWh(evid, &(ind->wh), &(ind->cmt), &(ind->wh100), &(ind->whI), &(ind->wh0));
  if (ind->wh0 == EVID0_SSINF){
  } else {
    // yp should be the current solve values
    //
    // Before solving the solve will be zero
    // After solving the yp will contain the solved values
    //
    if (update == 0) {
      if (ind->whI == EVIDF_INF_RATE) {
        return handleInfusionItem(idx, rx, op, ind);
      }
    } else {
      return getTimeCalculateInfusionTimes(idx, rx, op, ind);
    }
  }
  return getLag(ind, ind->id, ind->cmt, getAllTimes(ind, idx));
}

static inline double getTime_(int idx, rx_solving_options_ind *ind) {
  return getTime__(idx, ind, 0);
}


#endif

extern "C" {
#endif

#ifndef _isrxode2parse_
  double getTime(int idx, rx_solving_options_ind *ind);
#endif
	
#if defined(__cplusplus)
}
#endif
#endif
