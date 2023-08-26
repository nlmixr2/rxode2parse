// -*- mode: c++; c-basic-offset: 2; tab-width: 2; indent-tabs-mode: nil; -*-
#ifndef __RXODE2PARSETIMEEVIDSHARED_H___
#define __RXODE2PARSETIMEEVIDSHARED_H___

#if defined(__cplusplus)
#include "timsort.h"
#define SORT gfx::timsort
static inline void sortExtraDose(rx_solving_options_ind *ind) {
  if (ind->extraDoseN[0] && ind->extraSorted == 0) {
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
}


extern "C" {
#endif

  static inline int isIgnoredDose(rx_solving_options_ind *ind, int ixds) {
    for (int i = 0; i < ind->ignoredDosesN[0]; ++i) {
      if (ind->idx < 0 ) {
        return 0;
      } else if (ind->idx >= 0 && ind->ignoredDoses[i] == ixds) {
        return 1;
      }
    }
    return 0;
  }



  static inline void handleInfusionGetEndOfInfusionIndex(int idx, int *infEixds,
                                                         rx_solve *rx, rx_solving_options *op,
                                                         rx_solving_options_ind *ind) {
    int extraDose = idx < 0;
#if defined(__cplusplus)
    if (extraDose) sortExtraDose(ind);
#endif
    int trueIdx = extraDose ? idx : ind->idose[idx];
	int curEvid = getEvid(ind, trueIdx);
	double curAmt = getDose(ind, trueIdx);
	int lastKnownOff = 0;
	*infEixds = NA_INTEGER;
    int N = extraDose ? ind->extraDoseN[0] : ind->ndoses;
	for (int j = 0; j < N; j++) {
      int testIdx = extraDose ? -1-ind->extraDoseTimeIdx[j] : ind->idose[j];
      if (curEvid == getEvid(ind, testIdx) &&
          curAmt == getDose(ind, testIdx)) {
        // get the first dose combination
        if (lastKnownOff == 0) {
          lastKnownOff=j+1;
        } else {
          lastKnownOff++;
        }
        for (int k = lastKnownOff; k < N; k++) {
          int testIdxK = extraDose ? -1-ind->extraDoseTimeIdx[k] : ind->idose[k];
          if (curEvid == getEvid(ind, testIdxK) &&
              curAmt == -getDose(ind, testIdxK)) {
            lastKnownOff = k;
            if (extraDose) {
              if (testIdx == idx) {
                *infEixds = testIdxK;
              }
            } else {
              if (j == idx) {
                *infEixds = k;
              }
            }
            k = N;
          }
        }
      }
      if (*infEixds != NA_INTEGER) break;
	}
    if (!extraDose && *infEixds == NA_INTEGER) *infEixds=-1;
  }


#if defined(__cplusplus)
}
#undef SORT
#else
#endif

#endif
