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

#include "../inst/include/rxode2parse.h"
//#include "../inst/include/rxode2parseHandleEvid.h"
#include "../inst/include/rxode2parseGetTime.h"
#include "solComp.h"

int comp1solve3(double *yp, // prior solving information, will be updated with new information (like lsoda and the like)
                double *xout, // time to solve to
                double *xp, // last time
                double *rate, // rate in central compartment
                double *ka,
                double *k10,
                double *k12,
                double *k21,
                double *k13,
                double *k31) {
  double L[3], C1[9], C2[9], C3[9], E[3], Ea[3], Xo[3], Rm[3];
  int hasDepot = (*ka) != 0.0;
  double dT = (*xout)-(*xp);
  if (solComp3C(k10, k12, k21, k13, k31, L, C1, C2, C3) == 0) {
    return 0;
  } 
  E[0] = Ea[0] = exp(-L[0]*dT);
  E[1] = Ea[1] = exp(-L[1]*dT);
  E[2] = Ea[2] = exp(-L[2]*dT);
  const double one = 1.0, zero = 0.0;
  const int ione = 1, itwo = 2, ithree=3;
  //Xo = Xo + pX[1 + j] * Co[, , j] %*% E # Bolus
  F77_CALL(dgemm)("N", "N", &ithree, &ione, &ithree, &(yp[hasDepot+1]), C1, &ithree,
                  E, &ithree, &zero, Xo, &ithree FCONE FCONE);
  F77_CALL(dgemm)("N", "N", &ithree, &ione, &ithree, &(yp[hasDepot+2]), C2, &ithree,
                  E, &ithree, &one, Xo, &itwo FCONE FCONE);
  F77_CALL(dgemm)("N", "N", &ithree, &ione, &ithree, &(yp[hasDepot+3]), C3, &ithree,
                  E, &ithree, &one, Xo, &itwo FCONE FCONE);
  if (!isSameTime(*rate, 0.0)) {
    // Xo = Xo + ((cR*Co[, , 1]) %*% ((1 - E)/L)) # Infusion
    Rm[0] = (1.0 - E[0])/L[0];
    Rm[1] = (1.0 - E[1])/L[1];
    Rm[2] = (1.0 - E[2])/L[2];
    F77_CALL(dgemm)("N", "N", &ithree, &ione, &ithree, rate, C1, &ithree,
                    Rm, &ithree, &one, Xo, &ithree FCONE FCONE);
  }
  if (hasDepot == 1 && yp[0] > 0.0) {
    // Xo = Xo + Ka*pX[1]*(Co[, , 1] %*% ((E - Ea)/(Ka - L)))
    double expa = exp(-(*ka)*dT);
    Ea[0] = (E[0]- expa)/((*ka) - L[0]);
    Ea[1] = (E[1]- expa)/((*ka) - L[1]);
    Ea[2] = (E[2]- expa)/((*ka) - L[2]);
    expa = (*ka)*yp[0];
    F77_CALL(dgemm)("N", "N", &ithree, &ione, &ithree, &expa, C1, &ithree,
                    Ea, &ithree, &one, Xo, &ithree FCONE FCONE);
    yp[0] *= expa;
  }
  yp[hasDepot]   = Xo[0];
  yp[hasDepot+1] = Xo[1];
  yp[hasDepot+2] = Xo[2];
  return 1;
}
