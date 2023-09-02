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
#include "../inst/include/rxode2parseHandleEvid.h"
#include "../inst/include/rxode2parseGetTime.h"
#include "solComp.h"


int comp1solve2(double *yp, // prior solving information, will be updated with new information (like lsoda and the like)
                 double *xout, // time to solve to
                 double *xp, // last time
                 double *rate, // rate in central compartment
                 double *ka,
                 double *k10,
                 double *k12,
                 double *k21) {
  double L[2], C1[4], C2[4], E[2], Ea[2], Xo[2], Rm[2];
  int hasDepot = (*ka) != 0.0;
  double dT = (*xout)-(*xp);
  if (solComp2C(k10, k12, k21, L, C1, C2) == 0) {
    return 0;
  } 
  E[0] = Ea[0] = exp(-L[0]*dT);
  E[1] = Ea[1] = exp(-L[1]*dT);
  const double one = 1.0, zero = 0.0;
  const int ione = 1, itwo = 2;
  //Xo = Xo + pX[1 + j] * Co[, , j] %*% E # Bolus
  F77_CALL(dgemm)("N", "N", &itwo, &ione, &itwo, &(yp[hasDepot+1]), C1, &itwo, E, &itwo, &zero, Xo, &itwo FCONE FCONE);
  F77_CALL(dgemm)("N", "N", &itwo, &ione, &itwo, &(yp[hasDepot+2]), C2, &itwo, E, &itwo, &one, Xo, &itwo FCONE FCONE);
  if (!isSameTime(*rate, 0.0)) {
    // Xo = Xo + ((cR*Co[, , 1]) %*% ((1 - E)/L)) # Infusion
    Rm[0] = (1.0 - E[0])/L[0];
    Rm[1] = (1.0 - E[1])/L[1];
    F77_CALL(dgemm)("N", "N", &itwo, &ione, &itwo, rate, C1, &itwo, Rm, &itwo, &one, Xo, &itwo FCONE FCONE);
  }
  if (hasDepot == 1 && yp[0] > 0.0) {
    // Xo = Xo + Ka*pX[1]*(Co[, , 1] %*% ((E - Ea)/(Ka - L)))
    double expa = exp(-(*ka)*dT);
    Ea[0] = (E[0]- expa)/((*ka)-L[0]);
    Ea[1] = (E[1]- expa)/((*ka)-L[1]);
    expa = (*ka)*yp[0];
    F77_CALL(dgemm)("N", "N", &itwo, &ione, &itwo, &expa, C1, &itwo, Ea, &itwo, &one, Xo, &itwo FCONE FCONE);
    yp[0] *= expa;
  }
  yp[hasDepot] = Xo[0];
  yp[hasDepot+1] = Xo[1];
  return 1;
}
