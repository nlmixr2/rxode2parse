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

// Infinite infusion
static inline void comp2ssInf8(double *yp, double *rate, double *ka, double *k10,
                               double *k12, double *k21) {
  int hasDepot = (*ka) != 0.0;
  
  double E1 = (*k10)+(*k12);
  double s = E1+(*k21);
  double sqr = sqrt(s*s-4*(E1*(*k21)-(*k12)*(*k21)));
  double lambda1 = 0.5*(s+sqr);
  double lambda2 = 0.5*(s-sqr);
  double l12 = 1.0/(lambda1*lambda2);
  yp[hasDepot]=(*rate)*(*k21)*l12;
  yp[hasDepot+1]=(*rate)*(*k12)*l12;
}

static inline void comp2ssInf(double *yp, double *dur, double *ii, double *rate,
                              double *ka, double *k10, double *k12, double *k21) {
  int hasDepot = (*ka) != 0.0;
  double E1 = (*k10)+(*k12);
  double E2 = (*k21);

  double s = E1+E2;
  double sqr = sqrt(s*s-4*(E1*E2-(*k12)*(*k21)));
  double L1 = 0.5*(s+sqr);
  double L2 = 0.5*(s-sqr);

  double eTi1 = exp(-(*dur)*L1);
  double eTi2 = exp(-(*dur)*L2);
  double eT1 =exp(-L1*((*ii)-(*dur)))/(1.0-exp(-(*ii)*L1));
  double eT2 =exp(-L2*((*ii)-(*dur)))/(1.0-exp(-(*ii)*L2));
  yp[hasDepot]=(eT1*(E2*((eTi1*(*rate) - eTi2*(*rate))/(-L1 + L2) + (*rate)*E2*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))) - L1*((eTi1*(*rate) - eTi2*(*rate))/(-L1 + L2) + (*rate)*E2*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))) + (*rate)*(*k12)*(*k21)*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))) - eT2*(E2*((eTi1*(*rate) - eTi2*(*rate))/(-L1 + L2) + (*rate)*E2*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))) - L2*((eTi1*(*rate) - eTi2*(*rate))/(-L1 + L2) + (*rate)*E2*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))) + (*rate)*(*k12)*(*k21)*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))))/(-L1 + L2);
  yp[hasDepot+1]=(eT1*((*k12)*((eTi1*(*rate) - eTi2*(*rate))/(-L1 + L2) + (*rate)*E2*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))) + (*rate)*E1*(*k12)*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2)) - (*rate)*(*k12)*L1*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))) - eT2*((*k12)*((eTi1*(*rate) - eTi2*(*rate))/(-L1 + L2) + (*rate)*E2*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))) + (*rate)*E1*(*k12)*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2)) - (*rate)*(*k12)*L2*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))))/(-L1 + L2);
}

// Steady state central dosing
static inline void comp2ssBolusCentral(double *yp, double *ii, double *dose,
                                       double *ka, double *k10, double *k12, double *k21) {
  int hasDepot = (*ka) != 0.0;
  double E2 = (*k21);

  double s = (*k12)+(*k21)+(*k10);
  double sqr = sqrt(s*s-4*(*k21)*(*k10));
  double L1 = 0.5*(s+sqr);
  double L2 = 0.5*(s-sqr);

  double eL1 = 1.0/(1.0-exp(-(*ii)*L1));
  double eL2 = 1.0/(1.0-exp(-(*ii)*L2));
  
  yp[hasDepot]=(eL1*((*dose)*E2 - (*dose)*L1) - eL2*((*dose)*E2 - (*dose)*L2))/(-L1 + L2);
  yp[hasDepot+1]=(eL1*(*dose)*(*k12) - eL2*(*dose)*(*k12))/(-L1 + L2);
}


// Steady state bolus dosing
static inline void comp2ssBolusDepot(double *yp, double *ii, double *dose, 
                                     double *ka, double *k10, double *k12, double *k21) {
  double E2 = (*k10)+(*k12);
  double E3 = (*k21);
  double e2e3 = E2+E3;
  double s = sqrt(e2e3*e2e3-4*(E2*E3-(*k12)*(*k21)));

  double L1 = 0.5*(e2e3+s);
  double L2 = 0.5*(e2e3-s);
  double eKa=1.0/(1.0-exp(-(*ii)*(*ka)));
  double eL1=1.0/(1.0-exp(-(*ii)*L1));
  double eL2=1.0/(1.0-exp(-(*ii)*L2));
  yp[0]=eKa*(*dose);
  yp[1]=(*ka)*(*dose)*(eL1*(E3 - L1)/((-L1 + L2)*((*ka) - L1)) + eL2*(E3 - L2)/((L1 - L2)*((*ka) - L2)) + eKa*(E3 - (*ka))/((-(*ka) + L2)*(-(*ka) + L1)));
  yp[2]=(*ka)*(*dose)*(*k12)*(eL1/((-L1 + L2)*((*ka) - L1)) + eL2/((L1 - L2)*((*ka) - L2)) + eKa/((-(*ka) + L2)*(-(*ka) + L1)));
}

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
