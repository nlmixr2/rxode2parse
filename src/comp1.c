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

// Handle single point solve
void comp1solve1(double *yp, // prior solving information, will be updated with new information (like lsoda and the like)
                 double *xout, // time to solve to
                 double *xp, // last time
                 double *rate, // rate in central compartment
                 double *ka,
                 double *kel) {
  double dt = (*xout)-(*xp);
  double E  = exp(-(*kel)*dt);
  double Ea = E;
  int hasDepot = (*ka) != 0.0;
  if (hasDepot == 1) {
    Ea = exp(-(*ka)*dt);
  }
  double pDepot = yp[0];
  if (isSameTime((*ka), (*kel))) {
    yp[hasDepot] = yp[hasDepot]*E + (*rate)*(1.0-E)/(*kel) + pDepot*(*kel)*dt*E;
  } else {
    yp[hasDepot] = E*yp[hasDepot] + (*rate)*(1.0-E)/(*kel) + pDepot*(*ka)*(E-Ea)/((*ka)-(*kel));
  }
  if (hasDepot) {
    yp[0] *= Ea;
  }
}
