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

#include "compSSc.h"

int comp1solve1(double *yp, // prior solving information, will be updated with new information (like lsoda and the like)
                double *xout, // time to solve to
                double *xp, // last time
                double *rate, // rate in central compartment
                double *ka,
                double *kel);
int comp1solve2(double *yp, // prior solving information, will be updated with new information (like lsoda and the like)
                double *xout, // time to solve to
                double *xp, // last time
                double *rate, // rate in central compartment
                double *ka,
                double *k10,
                double *k12,
                double *k21);
int comp1solve3(double *yp, // prior solving information, will be updated with new information (like lsoda and the like)
                double *xout, // time to solve to
                double *xp, // last time
                double *rate, // rate in central compartment
                double *ka,
                double *k10,
                double *k12,
                double *k21,
                double *k13,
                double *k31);

