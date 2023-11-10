typedef struct lin_context_c_t {
  double *rate;
  int oral0;
  double ka;
  double k10;
  double k12;
  double k21;
  double k13;
  double k31;
  double v;
} lin_context_c_t;

static inline void printSqMat(double *in, int d) {
  // debugging function to print a square matrix in R.
  int pro=0;
  SEXP mat = PROTECT(Rf_allocVector(REALSXP, d*d)); pro++;
  SEXP dm  = PROTECT(Rf_allocVector(INTSXP, 2)); pro++;
  int *dmi = INTEGER(dm);
  dmi[0] =dmi[1] = d;
  double *matd = REAL(mat);
  for (int i = d*d; i--;)  {
    matd[i] = in[i];
  }
  Rf_setAttrib(mat, R_DimSymbol, dm);
  Rf_PrintValue(mat);
  UNPROTECT(pro);
}

static inline void printVec(double *in, int d) {
  int pro=0;
  SEXP vec = PROTECT(Rf_allocVector(REALSXP, d)); pro++;
  double *vecd =REAL(vec);
  for (int i = d; i--;) {
    vecd[i] = in[i];
  }
  Rf_PrintValue(vec);
  UNPROTECT(pro);
}

// Handle single point solve
static inline int comp1solve1(double *yp, // prior solving information, will be updated with new information (like lsoda and the like)
                              double *xout, // time to solve to
                              double *xp, // last time
                              double *rate, // rate in central compartment
                              double *ka,
                              double *kel) {
  double dt = (*xout)-(*xp);
  double E  = exp(-(*kel)*dt);
  double Ea = E;
  rx_solve *rx=(&rx_global);
  int hasDepot = rx->linKa;
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
  return 1;
}

static inline int comp1solve2(double *yp, // prior solving information, will be updated with new information (like lsoda and the like)
                 double *xout, // time to solve to
                 double *xp, // last time
                 double *rate, // rate in central compartment
                 double *ka,
                 double *k10,
                 double *k12,
                 double *k21) {
  double L[2], C1[4], C2[4], E[2], Ea[2], Xo[2], Rm[2];
  rx_solve *rx=(&rx_global);
  int hasDepot = rx->linKa;
  double dT = (*xout)-(*xp);
  if (solComp2C(k10, k12, k21, L, C1, C2) == 0) {
    return 0;
  }
  Xo[0] = Xo[1] = 0.0;
  E[0] = Ea[0] = exp(-L[0]*dT);
  E[1] = Ea[1] = exp(-L[1]*dT);
  const double one = 1.0, zero = 0.0;
  const int ione = 1, itwo = 2;
  //Xo = Xo + pX[1 + j] * Co[, , j] %*% E # Bolus
  F77_CALL(dgemm)("N", "N", &itwo, &ione, &itwo, &(yp[hasDepot]), C1, &itwo, E,
                  &itwo, &one, Xo, &itwo FCONE FCONE);

  F77_CALL(dgemm)("N", "N", &itwo, &ione, &itwo, &(yp[hasDepot+1]), C2, &itwo, E,
                  &itwo, &one, Xo, &itwo FCONE FCONE);
  if (!isSameTime(rate[hasDepot], 0.0)) {
    // Xo = Xo + ((cR*Co[, , 1]) %*% ((1 - E)/L)) # Infusion
    Rm[0] = (1.0 - E[0])/L[0];
    Rm[1] = (1.0 - E[1])/L[1];
    F77_CALL(dgemm)("N", "N", &itwo, &ione, &itwo, &(rate[hasDepot]),
                    C1, &itwo, Rm, &itwo, &one, Xo, &itwo FCONE FCONE);
  }
  if (hasDepot == 1 && yp[0] > 0.0) {
    // Xo = Xo + Ka*pX[1]*(Co[, , 1] %*% ((E - Ea)/(Ka - L)))
    double expa = exp(-(*ka)*dT);
    Ea[0] = (E[0]- expa)/((*ka)-L[0]);
    Ea[1] = (E[1]- expa)/((*ka)-L[1]);
    double cf = (*ka)*yp[0];
    F77_CALL(dgemm)("N", "N", &itwo, &ione, &itwo, &cf, C1, &itwo, Ea, &itwo, &one, Xo, &itwo FCONE FCONE);
    yp[0] *= expa;
  }
  yp[hasDepot] = Xo[0];
  yp[hasDepot+1] = Xo[1];
  return 1;
}


static inline int comp1solve3(double *yp, // prior solving information, will be updated with new information (like lsoda and the like)
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
  rx_solve *rx=(&rx_global);
  int hasDepot = rx->linKa;
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
  Xo[0]=Xo[1]=Xo[2]=0.0;
  F77_CALL(dgemm)("N", "N", &ithree, &ione, &ithree, &(yp[hasDepot]), C1, &ithree,
                  E, &ithree, &one, Xo, &ithree FCONE FCONE);
  F77_CALL(dgemm)("N", "N", &ithree, &ione, &ithree, &(yp[hasDepot+1]), C2, &ithree,
                  E, &ithree, &one, Xo, &ithree FCONE FCONE);
  F77_CALL(dgemm)("N", "N", &ithree, &ione, &ithree, &(yp[hasDepot+2]), C3, &ithree,
                  E, &ithree, &one, Xo, &ithree FCONE FCONE);
  if (!isSameTime(rate[hasDepot], 0.0)) {
    // Xo = Xo + ((cR*Co[, , 1]) %*% ((1 - E)/L)) # Infusion
    Rm[0] = (1.0 - E[0])/L[0];
    Rm[1] = (1.0 - E[1])/L[1];
    Rm[2] = (1.0 - E[2])/L[2];
    F77_CALL(dgemm)("N", "N", &ithree, &ione, &ithree, &(rate[hasDepot]), C1, &ithree,
                    Rm, &ithree, &one, Xo, &ithree FCONE FCONE);
  }
  if (hasDepot == 1 && yp[0] > 0.0) {
    // Xo = Xo + Ka*pX[1]*(Co[, , 1] %*% ((E - Ea)/(Ka - L)))
    double expa = exp(-(*ka)*dT);
    Ea[0] = (E[0]- expa)/((*ka) - L[0]);
    Ea[1] = (E[1]- expa)/((*ka) - L[1]);
    Ea[2] = (E[2]- expa)/((*ka) - L[2]);
    double expa2 = (*ka)*yp[0];
    F77_CALL(dgemm)("N", "N", &ithree, &ione, &ithree, &expa2, C1, &ithree,
                    Ea, &ithree, &one, Xo, &ithree FCONE FCONE);
    yp[0] = expa*yp[0];
  }
  yp[hasDepot]   = Xo[0];
  yp[hasDepot+1] = Xo[1];
  yp[hasDepot+2] = Xo[2];
  return 1;
}
