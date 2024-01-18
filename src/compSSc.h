//////////////////////////////////////////////////
// 1cmt
//////////////////////////////////////////////////
// Infinite infusion
static inline void comp1ssInf8(double *yp, double *rate, double *ka, double *kel) {
  rx_solve *rx=(&rx_global);
  int central=0;
  if (rx->linKa) {
    // has depot
    if (isSameTime(rate[0], 0.0)) {
      // central infusion
      central=1;
      yp[1] = rate[1]/(*kel);
    } else {
      // depot infusion
      yp[0] = rate[0]/(*kel);
      yp[1] = rate[0]/(*kel);
      return;
    }
  }
  yp[central] = rate[central]/(*kel);
}

static inline void comp1ssInf(rx_solving_options_ind *ind, double *yp, double *dur, double *ii, double *rate,
                              double *ka, double *kel) {
  rx_solve *rx=(&rx_global);
  int central = 0;
  if (rx->linKa) {
    // has depot
    if (isSameTime(rate[0], 0.0)) {
      // rate[0]; assume central infusion
      central = 1;
      yp[0]   = 0.0;
      rate[0] = 0.0;
    } else {
      // depot infusion
      double eKa = exp(-(*ka)*((*ii)-(*dur)))/(1.0-exp(-(*ii)*(*ka)));
      double eiKa = exp(-(*ka)*(*dur));
      double eiK = exp(-(*kel)*(*dur));
      double eK = exp(-(*kel)*((*ii)-(*dur)))/(1.0-exp(-(*ii)*(*kel)));
      double rDepot = rate[0];
      yp[0] = eKa*((rDepot)/(*ka) - eiKa*(rDepot)/(*ka));
      yp[1] = eK*((rDepot)/(*kel) +
                  eiKa*(rDepot)/(-(*kel) + (*ka)) -
                  eiK*(rDepot)*(*ka)/((*ka)*(*kel) - (*kel)*(*kel))) +
        (*ka)*(eK - eKa)*((rDepot)/(*ka) -
                          eiKa*(rDepot)/(*ka))/(-(*kel) + (*ka));
      double lag = ind->linCmtLag[0];
      if (lag + *dur < *ii) {
        rate[0] = 0.0;
      }
      rate[1] = 0.0;
      return;
    }
  }
  // central
  double eiK = exp(-(*kel)*(*dur));
  double eK = exp(-(*kel)*((*ii)-(*dur)))/(1.0-exp(-(*kel)*(*ii)));
  double rCentral = rate[central];
  yp[central]=eK*((rCentral)/(*kel) - eiK*(rCentral)*(-(*kel) + (*ka))/((*ka)*(*kel) - (*kel)*(*kel)));
  double lag = ind->linCmtLag[central];
  // lag + dur < ii
  if (lag + *dur < *ii) {
    // should be off
    rate[central] = 0.0;
  }
}

static inline void comp1ssBolus(int *cmtOff, double *yp, double *ii, double *dose, double *ka, double *kel) {
  rx_solve *rx=(&rx_global);
  int central = 0;
  if (rx->linKa) {
    // has depot
    if (*cmtOff == 0) {
      // dosing to depot
      double eKa = 1.0/(1.0-exp(-(*ii)*(*ka)));
      double eK =  1.0/(1.0-exp(-(*ii)*(*kel)));
      yp[0]=eKa*(*dose);
      yp[1]=(*ka)*(*dose)*(eK - eKa)/((*ka)-(*kel));
      return;
    } else {
      // dosing to central
      yp[0] = 0.0;
      central = 1;
    }
  }
  // dossing to central
  double eT = 1.0/(1.0-exp(-(*kel)*(*ii)));
  yp[central] = (*dose)*eT;
}

//////////////////////////////////////////////////
// 2cmt
//////////////////////////////////////////////////
// Infinite infusion
static inline void comp2ssInf8(double *yp, double *rate, double *ka, double *k10,
                               double *k12, double *k21) {
  rx_solve *rx=(&rx_global);
  int central = 0;
  if (rx->linKa) {
    // has depot
    if (isSameTime(rate[0], 0.0)) {
      // central infusion
      central = 1;
      yp[0] = 0.0;
      rate[0] = 0.0;
    } else {
      // depot infusion
      double s = (*k12)+(*k21)+(*k10);
      double beta  = 0.5*(s - sqrt(s*s - 4*(*k21)*(*k10)));
      double alpha = (*k21)*(*k10)/beta;
      yp[0]=(rate[0])/(*ka);
      yp[1]=(rate[0])*(*k21)/(beta*alpha);
      yp[2]=(rate[0])*(*k12)/(beta*alpha);
      rate[0] = rate[1] = 0.0;
      return;
    }
  }
  // central only
  double E1 = (*k10)+(*k12);
  double s = E1+(*k21);
  double sqr = sqrt(s*s-4*(E1*(*k21)-(*k12)*(*k21)));
  double lambda1 = 0.5*(s+sqr);
  double lambda2 = 0.5*(s-sqr);
  double l12 = 1.0/(lambda1*lambda2);
  yp[central]=(rate[central])*(*k21)*l12; // central
  yp[central+1]=(rate[central])*(*k12)*l12; // periph
  rate[central] = 0.0;
}

static inline void comp2ssInf(rx_solving_options_ind *ind, double *yp, double *dur, double *ii, double *rate,
                              double *ka, double *k10, double *k12, double *k21) {
  rx_solve *rx=(&rx_global);
  int central = 0;
  if (rx->linKa) {
    if (isSameTime(rate[0], 0.0)) {
      // central infusion
      central = 1;
      yp[0] = 0.0;
      rate[0] = 0.0;
    } else {
      // depot infusion
      double s = (*k12)+(*k21)+(*k10);
      double beta  = 0.5*(s - sqrt(s*s - 4*(*k21)*(*k10)));
      double alpha = (*k21)*(*k10)/beta;

      double eA = exp(-alpha*((*ii)-(*dur)))/(1.0-exp(-alpha*(*ii)));
      double eB = exp(-beta*((*ii)-(*dur)))/(1.0-exp(-beta*(*ii)));

      double eiA = exp(-alpha*(*dur));
      double eiB = exp(-beta*(*dur));

      double alpha2 = alpha*alpha;
      double alpha3 = alpha2*alpha;

      double beta2 = beta*beta;
      double beta3 = beta2*beta;

      double ka2 = (*ka)*(*ka);

      double eKa = exp(-(*ka)*((*ii)-(*dur)))/(1.0-exp(-(*ka)*(*ii)));
      double eiKa = exp(-(*ka)*(*dur));

      yp[0]=eKa*((rate[0])/(*ka) - eiKa*(rate[0])/(*ka));
      yp[1]=(eA*(-alpha*((rate[0])*(*k21)/(beta*alpha) + eiKa*(rate[0])*(-(*k21) + (*ka))/(beta*alpha + (*ka)*(-alpha - beta) + ka2) - eiA*(rate[0])*(*ka)*(-alpha + (*k21))/(-beta*alpha2 + (*ka)*(beta*alpha - alpha2) + alpha3) + eiB*(rate[0])*(*ka)*(-beta + (*k21))/(beta2*alpha + (*ka)*(-beta*alpha + beta2) - beta3)) + (*k21)*((rate[0])*(*k12)/(beta*alpha) - eiKa*(rate[0])*(*k12)/(beta*alpha + (*ka)*(-alpha - beta) + ka2) - eiA*(rate[0])*(*ka)*(*k12)/(-beta*alpha2 + (*ka)*(beta*alpha - alpha2) + alpha3) + eiB*(rate[0])*(*ka)*(*k12)/(beta2*alpha + (*ka)*(-beta*alpha + beta2) - beta3)) + (*k21)*((rate[0])*(*k21)/(beta*alpha) + eiKa*(rate[0])*(-(*k21) + (*ka))/(beta*alpha + (*ka)*(-alpha - beta) + ka2) - eiA*(rate[0])*(*ka)*(-alpha + (*k21))/(-beta*alpha2 + (*ka)*(beta*alpha - alpha2) + alpha3) + eiB*(rate[0])*(*ka)*(-beta + (*k21))/(beta2*alpha + (*ka)*(-beta*alpha + beta2) - beta3))) - eB*(-beta*((rate[0])*(*k21)/(beta*alpha) + eiKa*(rate[0])*(-(*k21) + (*ka))/(beta*alpha + (*ka)*(-alpha - beta) + ka2) - eiA*(rate[0])*(*ka)*(-alpha + (*k21))/(-beta*alpha2 + (*ka)*(beta*alpha - alpha2) + alpha3) + eiB*(rate[0])*(*ka)*(-beta + (*k21))/(beta2*alpha + (*ka)*(-beta*alpha + beta2) - beta3)) + (*k21)*((rate[0])*(*k12)/(beta*alpha) - eiKa*(rate[0])*(*k12)/(beta*alpha + (*ka)*(-alpha - beta) + ka2) - eiA*(rate[0])*(*ka)*(*k12)/(-beta*alpha2 + (*ka)*(beta*alpha - alpha2) + alpha3) + eiB*(rate[0])*(*ka)*(*k12)/(beta2*alpha + (*ka)*(-beta*alpha + beta2) - beta3)) + (*k21)*((rate[0])*(*k21)/(beta*alpha) + eiKa*(rate[0])*(-(*k21) + (*ka))/(beta*alpha + (*ka)*(-alpha - beta) + ka2) - eiA*(rate[0])*(*ka)*(-alpha + (*k21))/(-beta*alpha2 + (*ka)*(beta*alpha - alpha2) + alpha3) + eiB*(rate[0])*(*ka)*(-beta + (*k21))/(beta2*alpha + (*ka)*(-beta*alpha + beta2) - beta3))))/(-alpha + beta) + (*ka)*(eA*(-alpha + (*k21))/((-alpha + beta)*(-alpha + (*ka))) + eB*(-beta + (*k21))/((-beta + (*ka))*(alpha - beta)) + eKa*((*k21) - (*ka))/((beta - (*ka))*(alpha - (*ka))))*((rate[0])/(*ka) - eiKa*(rate[0])/(*ka));
      yp[2]=(eA*(-alpha*((rate[0])*(*k12)/(beta*alpha) - eiKa*(rate[0])*(*k12)/(beta*alpha + (*ka)*(-alpha - beta) + ka2) - eiA*(rate[0])*(*ka)*(*k12)/(-beta*alpha2 + (*ka)*(beta*alpha - alpha2) + alpha3) + eiB*(rate[0])*(*ka)*(*k12)/(beta2*alpha + (*ka)*(-beta*alpha + beta2) - beta3)) + (*k12)*((rate[0])*(*k21)/(beta*alpha) + eiKa*(rate[0])*(-(*k21) + (*ka))/(beta*alpha + (*ka)*(-alpha - beta) + ka2) - eiA*(rate[0])*(*ka)*(-alpha + (*k21))/(-beta*alpha2 + (*ka)*(beta*alpha - alpha2) + alpha3) + eiB*(rate[0])*(*ka)*(-beta + (*k21))/(beta2*alpha + (*ka)*(-beta*alpha + beta2) - beta3)) + ((*k10) + (*k12))*((rate[0])*(*k12)/(beta*alpha) - eiKa*(rate[0])*(*k12)/(beta*alpha + (*ka)*(-alpha - beta) + ka2) - eiA*(rate[0])*(*ka)*(*k12)/(-beta*alpha2 + (*ka)*(beta*alpha - alpha2) + alpha3) + eiB*(rate[0])*(*ka)*(*k12)/(beta2*alpha + (*ka)*(-beta*alpha + beta2) - beta3))) - eB*(-beta*((rate[0])*(*k12)/(beta*alpha) - eiKa*(rate[0])*(*k12)/(beta*alpha + (*ka)*(-alpha - beta) + ka2) - eiA*(rate[0])*(*ka)*(*k12)/(-beta*alpha2 + (*ka)*(beta*alpha - alpha2) + alpha3) + eiB*(rate[0])*(*ka)*(*k12)/(beta2*alpha + (*ka)*(-beta*alpha + beta2) - beta3)) + (*k12)*((rate[0])*(*k21)/(beta*alpha) + eiKa*(rate[0])*(-(*k21) + (*ka))/(beta*alpha + (*ka)*(-alpha - beta) + ka2) - eiA*(rate[0])*(*ka)*(-alpha + (*k21))/(-beta*alpha2 + (*ka)*(beta*alpha - alpha2) + alpha3) + eiB*(rate[0])*(*ka)*(-beta + (*k21))/(beta2*alpha + (*ka)*(-beta*alpha + beta2) - beta3)) + ((*k10) + (*k12))*((rate[0])*(*k12)/(beta*alpha) - eiKa*(rate[0])*(*k12)/(beta*alpha + (*ka)*(-alpha - beta) + ka2) - eiA*(rate[0])*(*ka)*(*k12)/(-beta*alpha2 + (*ka)*(beta*alpha - alpha2) + alpha3) + eiB*(rate[0])*(*ka)*(*k12)/(beta2*alpha + (*ka)*(-beta*alpha + beta2) - beta3))))/(-alpha + beta) + (*ka)*(*k12)*(eA/((-alpha + beta)*(-alpha + (*ka))) + eB/((-beta + (*ka))*(alpha - beta)) + eKa/((beta - (*ka))*(alpha - (*ka))))*((rate[0])/(*ka) - eiKa*(rate[0])/(*ka));
      rate[1] = 0.0;
      double lag = ind->linCmtLag[0];
      if (lag + *dur < *ii) {
        rate[0] = 0.0;
      }
      return;
    }
  }
  // central infusion
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
  yp[central]=(eT1*(E2*((eTi1*(*rate) - eTi2*(*rate))/(-L1 + L2) + (*rate)*E2*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))) - L1*((eTi1*(*rate) - eTi2*(*rate))/(-L1 + L2) + (*rate)*E2*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))) + (*rate)*(*k12)*(*k21)*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))) - eT2*(E2*((eTi1*(*rate) - eTi2*(*rate))/(-L1 + L2) + (*rate)*E2*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))) - L2*((eTi1*(*rate) - eTi2*(*rate))/(-L1 + L2) + (*rate)*E2*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))) + (*rate)*(*k12)*(*k21)*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))))/(-L1 + L2);
  yp[central+1]=(eT1*((*k12)*((eTi1*(*rate) - eTi2*(*rate))/(-L1 + L2) + (*rate)*E2*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))) + (*rate)*E1*(*k12)*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2)) - (*rate)*(*k12)*L1*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))) - eT2*((*k12)*((eTi1*(*rate) - eTi2*(*rate))/(-L1 + L2) + (*rate)*E2*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))) + (*rate)*E1*(*k12)*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2)) - (*rate)*(*k12)*L2*(1.0/(L1*L2) + eTi1/((L1 - L2)*L1) - eTi2/((L1 - L2)*L2))))/(-L1 + L2);
  double lag = ind->linCmtLag[central];
  if (lag + *dur < *ii) {
    rate[central] = 0.0;
  }
}

// Steady state central dosing
static inline void comp2ssBolus(int *cmtOff, double *yp, double *ii, double *dose,
                                double *ka, double *k10, double *k12, double *k21) {
  rx_solve *rx=(&rx_global);
  int central = 0;
  if (rx->linKa) {
    // has depot
    if (*cmtOff == 0) {
      // dosing to depot
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
      return;
    } else {
      central = 1;
      yp[0] = 0.0;
    }
  }
  double E2 = (*k21);

  double s = (*k12)+(*k21)+(*k10);
  double sqr = sqrt(s*s-4*(*k21)*(*k10));
  double L1 = 0.5*(s+sqr);
  double L2 = 0.5*(s-sqr);

  double eL1 = 1.0/(1.0-exp(-(*ii)*L1));
  double eL2 = 1.0/(1.0-exp(-(*ii)*L2));

  yp[central]=(eL1*((*dose)*E2 - (*dose)*L1) - eL2*((*dose)*E2 - (*dose)*L2))/(-L1 + L2);
  yp[central+1]=(eL1*(*dose)*(*k12) - eL2*(*dose)*(*k12))/(-L1 + L2);
}
//////////////////////////////////////////////////
// 3cmt
//////////////////////////////////////////////////
static inline void comp3ssInf8(double *yp, double *rate,
                               double *ka, double *k10,
                               double *k12, double *k21,
                               double *k13, double *k31) {
  rx_solve *rx=(&rx_global);
  int central = 0;
  if (rx->linKa) {
    if (isSameTime(rate[0], 0.0)) {
      // central infusion
      central=1;
      yp[0]=0.0;
    } else {
      // depot infusion
      double j = (*k12)+(*k10)+(*k21)+(*k31)+(*k13);
      double k = (*k12)*(*k31)+(*k10)*(*k21)+(*k10)*(*k31)+(*k21)*(*k31)+(*k13)*(*k21);
      double l = (*k10)*(*k21)*(*k31);

      double m = 0.3333333333333333*(3.0*k- j*j);
      double n = 0.03703703703703703*(2.0*j*j*j - 9.0*j*k + 27.0*l);
      double Q = 0.25*n*n + 0.03703703703703703*m*m*m;

      double alpha = sqrt(-Q);
      double beta = -0.5*n;
      double rho=sqrt(beta*beta+alpha*alpha);
      double theta = atan2(alpha,beta);
      double ct3 = cos(0.3333333333333333*theta);
      double rho3 = R_pow(rho,0.3333333333333333);
      double st3 = 1.732050807568877193177*sin(0.3333333333333333*theta);
      double j3 = 0.3333333333333333*j;
      double lam1 = j3  + rho3*(ct3 + st3);
      double lam2 = j3 + rho3*(ct3 - st3);
      double lam3 = j3 -(2.0*rho3*ct3);
      double l123 = 1.0/(lam1*lam2*lam3);
      yp[0]=(rate[0])/(*ka);
      yp[1]=(rate[0])*(*k31)*(*k21)*l123;
      yp[2]=(rate[0])*(*k31)*(*k12)*l123;
      yp[3]=(rate[0])*(*k13)*(*k21)*l123;
      return;
    }
  }
  double E1 = (*k10)+(*k12)+(*k13);
  double E2 = (*k21);
  double E3 = (*k31);

  double a = E1+E2+E3;
  double b = E1*E2+E3*(E1+E2)-(*k12)*(*k21)-(*k13)*(*k31);
  double c = E1*E2*E3-E3*(*k12)*(*k21)-E2*(*k13)*(*k31);

  double a2 = a*a;
  double m = 0.333333333333333*(3.0*b - a2);
  double n = 0.03703703703703703*(2.0*a2*a - 9.0*a*b + 27.0*c);
  double Q = 0.25*(n*n) + 0.03703703703703703*(m*m*m);

  double alpha = sqrt(-Q);
  double beta = -0.5*n;
  double gamma = sqrt(beta*beta+alpha*alpha);
  double theta = atan2(alpha,beta);
  double theta3 = 0.333333333333333*theta;
  double ctheta3 = cos(theta3);
  double stheta3 = 1.7320508075688771932*sin(theta3);
  double gamma3 = R_pow(gamma,0.333333333333333);

  double L1 = 0.333333333333333*a + gamma3*(ctheta3 + stheta3);
  double L2 = 0.333333333333333*a + gamma3*(ctheta3 -stheta3);
  double L3 = 0.333333333333333*a -(2.0*gamma3*ctheta3);
  double l123 = 1.0/(L1*L2*L3);
  yp[central]   = (rate[central])*E2*E3*l123;
  yp[central+1] = (rate[central])*E3*(*k12)*l123;
  yp[central+2] = (rate[central])*E2*(*k13)*l123;
}

static inline void comp3ssInf(rx_solving_options_ind *ind, double *yp, double *dur, double *ii, double *rate,
                              double *ka, double *k10,
                              double *k12, double *k21,
                              double *k13, double *k31) {
  rx_solve *rx=(&rx_global);
  int central = 0;
  if (rx->linKa) {
    if (isSameTime(rate[0], 0.0)) {
      // central infusion
      central = 1;
      yp[0] = 0.0;
    } else {
      // depot infusion
      double E2 = (*k10)+ (*k12) + (*k13);
      double E3 = (*k21);
      double E4 = (*k31);
      double j = (*k12)+(*k10)+(*k21)+(*k31)+(*k13);
      double k = (*k12)*(*k31)+(*k10)*(*k21)+(*k10)*(*k31)+(*k21)*(*k31)+(*k13)*(*k21);
      double l = (*k10)*(*k21)*(*k31);

      double m = 0.3333333333333333*(3.0*k - j*j);
      double n = 0.03703703703703703*(2.0*j*j*j - 9.0*j*k + 27.0*l);
      double Q = 0.25*n*n + 0.03703703703703703*m*m*m;

      double alpha = sqrt(-Q);
      double beta = -0.5*n;
      double rho=sqrt(beta*beta+alpha*alpha);
      double theta = atan2(alpha,beta);
      double ct3 = cos(0.3333333333333333*theta);
      double rho3 = R_pow(rho,0.3333333333333333);
      double st3 = 1.732050807568877193177*sin(0.3333333333333333*theta);
      double j3 = 0.3333333333333333*j;
      double lam1 = j3  + rho3*(ct3 + st3);
      double lam2 = j3 + rho3*(ct3 - st3);
      double lam3 = j3 -(2.0*rho3*ct3);

      double eKa = exp(-(*ka)*((*ii)-(*dur)))/(1.0-exp(-(*ka)*(*ii)));
      double eiKa = exp(-(*ka)*(*dur));

      double eL1 = exp(-lam1*((*ii)-(*dur)))/(1.0-exp(-lam1*(*ii)));
      double eiL1 = exp(-lam1*(*dur));

      double eL2 = exp(-lam2*((*ii)-(*dur)))/(1.0-exp(-lam2*(*ii)));
      double eiL2 = exp(-lam2*(*dur));

      double eL3 = exp(-lam3*((*ii)-(*dur)))/(1.0-exp(-lam3*(*ii)));
      double eiL3 = exp(-lam3*(*dur));

      double ka2 = (*ka)*(*ka);
      double ka3 = ka2*(*ka);

      double lam12 = lam1*lam1;
      double lam13 = lam12*lam1;
      double lam14 = lam13*lam1;

      double lam22 = lam2*lam2;
      double lam23 = lam22*lam2;
      double lam24 = lam23*lam2;

      double lam32 = lam3*lam3;
      double lam33 = lam32*lam3;
      double lam34 = lam33*lam3;
      yp[0]=eKa*((rate[0])/(*ka) - eiKa*(rate[0])/(*ka));
      yp[1]=(eL1*(E4 - lam1)*(E3 - lam1)/((-lam1 + lam3)*(-lam1 + lam2)) + eL2*(E4 - lam2)*(E3 - lam2)/((-lam2 + lam3)*(lam1 - lam2)) + eL3*(E3 - lam3)*(E4 - lam3)/((lam2 - lam3)*(lam1 - lam3)))*(-eiKa*(rate[0])*((*k31)*(*k21) + (-(*k21) - (*k31))*(*ka) + ka2)/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) - eiL1*(rate[0])*(-(*ka)*lam12 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) + eiL2*(rate[0])*(-(*ka)*lam22 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) - eiL3*(rate[0])*(-(*ka)*lam32 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k21)/(lam2*lam1*lam3)) + (*ka)*(eL1*(E4 - lam1)*(E3 - lam1)/((-lam1 + lam3)*(-lam1 + lam2)*((*ka) - lam1)) + eL2*(E4 - lam2)*(E3 - lam2)/((-lam2 + lam3)*(lam1 - lam2)*((*ka) - lam2)) + eL3*(E3 - lam3)*(E4 - lam3)/((lam2 - lam3)*(lam1 - lam3)*((*ka) - lam3)) + eKa*(E3 - (*ka))*(E4 - (*ka))/((-(*ka) + lam1)*(-(*ka) + lam2)*(-(*ka) + lam3)))*((rate[0])/(*ka) - eiKa*(rate[0])/(*ka)) + eL1*(-lam1*((*k21)*(eiKa*(rate[0])*(-(*k31)*(*k12) + (*ka)*(*k12))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*(-(*ka)*(*k12)*lam1 + (*ka)*(*k31)*(*k12))/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*(-(*ka)*(*k12)*lam2 + (*ka)*(*k31)*(*k12))/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*(-(*ka)*(*k12)*lam3 + (*ka)*(*k31)*(*k12))/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k12)/(lam2*lam1*lam3)) + (*k31)*(eiKa*(rate[0])*(-(*k13)*(*k21) + (*ka)*(*k13))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k13)*(*k21)/(lam2*lam1*lam3))) + E3*(*k31)*(eiKa*(rate[0])*(-(*k13)*(*k21) + (*ka)*(*k13))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k13)*(*k21)/(lam2*lam1*lam3)) + E4*(*k21)*(eiKa*(rate[0])*(-(*k31)*(*k12) + (*ka)*(*k12))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*(-(*ka)*(*k12)*lam1 + (*ka)*(*k31)*(*k12))/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*(-(*ka)*(*k12)*lam2 + (*ka)*(*k31)*(*k12))/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*(-(*ka)*(*k12)*lam3 + (*ka)*(*k31)*(*k12))/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k12)/(lam2*lam1*lam3)))/((lam1 - lam3)*(lam1 - lam2)) + eL3*(lam3*((*k21)*(eiKa*(rate[0])*(-(*k31)*(*k12) + (*ka)*(*k12))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*(-(*ka)*(*k12)*lam1 + (*ka)*(*k31)*(*k12))/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*(-(*ka)*(*k12)*lam2 + (*ka)*(*k31)*(*k12))/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*(-(*ka)*(*k12)*lam3 + (*ka)*(*k31)*(*k12))/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k12)/(lam2*lam1*lam3)) + (*k31)*(eiKa*(rate[0])*(-(*k13)*(*k21) + (*ka)*(*k13))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k13)*(*k21)/(lam2*lam1*lam3))) - (E3*(*k31)*(eiKa*(rate[0])*(-(*k13)*(*k21) + (*ka)*(*k13))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k13)*(*k21)/(lam2*lam1*lam3)) + E4*(*k21)*(eiKa*(rate[0])*(-(*k31)*(*k12) + (*ka)*(*k12))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*(-(*ka)*(*k12)*lam1 + (*ka)*(*k31)*(*k12))/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*(-(*ka)*(*k12)*lam2 + (*ka)*(*k31)*(*k12))/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*(-(*ka)*(*k12)*lam3 + (*ka)*(*k31)*(*k12))/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k12)/(lam2*lam1*lam3))))/((-lam2 + lam3)*(lam1 - lam3)) + eL2*(lam2*((*k21)*(eiKa*(rate[0])*(-(*k31)*(*k12) + (*ka)*(*k12))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*(-(*ka)*(*k12)*lam1 + (*ka)*(*k31)*(*k12))/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*(-(*ka)*(*k12)*lam2 + (*ka)*(*k31)*(*k12))/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*(-(*ka)*(*k12)*lam3 + (*ka)*(*k31)*(*k12))/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k12)/(lam2*lam1*lam3)) + (*k31)*(eiKa*(rate[0])*(-(*k13)*(*k21) + (*ka)*(*k13))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k13)*(*k21)/(lam2*lam1*lam3))) - (E3*(*k31)*(eiKa*(rate[0])*(-(*k13)*(*k21) + (*ka)*(*k13))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k13)*(*k21)/(lam2*lam1*lam3)) + E4*(*k21)*(eiKa*(rate[0])*(-(*k31)*(*k12) + (*ka)*(*k12))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*(-(*ka)*(*k12)*lam1 + (*ka)*(*k31)*(*k12))/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*(-(*ka)*(*k12)*lam2 + (*ka)*(*k31)*(*k12))/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*(-(*ka)*(*k12)*lam3 + (*ka)*(*k31)*(*k12))/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k12)/(lam2*lam1*lam3))))/((lam2 - lam3)*(lam1 - lam2));
      yp[2]=(eL1*(E4 - lam1)*(E2 - lam1)/((-lam1 + lam3)*(-lam1 + lam2)) + eL2*(E4 - lam2)*(E2 - lam2)/((-lam2 + lam3)*(lam1 - lam2)) + eL3*(E2 - lam3)*(E4 - lam3)/((lam2 - lam3)*(lam1 - lam3)))*(eiKa*(rate[0])*(-(*k31)*(*k12) + (*ka)*(*k12))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*(-(*ka)*(*k12)*lam1 + (*ka)*(*k31)*(*k12))/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*(-(*ka)*(*k12)*lam2 + (*ka)*(*k31)*(*k12))/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*(-(*ka)*(*k12)*lam3 + (*ka)*(*k31)*(*k12))/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k12)/(lam2*lam1*lam3)) + (*ka)*(*k12)*(eL1*(E4 - lam1)/((-lam1 + lam3)*(-lam1 + lam2)*((*ka) - lam1)) + eL2*(E4 - lam2)/((-lam2 + lam3)*(lam1 - lam2)*((*ka) - lam2)) + eL3*(E4 - lam3)/((lam2 - lam3)*(lam1 - lam3)*((*ka) - lam3)) + eKa*(E4 - (*ka))/((-(*ka) + lam1)*(-(*ka) + lam2)*(-(*ka) + lam3)))*((rate[0])/(*ka) - eiKa*(rate[0])/(*ka)) + eL1*(E4*(*k12)*(-eiKa*(rate[0])*((*k31)*(*k21) + (-(*k21) - (*k31))*(*ka) + ka2)/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) - eiL1*(rate[0])*(-(*ka)*lam12 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) + eiL2*(rate[0])*(-(*ka)*lam22 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) - eiL3*(rate[0])*(-(*ka)*lam32 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k21)/(lam2*lam1*lam3)) - (*k12)*lam1*(-eiKa*(rate[0])*((*k31)*(*k21) + (-(*k21) - (*k31))*(*ka) + ka2)/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) - eiL1*(rate[0])*(-(*ka)*lam12 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) + eiL2*(rate[0])*(-(*ka)*lam22 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) - eiL3*(rate[0])*(-(*ka)*lam32 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k21)/(lam2*lam1*lam3)) + (*k31)*(*k12)*(eiKa*(rate[0])*(-(*k13)*(*k21) + (*ka)*(*k13))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k13)*(*k21)/(lam2*lam1*lam3)) - (*k31)*(*k13)*(eiKa*(rate[0])*(-(*k31)*(*k12) + (*ka)*(*k12))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*(-(*ka)*(*k12)*lam1 + (*ka)*(*k31)*(*k12))/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*(-(*ka)*(*k12)*lam2 + (*ka)*(*k31)*(*k12))/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*(-(*ka)*(*k12)*lam3 + (*ka)*(*k31)*(*k12))/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k12)/(lam2*lam1*lam3)))/((lam1 - lam3)*(lam1 - lam2)) + eL3*((*k12)*lam3*(-eiKa*(rate[0])*((*k31)*(*k21) + (-(*k21) - (*k31))*(*ka) + ka2)/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) - eiL1*(rate[0])*(-(*ka)*lam12 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) + eiL2*(rate[0])*(-(*ka)*lam22 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) - eiL3*(rate[0])*(-(*ka)*lam32 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k21)/(lam2*lam1*lam3)) - (E4*(*k12)*(-eiKa*(rate[0])*((*k31)*(*k21) + (-(*k21) - (*k31))*(*ka) + ka2)/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) - eiL1*(rate[0])*(-(*ka)*lam12 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) + eiL2*(rate[0])*(-(*ka)*lam22 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) - eiL3*(rate[0])*(-(*ka)*lam32 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k21)/(lam2*lam1*lam3)) + (*k31)*(*k12)*(eiKa*(rate[0])*(-(*k13)*(*k21) + (*ka)*(*k13))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k13)*(*k21)/(lam2*lam1*lam3)) - (*k31)*(*k13)*(eiKa*(rate[0])*(-(*k31)*(*k12) + (*ka)*(*k12))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*(-(*ka)*(*k12)*lam1 + (*ka)*(*k31)*(*k12))/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*(-(*ka)*(*k12)*lam2 + (*ka)*(*k31)*(*k12))/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*(-(*ka)*(*k12)*lam3 + (*ka)*(*k31)*(*k12))/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k12)/(lam2*lam1*lam3))))/((-lam2 + lam3)*(lam1 - lam3)) + eL2*((*k12)*lam2*(-eiKa*(rate[0])*((*k31)*(*k21) + (-(*k21) - (*k31))*(*ka) + ka2)/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) - eiL1*(rate[0])*(-(*ka)*lam12 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) + eiL2*(rate[0])*(-(*ka)*lam22 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) - eiL3*(rate[0])*(-(*ka)*lam32 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k21)/(lam2*lam1*lam3)) - (E4*(*k12)*(-eiKa*(rate[0])*((*k31)*(*k21) + (-(*k21) - (*k31))*(*ka) + ka2)/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) - eiL1*(rate[0])*(-(*ka)*lam12 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) + eiL2*(rate[0])*(-(*ka)*lam22 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) - eiL3*(rate[0])*(-(*ka)*lam32 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k21)/(lam2*lam1*lam3)) + (*k31)*(*k12)*(eiKa*(rate[0])*(-(*k13)*(*k21) + (*ka)*(*k13))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k13)*(*k21)/(lam2*lam1*lam3)) - (*k31)*(*k13)*(eiKa*(rate[0])*(-(*k31)*(*k12) + (*ka)*(*k12))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*(-(*ka)*(*k12)*lam1 + (*ka)*(*k31)*(*k12))/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*(-(*ka)*(*k12)*lam2 + (*ka)*(*k31)*(*k12))/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*(-(*ka)*(*k12)*lam3 + (*ka)*(*k31)*(*k12))/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k12)/(lam2*lam1*lam3))))/((lam2 - lam3)*(lam1 - lam2));
      yp[3]=(eL1*(E3 - lam1)*(E2 - lam1)/((-lam1 + lam3)*(-lam1 + lam2)) + eL2*(E2 - lam2)*(E3 - lam2)/((-lam2 + lam3)*(lam1 - lam2)) + eL3*(E2 - lam3)*(E3 - lam3)/((lam2 - lam3)*(lam1 - lam3)))*(eiKa*(rate[0])*(-(*k13)*(*k21) + (*ka)*(*k13))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k13)*(*k21)/(lam2*lam1*lam3)) + (*ka)*(*k13)*(eL1*(E3 - lam1)/((-lam1 + lam3)*(-lam1 + lam2)*((*ka) - lam1)) + eL2*(E3 - lam2)/((-lam2 + lam3)*(lam1 - lam2)*((*ka) - lam2)) + eL3*(E3 - lam3)/((lam2 - lam3)*(lam1 - lam3)*((*ka) - lam3)) + eKa*(E3 - (*ka))/((-(*ka) + lam1)*(-(*ka) + lam2)*(-(*ka) + lam3)))*((rate[0])/(*ka) - eiKa*(rate[0])/(*ka)) + eL1*(E3*(*k13)*(-eiKa*(rate[0])*((*k31)*(*k21) + (-(*k21) - (*k31))*(*ka) + ka2)/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) - eiL1*(rate[0])*(-(*ka)*lam12 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) + eiL2*(rate[0])*(-(*ka)*lam22 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) - eiL3*(rate[0])*(-(*ka)*lam32 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k21)/(lam2*lam1*lam3)) - (*k12)*(*k21)*(eiKa*(rate[0])*(-(*k13)*(*k21) + (*ka)*(*k13))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k13)*(*k21)/(lam2*lam1*lam3)) + (*k13)*(*k21)*(eiKa*(rate[0])*(-(*k31)*(*k12) + (*ka)*(*k12))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*(-(*ka)*(*k12)*lam1 + (*ka)*(*k31)*(*k12))/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*(-(*ka)*(*k12)*lam2 + (*ka)*(*k31)*(*k12))/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*(-(*ka)*(*k12)*lam3 + (*ka)*(*k31)*(*k12))/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k12)/(lam2*lam1*lam3)) - (*k13)*lam1*(-eiKa*(rate[0])*((*k31)*(*k21) + (-(*k21) - (*k31))*(*ka) + ka2)/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) - eiL1*(rate[0])*(-(*ka)*lam12 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) + eiL2*(rate[0])*(-(*ka)*lam22 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) - eiL3*(rate[0])*(-(*ka)*lam32 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k21)/(lam2*lam1*lam3)))/((lam1 - lam3)*(lam1 - lam2)) + eL3*((*k13)*lam3*(-eiKa*(rate[0])*((*k31)*(*k21) + (-(*k21) - (*k31))*(*ka) + ka2)/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) - eiL1*(rate[0])*(-(*ka)*lam12 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) + eiL2*(rate[0])*(-(*ka)*lam22 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) - eiL3*(rate[0])*(-(*ka)*lam32 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k21)/(lam2*lam1*lam3)) - (E3*(*k13)*(-eiKa*(rate[0])*((*k31)*(*k21) + (-(*k21) - (*k31))*(*ka) + ka2)/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) - eiL1*(rate[0])*(-(*ka)*lam12 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) + eiL2*(rate[0])*(-(*ka)*lam22 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) - eiL3*(rate[0])*(-(*ka)*lam32 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k21)/(lam2*lam1*lam3)) - (*k12)*(*k21)*(eiKa*(rate[0])*(-(*k13)*(*k21) + (*ka)*(*k13))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k13)*(*k21)/(lam2*lam1*lam3)) + (*k13)*(*k21)*(eiKa*(rate[0])*(-(*k31)*(*k12) + (*ka)*(*k12))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*(-(*ka)*(*k12)*lam1 + (*ka)*(*k31)*(*k12))/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*(-(*ka)*(*k12)*lam2 + (*ka)*(*k31)*(*k12))/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*(-(*ka)*(*k12)*lam3 + (*ka)*(*k31)*(*k12))/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k12)/(lam2*lam1*lam3))))/((-lam2 + lam3)*(lam1 - lam3)) + eL2*((*k13)*lam2*(-eiKa*(rate[0])*((*k31)*(*k21) + (-(*k21) - (*k31))*(*ka) + ka2)/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) - eiL1*(rate[0])*(-(*ka)*lam12 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) + eiL2*(rate[0])*(-(*ka)*lam22 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) - eiL3*(rate[0])*(-(*ka)*lam32 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k21)/(lam2*lam1*lam3)) - (E3*(*k13)*(-eiKa*(rate[0])*((*k31)*(*k21) + (-(*k21) - (*k31))*(*ka) + ka2)/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) - eiL1*(rate[0])*(-(*ka)*lam12 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) + eiL2*(rate[0])*(-(*ka)*lam22 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) - eiL3*(rate[0])*(-(*ka)*lam32 - (*ka)*(*k31)*(*k21) + ((*k21) + (*k31))*(*ka)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k21)/(lam2*lam1*lam3)) - (*k12)*(*k21)*(eiKa*(rate[0])*(-(*k13)*(*k21) + (*ka)*(*k13))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam1)/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam2)/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*((*ka)*(*k13)*(*k21) - (*ka)*(*k13)*lam3)/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k13)*(*k21)/(lam2*lam1*lam3)) + (*k13)*(*k21)*(eiKa*(rate[0])*(-(*k31)*(*k12) + (*ka)*(*k12))/(ka2*lam1 + lam2*(-(*ka)*lam1 + ka2) + lam3*(-(*ka)*lam1 + lam2*(-(*ka) + lam1) + ka2) - ka3) + eiL1*(rate[0])*(-(*ka)*(*k12)*lam1 + (*ka)*(*k31)*(*k12))/(-(*ka)*lam13 + lam2*((*ka)*lam12 - lam13) + lam3*((*ka)*lam12 + lam2*(-(*ka)*lam1 + lam12) - lam13) + lam14) - eiL2*(rate[0])*(-(*ka)*(*k12)*lam2 + (*ka)*(*k31)*(*k12))/(lam23*((*ka) + lam1) + lam3*(lam22*(-(*ka) - lam1) + (*ka)*lam2*lam1 + lam23) - (*ka)*lam22*lam1 - lam24) + eiL3*(rate[0])*(-(*ka)*(*k12)*lam3 + (*ka)*(*k31)*(*k12))/(lam32*((*ka)*lam1 + lam2*((*ka) + lam1)) + (-(*ka) - lam1 - lam2)*lam33 - (*ka)*lam2*lam1*lam3 + lam34) + (rate[0])*(*k31)*(*k12)/(lam2*lam1*lam3))))/((lam2 - lam3)*(lam1 - lam2));
      return;
    }
  }
  // central infusion
  double E1 = (*k10)+(*k12)+(*k13);
  double E2 = (*k21);
  double E3 = (*k31);

  double a = E1+E2+E3;
  double b = E1*E2+E3*(E1+E2)-(*k12)*(*k21)-(*k13)*(*k31);
  double c = E1*E2*E3-E3*(*k12)*(*k21)-E2*(*k13)*(*k31);

  double a2 = a*a;
  double m = 0.333333333333333*(3.0*b - a2);
  double n = 0.03703703703703703*(2.0*a2*a - 9.0*a*b + 27.0*c);
  double Q = 0.25*(n*n) + 0.03703703703703703*(m*m*m);

  double alpha = sqrt(-Q);
  double beta = -0.5*n;
  double gamma = sqrt(beta*beta+alpha*alpha);
  double theta = atan2(alpha,beta);
  double theta3 = 0.333333333333333*theta;
  double ctheta3 = cos(theta3);
  double stheta3 = 1.7320508075688771932*sin(theta3);
  double gamma3 = R_pow(gamma,0.333333333333333);

  double L1 = 0.333333333333333*a + gamma3*(ctheta3 + stheta3);
  double L2 = 0.333333333333333*a + gamma3*(ctheta3 -stheta3);
  double L3 = 0.333333333333333*a -(2.0*gamma3*ctheta3);
  double eTi1 = exp(-(*dur)*L1);
  double eTi2 = exp(-(*dur)*L2);
  double eTi3 = exp(-(*dur)*L3);
  double eT1 = exp(-L1*((*ii)-(*dur)))/(1.0-exp(-(*ii)*L1));
  double eT2 = exp(-L2*((*ii)-(*dur)))/(1.0-exp(-(*ii)*L2));
  double eT3 = exp(-L3*((*ii)-(*dur)))/(1.0-exp(-(*ii)*L3));
  yp[central] = (rate[central])*(eT1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)) + eT2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)) + eT3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)))*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) + eT2*(L2*((rate[central])*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (rate[central])*(*k13)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3))) - ((rate[central])*E2*(*k13)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (rate[central])*E3*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))))/((L1 - L2)*(L2 - L3)) + eT1*(-L1*((rate[central])*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (rate[central])*(*k13)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3))) + (rate[central])*E2*(*k13)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (rate[central])*E3*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)))/((L1 - L3)*(L1 - L2)) + eT3*(L3*((rate[central])*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (rate[central])*(*k13)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3))) - ((rate[central])*E2*(*k13)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (rate[central])*E3*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))))/((L1 - L3)*(-L2 + L3));
  yp[central+1]=(rate[central])*(*k12)*(eT1*(E3 - L1)*(E1 - L1)/((-L1 + L3)*(-L1 + L2)) + eT2*(E1 - L2)*(E3 - L2)/((L1 - L2)*(-L2 + L3)) + eT3*(E1 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)))*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) + eT2*((rate[central])*(*k12)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))*L2 - ((rate[central])*E3*(*k12)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (rate[central])*(*k13)*(*k12)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) - (rate[central])*(*k13)*(*k12)*(*k31)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))))/((L1 - L2)*(L2 - L3)) + eT1*((rate[central])*E3*(*k12)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) - (rate[central])*(*k12)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))*L1 + (rate[central])*(*k13)*(*k12)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) - (rate[central])*(*k13)*(*k12)*(*k31)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)))/((L1 - L3)*(L1 - L2)) + eT3*((rate[central])*(*k12)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))*L3 - ((rate[central])*E3*(*k12)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (rate[central])*(*k13)*(*k12)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) - (rate[central])*(*k13)*(*k12)*(*k31)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))))/((L1 - L3)*(-L2 + L3));
  yp[central+2] =(rate[central])*(*k13)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3))*(eT1*(E2 - L1)*(E1 - L1)/((-L1 + L3)*(-L1 + L2)) + eT2*(E1 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)) + eT3*(E1 - L3)*(E2 - L3)/((L1 - L3)*(L2 - L3))) + eT2*((rate[central])*(*k13)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))*L2 - ((rate[central])*E2*(*k13)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) - (rate[central])*(*k13)*(*k12)*(*k21)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (rate[central])*(*k13)*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))))/((L1 - L2)*(L2 - L3)) + eT1*((rate[central])*E2*(*k13)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) - (rate[central])*(*k13)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))*L1 - (rate[central])*(*k13)*(*k12)*(*k21)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (rate[central])*(*k13)*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)))/((L1 - L3)*(L1 - L2)) + eT3*((rate[central])*(*k13)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))*L3 - ((rate[central])*E2*(*k13)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) - (rate[central])*(*k13)*(*k12)*(*k21)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (rate[central])*(*k13)*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))))/((L1 - L3)*(-L2 + L3));
}

// Steady state central dosing
static inline void comp3ssBolus(int *cmtOff, double *yp, double *ii, double *dose,
                                double *ka, double *k10,
                                double *k12, double *k21,
                                double *k13, double *k31) {
  rx_solve *rx=(&rx_global);
  int central = 0;
  if (rx->linKa) {
    // has depot
    if (*cmtOff == 0) {
      // dosing to depot
      double E2 = (*k10)+(*k12)+(*k13);
      double E3 = (*k21);
      double E4 = (*k31);

      double a = E2+E3+E4;
      double b = E2*E3+E4*(E2+E3)-(*k12)*(*k21)-(*k13)*(*k31);
      double c = E2*E3*E4-E4*(*k12)*(*k21)-E3*(*k13)*(*k31);

      double a2 = a*a;
      double m = 0.333333333333333*(3.0*b - a2);
      double n = 0.03703703703703703*(2.0*a2*a - 9.0*a*b + 27.0*c);
      double Q = 0.25*(n*n) + 0.03703703703703703*(m*m*m);

      double alpha = sqrt(-Q);
      double beta = -0.5*n;
      double gamma = sqrt(beta*beta+alpha*alpha);
      double theta = atan2(alpha,beta);
      double theta3 = 0.333333333333333*theta;
      double ctheta3 = cos(theta3);
      double stheta3 = 1.7320508075688771932*sin(theta3);
      double gamma3 = R_pow(gamma,0.333333333333333);

      double L1= 0.333333333333333*a + gamma3*(ctheta3 + stheta3);
      double L2 = 0.333333333333333*a + gamma3*(ctheta3 -stheta3);
      double L3 = 0.333333333333333*a -(2.0*gamma3*ctheta3);

      double eKa = 1.0/(1.0-exp(-(*ii)*(*ka)));
      double eL1 = 1.0/(1.0-exp(-(*ii)*L1));
      double eL2 = 1.0/(1.0-exp(-(*ii)*L2));
      double eL3 = 1.0/(1.0-exp(-(*ii)*L3));

      yp[0]=eKa*(*dose);
      yp[1]=(*ka)*(*dose)*(eL1*(E3 - L1)*(E4 - L1)/((-L1 + L3)*(-L1 + L2)*((*ka) - L1)) + eL2*(E3 - L2)*(E4 - L2)/((L1 - L2)*(-L2 + L3)*((*ka) - L2)) + eL3*(E3 - L3)*(E4 - L3)/((L1 - L3)*(L2 - L3)*((*ka) - L3)) + eKa*(E3 - (*ka))*(E4 - (*ka))/((-(*ka) + L1)*(-(*ka) + L3)*(-(*ka) + L2)));
      yp[2]=(*ka)*(*dose)*(*k12)*(eL1*(E4 - L1)/((-L1 + L3)*(-L1 + L2)*((*ka) - L1)) + eL2*(E4 - L2)/((L1 - L2)*(-L2 + L3)*((*ka) - L2)) + eL3*(E4 - L3)/((L1 - L3)*(L2 - L3)*((*ka) - L3)) + eKa*(E4 - (*ka))/((-(*ka) + L1)*(-(*ka) + L3)*(-(*ka) + L2)));
      yp[3]=(*ka)*(*dose)*(*k13)*(eL1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*((*ka) - L1)) + eL2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*((*ka) - L2)) + eL3*(E3 - L3)/((L1 - L3)*(L2 - L3)*((*ka) - L3)) + eKa*(E3 - (*ka))/((-(*ka) + L1)*(-(*ka) + L3)*(-(*ka) + L2)));
      return;
    } else {
      // dosing to central
      central=1;
      yp[0] = 0.0;
    }
  }
  double E2 = (*k10)+(*k12)+(*k13);
  double E3 = (*k21);
  double E4 = (*k31);

  double a = E2+E3+E4;
  double b = E2*E3+E4*(E2+E3)-(*k12)*(*k21)-(*k13)*(*k31);
  double c = E2*E3*E4-E4*(*k12)*(*k21)-E3*(*k13)*(*k31);

  double a2 = a*a;
  double m = 0.333333333333333*(3.0*b - a2);
  double n = 0.03703703703703703*(2.0*a2*a - 9.0*a*b + 27.0*c);
  double Q = 0.25*(n*n) + 0.03703703703703703*(m*m*m);

  double alpha = sqrt(-Q);
  double beta = -0.5*n;
  double gamma = sqrt(beta*beta+alpha*alpha);
  double theta = atan2(alpha,beta);
  double theta3 = 0.333333333333333*theta;
  double ctheta3 = cos(theta3);
  double stheta3 = 1.7320508075688771932*sin(theta3);
  double gamma3 = R_pow(gamma,0.333333333333333);

  double L1 = 0.333333333333333*a + gamma3*(ctheta3 + stheta3);
  double L2 = 0.333333333333333*a + gamma3*(ctheta3 -stheta3);
  double L3 = 0.333333333333333*a -(2.0*gamma3*ctheta3);

  double eL1 = 1.0/(1.0-exp(-(*ii)*L1));
  double eL2 = 1.0/(1.0-exp(-(*ii)*L2));
  double eL3 = 1.0/(1.0-exp(-(*ii)*L3));

  yp[central]=(*dose)*(eL1*(E3 - L1)*(E4 - L1)/((-L1 + L3)*(-L1 + L2)) + eL2*(E3 - L2)*(E4 - L2)/((L1 - L2)*(-L2 + L3)) + eL3*(E3 - L3)*(E4 - L3)/((L1 - L3)*(L2 - L3)));
  yp[central+1]=eL2*(-(*dose)*E4*(*k12) + (*dose)*(*k12)*L2)/((L1 - L2)*(L2 - L3)) + eL1*((*dose)*E4*(*k12) - (*dose)*(*k12)*L1)/((L1 - L3)*(L1 - L2)) + eL3*(-(*dose)*E4*(*k12) + (*dose)*(*k12)*L3)/((L1 - L3)*(-L2 + L3));
  yp[central+2]=eL2*(-(*dose)*E3*(*k13) + (*dose)*(*k13)*L2)/((L1 - L2)*(L2 - L3)) + eL1*((*dose)*E3*(*k13) - (*dose)*(*k13)*L1)/((L1 - L3)*(L1 - L2)) + eL3*(-(*dose)*E3*(*k13) + (*dose)*(*k13)*L3)/((L1 - L3)*(-L2 + L3));
}
