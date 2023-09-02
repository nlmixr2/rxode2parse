//////////////////////////////////////////////////
// 1cmt
//////////////////////////////////////////////////
// Infinite infusion
static inline void comp1ssInf8(double *yp, double *rate, double *ka, double *kel) {
  int hasDepot = (*ka) != 0.0;
  yp[hasDepot] = (*rate)/(*kel);
}

static inline void comp1ssInf(double *yp, double *dur, double *ii, double *rate,
                              double *ka, double *kel) {
  double eiK = exp(-(*kel)*(*dur));
  double eK = exp(-(*kel)*((*ii)-(*dur)))/(1.0-exp(-(*kel)*(*ii)));
  int hasDepot = (*ka) != 0.0;
  yp[hasDepot]=eK*((*rate)/(*kel) - eiK*(*rate)*(-(*kel) + (*ka))/((*ka)*(*kel) - (*kel)*(*kel)));
}

// Steady state central dosing
static inline void comp1ssBolusCentral(double *yp, double *ii, double *dose, double *ka, double *kel) {
  int hasDepot = (*ka) != 0.0;
  double eT = 1.0/(1.0-exp(-(*kel)*(*ii)));
  yp[hasDepot] = (*dose)*eT;
}

// Steady state bolus dosing
static inline void comp1ssBolusDepot(double *yp, double *ii, double *dose, double *ka, double *kel) {
  double eKa = 1.0/(1.0-exp(-(*ii)*(*ka)));
  double eK =  1.0/(1.0-exp(-(*ii)*(*kel)));
  yp[0]=eKa*(*dose);
  yp[1]=(*ka)*(*dose)*(eK - eKa)/((*ka)-(*kel));
}

//////////////////////////////////////////////////
// 2cmt
//////////////////////////////////////////////////
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

//////////////////////////////////////////////////
// 3cmt
//////////////////////////////////////////////////
static inline void comp3ssInf8(double *yp, double *rate,
                               double *ka, double *k10,
                               double *k12, double *k21,
                               double *k13, double *k31) {
  int hasDepot = (*ka) != 0.0;
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
  yp[hasDepot]   = (*rate)*E2*E3*l123;
  yp[hasDepot+1] = (*rate)*E3*(*k12)*l123;
  yp[hasDepot+2] = (*rate)*E2*(*k13)*l123;
}

static inline void comp3ssInf(double *yp, double *dur, double *ii, double *rate,
                              double *ka, double *k10,
                              double *k12, double *k21,
                              double *k13, double *k31) {
  int hasDepot = (*ka) != 0.0;
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
  yp[hasDepot] = (*rate)*(eT1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)) + eT2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)) + eT3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)))*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) + eT2*(L2*((*rate)*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (*rate)*(*k13)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3))) - ((*rate)*E2*(*k13)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (*rate)*E3*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))))/((L1 - L2)*(L2 - L3)) + eT1*(-L1*((*rate)*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (*rate)*(*k13)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3))) + (*rate)*E2*(*k13)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (*rate)*E3*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)))/((L1 - L3)*(L1 - L2)) + eT3*(L3*((*rate)*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (*rate)*(*k13)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3))) - ((*rate)*E2*(*k13)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (*rate)*E3*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))))/((L1 - L3)*(-L2 + L3));
  yp[hasDepot+1]=(*rate)*(*k12)*(eT1*(E3 - L1)*(E1 - L1)/((-L1 + L3)*(-L1 + L2)) + eT2*(E1 - L2)*(E3 - L2)/((L1 - L2)*(-L2 + L3)) + eT3*(E1 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)))*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) + eT2*((*rate)*(*k12)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))*L2 - ((*rate)*E3*(*k12)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (*rate)*(*k13)*(*k12)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) - (*rate)*(*k13)*(*k12)*(*k31)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))))/((L1 - L2)*(L2 - L3)) + eT1*((*rate)*E3*(*k12)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) - (*rate)*(*k12)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))*L1 + (*rate)*(*k13)*(*k12)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) - (*rate)*(*k13)*(*k12)*(*k31)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)))/((L1 - L3)*(L1 - L2)) + eT3*((*rate)*(*k12)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))*L3 - ((*rate)*E3*(*k12)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (*rate)*(*k13)*(*k12)*(*k31)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) - (*rate)*(*k13)*(*k12)*(*k31)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))))/((L1 - L3)*(-L2 + L3));
  yp[hasDepot+2] =(*rate)*(*k13)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3))*(eT1*(E2 - L1)*(E1 - L1)/((-L1 + L3)*(-L1 + L2)) + eT2*(E1 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)) + eT3*(E1 - L3)*(E2 - L3)/((L1 - L3)*(L2 - L3))) + eT2*((*rate)*(*k13)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))*L2 - ((*rate)*E2*(*k13)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) - (*rate)*(*k13)*(*k12)*(*k21)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (*rate)*(*k13)*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))))/((L1 - L2)*(L2 - L3)) + eT1*((*rate)*E2*(*k13)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) - (*rate)*(*k13)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))*L1 - (*rate)*(*k13)*(*k12)*(*k21)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (*rate)*(*k13)*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)))/((L1 - L3)*(L1 - L2)) + eT3*((*rate)*(*k13)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))*L3 - ((*rate)*E2*(*k13)*(E2*E3/(L1*L2*L3) - eTi1*(E2 - L1)*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3)) - (*rate)*(*k13)*(*k12)*(*k21)*(E2/(L1*L2*L3) - eTi1*(E2 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E2 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E2 - L3)/((L1 - L3)*(L2 - L3)*L3)) + (*rate)*(*k13)*(*k12)*(*k21)*(E3/(L1*L2*L3) - eTi1*(E3 - L1)/((-L1 + L3)*(-L1 + L2)*L1) - eTi2*(E3 - L2)/((L1 - L2)*(-L2 + L3)*L2) - eTi3*(E3 - L3)/((L1 - L3)*(L2 - L3)*L3))))/((L1 - L3)*(-L2 + L3));
}

// Steady state central dosing
static inline void comp3ssBolusCentral(double *yp, double *ii, double *dose,
                                       double *ka, double *k10,
                                       double *k12, double *k21,
                                       double *k13, double *k31) {
  int hasDepot = (*ka) != 0.0;
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

  yp[hasDepot]=(*dose)*(eL1*(E3 - L1)*(E4 - L1)/((-L1 + L3)*(-L1 + L2)) + eL2*(E3 - L2)*(E4 - L2)/((L1 - L2)*(-L2 + L3)) + eL3*(E3 - L3)*(E4 - L3)/((L1 - L3)*(L2 - L3)));
  yp[hasDepot+1]=eL2*(-(*dose)*E4*(*k12) + (*dose)*(*k12)*L2)/((L1 - L2)*(L2 - L3)) + eL1*((*dose)*E4*(*k12) - (*dose)*(*k12)*L1)/((L1 - L3)*(L1 - L2)) + eL3*(-(*dose)*E4*(*k12) + (*dose)*(*k12)*L3)/((L1 - L3)*(-L2 + L3));
  yp[hasDepot+2]=eL2*(-(*dose)*E3*(*k13) + (*dose)*(*k13)*L2)/((L1 - L2)*(L2 - L3)) + eL1*((*dose)*E3*(*k13) - (*dose)*(*k13)*L1)/((L1 - L3)*(L1 - L2)) + eL3*(-(*dose)*E3*(*k13) + (*dose)*(*k13)*L3)/((L1 - L3)*(-L2 + L3));
}

// Steady state bolus dosing
static inline void comp3ssBolusDepot(double *yp, double *ii, double *dose,
                                     double *ka, double *k10,
                                     double *k12, double *k21,
                                     double *k13, double *k31) {
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
}
