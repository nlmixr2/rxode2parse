#define USE_FC_LEN_T
#define STRICT_R_HEADERS
#include <stdio.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

#define op_global _rxode2parse_op_global
#define rx_global _rxode2parse_rx_global
#define AMT _rxode2parse_AMT
#define LAG _rxode2parse_LAG
#define RATE _rxode2parse_RATE
#define DUR _rxode2parse_DUR
#define calc_mtime _rxode2parse_calc_mtime
#define getTime_ _rxode2parse_getTime_
#define getTime _rxode2parse_getTime
#define _locateTimeIndex _rxode2parse_locateTimeIndex
#define handle_evidL _rxode2parse_handle_evidL

#include "../inst/include/rxode2parse.h"
#define _calcDerived _rxode2parse_calcDerived

extern rx_solving_options _rxode2parse_op_global;
extern rx_solve _rxode2parse_rx_global;
extern t_handle_evidL _rxode2parse_handle_evidL;
extern t_getDur _rxode2parse_getDur;
#define _getDur _rxode2parse_getDur

#include "../inst/include/rxode2parse.h"
#include "../inst/include/rxode2parseHandleEvid.h"
#include "../inst/include/rxode2parseGetTime.h"

extern t_getTime _rxode2parse_getTime;

extern t_locateTimeIndex _rxode2parse_locateTimeIndex;

#define safe_zero(a) ((a) == 0 ? DBL_EPSILON : (a))
#define _as_zero(a) (fabs(a) < sqrt(DBL_EPSILON) ? 0.0 : a)
#define _as_dbleps(a) (fabs(a) < sqrt(DBL_EPSILON) ? ((a) < 0 ? -sqrt(DBL_EPSILON)  : sqrt(DBL_EPSILON)) : a)

void _rxode2parse_unprotect();


#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("rxode2parse", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

static inline int parTrans(int *trans, 
                           double *p1, double *v1,
                           double *p2, double *p3,
                           double *p4, double *p5,
                           unsigned int *ncmt, double *rx_k, double *rx_v, double *rx_k12,
                           double *rx_k21, double *rx_k13, double *rx_k31){
  double btemp, ctemp, dtemp;
  if ((*p5) > 0.) {
    (*ncmt) = 3;
    switch (*trans) {
    case 1: // cl v q vp
      (*rx_k) = (*p1)/(*v1); // k = CL/V
      (*rx_v) = (*v1);
      (*rx_k12) = (*p2)/(*v1); // k12 = Q/V
      (*rx_k21) = (*p2)/(*p3); // k21 = Q/Vp
      (*rx_k13) = (*p4)/(*v1); // k31 = Q2/V
      (*rx_k31) = (*p4)/(*p5); // k31 = Q2/Vp2
      break;
    case 2: // k=(*p1) v=(*v1) k12=(*p2) k21=(*p3) k13=(*p4) k31=(*p5)
      (*rx_k) = (*p1);
      (*rx_v) = (*v1);
      (*rx_k12) = (*p2);
      (*rx_k21) = (*p3);
      (*rx_k13) = (*p4);
      (*rx_k31) = (*p5);
      break;
    case 11:
#undef beta
#define A (1/(*v1))
#define B (*p3)
#define C (*p5)
#define alpha (*p1)
#define beta (*p2)
#define gamma (*p4)
      (*ncmt)=3;
      (*rx_v)=1/(A+B+C);
      btemp = -(alpha*C + alpha*B + gamma*A + gamma*B + beta*A + beta*C)*(*rx_v);
      ctemp = (alpha*beta*C + alpha*gamma*B + beta*gamma*A)*(*rx_v);
      dtemp = sqrt(btemp*btemp-4*ctemp);
      (*rx_k21) = 0.5*(-btemp+dtemp);
      (*rx_k31) = 0.5*(-btemp-dtemp);
      (*rx_k)   = alpha*beta*gamma/(*rx_k21)/(*rx_k31);
      (*rx_k12) = ((beta*gamma + alpha*beta + alpha*gamma) -
                   (*rx_k21)*(alpha+beta+gamma) - (*rx_k) * (*rx_k31) + (*rx_k21)*(*rx_k21))/((*rx_k31) - (*rx_k21));
      (*rx_k13) = alpha + beta + gamma - ((*rx_k) + (*rx_k12) + (*rx_k21) + (*rx_k31));
      break;
    case 10:
#undef A
#define A (*v1)
      (*ncmt)=3;
      (*rx_v)=1/(A+B+C);
      btemp = -(alpha*C + alpha*B + gamma*A + gamma*B + beta*A + beta*C)*(*rx_v);
      ctemp = (alpha*beta*C + alpha*gamma*B + beta*gamma*A)*(*rx_v);
      dtemp = sqrt(btemp*btemp-4*ctemp);
      (*rx_k21) = 0.5*(-btemp+dtemp);
      (*rx_k31) = 0.5*(-btemp-dtemp);
      (*rx_k)   = alpha*beta*gamma/(*rx_k21)/(*rx_k31);
      (*rx_k12) = ((beta*gamma + alpha*beta + alpha*gamma) -
                   (*rx_k21)*(alpha+beta+gamma) - (*rx_k) * (*rx_k31) + (*rx_k21)*(*rx_k21))/((*rx_k31) - (*rx_k21));
      (*rx_k13) = alpha + beta + gamma - ((*rx_k) + (*rx_k12) + (*rx_k21) + (*rx_k31));
#undef A
#undef B
#undef C
#undef alpha
#undef beta
#undef gamma
#define beta Rf_beta
      break;
    default:
      return NA_REAL;
    }
  } else if ((*p3) > 0.) {
    (*ncmt) = 2;
    switch (*trans){
    case 1: // cl=(*p1) v=(*v1) q=(*p2) vp=(*p3)
      (*rx_k) = (*p1)/(*v1); // k = CL/V
      (*rx_v) = (*v1);
      (*rx_k12) = (*p2)/(*v1); // k12 = Q/V
      (*rx_k21) = (*p2)/(*p3); // k21 = Q/Vp
      break;
    case 2: // k=(*p1), (*v1)=v k12=(*p2) k21=(*p3)
      (*rx_k) = (*p1);
      (*rx_v) = (*v1);
      (*rx_k12) = (*p2);
      (*rx_k21) = (*p3);
      break;
    case 3: // cl=(*p1) v=(*v1) q=(*p2) vss=(*p3)
      (*rx_k) = (*p1)/(*v1); // k = CL/V
      (*rx_v) = (*v1);
      (*rx_k12) = (*p2)/(*v1); // k12 = Q/V
      (*rx_k21) = (*p2)/((*p3)-(*v1)); // k21 = Q/(Vss-V)
      break;
    case 4: // alpha=(*p1) beta=(*p2) k21=(*p3)
      (*rx_v) = (*v1);
      (*rx_k21) = (*p3);
      (*rx_k) = (*p1)*(*p2)/(*rx_k21); // (*p1) = alpha (*p2) = beta
      (*rx_k12) = (*p1) + (*p2) - (*rx_k21) - (*rx_k);
      break;
    case 5: // alpha=(*p1) beta=(*p2) aob=(*p3)
      (*rx_v)=(*v1);
      (*rx_k21) = ((*p3)*(*p2)+(*p1))/((*p3)+1.0);
      (*rx_k) = ((*p1)*(*p2))/(*rx_k21);
      (*rx_k12) = (*p1) + (*p2) - (*rx_k21) - (*rx_k);
      break;
    case 11: // A2 V, alpha=(*p1), beta=(*p2), k21
#undef beta
#define A (1/(*v1))
#define B (*p3)
#define alpha (*p1)
#define beta (*p2)
      (*ncmt)=2;
      (*rx_v)   = 1/(A+B);
      (*rx_k21) = (A*beta + B*alpha)*(*rx_v);
      (*rx_k)   = alpha*beta/(*rx_k21);
      (*rx_k12) = alpha+beta-(*rx_k21)-(*rx_k);
      break;
    case 10: // A=(*v1), alpha=(*p1), beta=(*p2), B=(*p3)
      // Convert to A (right now A=(*v1) or A=1/(*v1))
#undef A
#define A (*v1)
      (*ncmt)=2;
      (*rx_v)   = 1/(A + B);
      (*rx_k21) = (A*beta + B*alpha)*(*rx_v);
      (*rx_k)   = alpha*beta/(*rx_k21);
      (*rx_k12) = alpha + beta - (*rx_k21) - (*rx_k);
#undef A
#undef B
#undef alpha
#undef beta
#define beta Rf_beta
      break;
    default:
      return NA_REAL;
    }
  } else if ((*p1) > 0.) {
    (*ncmt) = 1;
    switch(*trans){
    case 1: // cl v
      (*rx_k) = (*p1)/(*v1); // k = CL/V
      (*rx_v) = (*v1);
      break;
    case 2: // k V
      (*rx_k) = (*p1);
      (*rx_v) = (*v1);
      break;
    case 11: // alpha V
      (*rx_k) = (*p1);
      (*rx_v) = (*v1);
      break;
    case 10: // alpha A
      (*rx_k) = (*p1);
      (*rx_v) = 1/(*v1);
      break;
    default:
      return 0;
    }
  } else {
    return 0;
  }
  return 1;
}

void linCmtPar1(double *v, double *k, 
                double *vss,
                double *cl,
                double *A,
                double *Af,
                double *alpha,
                double *t12alpha) {
  *vss = *v;
  *cl = (*v)*(*k);
  *A = 1/(*v);
  *alpha = (*k);
  *t12alpha = M_LN2/(*k);
  *Af = (*A)*(*v); // Always 1.
}

void linCmtPar2(double *v, double *k,
                double *k12, double *k21,
                double *vp, double *vss,
                double *cl, double *q,
                double *A, double *B,
                double *Af, double *Bf,
                double *alpha, double *beta,
                double *t12alpha, double *t12beta){
  *vp = (*v)*(*k12)/(*k21);
  *vss = (*v)+(*vp);
  *cl = (*v)*(*k);
  *q  = (*v)*(*k12);
  double a0 = (*k) * (*k21);
  double a1 = -((*k) + (*k12) + (*k21));
  double sq = sqrt(a1*a1-4*a0);
  *alpha = 0.5*(-a1+sq);
  *beta = 0.5*(-a1-sq);
  *A = ((*k21)-(*alpha))/((*beta)-(*alpha))/(*v);
  *B = ((*k21)-(*beta))/((*alpha)-(*beta))/(*v);
  *Af = (*A)*(*v);
  *Bf = (*B)*(*v);
  *t12alpha = M_LN2/(*alpha);
  *t12beta = M_LN2/(*beta);
}

void linCmtPar3(double *v, double *k10,
                double *k12, double *k21, double *k13, double *k31,
                double *vp, double *vp2, double *vss,
                double *cl, double *q, double *q2,
                double *A, double *B, double *C,
                double *Af, double *Bf, double *Cf,
                double *alpha, double *beta, double *gamma,
                double *t12alpha, double *t12beta, double *t12gamma) {
  double a0 = (*k10) * (*k21) * (*k31);
  double a1 = ((*k10) * (*k31)) + ((*k21) * (*k31)) + ((*k21) * (*k13)) + ((*k10) * (*k21)) + ((*k31) * (*k12));
  double a2 = (*k10) + (*k12) + (*k13) + (*k21) + (*k31);
  double p   = a1 - (a2 * a2 / 3.0);
  double qq   = (2.0 * a2 * a2 * a2 / 27.0) - (a1 * a2 / 3.0) + a0;
  double r1  = sqrt(-(p * p * p)/27.0);
  double phi = acos((-qq/2)/r1)/3.0;
  double r2  = 2.0 * exp(log(r1)/3.0);
  *alpha = -(cos(phi) * r2 - a2/3.0);
  *beta = -(cos(phi + 2.0 * M_PI/3.0) * r2 - a2/3.0);
  *gamma = -(cos(phi + 4.0 * M_PI/3.0) * r2 - a2/3.0);
  double a;
  if ((*alpha) < (*beta)) {
    a      = *beta;
    *beta  = *alpha;
    *alpha = a;
  } // now alpha >= beta
  if ((*beta) < (*gamma)) {
    a      = *beta;
    *beta  = *gamma;
    *gamma = a;
  } // now beta >= gamma
  if ((*alpha) < (*beta)) {
    a      = *alpha;
    *alpha = *beta;
    *beta  = a;
  } // now alpha >= beta >= gamma
  *A = ((*k21) - (*alpha)) * ((*k31) - (*alpha)) / ((*alpha) - (*beta)) / ((*alpha) - (*gamma))/(*v);
  *B = ((*k21) - (*beta)) * ((*k31) - (*beta)) / ((*beta) - (*alpha)) / ((*beta) - (*gamma))/(*v);
  *C = ((*k21) - (*gamma)) * ((*k31) - (*gamma)) / ((*gamma) - (*beta)) / ((*gamma) - (*alpha))/(*v);
  *vp  = (*v) * (*k12)/(*k21);
  *vp2 = (*v) * (*k13)/(*k31);
  *vss = (*v) + (*vp) + (*vp2);
  *cl  = (*v) * (*k10);
  *q   = (*v) * (*k12);
  *q2  = (*v) * (*k13);
  *Af  = (*A) * (*v);
  *Bf  = (*B) * (*v);
  *Cf  = (*C) * (*v);
  *t12alpha = M_LN2/(*alpha);
  *t12beta  = M_LN2/(*beta);
  *t12gamma = M_LN2/(*gamma);
}

SEXP toReal(SEXP in){
  int type = TYPEOF(in);
  if (type == REALSXP) return in;
  if (type == INTSXP) {
    SEXP ret = PROTECT(Rf_allocVector(REALSXP, Rf_length(in)));
    int *inI = INTEGER(in);
    double *retR = REAL(ret);
    for (int i = Rf_length(in); i--;){
      retR[i] = (double)(inI[i]);
    }
    UNPROTECT(1);
    return ret;
  }
  _rxode2parse_unprotect();
  Rf_errorcall(R_NilValue, _("not an integer/real"));
  return R_NilValue;
}

SEXP derived1(int trans, SEXP inp, double dig) {
  double zer = 0;
  int lenP = Rf_length(VECTOR_ELT(inp, 0));
  int pro=0;
  SEXP tmp = PROTECT(toReal(VECTOR_ELT(inp, 0))); pro++;
  double *p1 = REAL(tmp);
  int lenV = Rf_length(VECTOR_ELT(inp, 1));
  tmp = PROTECT(toReal(VECTOR_ELT(inp, 1))); pro++;
  double *v1 = REAL(tmp);
  int lenOut = lenP;
  if (lenV != lenP){
    if (lenP == 1){
      lenOut = lenV;
    } else if (lenV != 1){
      _rxode2parse_unprotect();
      Rf_errorcall(R_NilValue, _("The dimensions of the parameters must match"));
    }
  }
  // vc, kel, vss, cl, thalf, alpha, A, fracA
  SEXP ret  = PROTECT(allocVector(VECSXP, 8)); pro++;
  SEXP retN = PROTECT(allocVector(STRSXP, 8)); pro++;

  SET_STRING_ELT(retN,0,mkChar("vc"));
  SEXP vcS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *vc = REAL(vcS);
  SET_VECTOR_ELT(ret, 0, vcS);
  
  SET_STRING_ELT(retN,1,mkChar("kel"));
  SEXP kelS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *kel = REAL(kelS);
  SET_VECTOR_ELT(ret, 1, kelS);
  
  SET_STRING_ELT(retN,2,mkChar("vss"));
  SEXP vssS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *vss = REAL(vssS);
  SET_VECTOR_ELT(ret, 2, vssS);
  
  SET_STRING_ELT(retN,3,mkChar("cl"));
  SEXP clS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *cl = REAL(clS);
  SET_VECTOR_ELT(ret, 3, clS);
  
  SET_STRING_ELT(retN,4,mkChar("t12alpha"));
  SEXP thalfS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *thalf = REAL(thalfS);
  SET_VECTOR_ELT(ret, 4, thalfS);
  
  SET_STRING_ELT(retN,5,mkChar("alpha"));
  SEXP alphaS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *alpha = REAL(alphaS);
  SET_VECTOR_ELT(ret, 5, alphaS);
  
  SET_STRING_ELT(retN,6,mkChar("A"));
  SEXP AS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *A = REAL(AS);
  SET_VECTOR_ELT(ret, 6, AS);
  
  SET_STRING_ELT(retN,7,mkChar("fracA"));
  SEXP fracAS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *fracA = REAL(fracAS);
  SET_VECTOR_ELT(ret, 7, fracAS);

  SEXP sexp_class = PROTECT(allocVector(STRSXP, 1)); pro++;
  SET_STRING_ELT(sexp_class,0,mkChar("data.frame"));
  setAttrib(ret, R_ClassSymbol, sexp_class);

  SEXP sexp_rownames = PROTECT(allocVector(INTSXP,2)); pro++;
  INTEGER(sexp_rownames)[0] = NA_INTEGER;
  INTEGER(sexp_rownames)[1] = -lenOut;
  setAttrib(ret, R_RowNamesSymbol, sexp_rownames);

  setAttrib(ret, R_NamesSymbol, retN);

  unsigned int ncmta=0;

  for (int i = lenOut; i--;){
    parTrans(&trans, ((lenP == 1) ? p1 : p1++),
             ((lenV == 1) ? v1 : v1++), &zer, &zer, &zer, &zer,
             &ncmta, kel, vc, &zer, &zer, &zer, &zer);
    linCmtPar1(vc, kel, vss, cl, A, fracA, alpha, thalf);
    if (dig > 0){
      (*vc) = fprec((*vc), dig);
      (*kel) = fprec((*kel), dig);
      (*vss) = fprec((*vss), dig);
      (*cl) = fprec((*cl), dig);
      (*A) = fprec((*A), dig);
      (*alpha) = fprec((*alpha), dig);
      (*thalf) = fprec((*thalf), dig);
    }
    vc++; kel++; vss++; cl++; A++; fracA++; alpha++; thalf++;

  }
  UNPROTECT(pro);
  return ret;
}

SEXP derived2(int trans, SEXP inp, double dig) {
  double zer = 0;
  int pro=0;

  SEXP tmp = PROTECT(toReal(VECTOR_ELT(inp, 0))); pro++;
  int lenP1 = Rf_length(tmp);
  double *p1 = REAL(tmp);

  tmp = PROTECT(toReal(VECTOR_ELT(inp, 1))); pro++;
  int lenV = Rf_length(tmp);
  double *v1 = REAL(tmp);

  tmp = PROTECT(toReal(VECTOR_ELT(inp, 2))); pro++;
  int lenP2 = Rf_length(tmp);
  double *p2 = REAL(tmp);

  tmp = PROTECT(toReal(VECTOR_ELT(inp, 3))); pro++;
  int lenP3 = Rf_length(tmp);
  double *p3 = REAL(tmp);
  
  int lenOut = max2(lenV, lenP1);
  lenOut = max2(lenOut, lenP2);
  lenOut = max2(lenOut, lenP3);
  lenOut = max2(lenOut, lenV);
  if (lenOut != 1) {
    if ((lenP1 != 1 && lenP1 != lenOut) ||
        (lenP2 != 1 && lenP2 != lenOut) ||
        (lenP3 != 1 && lenP3 != lenOut) ||
        (lenV != 1  && lenV != lenOut)) {
      _rxode2parse_unprotect();
      Rf_errorcall(R_NilValue, _("The dimensions of the parameters must match"));
    }
  }
  // vc, kel, k12, k21, vp, vss, cl, q, thalfAlpha, thalfBeta,
  // alpha, beta, A, B, fracA, fracB
  SEXP ret  = PROTECT(allocVector(VECSXP, 16)); pro++;
  SEXP retN = PROTECT(allocVector(STRSXP, 16)); pro++;

  SET_STRING_ELT(retN,0,mkChar("vc"));
  SEXP vcS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *vc = REAL(vcS);
  SET_VECTOR_ELT(ret, 0, vcS);
  
  SET_STRING_ELT(retN,1,mkChar("kel"));
  SEXP kelS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *kel = REAL(kelS);
  SET_VECTOR_ELT(ret, 1, kelS);

  SET_STRING_ELT(retN,2,mkChar("k12"));
  SEXP k12S = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *k12 = REAL(k12S);
  SET_VECTOR_ELT(ret, 2, k12S);

  SET_STRING_ELT(retN,3,mkChar("k21"));
  SEXP k21S = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *k21 = REAL(k21S);
  SET_VECTOR_ELT(ret, 3, k21S);

  SET_STRING_ELT(retN,4,mkChar("vp"));
  SEXP vpS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *vp = REAL(vpS);
  SET_VECTOR_ELT(ret, 4, vpS);

  SET_STRING_ELT(retN,5,mkChar("vss"));
  SEXP vssS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *vss = REAL(vssS);
  SET_VECTOR_ELT(ret, 5, vssS);

  SET_STRING_ELT(retN,6,mkChar("cl"));
  SEXP clS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *cl = REAL(clS);
  SET_VECTOR_ELT(ret, 6, clS);

  SET_STRING_ELT(retN,7,mkChar("q"));
  SEXP qS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *q = REAL(qS);
  SET_VECTOR_ELT(ret, 7, qS);

  SET_STRING_ELT(retN,8,mkChar("t12alpha"));
  SEXP thalfAlphaS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *thalfAlpha = REAL(thalfAlphaS);
  SET_VECTOR_ELT(ret, 8, thalfAlphaS);

  SET_STRING_ELT(retN,9,mkChar("t12beta"));
  SEXP thalfBetaS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *thalfBeta = REAL(thalfBetaS);
  SET_VECTOR_ELT(ret, 9, thalfBetaS);

  SET_STRING_ELT(retN,10,mkChar("alpha"));
  SEXP alphaS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *alpha = REAL(alphaS);
  SET_VECTOR_ELT(ret, 10, alphaS);

  SET_STRING_ELT(retN,11,mkChar("beta"));
  SEXP betaS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *beta = REAL(betaS);
  SET_VECTOR_ELT(ret, 11, betaS);

  SET_STRING_ELT(retN,12,mkChar("A"));
  SEXP AS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *A = REAL(AS);
  SET_VECTOR_ELT(ret, 12, AS);

  SET_STRING_ELT(retN,13,mkChar("B"));
  SEXP BS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *B = REAL(BS);
  SET_VECTOR_ELT(ret, 13, BS);

  SET_STRING_ELT(retN,14,mkChar("fracA"));
  SEXP fracAS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *fracA = REAL(fracAS);
  SET_VECTOR_ELT(ret, 14, fracAS);

  SET_STRING_ELT(retN,15,mkChar("fracB"));
  SEXP fracBS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *fracB = REAL(fracBS);
  SET_VECTOR_ELT(ret, 15, fracBS);

  SEXP sexp_class = PROTECT(allocVector(STRSXP, 1)); pro++;
  SET_STRING_ELT(sexp_class,0,mkChar("data.frame"));
  setAttrib(ret, R_ClassSymbol, sexp_class);

  SEXP sexp_rownames = PROTECT(allocVector(INTSXP,2)); pro++;
  INTEGER(sexp_rownames)[0] = NA_INTEGER;
  INTEGER(sexp_rownames)[1] = -lenOut;
  setAttrib(ret, R_RowNamesSymbol, sexp_rownames);

  setAttrib(ret, R_NamesSymbol, retN);

  unsigned int ncmta=0;

  for (int i = lenOut; i--;){
    parTrans(&trans, ((lenP1 == 1) ? p1 : p1++), ((lenV == 1) ? v1 : v1++),
             ((lenP2 == 1) ? p2 : p2++), ((lenP3 == 1) ? p3 : p3++), &zer, &zer,
             &ncmta, kel, vc, k12, k21, &zer, &zer);
    linCmtPar2(vc, kel, k12, k21, vp, vss, cl, q, A, B, fracA, fracB,
               alpha, beta, thalfAlpha, thalfBeta);
    if (dig > 0){
      (*vc) = fprec((*vc), dig);
      (*kel) = fprec((*kel), dig);
      (*k12) = fprec((*k12), dig);
      (*k21) = fprec((*k21), dig);
      (*vp) = fprec((*vp), dig);
      (*vss) = fprec((*vss), dig);
      (*cl) = fprec((*cl), dig);
      (*q) = fprec((*q), dig);
      (*A) = fprec((*A), dig);
      (*B) = fprec((*B), dig);
      (*fracA) = fprec((*fracA), dig);
      (*fracB) = fprec((*fracB), dig);
      (*alpha) = fprec((*alpha), dig);
      (*beta) = fprec((*beta), dig);
      (*thalfAlpha) = fprec((*thalfAlpha), dig);
      (*thalfBeta) = fprec((*thalfBeta), dig);
    }
    vc++; kel++; k12++; k21++; vp++; vss++; cl++; q++;
    A++; B++; fracA++; fracB++; alpha++; beta++;
    thalfAlpha++; thalfBeta++;
  }
  UNPROTECT(pro);
  return ret;
}

SEXP derived3(int trans, SEXP inp, double dig) {
  int pro = 0;
  SEXP tmp = PROTECT(toReal(VECTOR_ELT(inp, 0))); pro++;
  int lenP1 = Rf_length(tmp);
  double *p1 = REAL(tmp);

  tmp = PROTECT(toReal(VECTOR_ELT(inp, 1))); pro++;
  int lenV = Rf_length(tmp);
  double *v1 = REAL(tmp);

  tmp = PROTECT(toReal(VECTOR_ELT(inp, 2))); pro++;
  int lenP2 = Rf_length(tmp);
  double *p2 = REAL(tmp);

  tmp = PROTECT(toReal(VECTOR_ELT(inp, 3))); pro++;
  int lenP3 = Rf_length(tmp);
  double *p3 = REAL(tmp);

  tmp = PROTECT(toReal(VECTOR_ELT(inp, 4))); pro++;
  int lenP4 = Rf_length(tmp);
  double *p4 = REAL(tmp);

  tmp = PROTECT(toReal(VECTOR_ELT(inp, 5))); pro++;
  int lenP5 = Rf_length(tmp);
  double *p5 = REAL(tmp);
  
  int lenOut = max2(lenV, lenP1);
  lenOut = max2(lenOut, lenP2);
  lenOut = max2(lenOut, lenP3);
  lenOut = max2(lenOut, lenV);
  lenOut = max2(lenOut, lenP4);
  lenOut = max2(lenOut, lenP5);
  if (lenOut != 1) {
    if ((lenP1 != 1 && lenP1 != lenOut) ||
        (lenP2 != 1 && lenP2 != lenOut) ||
        (lenP3 != 1 && lenP3 != lenOut) ||
        (lenP4 != 1 && lenP4 != lenOut) ||
        (lenP5 != 1 && lenP5 != lenOut) ||
        (lenV != 1  && lenV != lenOut)) {
      _rxode2parse_unprotect();
      Rf_errorcall(R_NilValue, _("The dimensions of the parameters must match"));
    }
  }
  // vc, kel, k12, k21, vp, vss, cl, q, thalfAlpha, thalfBeta,
  // alpha, beta, A, B, fracA, fracB
  SEXP ret  = PROTECT(allocVector(VECSXP, 24)); pro++;
  SEXP retN = PROTECT(allocVector(STRSXP, 24)); pro++;

  SET_STRING_ELT(retN,0,mkChar("vc"));
  SEXP vcS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *vc = REAL(vcS);
  SET_VECTOR_ELT(ret, 0, vcS);
  
  SET_STRING_ELT(retN,1,mkChar("kel"));
  SEXP kelS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *kel = REAL(kelS);
  SET_VECTOR_ELT(ret, 1, kelS);

  SET_STRING_ELT(retN,2,mkChar("k12"));
  SEXP k12S = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *k12 = REAL(k12S);
  SET_VECTOR_ELT(ret, 2, k12S);

  SET_STRING_ELT(retN,3,mkChar("k21"));
  SEXP k21S = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *k21 = REAL(k21S);
  SET_VECTOR_ELT(ret, 3, k21S);

  SET_STRING_ELT(retN,4,mkChar("k13"));
  SEXP k13S = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *k13 = REAL(k13S);
  SET_VECTOR_ELT(ret, 4, k13S);

  SET_STRING_ELT(retN,5,mkChar("k31"));
  SEXP k31S = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *k31 = REAL(k31S);
  SET_VECTOR_ELT(ret, 5, k31S);

  SET_STRING_ELT(retN,6,mkChar("vp"));
  SEXP vpS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *vp = REAL(vpS);
  SET_VECTOR_ELT(ret, 6, vpS);

  SET_STRING_ELT(retN,7,mkChar("vp2"));
  SEXP vp2S = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *vp2 = REAL(vp2S);
  SET_VECTOR_ELT(ret, 7, vp2S);

  SET_STRING_ELT(retN,8,mkChar("vss"));
  SEXP vssS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *vss = REAL(vssS);
  SET_VECTOR_ELT(ret, 8, vssS);

  SET_STRING_ELT(retN,9,mkChar("cl"));
  SEXP clS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *cl = REAL(clS);
  SET_VECTOR_ELT(ret, 9, clS);

  SET_STRING_ELT(retN,10,mkChar("q"));
  SEXP qS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *q = REAL(qS);
  SET_VECTOR_ELT(ret, 10, qS);

  SET_STRING_ELT(retN,11,mkChar("q2"));
  SEXP q2S = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *q2 = REAL(q2S);
  SET_VECTOR_ELT(ret, 11, q2S);

  SET_STRING_ELT(retN,12,mkChar("t12alpha"));
  SEXP thalfAlphaS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *thalfAlpha = REAL(thalfAlphaS);
  SET_VECTOR_ELT(ret, 12, thalfAlphaS);

  SET_STRING_ELT(retN,13,mkChar("t12beta"));
  SEXP thalfBetaS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *thalfBeta = REAL(thalfBetaS);
  SET_VECTOR_ELT(ret, 13, thalfBetaS);

  SET_STRING_ELT(retN,14,mkChar("t12gamma"));
  SEXP thalfGammaS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *thalfGamma = REAL(thalfGammaS);
  SET_VECTOR_ELT(ret, 14, thalfGammaS);

  SET_STRING_ELT(retN,15,mkChar("alpha"));
  SEXP alphaS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *alpha = REAL(alphaS);
  SET_VECTOR_ELT(ret, 15, alphaS);

  SET_STRING_ELT(retN,16,mkChar("beta"));
  SEXP betaS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *beta = REAL(betaS);
  SET_VECTOR_ELT(ret, 16, betaS);

  SET_STRING_ELT(retN,17,mkChar("gamma"));
  SEXP gammaS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *gamma = REAL(gammaS);
  SET_VECTOR_ELT(ret, 17, gammaS);

  SET_STRING_ELT(retN,18,mkChar("A"));
  SEXP AS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *A = REAL(AS);
  SET_VECTOR_ELT(ret, 18, AS);


  SET_STRING_ELT(retN,19,mkChar("B"));
  SEXP BS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *B = REAL(BS);
  SET_VECTOR_ELT(ret, 19, BS);

  SET_STRING_ELT(retN,20,mkChar("C"));
  SEXP CS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *C = REAL(CS);
  SET_VECTOR_ELT(ret, 20, CS);

  SET_STRING_ELT(retN,21,mkChar("fracA"));
  SEXP fracAS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *fracA = REAL(fracAS);
  SET_VECTOR_ELT(ret, 21, fracAS);

  SET_STRING_ELT(retN,22,mkChar("fracB"));
  SEXP fracBS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *fracB = REAL(fracBS);
  SET_VECTOR_ELT(ret, 22, fracBS);

  SET_STRING_ELT(retN,23,mkChar("fracC"));
  SEXP fracCS = PROTECT(allocVector(REALSXP, lenOut)); pro++;
  double *fracC = REAL(fracCS);
  SET_VECTOR_ELT(ret, 23, fracCS);

  SEXP sexp_class = PROTECT(allocVector(STRSXP, 1)); pro++;
  SET_STRING_ELT(sexp_class,0,mkChar("data.frame"));
  setAttrib(ret, R_ClassSymbol, sexp_class);

  SEXP sexp_rownames = PROTECT(allocVector(INTSXP,2)); pro++;
  INTEGER(sexp_rownames)[0] = NA_INTEGER;
  INTEGER(sexp_rownames)[1] = -lenOut;
  setAttrib(ret, R_RowNamesSymbol, sexp_rownames);

  setAttrib(ret, R_NamesSymbol, retN);

  unsigned int ncmta=0;

  for (int i = lenOut; i--;){
    parTrans(&trans, ((lenP1 == 1) ? p1 : p1++), ((lenV == 1) ? v1 : v1++),
             ((lenP2 == 1) ? p2 : p2++), ((lenP3 == 1) ? p3 : p3++),
             ((lenP4 == 1) ? p4 : p4++), ((lenP5 == 1) ? p5 : p5++),
             &ncmta, kel, vc, k12, k21, k13, k31);
    linCmtPar3(vc, kel, k12, k21, k13, k31,
               vp, vp2, vss, cl, q, q2,
               A, B, C, fracA, fracB, fracC, alpha, beta, gamma,
               thalfAlpha, thalfBeta, thalfGamma);
    if (dig > 0) {
      (*vc)  = fprec((*vc), dig);
      (*kel) = fprec((*kel), dig);
      (*k12) = fprec((*k12), dig);
      (*k21) = fprec((*k21), dig);
      (*k13) = fprec((*k13), dig);
      (*k31) = fprec((*k31), dig);
      (*vp)  = fprec((*vp), dig);
      (*vss) = fprec((*vss), dig);
      (*vp2) = fprec((*vp2), dig);
      (*cl)  = fprec((*cl), dig);
      (*q)   = fprec((*q), dig);
      (*q2)  = fprec((*q2), dig);
      (*A)   = fprec((*A), dig);
      (*B)   = fprec((*B), dig);
      (*C)   = fprec((*C), dig);
      (*fracA)=fprec((*fracA), dig);
      (*fracB)=fprec((*fracB), dig);
      (*fracC)=fprec((*fracC), dig);
      (*alpha)=fprec((*alpha), dig);
      (*beta) =fprec((*beta), dig);
      (*gamma)=fprec((*gamma), dig);
      (*thalfAlpha)=fprec((*thalfAlpha), dig);
      (*thalfBeta)=fprec((*thalfBeta), dig);
      (*thalfGamma)=fprec((*thalfGamma), dig);
    }
    vc++; kel++; k12++; k21++; k13++; k31++;
    vp++; vp2++; vss++; cl++; q++; q2++;
    A++; B++; C++; fracA++; fracB++; fracC++; alpha++; beta++; gamma++;
    thalfAlpha++; thalfBeta++; thalfGamma++;
  }
  UNPROTECT(pro);
  return ret;
}

SEXP _calcDerived(SEXP ncmtSXP, SEXP transSXP, SEXP inp, SEXP sigdigSXP) {
  int tInp = TYPEOF(inp);
  int trans=-1;
  if (TYPEOF(transSXP) == REALSXP){
    trans = (int)(REAL(transSXP)[0]);
  }
  int ncmt=-1;
  if (TYPEOF(ncmtSXP) == REALSXP) {
    ncmt = (int)(REAL(ncmtSXP)[0]);
  }
  double dig=0.0;
  int tDig = TYPEOF(sigdigSXP);
  if (tDig == INTSXP) {
    dig = (double)(INTEGER(sigdigSXP)[0]);
  } else if (tDig == REALSXP) {
    dig = REAL(sigdigSXP)[0];
  }
  if (tInp == VECSXP){
    switch (ncmt){
    case 1:
      return derived1(trans, inp, dig);
      break;
    case 2:
      return derived2(trans, inp, dig);
      break;
    case 3:
      return derived3(trans, inp, dig);
      break;
    default:
      _rxode2parse_unprotect();
      Rf_errorcall(R_NilValue, _("'ncmt' needs to be 1-3"));
    }
  } else {
    _rxode2parse_unprotect();
    Rf_errorcall(R_NilValue, _("'inp' needs to be list/data frame"));
  }
  return R_NilValue;
}

double linCmtA(rx_solve *rx, unsigned int id, double _t, int linCmt,
               int i_cmt, int trans,
               double p1, double v1,
               double p2, double p3,
               double p4, double p5,
               double d_tlag, double d_F, double d_rate1, double d_dur1,
               // Oral parameters
               double d_ka, double d_tlag2, double d_F2,  double d_rate2, double d_dur2) {
  rx_solving_options_ind *ind = &(rx->subjects[id]);
  double t = _t - ind->curShift;
  int evid = 0;
  int idx = ind->idx;
  double Alast0[4] = {0, 0, 0, 0};
  rx_solving_options *op = rx->op;
  int oral0;
  oral0 = (d_ka > 0) ? 1 : 0;
  double *A;
  double *Alast;
  /* A = Alast0; Alast=Alast0; */
  double tlast;
  unsigned int ncmt = 1;
  double rx_k=0, rx_v=0;
  double rx_k12=0;
  double rx_k21=0;
  double rx_k13=0;
  double rx_k31=0;
  double b1=0, b2=0, r1 = 0, r2 = 0;
  double curTime = getTime(ind->ix[idx], ind);
  int sameTime = isSameTime(t, curTime);
  if (!parTrans(&trans, &p1, &v1, &p2, &p3, &p4, &p5,
                &ncmt, &rx_k, &rx_v, &rx_k12,
                &rx_k21, &rx_k13, &rx_k31)){
    return NA_REAL;
  }
  return 0.0;
}
