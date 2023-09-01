static inline int solComp2C(double *k10, double *k12, double *k21,
                            double *L, double *C1, double *C2) {
  // blas and R uses row major matrices
  double sum = (*k10) + (*k12) + (*k21);
  double disc= sqrt(sum*sum - 4.0* (*k10)*(*k21));
  double div[2];
  double tmp;
  L[0] = 0.5*(sum + disc);
  L[1] = 0.5*(sum - disc);
  div[0] = L[1] - L[0];
  div[1] = L[0] - L[1];
  if (div[0]*div[1] == 0) return 0;
  C1[0] = *k21 - L[0];
  C1[2] = *k21 - L[1];
  C2[0] = C2[2] = *k21;
  C1[1] = C1[3] = *k12;
  tmp = *k10 + *k12;
  C2[1] = tmp - L[0];
  C2[3] = tmp - L[1];
  C1[0] = C1[0]/div[0];
  C1[1] = C1[1]/div[0];
  C2[0] = C2[0]/div[0];
  C2[1] = C2[1]/div[0];
  C1[2] = C1[2]/div[1];
  C1[3] = C1[3]/div[1];
  C2[2] = C2[2]/div[1];
  C2[3] = C2[3]/div[1];
  return 1;
}
