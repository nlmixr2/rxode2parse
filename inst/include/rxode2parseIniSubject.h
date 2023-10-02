#if defined(__cplusplus)
extern "C" {
#endif

#ifdef _isrxode2parse_
#define setIndPointersByThread rxode2parse_setIndPointersByThread
  void setIndPointersByThread(rx_solving_options_ind *ind);
#define sortInd rxode2parse_sortInd0
  void sortInd(rx_solving_options_ind *ind);
#endif

  static inline int iniSubject(int solveid, int inLhs, rx_solving_options_ind *ind,
                               rx_solving_options *op, rx_solve *rx,
                               t_update_inis u_inis) {
    ind->_rxFlag=1;
    setIndPointersByThread(ind);
    for (int i=rxLlikSaveSize*op->nLlik; i--;) {
      ind->llikSave[i] = 0.0;
    }
    ind->ixds = ind->idx = ind->_update_par_ptr_in = 0; // reset dosing
    ind->id=solveid;
    ind->cacheME=0;
    ind->curShift=0.0;
    ind->doSS = 0;
    ind->lastIsSs2 = false;
    // neq[0] = op->neq
    int ncmt = (op->neq + op->extraCmt);
    for (int j = ncmt; j--;) {
      ind->InfusionRate[j] = 0;
      ind->on[j] = 1;
      ind->tlastS[j] = NA_REAL;
      ind->tfirstS[j] = NA_REAL;
      ind->curDoseS[j] = NA_REAL;
    }
    ind->inLhs = inLhs;
    if (rx->nMtime) calc_mtime(solveid, ind->mtime);
    for (int j = op->nlhs; j--;) ind->lhs[j] = NA_REAL;
    if ((inLhs == 0 && op->neq > 0) ||
        (inLhs == 1 && op->neq == 0 && (rx->nIndSim > 0 || (rx->simflg & 1) != 0 ))) {
      ind->isIni = 1;
      // Also can update individual random variables (if needed)
      if (inLhs == 0) memcpy(ind->solve, op->inits, op->neq*sizeof(double));
      u_inis(solveid, ind->solve); // Update initial conditions @ current time
      ind->isIni = 0;
    }
    ind->_newind = 1;
    ind->dosenum = 0;
    ind->tlast = NA_REAL;
    ind->tfirst = NA_REAL;
    ind->curDose = NA_REAL;
    if (inLhs == 0 || (inLhs == 1 && op->neq==0)) {
      ind->solved = -1;
    }
    sortInd(ind);
    if (op->badSolve) return 0;
    ind->ixds=ind->idx=0;
    if (ncmt) ind->pendingDosesN[0] = 0;
    return 1;
  }

#if defined(__cplusplus)
}
#endif
