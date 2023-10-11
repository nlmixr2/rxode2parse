static inline void rxode2parse_sortInd(rx_solving_options_ind *ind) {
// #ifdef _OPENMP
//   int core = omp_get_thread_num();
// #else
//   int core = 0;
// #endif
  rx_solve *rx = &rx_global;
  rx_solving_options *op = &op_global;
  // Reset times for infusion
  int doSort = 1;
  double *time = ind->timeThread;
  ind->ixds = 0;
  ind->curShift = 0;
  for (int i = 0; i < ind->n_all_times; i++) {
    ind->ix[i] = i;
    ind->idx = i;
    int evid = getEvid(ind, i);
    if (!isObs(evid)) {
      time[i] = getTime__(ind->ix[i], ind, 1);
      ind->ixds++;
    } else {
      if (evid == 3) {
        ind->curShift -= rx->maxShift;
      }
      time[i] = getTime__(ind->ix[i], ind, 1);
    }
    /* REprintf("i: %d; evid: %d; %f -> %f\n", i, evid, getAllTimes(ind, i), time[i]); */
    if (op->naTime == 1){
      doSort=0;
      break;
    }
  }
  if (doSort) {
    gfx::timsort(ind->ix, ind->ix + ind->n_all_times,
         [ind, time](int a, int b){
           double timea = time[a],
             timeb = time[b];
           if (timea == timeb) {
             return a < b;
           }
           return timea < timeb;
         });
  }
}
