static inline void rxode2parse_sortRest(rx_solving_options_ind *ind, int i0) {
  // Here the ix has been calculated at least once. Sort it if it changed
  rx_solve *rx = &rx_global;
  rx_solving_options *op = &op_global;
  int evid0 = getEvid(ind, i0);
  if (evid0 == 0 || evid0 == 2) return;
  // Reset times for infusion
  double time0 = getTime__(ind->ix[i0], ind, 1);
  int doSort = 0;
  double *time = ind->timeThread;
  int idx0 = ind->idx;
  double curTime = time0;
  for (int i = i0; i < ind->n_all_times; i++) {
    // no need to assign     ind->ix[i] = i;
    ind->idx = ind->ix[i];
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
    if (curTime > time[i]) {
      // 0.4 0.1
      if (time0 > curTime) {
        // can't go and change realized times
        // error here
      } else {
        doSort = 1;
      }
    }
    curTime = time[i];
  }
  if (doSort) {
    gfx::timsort(ind->ix + i0, ind->ix + ind->n_all_times,
                 [ind, time](int a, int b){
                   double timea = time[a],
                     timeb = time[b];
                   if (timea == timeb) {
                     int evida = -getEvid(ind, a),
                       evidb = -getEvid(ind, b);
                     if (evida == evidb){
                       return a < b;
                     }
                     return evida < evidb;
                   }
                   return timea < timeb;
                 });
    ind->ixds = 0;
    for (int i = 0; i < i0; i++) {
      int evid = getEvid(ind, i);
      if (!isObs(evid)) {
        ind->ixds++;
      }
    }
  }
  ind->idx = idx0;
}

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
