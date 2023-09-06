#define USE_FC_LEN_T
#define STRICT_R_HEADERS

#include <Rcpp.h>
#include <algorithm>
#include "../inst/include/rxode2parse.h"
#include "../inst/include/timsort.h"
#include "needSortDefines.h"


#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("rxode2parse", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

#ifndef max2
#define max2( a , b )  ( (a) > (b) ? (a) : (b) )
#endif

using namespace Rcpp;

extern "C" {
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

  extern rx_solving_options _rxode2parse_op_global;
  extern rx_solve _rxode2parse_rx_global;
  extern t_handle_evidL _rxode2parse_handle_evidL;
  extern t_getDur _rxode2parse_getDur;
  extern t_getTime _rxode2parse_getTime;
#define _getDur _rxode2parse_getDur
}

#include "../inst/include/rxode2parseHandleEvid.h"
#include "../inst/include/rxode2parseGetTime.h"
#include "../inst/include/rxode2parseHandleSs.h"
#include "../inst/include/rxode2parseSortInd.h"

extern "C" void rxode2parse_sortInd0(rx_solving_options_ind *ind) {
  rxode2parse_sortInd(ind);
}
