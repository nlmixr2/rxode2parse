#define USE_FC_LEN_T
#define STRICT_R_HEADERS

#include <Rcpp.h>
#include <algorithm>
#include "../inst/include/timsort.h"
#include "needSortDefines.h"
#include "../inst/include/rxode2parseHandleSs.h"
#include "../inst/include/rxode2parseSortInd.h"

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


extern "C" void rxode2parse_sortInd0(rx_solving_options_ind *ind) {
  rxode2parse_sortInd(ind);
}
