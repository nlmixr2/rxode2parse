#define USE_FC_LEN_T
#define STRICT_R_HEADERS
#include <sys/stat.h> 
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>   /* dj: import intptr_t */
#include <errno.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("rxode2parse", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif
#include "../inst/include/rxode2parse.h"

int _setSilentErr=0, _isRstudio2=0;
extern void setSilentErr(int silent){
  _setSilentErr = silent;
}

extern void setRstudioPrint(int rstudio){
  _isRstudio2=rstudio;
}


extern int getSilentErr(){return _setSilentErr;}

extern int getRstudioPrint(){return _isRstudio2;}

extern void RSprintf(const char *format, ...) {
  if (_setSilentErr == 0) {
    if(_isRstudio2){
      va_list args;
      va_start(args, format);
      REvprintf(format, args);
      va_end(args);
    } else{
      va_list args;
      va_start(args, format);
      Rvprintf(format, args);
      va_end(args);
    } 
  }
}

