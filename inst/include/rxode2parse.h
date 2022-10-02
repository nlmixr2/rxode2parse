#pragma once
#define STRICT_R_HEADERS
#ifndef __rxode2parse_H__
#define __rxode2parse_H__
#define rxLlikSaveSize 9


#include <R.h>
#include <stdbool.h>

#include <float.h>
#include <stdio.h>
#include <stdarg.h>

#include "rxode2parse_control.h"
#include <stdint.h>    // for uint64_t rather than unsigned long long

#ifdef _isrxode2parse_
#define _linCmtParse _rxode2parse_linCmtParse
#define _rxode2_linCmtGen _rxode2parse_linCmtGen
#define rc_buf_read _rxode2parse_rc_buf_read
#define sIniTo _rxode2parse_sIniTo
#define sFree _rxode2parse_sFree
#define sFreeIni _rxode2parse_sFreeIni
#define sAppendN _rxode2parse_sAppendN
#define sAppend _rxode2parse_sAppend
#define sPrint _rxode2parse_sPrint
#define lineIni _rxode2parse_lineIni
#define lineFree _rxode2parse_lineFree
#define addLine _rxode2parse_addLine
#define curLineProp _rxode2parse_curLineProp
#define curLineType _rxode2parse_curLineType
#define doDot _rxode2parse_doDot
#define doDot2 _rxode2parse_doDot2
#define _setSilentErr _rxode2parse__setSilentErr
#define _isRstudio2 _rxode2parse_isRstudio2
#define setSilentErr _rxode2parse_setSilentErr
#define setRstudioPrint _rxode2parse_setRstudioPrint
#define getSilentErr _rxode2parse_getSilentErr
#define getRstudioPrint _rxode2parse_getRstudioPrint
#define RSprintf _rxode2parse_RSprintf
#endif


typedef struct sbuf {
  char *s;        /* curr print buffer */
  int sN;
  int o;                        /* offset of print buffer */
} sbuf;
  
typedef struct vLines {
  char *s;
  int sN;
  int o;
  int n;
  int nL;
  char **line;
  int *lProp;
  int *lType;
  int *os;
} vLines;

static inline void sNull(sbuf *sbb) {
  sbb->s = NULL;
  sbb->sN=0;
  sbb->o=0;
}

static inline void lineNull(vLines *sbb) {
  sbb->s = NULL;
  sbb->lProp = NULL;
  sbb->lType = NULL;
  sbb->line = NULL;
  sbb->os = NULL;
  sbb->sN = 0;
  sbb->nL = 0;
  sbb->n  = 0;
  sbb->o  = 0;
}



#endif
