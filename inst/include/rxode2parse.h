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
