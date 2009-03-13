#ifndef _INC_fl__r
#define _INC_fl__r

#include <stdio.h>
#include <stdlib.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "fl__types.hpp"
#include "flalloc.hpp"

int  NElemList(SEXP);
bool InputRange(SEXP, int *, int *, int *, int *, int *);

bool InputFLQuant(SEXP, double ***, int,   int,   int,   int);
bool InputFLQuant(SEXP, double ***, short *, short *, short *, short *);

SEXP CreateArray(  double ***, int, int, int, int);
SEXP CreateFLQuant(double ***, int, int, int, int);

bool GetFLQuantDims(SEXP, int (*)[5], int (*)[5]);

SEXP ReturnError(char *);

#endif /* _INC_fl__r  */
