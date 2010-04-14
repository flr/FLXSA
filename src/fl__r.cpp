#include <math.h>
#include "fl__r.hpp" 
//#include <float.h>

SEXP ReturnError(char *ErrStr) 
    {
    SEXP Err;
    
    PROTECT(Err = allocVector(STRSXP, 1));
    SET_STRING_ELT(Err, 0, mkChar(ErrStr));
    
    UNPROTECT(1);

    return Err;
    }

int NElemList(SEXP x)
   {
   //Check that it is a list
   if (!IS_LIST(x) || TYPEOF(x) != VECSXP) 
      return 0;
   else
      return length(x);
  }

bool InputRange(SEXP _range, int *pMinAge, int *pMaxAge, int *pPlusgroup, int *pMinYear, int *pMaxYear)
    {
    int n;
    double *x;


    if (!isVector(_range) || !isNumeric(_range)) 
         return false;

    PROTECT(_range = AS_NUMERIC(_range));
    n = LENGTH(_range);

    x = NUMERIC_POINTER(_range); 
    
    switch (n)
      {
      case 5:
         *pMaxYear   = (int)x[4];
      case 4:
         *pMinYear   = (int)x[3];
      case 3:
         *pPlusgroup = (int)x[2];
      case 2:
         *pMaxAge    = (int)x[1];
      case 1:
         *pMinAge    = (int)x[0];
      }

    UNPROTECT(1);

    return true;
    }

bool InputFLQuant(SEXP FLQuant, double ***pD, int MinAge, int MaxAge, int MinYear, int MaxYear)
    {
    SEXP Quant    = PROTECT(duplicate(GET_SLOT(FLQuant, install(".Data")))),
         dims     = GET_DIM(Quant),
         dimnames = GET_DIMNAMES(Quant);

    double *Q     = NUMERIC_POINTER(Quant);

  	 int dim[6], n = length(dims);
    
    dim[0] = INTEGER(dims)[0];
    dim[1] = INTEGER(dims)[1];
    dim[2] = INTEGER(dims)[2];
    dim[3] = INTEGER(dims)[3];
    dim[4] = INTEGER(dims)[4];
    dim[5] = INTEGER(dims)[5];

    if (((int)dim[0]) <  1 || ((int)dim[1]) <  1 || 
        ((int)dim[2]) != 1 || ((int)dim[3]) != 1 || ((int)dim[4]) != 1)
      {
      UNPROTECT(1);

      return FALSE;
      }

    int   _MinAge  = 1,
          _MinYear = 1,
          _MaxAge  = (short)dim[0],
          _MaxYear = (short)dim[1];
    	   
    if (dimnames != R_NilValue) 
      if (TYPEOF(dimnames) == VECSXP) 
         {
         int t;
         
         if (n >= 1 && INTEGER(dims)[0] >= 1) 
               {
	            t = atoi(CHAR(STRING_ELT(VECTOR_ELT(dimnames, 0), 0))) - 1;
               _MinAge += t;
               _MaxAge += t;
  	            }
		   if (n >= 2 && INTEGER(dims)[1] >= 1) 
               {
	            t = atoi(CHAR(STRING_ELT(VECTOR_ELT(dimnames, 1), 0))) - 1;
               _MinYear += t;
               _MaxYear += t;
    	         }
		   }

    if (MinAge  < _MinAge  || MaxAge  > _MaxAge || 
	     MinYear < _MinYear || MaxYear > _MaxYear )
      {
      UNPROTECT(1);
   
      return FALSE;
      }

    int iAge, iYear,
        i,    j;

    for (iAge = MinAge, i = 0; iAge <= MaxAge; iAge++, i++)
       for (iYear = MinYear, j = 0; iYear <= MaxYear; iYear++, j++)
          (*pD)[iAge][iYear] = (Q)[i + j*(_MaxAge-_MinAge+1)];       
      
    UNPROTECT(1);
  
    return TRUE;
    }

bool InputFLQuant(SEXP FLQuant, double ***pD, short *pMin1, short *pMax1, short *pMin2, short *pMax2)
    {
    SEXP v        = PROTECT(duplicate(GET_SLOT(FLQuant, install(".Data")))),
         dims     = GET_DIM(v),
         dimnames = GET_DIMNAMES(v);

    double *a     = NUMERIC_POINTER(v);

  	 short dim[2], n = length(dims);

    if (n < 2)
       {
       UNPROTECT(1);
  
       return FALSE;
       }

    dim[0] = INTEGER(dims)[0];
    dim[1] = INTEGER(dims)[1];
    
    *pMin1  = 1;
    *pMin2  = 1;
    *pMax1  = (short)dim[0];
    *pMax2  = (short)dim[1];
    	   
    if (dimnames != R_NilValue) 
      if (TYPEOF(dimnames) == VECSXP) 
         {
         short t;

         if (n >= 1 && INTEGER(dims)[0] >= 1) 
               {
	            t = atoi(CHAR(STRING_ELT(VECTOR_ELT(dimnames, 0), 0))) - 1;
               *pMin1 += t;
               *pMax1 += t;
  	            }
		   if (n >= 2 && INTEGER(dims)[1] >= 1) 
               {
	            t = atoi(CHAR(STRING_ELT(VECTOR_ELT(dimnames, 1), 0))) - 1;
               *pMin2 += t;
               *pMax2 += t;
    	         }
		   }
		 
    flallocArray(pD, *pMin1, *pMax1, *pMin2, *pMax2);

    short i1, i2,
        i,    j;

    for (i1 = *pMin1, i = 0; i1 <= *pMax1; i1++, i++)
       for (i2 = *pMin2, j = 0; i2 <= *pMax2; i2++, j++)
          (*pD)[i1][i2] = (a)[i + j*(*pMax1-*pMin1+1)];       
      
    UNPROTECT(1);
  
    return TRUE;
    }

SEXP CreateFLQuant(double ***pD, int MinAge, int MaxAge, int MinYear, int MaxYear)
    {
    int i, j, iAge, iYear;
    
    SEXP FLQuant, v, 
         d1, d2, d3, d4, d5, d6, 
         dim, dimnames, names;    

    //Create new S4 object    
    //PROTECT(v = NEW_OBJECT(MAKE_CLASS("array")));
    PROTECT(FLQuant = NEW_OBJECT(MAKE_CLASS("FLQuant")));

    //Create array for slot    
    //Set dimensions of array
    PROTECT(dim     = allocVector(INTSXP, 6));       
    INTEGER(dim)[0] = MaxAge -MinAge +1;
    INTEGER(dim)[1] = MaxYear-MinYear+1;
    INTEGER(dim)[2] = 1; 
    INTEGER(dim)[3] = 1; 
    INTEGER(dim)[4] = 1; 
    INTEGER(dim)[5] = 1; 
        
    //Allocate memory
    PROTECT(v = Rf_allocArray(REALSXP, dim)); 
    
    //Create dimension names
    PROTECT(dimnames = allocVector(VECSXP, 6));
    
    PROTECT(d1 = allocVector(INTSXP, MaxAge-MinAge +1));
    for (iAge=MinAge, i=0; iAge<=MaxAge; iAge++, i++)
        INTEGER(d1)[i] = iAge; 
    SET_VECTOR_ELT(dimnames, 0, d1);
    
    PROTECT(d2 = allocVector(INTSXP, MaxYear-MinYear+1));
    for (iYear=MinYear, i=0; iYear<=MaxYear; iYear++, i++)
        INTEGER(d2)[i] = iYear; 
    SET_VECTOR_ELT(dimnames, 1, d2);
     
    PROTECT(d3 = allocVector(STRSXP, 1));
    SET_STRING_ELT(d3, 0, mkChar("unique"));
    SET_VECTOR_ELT(dimnames, 2, d3);
    
    PROTECT(d4 = allocVector(STRSXP, 1));
    SET_STRING_ELT(d4, 0, mkChar("all"));
    SET_VECTOR_ELT(dimnames, 3, d4);
    
    PROTECT(d5 = allocVector(STRSXP, 1));
    SET_STRING_ELT(d5, 0, mkChar("unique"));
    SET_VECTOR_ELT(dimnames, 4, d5);
    
    PROTECT(d6 = allocVector(STRSXP, 1));
    SET_STRING_ELT(d6, 0, mkChar("1"));
    SET_VECTOR_ELT(dimnames, 5, d6);
    
    //Create names for dimensions
    PROTECT(names = allocVector(STRSXP, 6));
    SET_STRING_ELT(names, 0, mkChar("age"));
    SET_STRING_ELT(names, 1, mkChar("year"));
    SET_STRING_ELT(names, 2, mkChar("unit"));
    SET_STRING_ELT(names, 3, mkChar("season"));
    SET_STRING_ELT(names, 4, mkChar("area"));
    SET_STRING_ELT(names, 5, mkChar("iter"));
    setAttrib(dimnames, R_NamesSymbol, names);
    setAttrib(v, R_DimNamesSymbol, dimnames);
    //setAttrib(v, R_ClassSymbol, mkChar("FLQuant"));

  //SEXP klass;
  //PROTECT(klass = allocVector(STRSXP, 1));
  //SET_STRING_ELT(klass, 0, mkChar("FLQuant"));
  //setAttrib(v, R_ClassSymbol, klass);

    
    //Set data
    for (iAge=MinAge, i=0; iAge<=MaxAge; iAge++, i++)
       for (iYear=MinYear, j=0; iYear<=MaxYear; iYear++, j++)
          {
          //char msg[1024];
          //sprintf(msg, "%d, %d, %d, %e", i + j*(MaxAge-MinAge+1), iAge, iYear, (*pD)[iAge][iYear]);
          //R_ShowMessage(msg);
          REAL(v)[i + j*(MaxAge-MinAge+1)] =(*pD)[iAge][iYear];
    	    }
           
    //Set slot
    //FLQuant = SET_SLOT(FLQuant, install("v"), v);
    FLQuant = R_do_slot_assign(FLQuant, install(".Data"), v);

    UNPROTECT(11);
    
    return FLQuant;
    }

SEXP CreateArray(double ***pD, int MinAge, int MaxAge, int MinYear, int MaxYear)
    {
    int i, j, iAge, iYear;
    
    SEXP FLQuant, v, 
         d1, d2,  
         dim,   dimnames, names;    

    //Create new S4 object    
    //PROTECT(v = NEW_OBJECT(MAKE_CLASS("array")));
    PROTECT(FLQuant = NEW_OBJECT(MAKE_CLASS("FLQuant")));

    //Create array for slot    
    //Set dimensions of array
    PROTECT(dim     = allocVector(INTSXP, 2));       
    INTEGER(dim)[0] = MaxAge -MinAge +1;
    INTEGER(dim)[1] = MaxYear-MinYear+1; 
        
    //Allocate memory
    PROTECT(v = Rf_allocArray(REALSXP, dim)); 
    
    //Create dimension names
    PROTECT(dimnames = allocVector(VECSXP, 2));
    
    PROTECT(d1 = allocVector(INTSXP, MaxAge-MinAge +1));
    for (iAge=MinAge, i=0; iAge<=MaxAge; iAge++, i++)
        INTEGER(d1)[i] = iAge; 
    SET_VECTOR_ELT(dimnames, 0, d1);
    
    PROTECT(d2 = allocVector(INTSXP, MaxYear-MinYear+1));
    for (iYear=MinYear, i=0; iYear<=MaxYear; iYear++, i++)
        INTEGER(d2)[i] = iYear; 
    SET_VECTOR_ELT(dimnames, 1, d2);
     
    
    //Create names for dimensions
    PROTECT(names = allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("age"));
    setAttrib(dimnames, R_NamesSymbol, names);
    setAttrib(v, R_DimNamesSymbol, dimnames);
    
    //Set data
    for (iAge=MinAge, i=0; iAge<=MaxAge; iAge++, i++)
       for (iYear=MinYear, j=0; iYear<=MaxYear; iYear++, j++)
          {
          REAL(v)[i + j*(MaxAge-MinAge+1)] =(*pD)[iAge][iYear];
    	    }
           
    UNPROTECT(7);
    
    return v;
    }
