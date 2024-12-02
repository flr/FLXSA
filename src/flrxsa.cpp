//
// flxrsa.cpp
//
// Last Change: Wed Jul 10, 2013 at 02:21 PM +0200
// $Id: flrxsa.cpp,v 1.11 2007/10/22 09:25:31 ltkell Exp $

#include "flrxsa.h"

ExtendedSurvivorsAnalysisR::ExtendedSurvivorsAnalysisR(SEXP xControl, SEXP xM, CatchData *pCatchData, TuningFleets *pTuningFleets)
   {                   
   int iFleet;
   
   pCatch = pCatchData;
   pTune  = pTuningFleets;
   
   Error = !FiFiErrNull;

   //Controls
   Control.MSSTol        =     max(   REAL(GET_SLOT(xControl,Rf_install("tol"           )))[0],0.0);
   Control.MaxIters      = (short) INTEGER(GET_SLOT(xControl,Rf_install("maxit"         )))[0];
   Control.MinNSE        =     max(   REAL(GET_SLOT(xControl,Rf_install("min.nse"       )))[0],0.0);
   Control.FShkSE        =     max(   REAL(GET_SLOT(xControl,Rf_install("fse"           )))[0],0.0); 
   Control.LRcrtAge      = (short) INTEGER(GET_SLOT(xControl,Rf_install("rage"          )))[0];
   Control.ConQAge       = (short) INTEGER(GET_SLOT(xControl,Rf_install("qage"          )))[0];
   Control.Shk2N         =         LOGICAL(GET_SLOT(xControl,Rf_install("shk.n"         )))[0];
   Control.Shk2F         =         LOGICAL(GET_SLOT(xControl,Rf_install("shk.f"         )))[0];
   Control.Shk2FYr       = (short) INTEGER(GET_SLOT(xControl,Rf_install("shk.yrs"       )))[0];
   Control.Shk2FAge      = (short) INTEGER(GET_SLOT(xControl,Rf_install("shk.ages"      )))[0];
   Control.TSRange       = (short) INTEGER(GET_SLOT(xControl,Rf_install("tsrange"       )))[0];
   Control.TSPower       = (short) INTEGER(GET_SLOT(xControl,Rf_install("tspower"       )))[0];
   Control.PlusGroup     = true;
   Control.TuningWindow  = (short) INTEGER(GET_SLOT(xControl,Rf_install("window"        )))[0];
   Control.VPA           =         LOGICAL(GET_SLOT(xControl,Rf_install("vpa"           )))[0];
   
   MinAge  = pCatch->GetMinAge();             
   
   for (iFleet=1, MaxAge=MinAge; iFleet <= pTune->GetNFleet(); iFleet++)
      MaxAge = pTune->GetMaxAge(iFleet);

   if (Control.PlusGroup)
      Control.PlusGroup = pCatch->GetMaxAge();
   else
      Control.PlusGroup = pCatch->GetMaxAge()+1;
  
   MaxAge = min(pCatch->GetMaxAge(),Control.PlusGroup-1);
   Control.ConQAge = max(min(Control.ConQAge, MaxAge), MinAge);
   
   MinYear = pCatch->GetMinYear();
   MaxYear = pCatch->GetMaxYear();
      
   AllocArrays();

   if (!InputFLQuant(xM, &M, MinAge, MaxAge + 1, MinYear, MaxYear))
      return;

   Error = FiFiErrNull;
   }

SEXP ExtendedSurvivorsAnalysisR::Return(int _MaxYear)
   {
   if (_MaxYear > 0)
      MaxYear = min(MaxYear,_MaxYear);

   SEXP ReturnObject = R_NilValue;

   PROTECT(ReturnObject = NEW_OBJECT(MAKE_CLASS("FLXSA")));

   int _MaxAge;

   if (MaxAge == Control.PlusGroup-1) _MaxAge = MaxAge+1; else _MaxAge = MaxAge;

   //N
   SET_SLOT(ReturnObject, Rf_install("stock.n"),  CreateFLQuant(&N, MinAge, _MaxAge, MinYear, MaxYear));

   //Terminal N std errors
   SET_SLOT(ReturnObject, Rf_install("se.int"), ReturnInternalSE());
   SET_SLOT(ReturnObject, Rf_install("se.ext"), ReturnExternalSE());

   //F
   SET_SLOT(ReturnObject, Rf_install("harvest"),  CreateFLQuant(&F, MinAge, _MaxAge, MinYear, MaxYear));

   //Control
   SET_SLOT(ReturnObject, Rf_install("control"),  ReturnControl());

   //index hat
   SEXP indexhat;
   PROTECT(indexhat = Rf_allocVector(VECSXP,pTune->GetNFleet()));
   int i;
   for (i=1; i<=pTune->GetNFleet(); i++)
      SET_VECTOR_ELT(indexhat, i-1,  ReturnIndexHat(i));

   SET_SLOT(ReturnObject, Rf_install("index.hat"), indexhat);

   //index var
   SEXP indexvar;
   PROTECT(indexvar = Rf_allocVector(VECSXP,pTune->GetNFleet()));
   for (i=1; i<=pTune->GetNFleet(); i++)
      SET_VECTOR_ELT(indexvar, i-1,  ReturnIndexVar(i));

   SET_SLOT(ReturnObject, Rf_install("index.var"), indexvar);

   //q residuals
   SEXP qreslist;
   PROTECT(qreslist = Rf_allocVector(VECSXP,pTune->GetNFleet()));
   for (i=1; i<=pTune->GetNFleet(); i++)
      SET_VECTOR_ELT(qreslist, i-1,  ReturnQRes(i));

   SET_SLOT(ReturnObject, Rf_install("index.res"), qreslist);
   
   //q 
   SEXP qlist, q2list;
   PROTECT(qlist  = Rf_allocVector(VECSXP,pTune->GetNFleet()));
   PROTECT(q2list = Rf_allocVector(VECSXP,pTune->GetNFleet()));
   for (i=1; i<=pTune->GetNFleet(); i++)
      {
      SET_VECTOR_ELT(qlist,  i-1,  ReturnQ(1, i));
      SET_VECTOR_ELT(q2list, i-1,  ReturnQ(2, i));
      }

   SET_SLOT(ReturnObject, Rf_install("q.hat"),  qlist);
   SET_SLOT(ReturnObject, Rf_install("q2.hat"), q2list);
   
   //corrected cpue
   SEXP cpuelist;
   PROTECT(cpuelist = Rf_allocVector(VECSXP,pTune->GetNFleet()));
   for (i=1; i<=pTune->GetNFleet(); i++)
      SET_VECTOR_ELT(cpuelist, i-1,  ReturnCorrectedCPUE(i));

   SET_SLOT(ReturnObject, Rf_install("index"), cpuelist);
   
   double **t; int nr, nc;
   OutputTermPopWts(&t, &nr, &nc);

   SET_SLOT(ReturnObject, Rf_install("wts"),  CreateArray(&t, 1, nr, 1, nc));

   free_flallocArray(t, 1, nr, 1, nc);

   //Survivors
   CalcTermPops();
   SET_SLOT(ReturnObject, Rf_install("survivors"), CreateFLQuant(&N, MinAge, _MaxAge, MinYear, MaxYear+1));

   UNPROTECT(7);

   return ReturnObject;
   }

SEXP ExtendedSurvivorsAnalysisR::ReturnQ(short param, short iFleet)
    {
    SEXP v = R_NilValue;   

    if (param < 1 || param > 2)
        return v;
           
    int i, 
        iAge;
    
    SEXP d1,
         dim, dimnames, names;    

    //Create array    
    //Set dimensions of array
    PROTECT(dim     = Rf_allocVector(INTSXP, 1));       
    INTEGER(dim)[0] = MaxTuneAge[iFleet] -MinTuneAge[iFleet] +1;
        
    //Allocate memory
    PROTECT(v = Rf_allocArray(REALSXP, dim)); 
    
    //Create dimension names
    PROTECT(dimnames = Rf_allocVector(VECSXP, 1));
    
    PROTECT(d1 = Rf_allocVector(INTSXP, MaxTuneAge[iFleet]-MinTuneAge[iFleet] +1));
    for (iAge=MinTuneAge[iFleet], i=0; iAge<=MaxTuneAge[iFleet]; iAge++, i++)
        INTEGER(d1)[i] = iAge; 
    SET_VECTOR_ELT(dimnames, 0, d1);


    //Create names for dimensions
    PROTECT(names = Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(names, 0, Rf_mkChar("age"));
    Rf_setAttrib(dimnames, R_NamesSymbol, names);
    Rf_setAttrib(v, R_DimNamesSymbol, dimnames);
    
    //Set data
    for (iAge=MinTuneAge[iFleet], i=0; iAge<=MaxTuneAge[iFleet]; iAge++, i++)
       REAL(v)[i] = (param == 1 ? Getq(iFleet,iAge) : Getq_b(iFleet,iAge));

    UNPROTECT(5);

    return v;
    }

SEXP ExtendedSurvivorsAnalysisR::ReturnSimple(int _MaxYear)
   {
   if (_MaxYear > 0)
      MaxYear = min(MaxYear,_MaxYear);

   SEXP ReturnObject = R_NilValue;

   PROTECT(ReturnObject = NEW_OBJECT(MAKE_CLASS("FLXSA")));

   int _MaxAge;

   if (MaxAge == Control.PlusGroup-1) _MaxAge = MaxAge+1; else _MaxAge = MaxAge;

   //N
   SET_SLOT(ReturnObject, Rf_install("stock.n"),  CreateFLQuant(&N, MinAge, _MaxAge, MinYear, MaxYear));

   //F
   SET_SLOT(ReturnObject, Rf_install("harvest"),  CreateFLQuant(&F, MinAge, _MaxAge, MinYear, MaxYear));


   UNPROTECT(1);

   return ReturnObject;
   }

SEXP ExtendedSurvivorsAnalysisR::ReturnIndexHat(short iFleet)
   {
    int i,    j, 
        iAge, iYear;
    
    SEXP v,   d1,       d2,
         dim, dimnames, names;    

    //Create array    
    //Set dimensions of array
    PROTECT(dim     = Rf_allocVector(INTSXP, 2));       
    INTEGER(dim)[0] = MaxTuneAge[iFleet] -MinTuneAge[iFleet] +1;
    INTEGER(dim)[1] = MaxTuneYear[iFleet]-MinTuneYear[iFleet]+1;
        
    //Allocate memory
    PROTECT(v = Rf_allocArray(REALSXP, dim)); 
    
    //Create dimension names
    PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
    
    PROTECT(d1 = Rf_allocVector(INTSXP, MaxTuneAge[iFleet]-MinTuneAge[iFleet] +1));
    for (iAge=MinTuneAge[iFleet], i=0; iAge<=MaxTuneAge[iFleet]; iAge++, i++)
        INTEGER(d1)[i] = iAge; 
    SET_VECTOR_ELT(dimnames, 0, d1);
    
    PROTECT(d2 = Rf_allocVector(INTSXP, MaxTuneYear[iFleet]-MinTuneYear[iFleet]+1));
    for (iYear=MinTuneYear[iFleet], i=0; iYear<=MaxTuneYear[iFleet]; iYear++, i++)
        INTEGER(d2)[i] = iYear; 
    SET_VECTOR_ELT(dimnames, 1, d2);
     
    //Create names for dimensions
    PROTECT(names = Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, Rf_mkChar("age"));
    SET_STRING_ELT(names, 1, Rf_mkChar("year"));
    Rf_setAttrib(dimnames, R_NamesSymbol, names);
    Rf_setAttrib(v, R_DimNamesSymbol, dimnames);
    
    //Set data
    for (iAge=MinTuneAge[iFleet], i=0; iAge<=MaxTuneAge[iFleet]; iAge++, i++)
       for (iYear=MinTuneYear[iFleet], j=0; iYear<=MaxTuneYear[iFleet]; iYear++, j++)
          {
          if (Theta(iFleet, iAge, iYear))
               REAL(v)[i + j*(MaxTuneAge[iFleet]-MinTuneAge[iFleet]+1)] = exp(Getq(iFleet,iAge)+Getq_b(iFleet,iAge)*log(GetNHat(iFleet,iAge,iYear)));
            else
               REAL(v)[i + j*(MaxTuneAge[iFleet]-MinTuneAge[iFleet]+1)] = 0.0;
          }
           
    UNPROTECT(6);
    
    return v;
    }

SEXP ExtendedSurvivorsAnalysisR::ReturnIndexVar(short iFleet)
   {
    int i,    j, 
        iAge, iYear;
    
    SEXP v,   d1,       d2,
         dim, dimnames, names;    

    //Create array    
    //Set dimensions of array
    PROTECT(dim     = Rf_allocVector(INTSXP, 2));       
    INTEGER(dim)[0] = MaxTuneAge[iFleet] -MinTuneAge[iFleet] +1;
    INTEGER(dim)[1] = MaxTuneYear[iFleet]-MinTuneYear[iFleet]+1;
        
    //Allocate memory
    PROTECT(v = Rf_allocArray(REALSXP, dim)); 
    
    //Create dimension names
    PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
    
    PROTECT(d1 = Rf_allocVector(INTSXP, MaxTuneAge[iFleet]-MinTuneAge[iFleet] +1));
    for (iAge=MinTuneAge[iFleet], i=0; iAge<=MaxTuneAge[iFleet]; iAge++, i++)
        INTEGER(d1)[i] = iAge; 
    SET_VECTOR_ELT(dimnames, 0, d1);
    
    PROTECT(d2 = Rf_allocVector(INTSXP, MaxTuneYear[iFleet]-MinTuneYear[iFleet]+1));
    for (iYear=MinTuneYear[iFleet], i=0; iYear<=MaxTuneYear[iFleet]; iYear++, i++)
        INTEGER(d2)[i] = iYear; 
    SET_VECTOR_ELT(dimnames, 1, d2);
     
    //Create names for dimensions
    PROTECT(names = Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, Rf_mkChar("age"));
    SET_STRING_ELT(names, 1, Rf_mkChar("year"));
    Rf_setAttrib(dimnames, R_NamesSymbol, names);
    Rf_setAttrib(v, R_DimNamesSymbol, dimnames);
    
    //Set data
    for (iAge=MinTuneAge[iFleet], i=0; iAge<=MaxTuneAge[iFleet]; iAge++, i++)
       for (iYear=MinTuneYear[iFleet], j=0; iYear<=MaxTuneYear[iFleet]; iYear++, j++)
          {
          if (Theta(iFleet, iAge, iYear) && Getqn(iFleet, iAge) > 1)
               {
               double mn = exp(Getq(iFleet,iAge)+Getq_b(iFleet,iAge)*log(GetNHat(iFleet,iAge,iYear))),
                      cv = GetNSE(iFleet,iAge,iYear)/GetNHat(iFleet,iAge,iYear);   

               REAL(v)[i + j*(MaxTuneAge[iFleet]-MinTuneAge[iFleet]+1)] = SQR(cv*mn);
               }
           else
               REAL(v)[i + j*(MaxTuneAge[iFleet]-MinTuneAge[iFleet]+1)] = 0.0;
          }

    UNPROTECT(6);
    
    return v;
    }

SEXP ExtendedSurvivorsAnalysisR::ReturnQRes(short iFleet)
   {
    int i,    j, 
        iAge, iYear;
    
    SEXP v,   d1,       d2,
         dim, dimnames, names;    

    //Create array    
    //Set dimensions of array
    PROTECT(dim     = Rf_allocVector(INTSXP, 2));       
    INTEGER(dim)[0] = MaxTuneAge[iFleet] -MinTuneAge[iFleet] +1;
    INTEGER(dim)[1] = MaxTuneYear[iFleet]-MinTuneYear[iFleet]+1;
        
    //Allocate memory
    PROTECT(v = Rf_allocArray(REALSXP, dim)); 
    
    //Create dimension names
    PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
    
    PROTECT(d1 = Rf_allocVector(INTSXP, MaxTuneAge[iFleet]-MinTuneAge[iFleet] +1));
    for (iAge=MinTuneAge[iFleet], i=0; iAge<=MaxTuneAge[iFleet]; iAge++, i++)
        INTEGER(d1)[i] = iAge; 
    SET_VECTOR_ELT(dimnames, 0, d1);
    
    PROTECT(d2 = Rf_allocVector(INTSXP, MaxTuneYear[iFleet]-MinTuneYear[iFleet]+1));
    for (iYear=MinTuneYear[iFleet], i=0; iYear<=MaxTuneYear[iFleet]; iYear++, i++)
        INTEGER(d2)[i] = iYear; 
    SET_VECTOR_ELT(dimnames, 1, d2);
     
    //Create names for dimensions
    PROTECT(names = Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, Rf_mkChar("age"));
    SET_STRING_ELT(names, 1, Rf_mkChar("year"));
    Rf_setAttrib(dimnames, R_NamesSymbol, names);
    Rf_setAttrib(v, R_DimNamesSymbol, dimnames);
    
    //Set data
    for (iAge=MinTuneAge[iFleet], i=0; iAge<=MaxTuneAge[iFleet]; iAge++, i++)
       for (iYear=MinTuneYear[iFleet], j=0; iYear<=MaxTuneYear[iFleet]; iYear++, j++)
          {
          if (Theta(iFleet, iAge, iYear))
              REAL(v)[i + j*(MaxTuneAge[iFleet]-MinTuneAge[iFleet]+1)] = log(pow(GetNHat(iFleet,iAge,iYear)/N[iAge][iYear],Getq_b(iFleet,iAge)));
          }
           
    UNPROTECT(6);
    
    return v;
    }

SEXP ExtendedSurvivorsAnalysisR::ReturnCorrectedCPUE(short iFleet, short _MaxYear)
   {
    int i,    j, 
        iAge, iYear,
        MaxYear;
    
    if (_MaxYear > MinTuneYear[iFleet])
      MaxYear = min(_MaxYear,MaxTuneYear[iFleet]);
    else
      MaxYear = MaxTuneYear[iFleet];

    SEXP v,   d1,       d2,
         dim, dimnames, names;    

    //Create array    
    //Set dimensions of array
    PROTECT(dim     = Rf_allocVector(INTSXP, 2));       
    INTEGER(dim)[0] = MaxTuneAge[iFleet] -MinTuneAge[iFleet] +1;
    INTEGER(dim)[1] = MaxYear-MinTuneYear[iFleet]+1;
        
    //Allocate memory
    PROTECT(v = Rf_allocArray(REALSXP, dim)); 
    
    //Create dimension names
    PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
    
    PROTECT(d1 = Rf_allocVector(INTSXP, MaxTuneAge[iFleet]-MinTuneAge[iFleet] +1));
    for (iAge=MinTuneAge[iFleet], i=0; iAge<=MaxTuneAge[iFleet]; iAge++, i++)
        INTEGER(d1)[i] = iAge; 
    SET_VECTOR_ELT(dimnames, 0, d1);
    
    PROTECT(d2 = Rf_allocVector(INTSXP, MaxTuneYear[iFleet]-MinTuneYear[iFleet]+1));
    for (iYear=MinTuneYear[iFleet], i=0; iYear<=MaxYear; iYear++, i++)
        INTEGER(d2)[i] = iYear; 
    SET_VECTOR_ELT(dimnames, 1, d2);
     
    //Create names for dimensions
    PROTECT(names = Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, Rf_mkChar("age"));
    SET_STRING_ELT(names, 1, Rf_mkChar("year"));
    Rf_setAttrib(dimnames, R_NamesSymbol, names);
    Rf_setAttrib(v, R_DimNamesSymbol, dimnames);
    
    //Set data
    for (iAge=MinTuneAge[iFleet], i=0; iAge<=MaxTuneAge[iFleet]; iAge++, i++)
       for (iYear=MinTuneYear[iFleet], j=0; iYear<=MaxYear; iYear++, j++)
          {
          if (Theta(iFleet, iAge, iYear))
               REAL(v)[i + j*(MaxTuneAge[iFleet]-MinTuneAge[iFleet]+1)] = GetCorrectedCPUE(iFleet, iAge, iYear);
            else
               REAL(v)[i + j*(MaxTuneAge[iFleet]-MinTuneAge[iFleet]+1)] = 0.0;
          }
           
    UNPROTECT(6);
    
    return v;
    }

SEXP ExtendedSurvivorsAnalysisR::ReturnInternalSE(void)
    {
    int i, iAge;
    
    SEXP FLQuant, v, 
         d1, d2, d3, d4, d5, d6, 
         dim,   dimnames, names;    

    //Create new S4 object    
    //PROTECT(v = NEW_OBJECT(MAKE_CLASS("array")));
    PROTECT(FLQuant = NEW_OBJECT(MAKE_CLASS("FLQuant")));

    //Create array for slot    
    //Set dimensions of array
    PROTECT(dim     = Rf_allocVector(INTSXP, 6));       
    INTEGER(dim)[0] = MaxAge -MinAge +1;
    INTEGER(dim)[1] = 1;
    INTEGER(dim)[2] = 1;
    INTEGER(dim)[3] = 1;
    INTEGER(dim)[4] = 1;
    INTEGER(dim)[5] = 1;
        
    //Allocate memory
    PROTECT(v = Rf_allocArray(REALSXP, dim)); 
    
    //Create dimension names
    PROTECT(dimnames = Rf_allocVector(VECSXP, 6));
    
    PROTECT(d1 = Rf_allocVector(INTSXP, MaxAge-MinAge +1));
    for (iAge=MinAge, i=0; iAge<=MaxAge; iAge++, i++)
        INTEGER(d1)[i] = iAge; 
    SET_VECTOR_ELT(dimnames, 0, d1);
    
    PROTECT(d2 = Rf_allocVector(INTSXP, 1));
    INTEGER(d2)[0] = MaxYear+1;
    SET_VECTOR_ELT(dimnames, 1, d2);
     
    PROTECT(d3 = Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(d3, 0, Rf_mkChar("unique"));
    SET_VECTOR_ELT(dimnames, 2, d3);
    
    PROTECT(d4 = Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(d4, 0, Rf_mkChar("all"));
    SET_VECTOR_ELT(dimnames, 3, d4);

    PROTECT(d5 = Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(d5, 0, Rf_mkChar("unique"));
    SET_VECTOR_ELT(dimnames, 4, d5);

    PROTECT(d6 = Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(d6, 0, Rf_mkChar("1"));
    SET_VECTOR_ELT(dimnames, 5, d6);
    
    //Create names for dimensions
    PROTECT(names = Rf_allocVector(STRSXP, 6));
    SET_STRING_ELT(names, 0, Rf_mkChar("age"));
    SET_STRING_ELT(names, 1, Rf_mkChar("year"));
    SET_STRING_ELT(names, 2, Rf_mkChar("unit"));
    SET_STRING_ELT(names, 3, Rf_mkChar("season"));
    SET_STRING_ELT(names, 4, Rf_mkChar("area"));
    SET_STRING_ELT(names, 5, Rf_mkChar("iter"));
    Rf_setAttrib(dimnames, R_NamesSymbol, names);
    Rf_setAttrib(v, R_DimNamesSymbol, dimnames);

    //Set data
    for (iAge=MinAge, i=0; iAge<=MaxAge; iAge++, i++) 
       REAL(v)[i] = InternalSE[0][MaxYear-iAge];
        
    //Set slot
    //FLQuant = SET_SLOT(FLQuant, Rf_install("v"), v);
    FLQuant = R_do_slot_assign(FLQuant, Rf_install(".Data"), v);

    UNPROTECT(11);
    
    return FLQuant;
    }

SEXP ExtendedSurvivorsAnalysisR::ReturnExternalSE(void)
    {
    int i, iAge;
    
    SEXP FLQuant, v, 
         d1, d2, d3, d4, d5, d6, 
         dim,   dimnames, names;    

    //Create new S4 object    
    PROTECT(FLQuant = NEW_OBJECT(MAKE_CLASS("FLQuant")));

    //Create array for slot    
    //Set dimensions of array
    PROTECT(dim     = Rf_allocVector(INTSXP, 6));       
    INTEGER(dim)[0] = MaxAge -MinAge +1;
    INTEGER(dim)[1] = 1;
    INTEGER(dim)[2] = 1;
    INTEGER(dim)[3] = 1;
    INTEGER(dim)[4] = 1;
    INTEGER(dim)[5] = 1;
        
    //Allocate memory
    PROTECT(v = Rf_allocArray(REALSXP, dim)); 
    
    //Create dimension names
    PROTECT(dimnames = Rf_allocVector(VECSXP, 6));
    
    PROTECT(d1 = Rf_allocVector(INTSXP, MaxAge-MinAge +1));
    for (iAge=MinAge, i=0; iAge<=MaxAge; iAge++, i++)
        INTEGER(d1)[i] = iAge; 
    SET_VECTOR_ELT(dimnames, 0, d1);
    
    PROTECT(d2 = Rf_allocVector(INTSXP, 1));
    INTEGER(d2)[0] = MaxYear+1;
    SET_VECTOR_ELT(dimnames, 1, d2);
     
    PROTECT(d3 = Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(d3, 0, Rf_mkChar("unique"));
    SET_VECTOR_ELT(dimnames, 2, d3);
    
    PROTECT(d4 = Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(d4, 0, Rf_mkChar("all"));
    SET_VECTOR_ELT(dimnames, 3, d4);

    PROTECT(d5 = Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(d5, 0, Rf_mkChar("unique"));
    SET_VECTOR_ELT(dimnames, 4, d5);
    
    PROTECT(d6 = Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(d6, 0, Rf_mkChar("1"));
    SET_VECTOR_ELT(dimnames, 5, d6);
    
    //Create names for dimensions
    PROTECT(names = Rf_allocVector(STRSXP, 6));
    SET_STRING_ELT(names, 0, Rf_mkChar("age"));
    SET_STRING_ELT(names, 1, Rf_mkChar("year"));
    SET_STRING_ELT(names, 2, Rf_mkChar("unit"));
    SET_STRING_ELT(names, 3, Rf_mkChar("season"));
    SET_STRING_ELT(names, 4, Rf_mkChar("area"));
    SET_STRING_ELT(names, 5, Rf_mkChar("iter"));
    Rf_setAttrib(dimnames, R_NamesSymbol, names);
    Rf_setAttrib(v, R_DimNamesSymbol, dimnames);

    //Set data
    for (iAge=MinAge, i=0; iAge<=MaxAge; iAge++, i++) 
       REAL(v)[i] = ExternalSE[0][MaxYear-iAge];
        
    //Set slot
    //FLQuant = SET_SLOT(FLQuant, Rf_install("v"), v);
    FLQuant = R_do_slot_assign(FLQuant, Rf_install(".Data"), v);

    UNPROTECT(11);
    
    return FLQuant;
    }

SEXP ExtendedSurvivorsAnalysisR::ReturnControl(void)
   {
   SEXP xControl = R_NilValue;

   PROTECT(xControl = NEW_OBJECT(MAKE_CLASS("FLXSA.control")));

   SEXP MSStol       = R_NilValue,  
        MaxIters     = R_NilValue,  
        minNse       = R_NilValue,  
        FShkSE       = R_NilValue,  
        LRcrtAge     = R_NilValue,  
        ConQAge      = R_NilValue,  
        shkN         = R_NilValue,  
        shkF         = R_NilValue,  
        shkyrs       = R_NilValue,  
        shkages      = R_NilValue,  
        tsrange      = R_NilValue,  
        tspower      = R_NilValue,  
        plusgroup    = R_NilValue,
        TuningWindow = R_NilValue,
        VPA          = R_NilValue;

   MSStol        = PROTECT(NEW_NUMERIC(1));        
   MaxIters      = PROTECT(NEW_INTEGER(1));
   minNse        = PROTECT(NEW_NUMERIC(1));
   FShkSE        = PROTECT(NEW_NUMERIC(1));
   LRcrtAge      = PROTECT(NEW_INTEGER(1));
   ConQAge       = PROTECT(NEW_INTEGER(1));
   shkN          = PROTECT(NEW_INTEGER(1));
   shkF          = PROTECT(NEW_INTEGER(1));
   shkyrs        = PROTECT(NEW_INTEGER(1));
   shkages       = PROTECT(NEW_INTEGER(1));
   tsrange       = PROTECT(NEW_NUMERIC(1));
   tspower       = PROTECT(NEW_NUMERIC(1));
   plusgroup     = PROTECT(NEW_LOGICAL(1));
   TuningWindow  = PROTECT(NEW_INTEGER(1));
   VPA           = PROTECT(NEW_INTEGER(1));
   
   REAL(MSStol)[0]          = MSS;
   INTEGER(MaxIters)[0]     = Iters;   
   REAL(minNse)[0]          = Control.MinNSE;
   REAL(FShkSE)[0]          = Control.FShkSE;      
   INTEGER(LRcrtAge)[0]     = Control.LRcrtAge;  
   INTEGER(ConQAge)[0]      = Control.ConQAge;   
   INTEGER(shkN)[0]         = Control.Shk2N;     
   INTEGER(shkF)[0]         = Control.Shk2F;      
   INTEGER(shkyrs)[0]       = Control.Shk2FYr;    
   INTEGER(shkages)[0]      = Control.Shk2FAge;    
   REAL(tsrange)[0]         = Control.TSRange;     
   REAL(tspower)[0]         = Control.TSPower;
   INTEGER(TuningWindow)[0] = Control.TuningWindow;
   INTEGER(VPA)[0]          = Control.VPA;

   SET_SLOT(xControl,Rf_install("tol"          ), MSStol);
   SET_SLOT(xControl,Rf_install("maxit"        ), MaxIters);   
   SET_SLOT(xControl,Rf_install("min.nse"      ), minNse);
   SET_SLOT(xControl,Rf_install("fse"          ), FShkSE);      
   SET_SLOT(xControl,Rf_install("rage"         ), LRcrtAge);  
   SET_SLOT(xControl,Rf_install("qage"         ), ConQAge);   
   SET_SLOT(xControl,Rf_install("shk.n"        ), shkN);     
   SET_SLOT(xControl,Rf_install("shk.f"        ), shkF);      
   SET_SLOT(xControl,Rf_install("shk.yrs"      ), shkyrs);    
   SET_SLOT(xControl,Rf_install("shk.ages"     ), shkages);    
   SET_SLOT(xControl,Rf_install("tsrange"      ), tsrange);     
   SET_SLOT(xControl,Rf_install("tspower"      ), tspower);     
   SET_SLOT(xControl,Rf_install("window"       ), TuningWindow);
   SET_SLOT(xControl,Rf_install("vpa"          ), VPA);
   
   UNPROTECT(16);
   
   return xControl;
   }

SEXP ExtendedSurvivorsAnalysisR::dataframe(void)
   {
   SEXP data, tmp;
   
   LP2DOUBLE t;

   flallocArray(&t, 1, 10, 10, 19);
   for (int i=1; i<=10; i++)
      for (int j=10; j<=19; j++)
         t[i][j] = i*1000+j;

   data = CreateArray(&t, 1, 10, 10, 19);
   free_flallocArray(t, 1, 10, 10, 19);
      
   /* Assemble the base data frame. */
//   PROTECT(names = Rf_allocVector(STRSXP, 1));

   //SET_STRING_ELT(names, 0, Rf_mkChar("fleet"));

   //Rf_setAttrib(data, R_NamesSymbol, names);
   
   UNPROTECT(2);

   /* Turn the data "list" into a "data.frame" */
   /* so that subsetting methods will work. */
   /* To do this we must attach "class"  and */
   /* "row.names" attributes */

   PROTECT(tmp = Rf_mkString("data.frame"));
//   Rf_setAttrib(data, R_ClassSymbol, tmp);
   UNPROTECT(1);

   return data;
   }

bool ExtendedSurvivorsAnalysisR::OutputTermPopWts(LP2DOUBLE *pResults, int *NRow, int *NCol)
   {     
   short  iAge, iFleet, iYear = 0, iYearClass, RecruitAge;
   double ECF, ECZ, Wt;           

   //"Source, log(NHat*Wt), Wt, NHat, YearClass, Fleet, Age, Year";
   //"Wt, NHat, YearClass, Fleet, Age";

   //count up number of rows
   (*NCol) = 6;
   (*NRow) = 0;
   for (iYearClass=MinYear-MaxAge; iYearClass <= MaxYear-MinAge; iYearClass++)
      { 
      for (iAge = min(MaxAge, MaxYear - iYearClass); iAge >= max(MinAge, MinYear - iYearClass); iAge--)
         {
         iYear = iYearClass + iAge;
         for (iFleet=0, pTune->IncrFleet(&iFleet); iFleet<=NTuneFleet; pTune->IncrFleet(&iFleet))
            {
            if (GetMaxTuneAge(iFleet)  >= iAge  && GetMinTuneAge(iFleet)  <= iAge &&
                GetMaxTuneYear(iFleet) >= iYear && GetMinTuneYear(iFleet) <= iYear  )
                if (Theta(iFleet, iAge, iYear) && Getqn(iFleet, iAge) > 1)
                   {
                   (*NRow)++;
                   }
            }       
         }
      if (iYear >= MinYear && Control.Shk2F == TRUE && GetFShrinkMean(iYearClass) > 0 && Control.FShkSE > 0.0) 
         (*NRow)++;

      RecruitAge = MaxYear - iYearClass;                                                 

      if (iYear == MaxYear && Control.Shk2N == TRUE && RecruitAge <= Control.LRcrtAge && GetPopShrinkMean(RecruitAge) > 0.0 && PopShrinkSE[RecruitAge] > 0.0)             
         (*NRow)++;
      }

   flallocArray(pResults, 1, *NRow, 1, *NCol);

   int i=1;
   for (iYearClass=MinYear-MaxAge; iYearClass <= MaxYear-MinAge; iYearClass++)
      { 
      ECF            =                
      ECZ            = 1.0;
      for (iAge = min(MaxAge, MaxYear - iYearClass); iAge >= max(MinAge, MinYear - iYearClass); iAge--)
         {
         iYear = iYearClass + iAge; 
         ECF *= exp(F[iAge][iYear]);               
         ECZ *= exp(F[iAge][iYear] + GetM(iAge, iYear));
         for (iFleet=0, pTune->IncrFleet(&iFleet); iFleet<=NTuneFleet; pTune->IncrFleet(&iFleet))
            {
            if (GetMaxTuneAge(iFleet)  >= iAge  && GetMinTuneAge(iFleet)  <= iAge &&
                GetMaxTuneYear(iFleet) >= iYear && GetMinTuneYear(iFleet) <= iYear  )
                if (Theta(iFleet, iAge, iYear) && Getqn(iFleet, iAge) > 1)
                   {
                   if (Control.SmoothQ)
                      Wt = Weight(iFleet, iAge, iYear, iYearClass+MaxAge, ECF);                                                                          
                   else
                      Wt = Weight(iFleet, iAge, iYear, ECF);                                                                          
                          
                   (*pResults)[i][1]   = Wt;
                   (*pResults)[i][2]   = log(GetNHat(iFleet,iAge,iYear)) - log(ECZ);
                   (*pResults)[i][3]   = iYearClass;
                   (*pResults)[i][4]   = iFleet;
                   (*pResults)[i][5]   = iAge;
                   (*pResults)[i++][6] = iYear;
                   }
            }       
         }

      //"Source, log(NHat*Wt), Wt, NHat, YearClass, Fleet, Age, Year";
      //"Wt, NHat, YearClass, Fleet, Age, Year";

      if (iYear >= MinYear && Control.Shk2F == TRUE && GetFShrinkMean(iYearClass) > 0 && Control.FShkSE > 0.0) 
         {
         (*pResults)[i][1]   = 1/(SQR(Control.FShkSE));
         if (GetFShrinkMean(iYearClass)>0.0)
            (*pResults)[i][2]   = log(GetFShrinkMean(iYearClass));
         else
            (*pResults)[i][2]   = 10e-19;
         (*pResults)[i][3]   = iYearClass;
         (*pResults)[i][4]   = -1;
         (*pResults)[i][5]   = iYearClass <= (MaxYear - MaxAge) ? MaxAge : MaxYear-iYear+MinAge;
         (*pResults)[i++][6] = iYearClass <= (MaxYear - MaxAge) ? MaxAge+iYearClass : MaxYear;
         }

      iAge = min(MaxAge, MaxYear - iYearClass);
      RecruitAge = MaxYear - iYearClass;                                                 

      if (iYear == MaxYear && Control.Shk2N == TRUE && RecruitAge <= Control.LRcrtAge && GetPopShrinkMean(RecruitAge) > 0.0 && PopShrinkSE[RecruitAge] > 0.0)             
         {
         (*pResults)[i][1]   = 1/(SQR(GetPopShrinkSE(RecruitAge))) ;
         if (GetPopShrinkMean(RecruitAge)>0.0)
            (*pResults)[i][2]   = log(GetPopShrinkMean(RecruitAge));
         else
            (*pResults)[i][2]   = 10e-19;
         (*pResults)[i][3]   = iYearClass;
         (*pResults)[i][4]   = -2;
         (*pResults)[i][5]   = iAge;
         (*pResults)[i++][6] = MaxYear;
         }
      } 
      
   return true;
   }           

void ExtendedSurvivorsAnalysisR::InputNF(SEXP xStock)
   { 
   TestFlag = false;

   SEXP _N = PROTECT(Rf_duplicate(GET_SLOT(xStock, Rf_install("stock.n"))));
   SEXP _F = PROTECT(Rf_duplicate(GET_SLOT(xStock, Rf_install("harvest"))));
   
   if (InputFLQuant(_N, &N, MinAge, MaxAge, MinYear, MaxYear) &&
       InputFLQuant(_F, &F, MinAge, MaxAge, MinYear, MaxYear)    )
      TestFlag = true;
      
   UNPROTECT(2);      
   }

TuningFleetsR::TuningFleetsR(SEXP x, CatchData *pCatch, short TuningWindow)
   {
   Error = !FiFiErrNull;

   if (Catch == NULL) // array has not been initialized
      Error  = TuningDataErrInput;
   
   ArSmry.NFleet = NElemList(x);
   if (ArSmry.NFleet <= 0)
      return;

   MaxTuningAge  = pCatch->GetMaxAge();
   MinTuningAge  = pCatch->GetMinAge();
   MaxTuningYear = pCatch->GetMaxYear();
   MinTuningYear = pCatch->GetMinYear();
   if (TuningWindow>=2)
      MinTuningYear = max(pCatch->GetMinYear(),pCatch->GetMinYear()-TuningWindow+1);

   AllocArrays();
   
   //Get age & year ranges
   int iFleet;
   SEXP xIndex,
        vIndex;

   for (iFleet = 1; iFleet <= ArSmry.NFleet; iFleet++)
      {
      //int MinDimCatch[5],  MaxDimCatch[5],
      //    MinDimEffort[5], MaxDimEffort[5]; 

      xIndex = PROTECT(Rf_duplicate(GET_SLOT(VECTOR_ELT(x, iFleet-1), Rf_install("index"))));
      vIndex = PROTECT(Rf_duplicate(GET_SLOT(xIndex, Rf_install(".Data"))));
//      GetFLQuantDims(vIndex, &MinDimCatch, &MaxDimCatch);
      
//      SEXP xEffort = PROTECT(Rf_duplicate(GET_SLOT(VECTOR_ELT(x, iFleet-1), Rf_install("effort"))));     
//      SEXP vEffort = PROTECT(Rf_duplicate(GET_SLOT(xEffort, Rf_install(".Data"))));
//      GetFLQuantDims(vEffort, &MinDimEffort, &MaxDimEffort);
      
       //Get range
      SEXP range    = PROTECT(Rf_duplicate(GET_SLOT(VECTOR_ELT(x, iFleet-1), Rf_install("range"))));    

      alpha[iFleet] = REAL(range)[5];
      beta[ iFleet] = REAL(range)[6];
     
      ArSmry.MinAge[iFleet]  = (int)REAL(range)[0];
      ArSmry.MaxAge[iFleet]  = (int)REAL(range)[1];
      ArSmry.MinYear[iFleet] = (int)REAL(range)[3];
      ArSmry.MaxYear[iFleet] = (int)REAL(range)[4];
      
      UNPROTECT(3);
      }
      
   for (iFleet = 1; iFleet <= ArSmry.NFleet; iFleet++)
      {
      flallocArray(&(Catch[iFleet]),  ArSmry.MinAge[iFleet],  ArSmry.MaxAge[iFleet], ArSmry.MinYear[iFleet], ArSmry.MaxYear[iFleet]);
      flallocArray(&(Effort[iFleet]), ArSmry.MinYear[iFleet], ArSmry.MaxYear[iFleet]);
      }

   
   for (iFleet = 1; iFleet <= ArSmry.NFleet; iFleet++)
      {
      SEXP xIndex  = PROTECT(Rf_duplicate(GET_SLOT(VECTOR_ELT(x, iFleet-1), Rf_install("index"))));
      InputFLQuant(xIndex, &(Catch[iFleet]), ArSmry.MinAge[iFleet], ArSmry.MaxAge[iFleet], ArSmry.MinYear[iFleet], ArSmry.MaxYear[iFleet]);
      
      //SEXP xEffort = PROTECT(Rf_duplicate(GET_SLOT(VECTOR_ELT(x, iFleet-1), Rf_install("effort"))));     
      //InputFLQuant(xEffort, &(Effort[iFleet]), ArSmry.MinYear[iFleet], ArSmry.MaxYear[iFleet]);
      
      for (int iYear = ArSmry.MinYear[iFleet]; iYear <= ArSmry.MaxYear[iFleet]; iYear++)
         Effort[iFleet][iYear] = 1.0;
      
      UNPROTECT(1); 
      }  

   Error = FiFiErrNull;
   } 

CatchDataR::CatchDataR(SEXP x)
   {  
   ArraysSet = FALSE;
   Error     = FiFiErrNull;
   MinAge    = MaxAge  =   
   MinYear   = MaxYear = -9;
       
   if (InputFLQuant(x, &Catch, &MinAge, &MaxAge, &MinYear, &MaxYear))
      ArraysSet = TRUE;
   else
      Error = FiFiErrInput;
   }
