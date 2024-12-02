#include "flrxsa.h"

extern "C" SEXP runFLXSA(SEXP Stock, SEXP CPUE, SEXP Control, SEXP diags) 
   {
   int MinAge, MaxAge, Plusgroup, MinYear, MaxYear;

   char *err1 = const_cast<char*>("Error in FLXSA.control");
   char *err2 = const_cast<char*>("Error in FLIndices");
   char *err3 = const_cast<char*>("Error in FLXSA.control");
   char *err4 = const_cast<char*>("Error in running XSA");

   //Input Catch Data
   CatchDataR Catch(PROTECT(Rf_duplicate(GET_SLOT(Stock, Rf_install("catch.n"))))); 
   UNPROTECT(1);      
   if (Catch.Error != FiFiErrNull)
      return ReturnError(err1);

   //Input Tuning Data
   TuningFleetsR TuningData(CPUE, &Catch);       
   if (TuningData.Error != FiFiErrNull)
      return ReturnError(err2);
      
   //Input XSAControls
   ExtendedSurvivorsAnalysisR XSA(Control, PROTECT(Rf_duplicate(GET_SLOT(Stock, Rf_install("m")))), &Catch, &TuningData); 
   UNPROTECT(1);      
   if (XSA.Error != FiFiErrNull)
      return ReturnError(err3);

   InputRange(PROTECT(Rf_duplicate(GET_SLOT(Stock, Rf_install("range")))), &MinAge, &MaxAge, &Plusgroup, &MinYear, &MaxYear);
   UNPROTECT(1);      
    
   //Run XSA
   //if N & F supplied then do not iterate
  // XSA.InputNF(Stock);


    SEXP v = R_NilValue;   

   if (!XSA.Run()) 
       return ReturnError(err4);

   //Results  
   if (Rf_isLogical(diags) && !LOGICAL(diags)[0])
      return XSA.ReturnSimple(MaxYear); 
   else
         return XSA.Return(MaxYear); 
   }



