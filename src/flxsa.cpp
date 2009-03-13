#include <math.h>
#include <float.h>

#ifdef WIN32
   #include <windows.h>
   #define SEXPDLLExport __declspec(dllexport) SEXP __cdecl    
#else
   #define SEXPDLLExport SEXP    
#endif

#include "flrxsa.hpp"

extern "C" SEXPDLLExport FLXSA(SEXP Stock, SEXP CPUE, SEXP Control, SEXP diags) 
   {
   int MinAge, MaxAge, Plusgroup, MinYear, MaxYear;

   char* err1 = "Error in FLXSA.control";
   char* err2 = "Error in FLIndices";
   char* err3 = "Error in FLXSA.control";
   char* err4 = "Error in running XSA";

   //Input Catch Data
   CatchDataR Catch(PROTECT(duplicate(GET_SLOT(Stock, install("catch.n"))))); 
   UNPROTECT(1);      
   if (Catch.Error != FiFiErrNull)
      return ReturnError(err1);

   //Input Tuning Data
   TuningFleetsR TuningData(CPUE, &Catch);       
   if (TuningData.Error != FiFiErrNull)
      return ReturnError(err2);
      
   //Input XSAControls
   ExtendedSurvivorsAnalysisR XSA(Control, PROTECT(duplicate(GET_SLOT(Stock, install("m")))), &Catch, &TuningData); 
   UNPROTECT(1);      
   if (XSA.Error != FiFiErrNull)
      return ReturnError(err3);

   InputRange(PROTECT(duplicate(GET_SLOT(Stock, install("range")))), &MinAge, &MaxAge, &Plusgroup, &MinYear, &MaxYear);
   UNPROTECT(1);      
    
   //Run XSA
   //if N & F supplied then do not iterate
  // XSA.InputNF(Stock);


    SEXP v = R_NilValue;   

   if (!XSA.Run()) 
       return ReturnError(err4);

   //Results  
   if (isLogical(diags) && !LOGICAL(diags)[0])
      return XSA.ReturnSimple(MaxYear); 
   else
         return XSA.Return(MaxYear); 
   }



