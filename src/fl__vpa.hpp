#ifndef _INC_fl__vpa
#define _INC_fl__vpa
                           
#include "fl__types.hpp"
#include "fl__ctc.hpp"

#define VPA_TOL     1e-6
#define SEPVPA_TOL  1e-40
#define VPA_ITS     200     

class VirtualPopAnalysis
{
public:        
   short  Error;
   bool   ArraysSet;
   
   CatchData *pCatch;

   VirtualPopAnalysis(void);     
   VirtualPopAnalysis(LPSHORT,LPSHORT,LPDOUBLE,LPDOUBLE,LPSHORT,LPSHORT,LPDOUBLE,CatchData *);
   void Set(LPSHORT,LPSHORT,LPDOUBLE,LPDOUBLE,LPSHORT,LPSHORT,LPDOUBLE,CatchData *);
                                                           
   short   GetMinAge(void);
   short   GetMaxAge(void);   
   short   GetMinYear(void);  
   short   GetMaxYear(void);                     

   double  GetN(short, short);
   double  GetF(short, short); 
   double  GetM(short); 
   double  GetM(short, short); 

   bool    RunCohortAnalysis(void);
   bool    RunVPA(void);

   ~VirtualPopAnalysis(void);

protected: 
   short MinAge,
         MaxAge,
         MinYear,
         MaxYear; 

   bool FRatioFlag;
   
   double FRatio;
   
   LP2DOUBLE  M,
              N,
              F;   

   void SetRange(LPSHORT,LPSHORT,short);
   void InitialiseMatrix(LP2DOUBLE);
   void CheckRangeWithCatch(short);
                                             
   void SetM(short, double); 
   void SetM(short, short, double);
   void SetF(short, short, double);
   void SetTerminalNs(void);
    
   double Calcf(double, double, double, double);
   double Calcdfdx(double, double, double);
   double NewtonRhapson(double, double, double);

   inline double CalcY(double, double, double, double);
   inline double CalcdYdF(double, double, double, double);

   void CalcPlusGroupFN(void);

   bool AllocArrays(void);                
   void FreeAllocArrays(void);
};

#endif /* _INC_fl__vpa */
