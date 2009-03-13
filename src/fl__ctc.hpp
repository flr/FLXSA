#ifndef _INC_fl__ctc
#define _INC_fl__ctc

#include "fl__types.hpp"
#include "flalloc.hpp"

class CatchData
{
public:        
   short Error; 
   bool  ArraysSet;

   CatchData(void);           
   CatchData(LPSHORT,LPSHORT,LPDOUBLE,LPSHORT);
   
   bool Set(LPSHORT,LPSHORT,LPDOUBLE,LPSHORT);
                                       
   inline short  GetMinAge(void)  {return MinAge;}
   inline short  GetMaxAge(void)  {return MaxAge;}
   inline short  GetMinYear(void) {return MinYear;}
   inline short  GetMaxYear(void) {return MaxYear;}
                     
   double GetCatch(short, short);   
   
   void Check(void);
   
   
   ~CatchData();
 
protected:
   short MinAge,
         MaxAge,
         MinYear,
         MaxYear; 

   LP2DOUBLE Catch;                  

   void SetCatch(short, short, double); 
   void SetRange(LPSHORT,LPSHORT,short);
   void InitialiseMatrix(LP2DOUBLE);
                                       
   bool AllocArrays(void);
   void FreeAllocArrays(void);
};                  

#endif /* _INC_fl__ctc */

