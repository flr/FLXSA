#include "fl__ctc.hpp"

CatchData::CatchData(void)
   {
   Error     = !FiFiErrNull;
   ArraysSet = FALSE;
   }
   
CatchData::CatchData(LPSHORT pAge, LPSHORT pYear, LPDOUBLE pCatch, LPSHORT pn) 
   {
   Error     = !FiFiErrNull;
   ArraysSet = FALSE; 
   
   if (!Set(pAge, pYear, pCatch, pn))
      return;
   
   Error = FiFiErrNull;
   }

bool CatchData::Set(LPSHORT pAge, LPSHORT pYear, LPDOUBLE pCatch, LPSHORT pn) 
   {
   Error     = !FiFiErrNull;

   //This sets the Age and Year ranges
   SetRange(pAge,pYear,*pn);
   
   //This allocates the arrays
   if (!AllocArrays())
      return FALSE;           
      
   //Initialise Arrays in case missing values in input vectors
   InitialiseMatrix(Catch);
                        
   //Input Data
   for (short i = 0; i < *pn; i++)
      SetCatch((short)pAge[i],(short)pYear[i],pCatch[i]);

   Error = FiFiErrNull;   
   
   return TRUE;
   }

void CatchData::SetRange(LPSHORT Age,LPSHORT Year,short n)
   {   
   MinAge   = 
   MaxAge   = (short) Age[0];
   MinYear  = 
   MaxYear  = (short) Year[0];
   
   for (short i = 1; i < n; i++)
      {
      MinAge   = (short)min(MinAge,  Age[i]);
      MaxAge   = (short)max(MaxAge,  Age[i]);
      MinYear  = (short)min(MinYear, Year[i]);
      MaxYear  = (short)max(MaxYear, Year[i]);
      }
   }

void CatchData::InitialiseMatrix(LP2DOUBLE Value)
   {
   for (short iAge=MinAge; iAge<=MaxAge; iAge++)
      for (short iYear=MinYear; iYear<=MaxYear; iYear++)
//          if (!IsBadReadPtr(&Value[iAge][iYear], sizeof(short)))
              Value[iAge][iYear] = 0.0;
   }   

void CatchData::SetCatch(short Age, short Year, double Value)
   {
   if (Age < MinAge || Age > MaxAge || Year < MinYear || Year > MaxYear)
      return;
   else
     Catch[Age][Year] = (double)min(100e10,max(0.0f,Value));
   }   

double CatchData::GetCatch(short Age, short Year)  
   { 
   if (Age < MinAge || Age > MaxAge || Year < MinYear || Year > MaxYear)
      return 0.0;
   else
      return max(0.0,Catch[Age][Year]); 
   } 
               
bool CatchData::AllocArrays(void)
   {                           
   if (ArraysSet)
      FreeAllocArrays();  
   
   if (!flallocArray(&Catch, MinAge, MaxAge, MinYear, MaxYear))
      return FALSE;

   ArraysSet = TRUE;

   return TRUE;
   }      

void CatchData::FreeAllocArrays(void)
   {                                    
   if (!ArraysSet)
      return;
                                                              
   free_flallocArray(Catch, MinAge, MaxAge, MinYear, MaxYear);

   ArraysSet = FALSE;
   }                 
                
CatchData::~CatchData(void)
   {                    
   FreeAllocArrays();
   } 

void CatchData::Check(void)
   {
   int i;
   
   i=1;
   } 
