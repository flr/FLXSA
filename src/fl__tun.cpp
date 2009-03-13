#include "fl__tun.hpp"

TuningFleets::TuningFleets(void)
   {
   Error     = !FiFiErrNull;
   ArraysSet = FALSE;  

   Units     = FLcNumber;
   SelType   = FLcCPUEbyAge;
   } 

bool TuningFleets::SetCatch(short Fleet, short Age, short Year, double Value)
   {                
   if (Fleet < 1                     || Fleet > ArSmry.NFleet        ||   
       Age   < ArSmry.MinAge[Fleet]  || Age   > ArSmry.MaxAge[Fleet] ||
       Year  < ArSmry.MinYear[Fleet] || Year  > ArSmry.MaxYear[Fleet]  )
      return FALSE;

   Catch[Fleet][Age][Year] = (double)min(10e100,max(0.0f,Value));
   
   return TRUE;
   }   

short  TuningFleets::GetMinAge(short Fleet)    
   {               
   return max(ArSmry.MinAge[Fleet], MinTuningAge); 
   } 
   
short  TuningFleets::GetMaxAge(short Fleet)    
   { 
   return min(ArSmry.MaxAge[Fleet], MaxTuningAge); 
   }
   
short  TuningFleets::GetMinYear(short Fleet) 
   { 
   return max(ArSmry.MinYear[Fleet], MinTuningYear); 
   }
   
short  TuningFleets::GetMaxYear(short Fleet)  
   { 
   return min(ArSmry.MaxYear[Fleet], MaxTuningYear); 
   } 

double TuningFleets::GetEffort(short Fleet, short Year)  
   { 
   return Effort[max(min(Fleet,ArSmry.NFleet),0)][max(min(Year,MaxTuningYear),MinTuningYear)]; 
   }      

double TuningFleets::GetCatch(short Fleet, short Age, short Year)  
   { 
   if (Fleet < 1                     || Fleet > ArSmry.NFleet         ||   
       Age   < ArSmry.MinAge[Fleet]  || Age   > ArSmry.MaxAge[Fleet] ||
       Year  < ArSmry.MinYear[Fleet] || Year  > ArSmry.MaxYear[Fleet]  )
      return 0.0;

   return max(0.0,Catch[Fleet][Age][Year]); 
   }      
                       
short TuningFleets::GetNFleet()  
   { 
   return ArSmry.NFleet; 
   } 

short TuningFleets::AddSeason(short iFleet)
   {
   return Season[min(max(iFleet,1),ArSmry.NFleet)] - 1; 
   }

char *TuningFleets::GetFleetName(short iFleet)  
   { 
   return ArSmry.Fleet[iFleet]; 
   } 

void TuningFleets::IncrFleet(short * pFleet) 
   {
   do
      (*pFleet)++;
   while (*pFleet <= ArSmry.NFleet && TuningFleet[*pFleet] == false);
   }    

bool TuningFleets::Tuning(short iFleet) 
   {
   if (iFleet >= 1 && iFleet <= ArSmry.NFleet)
      return (TuningFleet[iFleet]==1);
   else
      return false;
   }    

struct FleetArrayRanges *TuningFleets::GetArSmry(void)  
   { 
   return &ArSmry; 
   } 

double TuningFleets::GetStartFishing(short Fleet)
   {
   return alpha[max(min(ArSmry.NFleet,Fleet),1)];
   }

double TuningFleets::GetEndFishing(short Fleet)
   {            
   return beta[max(min(ArSmry.NFleet,Fleet),1)];
   }  
   
double TuningFleets::GetCPUE(short Fleet, short Age, short Year)  
   {                
   if (Fleet < 1                     || Fleet > ArSmry.NFleet         ||
       Age   < ArSmry.MinAge[Fleet]  || Age   > ArSmry.MaxAge[Fleet]  ||
       Year  < ArSmry.MinYear[Fleet] || Year  > ArSmry.MaxYear[Fleet] ||
       Effort[Fleet][Year] <= 0.0)
      return 0.0f;
   else
      return Catch[Fleet][Age][Year]/Effort[Fleet][Year];
   }

bool TuningFleets::AllocArrays()    
   {    
   short iFleet;
                                                            
   if (ArraysSet)
      FreeAllocArrays();  
   
   if (!flallocArray(&alpha,       1, ArSmry.NFleet))
      return FALSE;
   if (!flallocArray(&beta,        1, ArSmry.NFleet))
      return FALSE;
   if (!flallocArray(&Season,      1, ArSmry.NFleet))
      return FALSE;
   if (!flallocArray(&TuningFleet, 1, ArSmry.NFleet))
      return FALSE;
   for (iFleet=1; iFleet<= ArSmry.NFleet; iFleet++)
      {
      alpha[iFleet]       = 0.0f; 
      beta[iFleet]        = 1.0f;
      Season[iFleet]      = 1;
      TuningFleet[iFleet] = true;
      }
   
   if (!flallocArray(&(ArSmry.MinAge),    1, ArSmry.NFleet))
      return FALSE;
   if (!flallocArray(&(ArSmry.MaxAge),    1, ArSmry.NFleet))
      return FALSE;
   if (!flallocArray(&(ArSmry.MinYear),   1, ArSmry.NFleet))
      return FALSE;
   if (!flallocArray(&(ArSmry.MaxYear),   1, ArSmry.NFleet))   
      return FALSE;
   for (iFleet=1; iFleet<= ArSmry.NFleet; iFleet++)
      {
      ArSmry.MinAge[iFleet]    =  
      ArSmry.MinYear[iFleet]   =  
      ArSmry.MaxAge[iFleet]    =    
      ArSmry.MaxYear[iFleet]   = 0;
      }
   
   if (!flallocArray(&Catch,  1, ArSmry.NFleet))
      return FALSE;
   if (!flallocArray(&Effort, 1, ArSmry.NFleet))
      return FALSE;

   if (!flallocArray(&ArSmry.Fleet, 1, ArSmry.NFleet)) 
      return FALSE;
   for (iFleet=1; iFleet<= ArSmry.NFleet; iFleet++)
      {
      if (!flallocArray(&ArSmry.Fleet[iFleet], 0, FLEETSTRLEN-1)) 
         return FALSE;
      for (short iChar=0; iChar < FLEETSTRLEN; iChar++)
         ArSmry.Fleet[iFleet][iChar] = 0;
      }

   ArraysSet = TRUE;  

   return TRUE;
   } 
     
void TuningFleets::FreeAllocArrays(void)
   {
   if (!ArraysSet)
      return;   
//return;

   for (short iFleet=1; iFleet<= ArSmry.NFleet; iFleet++)
      {
	   if (ArSmry.MinAge[iFleet] != 0 && ArSmry.MaxAge[iFleet] != 0 && ArSmry.MinYear[iFleet] != 0 && ArSmry.MaxYear[iFleet] != 0)
		{   	
		free_flallocArray(Catch[iFleet],  ArSmry.MinAge[iFleet],  ArSmry.MaxAge[iFleet], ArSmry.MinYear[iFleet], ArSmry.MaxYear[iFleet]);
		free_flallocArray(Effort[iFleet], ArSmry.MinYear[iFleet], ArSmry.MaxYear[iFleet]);
		}
	  }                                  
   
   if (ArSmry.NFleet>0)
      {
      free_flallocArray(Catch,          1, ArSmry.NFleet);
      free_flallocArray(Effort,         1, ArSmry.NFleet);
      free_flallocArray(ArSmry.MinAge,  1, ArSmry.NFleet);
      free_flallocArray(ArSmry.MaxAge,  1, ArSmry.NFleet);
      free_flallocArray(ArSmry.MinYear, 1, ArSmry.NFleet);
      free_flallocArray(ArSmry.MaxYear, 1, ArSmry.NFleet);
      free_flallocArray(alpha,          1, ArSmry.NFleet);
      free_flallocArray(beta,           1, ArSmry.NFleet);
      free_flallocArray(Season,         1, ArSmry.NFleet);
      free_flallocArray(TuningFleet,    1, ArSmry.NFleet);

      for (short iFleet=1; iFleet<= ArSmry.NFleet; iFleet++)
         free_flallocArray(ArSmry.Fleet[iFleet], 0, FLEETSTRLEN-1);
      free_flallocArray(ArSmry.Fleet, 1, ArSmry.NFleet); 
      }                                  
   }
               
TuningFleets::~TuningFleets(void)
   {
   FreeAllocArrays();
   }
    
