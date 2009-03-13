#ifndef _INC_fl__tun
#define _INC_fl__tun

#include "fl__types.hpp"
#include "flalloc.hpp"
                                    
#define MAXNFLEETS  16
#define FLEETSTRLEN 32       

#define CONST_TUNING_BIOMASS         1
#define CONST_TUNING_NUMBERS         2

#define CONST_TUNING_INDEXSEL_CPUE   1
#define CONST_TUNING_INDEXSEL_KNOWN  2
#define CONST_TUNING_INDEXSEL_CATCH  3

#define TuningDataErrInput  1
#define TuningDataErrInput2 2
#define TuningDataErrInput3 3 

struct FleetArrayRanges
   {
   short NFleet;
   
   LP2CHAR Fleet;
    
   LPSHORT MinAge,
           MaxAge,
           MinYear,
           MaxYear; 
   };

class TuningFleets
{  
public:  
   short Error; 
   bool  ArraysSet;

   TuningFleets(void);
   
   bool    SetCatch(short,short,short,double);
   
   inline void SetMinTuningYear(short _Year) {MinTuningYear = max(MinTuningYear,min(MaxTuningYear,MaxTuningYear - _Year+1));};
   
   short   GetMinAge(short);    
   short   GetMaxAge(short);    
   short   GetMinYear(short); 
   short   GetMaxYear(short);  
   double  GetCatch(short,short,short);
   double  GetEffort(short,short);  
   short   GetNFleet(void);  
   short   AddSeason(short);
   char   *GetFleetName(short);  
   struct  FleetArrayRanges *GetArSmry(void);  
   double  GetStartFishing(short);
   double  GetEndFishing(short);
   
   bool    Tuning(short);
   void    IncrFleet(short *); 

   double  GetCPUE(short,short,short);               

    ~TuningFleets();                                                       
                      
protected:                            
    short Units,
          SelType,
          MinTuningAge,
          MaxTuningAge,
          MinTuningYear,
          MaxTuningYear;
          
    LPSHORT Season,
            TuningFleet;        
        
    LPDOUBLE  alpha, beta;
    
    LP2DOUBLE Effort;
                                     
    LP3DOUBLE Catch ;                  
    
    struct FleetArrayRanges ArSmry;
     
    bool AllocArrays(void);
    void FreeAllocArrays(void);
};

#endif /* _INC_fl__tun */
