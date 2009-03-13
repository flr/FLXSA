#ifndef _INC_fl__xsa
#define _INC_fl__xsa

#include <math.h>
#include "fl__types.hpp"
#include "fl__ctc.hpp"
#include "fl__tun.hpp"
#include "fl__vpa.hpp"

#define NSE_WTS_AUTOMATIC 0
#define NSE_WTS_MANUAL    1
#define NSE_WTS_COMBINE   2


typedef enum tagConstXSAReturn 
	{
   FLcXSAReturnResidual    = 1,
   FLcXSAReturnNHat        = 2
   } FLConstXSAReturn;


//FLTypeXSAControl
struct FLTypeXSAControl
   {

   double MSSTol,
          MinNSE,
          FShkSE,
          FPreSpwn,
          MPreSpwn;      

   short TSRange,
         TSPower,
         LRcrtAge,
         ConQAge,
         Shk2N,
         Shk2F,
         Shk2FYr,         
         Shk2FAge,
         MaxIters,
         PlusGroup,          
         VPA,
         SmoothQ,
         Shk2FTarget,
         NSEWts,
         TuningWindow;
};


class ExtendedSurvivorsAnalysis : public VirtualPopAnalysis
{  
public:   
   bool OutputWts,
        FlagResetRecruitQ;


   TuningFleets *pTune;

   ExtendedSurvivorsAnalysis(void);

   bool SetControl(FLTypeXSAControl  *);

   bool  Run(void);
   bool  RunForInputParam(void);
   void  RunRetro(short,LP2DOUBLE,LP2SHORT,short,LP2SHORT,short,LP2DOUBLE,LP2DOUBLE,LP3DOUBLE,LP3DOUBLE);

   double GetNHat(short,short,short);  
   double GetNSE(short,short,short);  
   double GetManualNSE(short,short,short);  

   inline double GetResidual(short iFleet,short iAge ,short iYear) { return GetNHat(iFleet,iAge,iYear)/GetN(iAge,iYear); }; 
   inline double GetTermN(short iAge) {return GetN(iAge-1,MaxYear)*exp(-GetF(iAge-1,MaxYear)-GetM(iAge-1,MaxYear)); };
   
   short GetTSRange(void);
   short GetTSPower(void);
   short GetLRcrtAge(void);
   short GetConQAge(void);
   short GetShk2N(void);
   short GetShk2F(void);
   short GetShk2FYr(void);
   short GetShk2FAge(void);
   short GetMaxIters(void);
   short GetPlusGroup(void);
          
   double GetMSSTol(void);
   double GetMinNSE(void);
   double GetFShkSE(void);
   double GetFPreSpwn(void);
   double GetMPreSpwn(void);   

   short  GetNTuneFleet(void);    
   short  GetMinTuneAge(short);    
   short  GetMaxTuneAge(short);    
   short  GetMinTuneYear(short); 
   short  GetMaxTuneYear(short);

   double   AvCatch(short,short,double,double);
   void     CalcMeanF(short, short, short);
   double   GetMeanF(short); 
   short    GetMaxCatchAge(void);
   double   GetTermPopSE(short);  
   LPDOUBLE GetRecruits(void);                                
   double   Getq(short,short);
   double   Getq_b(short,short);
   short    Getqn(short,short);
   double   GetqSE(short,short);
   double   Getq(short,short,short);
   double   Getq_b(short,short,short);
   short    Getqn(short,short,short);
   double   GetqSE(short,short,short);
   double   GetCorrectedCPUE(short,short,short);  
   double   GetPopShrinkMean(short);
   void     SetPopShrinkMean(short,double);
   double   GetPopShrinkSE(short);
   void     SetPopShrinkSE(short,double);
   double   GetFShrinkMean(short); 
   double   GetCPUEHat(short, short, short);
   void     SetFShrinkMean(short, double); 
   void     SetNHat(short,short,short,double);
   void     SetNSE(short,short,short,double);
   void     SetN(short,short,double);
   void     Setq(short,short,double);
   void     Setq_b(short,short,double);
   void     Setqn(short,short,short);
   void     SetqSE(short,short,double);
   void     Setq(short,short,short,double);
   void     Setq_b(short,short,short,double);
   void     Setqn(short,short,short,short);
   void     SetqSE(short,short,short,double);
   void     GetPopShrinkSE(short,double);

   inline bool  UseManualNSEWts(void) { return (Control.NSEWts==1); };

   void     Project(double,double,bool,double,double,double,short,short,short,LP2DOUBLE,LP2DOUBLE,LP2DOUBLE,LPDOUBLE,LPDOUBLE,LPDOUBLE,LPDOUBLE,LPDOUBLE,LPDOUBLE,bool RCT3Flag=false);
   bool     CheckRange(short, short, short, short);

   ~ExtendedSurvivorsAnalysis();
protected:    
     short Iters,
           NTuneFleet,
           DegreesOfFreedom,
           NObs,
           NParam;
          
   bool   TestFlag;

   double MSS,
          LogLikelihood,
          CoeffDet,
          AIC;

   LPSHORT MinTuneAge,   MaxTuneAge,  
           MinTuneYear,  MaxTuneYear;
         
   struct XSAControl Control;      

   LPDOUBLE TermPopSE,
            VarNShrink,
            PopShrinkMean,
            PopShrinkSE,
            FShrinkMean,
            MeanF,
            GeoMeanN,
            GeoMeanNSE,
            TSWeight,
            LogLikelihoodRecord;

   LP2DOUBLE PrevF,
             InternalSE,
             ExternalSE,
             ExtraVariable;
                 
   LP3SHORT  qn;
                 
   LP3DOUBLE q,
             qSE,
             q_b;

   LP3DOUBLE NHat,
             NSE,
             ManualNSE,
             CorrectedCPUE;

   void   SetDefaultControls(void);
   bool   SetControls(LPDOUBLE,short);
   bool   SetControls(LP2DOUBLE,short,short);
   bool   SetControls(char *, char *);
   bool   SetCorrectedCPUE(void);
   bool   CalcRecruitQ(short Iters=6);
   bool   CalcQInsteadOfCalcRecruitQ(void);
   bool   XSARecruitsSE(void);
   bool   iFleetQ(void);
   bool   CalcQ(void);
   bool   XSANHatQSE(void);
   bool   Initialise(void);
   bool   RunVPAFromXSA(void);
   bool   XSATermPop(void);
   bool   CalcTermPops(void);
   double CalcTermPopByFleet(short,short);
   void   CalcTermPops(LP3DOUBLE,LP3DOUBLE,LPDOUBLE,LPDOUBLE,LPDOUBLE,LPDOUBLE);
   void   CalcTermPopsByMainEffects(short,LP3DOUBLE,LP3DOUBLE,LP2DOUBLE,LP2DOUBLE,LP2DOUBLE,LP2DOUBLE,LP2DOUBLE,LP2DOUBLE);
   void   ScaleTermPops(LP3DOUBLE,LPDOUBLE,LPDOUBLE,LPDOUBLE);
   double CalcTermPopNoShrink(short);

   bool   CalcTermPopSE(void);
   bool   OutputTermPopWts(void);
   bool   _OutputTermPopWts(void);
   bool   Shrink2PopMean();
   bool   Shrink2FMean(void);
   bool   FBarShrink(void);
   void   SetTau(void);
   double Tau(short);
   double Tau(short, short);
   bool   Theta(short,short,short);
   double Weight(short,short,short,double);
   double Weight(short,short,short,short,double);
   void   CalcPlusGroupNF(void);
   void   CalcInternalExternalSE(void);  
   void   CalcWeightedLogMeanN(short,double *,double *);

   double GetLogLikelihood(void);
   double GetRSquared(void);

   void   ResetRecruitQ(void);

   bool  AllocArrays(void);
   void  FreeAllocArrays(void); 
};
                                                       
#endif /* _INC_fl__xsa */
