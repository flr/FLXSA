#include "fl__xsa.hpp" 
#include <float.h>
#include <R_ext/Arith.h>

bool GM = FALSE;



void nrerror(char error_text[])
   { 
   short i;
   i = 0;
   }

double gammln(double xx)
   {      
   double x,tmp,ser;
   static double cof[6]={76.18009173,-86.50532033,24.01409822,
      -1.231739516,0.120858003e-2,-0.536382e-5};
   short j;

   x=xx-1.0f;
   tmp=x+5.5;
   tmp -= (x+0.5f)*log(tmp);
   ser=1.0f;
   for (j=0;j<=5;j++) {
      x += 1.0f;
      ser += cof[j]/x;
   }
   return (double) (-tmp+log(2.50662827465*ser));
}

#define ITMAX 100
#define EPS 3.0e-7

void gser(LPDOUBLE gamser,double a,double x,LPDOUBLE gln)
{
   short n;
   double sum,del,ap;

   *gln=gammln(a);
   if (x <= 0.0) {
      if (x < 0.0) return; //nrerror("x less than 0 in routine GSER");
      *gamser=0.0;
      return;
   } else {
      ap=a;
      del=sum=1.0f/a;
      for (n=1;n<=ITMAX;n++) {
         ap += 1.0f;
         del *= x/ap;
         sum += del;
         if (fabs(del) < fabs(sum)*EPS) {
            *gamser=sum*(double)exp(-x+a*log(x)-(*gln));
            return;
         }
      }
      //nrerror("a too large, ITMAX too small in routine GSER");
      return;
   }
}

void gcf(LPDOUBLE gammcf, double a, double x,LPDOUBLE gln)
{
   short n;
   double gold=0.0,g,fac=1.0f,b1=1.0f;
   double b0=0.0,anf,ana,an,a1,a0=1.0f;

   *gln=gammln(a);
   a1=x;
   for (n=1;n<=ITMAX;n++) {
      an=(double) n;
      ana=an-a;
      a0=(a1+a0*ana)*fac;
      b0=(b1+b0*ana)*fac;
      anf=an*fac;
      a1=x*a0+anf*a1;
      b1=x*b0+anf*b1;
      if (a1) {
         fac=1.0f/a1;
         g=b1*fac;
         if (fabs((g-gold)/g) < EPS) {
            *gammcf=(double)exp(-x+a*(double)log(x)-(*gln))*g;
            return;
         }
         gold=g;
      }
   }
   //nrerror("a too large, ITMAX too small in routine GCF");
}

#undef ITMAX
#undef EPS


double gammq(double a, double x)
{
   double gamser,gammcf,gln;

   if (x < 0.0 || a <= 0.0) return 0.0; //nrerror("Invalid arguments in routine GAMMQ");
   if (x < (a+1.0)) {
      gser(&gamser,a,x,&gln);
      return 1.0f-gamser;
   } else {
      gcf(&gammcf,a,x,&gln);
      return gammcf;
   }
}

void fit(double * x, double *  y, short ndata, short mwt, double *  sig, 
         double * a,  double * b, double * siga, double * sigb, 
         double * chi2, double * q)
{
   short i;
   double wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;
            
   *b=0.0;
   if (mwt) {
      ss=0.0;
      for (i=1;i<=ndata;i++) {
         wt= 1.0f/((sig[i])*(sig[i]));
         ss += wt;
         sx += (x[i]*wt);
         sy += (y[i]*wt);
      }
   } else {
      for (i=1;i<=ndata;i++) {
         sx += x[i];
         sy += y[i];
      }
      ss=ndata;
   }
   sxoss=(sx/ss);
   if (mwt) {
      for (i=1;i<=ndata;i++) {
         t=((x[i]-sxoss)/sig[i]);
         st2 += (t*t);
         *b += (t*y[i]/sig[i]);
      }
   } else {
      for (i=1;i<=ndata;i++) {
         t=(x[i]-sxoss);
         st2 += (t*t);
         *b += (t*y[i]);
      }
   }
   *b /= st2;
   *a=((sy-sx*(*b))/ss);
   *siga=sqrt((1.0f+sx*sx/(ss*st2))/ss);
   *sigb=sqrt(1.0f/st2);
   *chi2=0.0;
   if (mwt == 0) {
      for (i=1;i<=ndata;i++)
         *chi2 += (y[i]-(*a)-(*b)*x[i])*(y[i]-(*a)-(*b)*x[i]);
      *q=1.0f;
      sigdat=sqrt((*chi2)/(ndata-2));
      *siga *= sigdat;
      *sigb *= sigdat;
   } else {
      for (i=1;i<=ndata;i++)
         *chi2 += ((y[i]-(*a)-(*b)*x[i])/sig[i])*((y[i]-(*a)-(*b)*x[i])/sig[i]);
      *q=gammq(0.5f*(ndata-2),0.5f*(*chi2));
   }
}
/* (C) Copr. 1986-92 Numerical Recipes Software Y5j.B+J,. */                       

ExtendedSurvivorsAnalysis::ExtendedSurvivorsAnalysis(void)
   {
   SetDefaultControls();
   
   TestFlag = false;

   OutputWts = false;
   };


void ExtendedSurvivorsAnalysis::ResetRecruitQ()
   {
   if (!FlagResetRecruitQ)
      return;

   short iAge, iFleet;

   for (iFleet = 0, pTune->IncrFleet(&iFleet); iFleet <= NTuneFleet; pTune->IncrFleet(&iFleet))                                                       
      for (iAge = min(MaxTuneAge[iFleet],Control.LRcrtAge); iAge>=MinTuneAge[iFleet]; iAge--)
         if (Getq(iFleet,iAge) > fabs(100.0))
            Control.LRcrtAge = iAge-1;               
   }

void ExtendedSurvivorsAnalysis::SetDefaultControls(void)
   {
   Control.VPA             = FALSE;
   Control.SmoothQ         = FALSE;
// Control.TSRangeShk2N    = 
   Control.TSRange         = 20;
   Control.TSPower         = 3;
   Control.LRcrtAge        = 3;
   Control.ConQAge         = 6;
   Control.FShkSE          = 0.5f;
   Control.Shk2F           = TRUE;
   Control.Shk2N           = TRUE;
   Control.Shk2FYr         = 5; 
   Control.Shk2FAge        = 5;
   Control.MinNSE          = 0.3f;
   Control.MSSTol          = 1.0E-8f;                    
   Control.MaxIters        = 30;    
   Control.PlusGroup         = 50;   
   Control.Shk2FTarget     = 0;
   Control.NSEWts          = NSE_WTS_AUTOMATIC;
   Control.ReturnType      = FLcXSAReturnResidual;
   Control.TuningWindow    = 1000;
   Control.FPreSpwn        =
   Control.MPreSpwn        = 1.0;

   FlagResetRecruitQ = false;
   }


bool ExtendedSurvivorsAnalysis::Run(void)
    {        
   //Initialise Max row and column of pop matrix
   Initialise(); 
        
   Iters = 0;
   MSS = Control.MSSTol*2.0;   
       
   while (((++(Iters) <= Control.MaxIters) && (MSS > Control.MSSTol)) || Iters < 2)   
      {
      //Run VPA
      if (TestFlag)
         Iters = Control.MaxIters;
      else
         RunVPAFromXSA();
 
      //Corrected CPUE
      SetCorrectedCPUE();
   
      //estimate recruiting iYear classes
      CalcRecruitQ(Iters);
   
      //estimated fully recruited pop size
      CalcQ();

      //estimate terminal pop size (i.e. oldest iAge in each cohort)
      //if converged or last iter then estimate SE of terminal pop size
      XSATermPop();

      //if (TestFlag)
      //   RunVPAFromXSA();
 
      LogLikelihoodRecord[Iters] = GetLogLikelihood();

      //OutputTermPopWts();

      ResetRecruitQ();
      }                           

   Iters--;

   CalcInternalExternalSE();

   CalcPlusGroupNF();  

   return true;
   }
     
bool ExtendedSurvivorsAnalysis::RunForInputParam(void)
   {
   //Initialise Max row and column of pop matrix
   Initialise(); 
               
   //Run VPA
   RunVPAFromXSA();
 
   //Corrected CPUE
   SetCorrectedCPUE();
   
   //estimate recruiting iYear classes
   CalcRecruitQ();
   
   //estimated fully recruited pop size
   CalcQ();

   CalcInternalExternalSE();

   CalcPlusGroupNF();  

   return true;
   } 

void ExtendedSurvivorsAnalysis::CalcInternalExternalSE(void)
   {       
   short iAge, iFleet, iYear, iYearClass;
   
   double ECF, ECZ, Wt;

   LPSHORT  NPts;

   LPDOUBLE IntTop,
            ExtTop,
            Bottom;

   flallocArray(&NPts,   0, NTuneFleet);
   flallocArray(&IntTop, 0, NTuneFleet);
   flallocArray(&ExtTop, 0, NTuneFleet);
   flallocArray(&Bottom, 0, NTuneFleet);
   
   for (iYearClass=MinYear-MaxAge; iYearClass <= MaxYear-MinAge; iYearClass++)
      {                                                            
      ECF   = 
      ECZ   = 1.0;

      for (iFleet=0; iFleet<=NTuneFleet; iFleet++)
         IntTop[iFleet] = 
         ExtTop[iFleet] = 
         Bottom[iFleet] = 
         NPts[iFleet]     = 0;
      
      short  TermAge = min(MaxAge  + 1, MaxYear - iYearClass + 1);
      double LogTermN = log(N[TermAge][iYearClass+TermAge]);

      for (iAge  = min(MaxAge, MaxYear  - iYearClass); iAge >= max(MinAge, MinYear - iYearClass); iAge--)
         {
         iYear = iYearClass + iAge; 

         ECF *= exp(F[iAge][iYear]);               
         ECZ *= exp(F[iAge][iYear] + GetM(iAge,iYear));               
            
         for (iFleet=0, pTune->IncrFleet(&iFleet); iFleet<=NTuneFleet; pTune->IncrFleet(&iFleet))
            if (Theta(iFleet, iAge, iYear))
               {
               Wt              = Weight(iFleet, iAge, iYear, ECF);                                                                          
               
               IntTop[iFleet] += Wt/ECF;
               double ptest    = log(GetNHat(iFleet,iAge,iYear)/ECZ);
               double overmean = log(CalcTermPopByFleet(iFleet, iYearClass)); 

               ExtTop[iFleet] += Wt*SQR(log(GetNHat(iFleet,iAge,iYear)/ECZ) - log(CalcTermPopByFleet(iFleet, iYearClass))); 
               ExtTop[0]      += Wt*SQR(log(GetNHat(iFleet,iAge,iYear)/ECZ) - LogTermN); 
               Bottom[iFleet] += Wt; 
               NPts[iFleet]++;
               }

      
         if (Control.Shk2F == TRUE && GetFShrinkMean(iYearClass) > 0.0 && Control.FShkSE > 0.0) 
            {
            IntTop[0] += 1.0/SQR(Control.FShkSE);

            if (GetFShrinkMean(iYearClass) > 0.0)
               ExtTop[0] += SQR((log(GetFShrinkMean(iYearClass)) - LogTermN)/Control.FShkSE);
            
            Bottom[0] += 1.0/SQR(Control.FShkSE);
            }

         if ((Control.Shk2N == TRUE) && (TermAge < Control.LRcrtAge))             
            {
            IntTop[0] += 1.0/SQR(GetPopShrinkSE(TermAge));                        
         
            //weighted VPA Mean
            ExtTop[0] += SQR((log(GetPopShrinkMean(TermAge)) - LogTermN)/GetPopShrinkSE(TermAge));                        
            Bottom[0] += 1.0/SQR(GetPopShrinkSE(TermAge));
            } 
         }

      for (iFleet=0, pTune->IncrFleet(&iFleet); iFleet<=NTuneFleet; pTune->IncrFleet(&iFleet))
         {
         NPts[0]   += max(0,NPts[iFleet]);
         IntTop[0] += IntTop[iFleet]; 
         Bottom[0] += Bottom[iFleet]; 
         }       
                     
      
      //Overall
      if (Bottom[0] > 0.0)
         {
         InternalSE[0][iYearClass] = sqrt(IntTop[0]) / Bottom[0];

         if (NPts[0] > 0)
            ExternalSE[0][iYearClass] = 1.0/sqrt(NPts[0]-1)*sqrt(ExtTop[0]/Bottom[0]);
         else
            ExternalSE[0][iYearClass] = 0.0;
         }
      else
         InternalSE[0][iYearClass] = 
         ExternalSE[0][iYearClass] = 0.0;
      
      //Fleet
      for (iFleet=0, pTune->IncrFleet(&iFleet); iFleet<=NTuneFleet; pTune->IncrFleet(&iFleet))
         if (Bottom[iFleet] > 0.0 && pTune->Tuning(iFleet))
            {
            InternalSE[iFleet][iYearClass] = sqrt(IntTop[iFleet]) / Bottom[iFleet];

            if (NPts[iFleet] > 1 && iAge < Control.ConQAge)
               ExternalSE[iFleet][iYearClass] = 1.0/sqrt(NPts[iFleet]-1)*sqrt(ExtTop[iFleet]/Bottom[iFleet]);
            else
               ExternalSE[iFleet][iYearClass] = 0.0;
            }
         else
            InternalSE[iFleet][iYearClass] = 
            ExternalSE[iFleet][iYearClass] = 0.0;
      
         }
    
   //OutputWts = false;
                                        
   //if (OutputWts)
   //   OutputTermPopWts();
    
   free_flallocArray(NPts,   0, NTuneFleet);
   free_flallocArray(IntTop, 0, NTuneFleet);
   free_flallocArray(ExtTop, 0, NTuneFleet);
   free_flallocArray(Bottom, 0, NTuneFleet);
   }

double ExtendedSurvivorsAnalysis::CalcTermPopByFleet(short iFleet, short iYearClass)
   {
   short iAge, iYear;

   double ECF, 
          ECZ,
          Wt,
          SumTop, 
          SumBottom;

   ECF    = ECZ       = 1.0;
   SumTop = SumBottom = 0.0;
   
   for (iAge = min(min(GetPlusGroup()-1,MaxAge), MaxYear - iYearClass); iAge >= max(GetMinTuneAge(iFleet), GetMinTuneYear(iFleet) - iYearClass); iAge--)
      {
      iYear = iYearClass + iAge; 
      ECF *= exp(F[iAge][iYear]);               
      ECZ *= exp(F[iAge][iYear] + GetM(iAge, iYear));
      if (Theta(iFleet, iAge, iYear) && Getqn(iFleet, iAge) > 1)
         {
         if (Control.SmoothQ)
            Wt = Weight(iFleet, iAge, iYear, iYearClass+MaxAge, ECF);                                                                         
         else
            Wt = Weight(iFleet, iAge, iYear, ECF);                                                                          
               
            SumTop    += Wt*(log(GetNHat(iFleet,iAge,iYear) ) - log(ECZ)); 
            SumBottom += Wt;
            }
      }       
   
   double ReturnVal;
   
   if (SumBottom > 0)
      ReturnVal = exp(SumTop/SumBottom);
   else
      ReturnVal = 0.0;

   return ReturnVal;
   }

void ExtendedSurvivorsAnalysis::CalcTermPops(LP3DOUBLE NHat, LP3DOUBLE Wt, LPDOUBLE NHatShrinkN, LPDOUBLE WtShrinkN, LPDOUBLE NHatShrinkF, LPDOUBLE WtShrinkF)
   { 
   double ECF, ECZ;

   short iAge, iFleet, iYear, iYearClass, RecruitAge;

   for (iYearClass=MaxYear-MaxAge; iYearClass <= MaxYear-MinAge; iYearClass++)
      NHatShrinkN[iYearClass] =
      WtShrinkN[iYearClass]   = 0.0;
      
      
   for (iYearClass=MinYear-MaxAge; iYearClass <= MaxYear-MinAge; iYearClass++)
      { 
      NHatShrinkF[iYearClass] =
      WtShrinkF[iYearClass]   = 0.0;

      ECF                     =                
      ECZ                     = 1.0;

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
                     Wt[iFleet][iAge][iYear] = Weight(iFleet, iAge, iYear, iYearClass+MaxAge, ECF);                                                                          
                  else
                     Wt[iFleet][iAge][iYear] = Weight(iFleet, iAge, iYear, ECF)*ECZ;                                                                          
         
                  NHat[iFleet][iAge][iYear] = log(GetNHat(iFleet,iAge,iYear)) - log(ECZ)*Weight(iFleet, iAge, iYear, ECF); 
                  }
            }       
         }

      if (Control.Shk2F == TRUE && GetFShrinkMean(iYearClass) > 0 && Control.FShkSE > 0.0)
         {
         NHatShrinkF[iYearClass] = log(GetFShrinkMean(iYearClass));
         WtShrinkF[iYearClass]   = 1/(SQR(Control.FShkSE));
         }

      iAge = min(MaxAge, MaxYear - iYearClass);
      RecruitAge = MaxYear - iYearClass;                                                 

      if (Control.Shk2N == TRUE && RecruitAge <= Control.LRcrtAge && GetPopShrinkMean(RecruitAge) > 0.0 && PopShrinkSE[RecruitAge] > 0.0)             
        {
        NHatShrinkN[iYearClass] = log(GetPopShrinkMean(RecruitAge));                       
        WtShrinkN[iYearClass]   = 1/(SQR(GetPopShrinkSE(RecruitAge)));
        }
      } 
   }
 
void ExtendedSurvivorsAnalysis::ScaleTermPops(LP3DOUBLE Wt, LPDOUBLE WtShrinkN, LPDOUBLE WtShrinkF, LPDOUBLE SumWts)
   { 
   short iAge, iFleet, iYear, iYearClass, RecruitAge;


   for (iYearClass=MinYear-MaxAge; iYearClass <= MaxYear-MinAge; iYearClass++)
      { 
      SumWts[iYearClass] = 0.0;
   
      for (iAge = min(MaxAge, MaxYear - iYearClass); iAge >= max(MinAge, MinYear - iYearClass); iAge--)
         {
         iYear = iYearClass + iAge; 
         for (iFleet=0, pTune->IncrFleet(&iFleet); iFleet<=NTuneFleet; pTune->IncrFleet(&iFleet))
            {
            if (GetMaxTuneAge(iFleet)  >= iAge  && GetMinTuneAge(iFleet)  <= iAge &&
                GetMaxTuneYear(iFleet) >= iYear && GetMinTuneYear(iFleet) <= iYear  )
               SumWts[iYearClass] += Wt[iFleet][iAge][iYear];
            }       
         }

      if (Control.Shk2F == TRUE && GetFShrinkMean(iYearClass) > 0 && Control.FShkSE > 0.0)
         SumWts[iYearClass]  += WtShrinkF[iYearClass];

      iAge = min(MaxAge, MaxYear - iYearClass);
      RecruitAge = MaxYear - iYearClass;                                                 

      if (Control.Shk2N == TRUE && RecruitAge <= Control.LRcrtAge && GetPopShrinkMean(RecruitAge) > 0.0 && PopShrinkSE[RecruitAge] > 0.0)             
        SumWts[iYearClass]  += WtShrinkN[iYearClass];
      } 
   }

void ExtendedSurvivorsAnalysis::CalcTermPopsByMainEffects(short MaxTrueAge, LP3DOUBLE NHat, LP3DOUBLE Wt, LP2DOUBLE NHatByFleet, LP2DOUBLE WtByFleet, LP2DOUBLE NHatByAge, LP2DOUBLE WtByAge, LP2DOUBLE NHatByYear, LP2DOUBLE WtByYear)
   { 
   short iAge, iFleet, iYear, iYearClass;

   for (iYearClass = MinYear-MaxTrueAge; iYearClass <= MaxYear-MinAge; iYearClass++)
      {
      for (iAge = MinAge; iAge <= MaxTrueAge; iAge++)
         NHatByAge[iYearClass][iAge] = 
         WtByAge[iYearClass][iAge]   = 0.0;

      for (iFleet=0, pTune->IncrFleet(&iFleet); iFleet<=NTuneFleet; pTune->IncrFleet(&iFleet))
         NHatByFleet[iYearClass][iFleet] =
         WtByFleet[iYearClass][iFleet]   = 0.0;
      
      for (iYear = MinYear; iYear <= MaxYear; iYear++)
         NHatByYear[iYearClass][iYear] = 
         WtByYear[iYearClass][iYear]   = 0.0;
   
      for (iAge = min(MaxTrueAge, MaxYear - iYearClass); iAge >= max(MinAge, MinYear - iYearClass); iAge--)
         {
         iYear = iYearClass + iAge;
      
         for (iFleet=0, pTune->IncrFleet(&iFleet); iFleet<=NTuneFleet; pTune->IncrFleet(&iFleet))
            {
            if (GetMaxTuneAge(iFleet)  >= iAge  && GetMinTuneAge(iFleet)  <= iAge &&
                GetMaxTuneYear(iFleet) >= iYear && GetMinTuneYear(iFleet) <= iYear  )
               {
               NHatByFleet[iYearClass][iFleet] += NHat[iFleet][iAge][iYear]*Wt[iFleet][iAge][iYear];
               NHatByAge[iYearClass][iAge]     += NHat[iFleet][iAge][iYear]*Wt[iFleet][iAge][iYear];
               NHatByYear[iYearClass][iYear]   += NHat[iFleet][iAge][iYear]*Wt[iFleet][iAge][iYear];
                  
               WtByFleet[iYearClass][iFleet]   += Wt[iFleet][iAge][iYear];
               WtByAge[iYearClass][iAge]       += Wt[iFleet][iAge][iYear];
               WtByYear[iYearClass][iYear]     += Wt[iFleet][iAge][iYear];
               }
            }       
         }
       
      for (iAge = min(MaxTrueAge, MaxYear - iYearClass); iAge >= max(MinAge, MinYear - iYearClass); iAge--)
         if (WtByAge[iYearClass][iAge] > 0.0) 
            NHatByAge[iYearClass][iAge]     /=  WtByAge[iYearClass][iAge];

      for (iFleet=0, pTune->IncrFleet(&iFleet); iFleet<=NTuneFleet; pTune->IncrFleet(&iFleet))
         if (WtByFleet[iYearClass][iFleet] > 0.0) 
            NHatByFleet[iYearClass][iFleet] /=  WtByFleet[iYearClass][iFleet];

      for (iYear = MinYear; iYear <= MaxYear; iYear++)
         if (WtByYear[iYearClass][iYear] > 0.0) 
            NHatByYear[iYearClass][iYear]   /=  WtByYear[iYearClass][iYear];
      }
   
   }

double ExtendedSurvivorsAnalysis::CalcTermPopNoShrink(short iYearClass)
   {
   short iAge, iYear, iFleet;

   double ECF, 
          ECZ,
          Wt,
          SumTop, 
          SumBottom;

   ECF    = ECZ       = 1.0;
   SumTop = SumBottom = 0.0;

   for (iFleet=0, pTune->IncrFleet(&iFleet); iFleet<=NTuneFleet; pTune->IncrFleet(&iFleet))
      for (iAge = min(min(GetPlusGroup()-1,MaxAge), MaxYear - iYearClass); iAge >= max(GetMinTuneAge(iFleet), GetMinTuneYear(iFleet) - iYearClass); iAge--)
         {
         iYear = iYearClass + iAge; 
         ECF *= exp(F[iAge][iYear]);               
         ECZ *= exp(F[iAge][iYear] + GetM(iAge, iYear));
   
         if (Theta(iFleet, iAge, iYear) && Getqn(iFleet, iAge) > 1)
            {
            if (Control.SmoothQ)
               Wt = Weight(iFleet, iAge, iYear, iYearClass+MaxAge, ECF);                                                                         
            else
               Wt = Weight(iFleet, iAge, iYear, ECF);                                                                          
               
            SumTop    += Wt*(log(GetNHat(iFleet,iAge,iYear) ) - log(ECZ)); 
            SumBottom += Wt;
            }
         }       
   
   double ReturnVal;
   
   if (SumBottom > 0.0)
      ReturnVal = exp(SumTop/SumBottom);
   else
      ReturnVal = 0.0;

   return ReturnVal;
   }

double ExtendedSurvivorsAnalysis::AvCatch(short iAge, short iYear, double StartFishing, double EndFishing)
   {
   //Equation (2) VPA User Guide
   double Z;
                                          
   Z = (GetM(iAge, iYear)+F[iAge][iYear]);                   
   
   if (StartFishing == 0.0 && EndFishing == 0.0)
      return 1.0;
   else if (StartFishing == EndFishing)
      EndFishing += 0.00001;

   return (exp(-StartFishing*Z) - exp(-EndFishing*Z))/((EndFishing-StartFishing)*Z); 
   }

bool ExtendedSurvivorsAnalysis::SetCorrectedCPUE(void)
   { 
   //Equation (3) VPA User Guide
   short iFleet;

   for (iFleet = 0, pTune->IncrFleet(&iFleet); iFleet <= NTuneFleet; pTune->IncrFleet(&iFleet))
      for (short iAge = MinTuneAge[iFleet]; iAge <= MaxTuneAge[iFleet]; iAge++)
         for (short iYear = MinTuneYear[iFleet]; iYear <= MaxTuneYear[iFleet]; iYear++)
            CorrectedCPUE[iFleet][iAge][iYear] = 
                  pTune->GetCPUE(iFleet, iAge, iYear)
                 /AvCatch(iAge,iYear,pTune->GetStartFishing(iFleet), pTune->GetEndFishing(iFleet));

   return true;
   }   

//Regression to estimate recuiting pops 
bool ExtendedSurvivorsAnalysis::CalcRecruitQ(short Iters)
   {                                                                      
   short   iAge, iYear, iFleet, ndata;                             
   LPDOUBLE _U, _N, sig;
   double a,  b,  siga, sigb, chi2, qval;  

   if (Iters >= 6)
      {
      for (iFleet = 0, pTune->IncrFleet(&iFleet); iFleet <= NTuneFleet; pTune->IncrFleet(&iFleet))                                                       
         {                                                                  
         flallocArray(&_U,  1, GetMaxTuneYear(iFleet) - GetMinTuneYear(iFleet) + 1);
         flallocArray(&_N,  1, GetMaxTuneYear(iFleet) - GetMinTuneYear(iFleet) + 1);
         flallocArray(&sig,1, GetMaxTuneYear(iFleet) - GetMinTuneYear(iFleet) + 1);
     
         for (iAge = MinTuneAge[iFleet]; iAge <= min(MaxTuneAge[iFleet],Control.LRcrtAge); iAge++)
            {
            for (iYear = MinTuneYear[iFleet], ndata=0; iYear <= MaxTuneYear[iFleet]; iYear++)                                                       
               {                                           
               if (Theta(iFleet, iAge, iYear) && Tau(iYear) > 0.0)
                  {
                  ++ndata;
                  _U[ndata]   = log(GetCorrectedCPUE(iFleet,iAge,iYear));
                  _N[ndata]   = log(N[iAge][iYear]);
                  
                  sig[ndata] = Theta(iFleet, iAge, iYear)*Tau(iYear);
                  sig[ndata] = 1/sqrt(sig[ndata]);                    
                  }
                }  
            Setqn(iFleet,iAge, ndata);
            if (ndata >= 5)
               {
               //Equation (8) VPA User Guide
               fit(_N, _U, ndata, 1, sig, &a, &b, &siga, &sigb, &chi2, &qval);
               b = 1/b;
               a *= -b; 
               
               Setq(  iFleet,iAge, a);
               Setq_b(iFleet,iAge, b);
               SetqSE(iFleet,iAge, sqrt(chi2/(ndata-2)));
               //SetqSE(iFleet,iAge, pow(1.0/exp(sqrt(chi2/(ndata-2))), 1/b));
               
               //Equation (9) VPA User Guide
               //fleet estimates 
               for (iYear = MinTuneYear[iFleet]; iYear <= MaxTuneYear[iFleet]; iYear++)                                                       
                  if (Theta(iFleet, iAge, iYear))
                     SetNHat(iFleet,iAge,iYear, exp(a + b*log(GetCorrectedCPUE(iFleet,iAge,iYear))));
                  else
                     SetNHat(iFleet,iAge,iYear, 0.0);
               }
            }
                      
         free_flallocArray(_U,   1, GetMaxTuneYear(iFleet) - GetMinTuneYear(iFleet) + 1);
         free_flallocArray(_N,   1, GetMaxTuneYear(iFleet) - GetMinTuneYear(iFleet) + 1);
         free_flallocArray(sig,  1, GetMaxTuneYear(iFleet) - GetMinTuneYear(iFleet) + 1);
         }
     
     //Equation (10) VPA User Guide
     XSARecruitsSE();
     }
   else 
      CalcQInsteadOfCalcRecruitQ();
         
   return true;          
   }                     
   
bool ExtendedSurvivorsAnalysis::CalcQInsteadOfCalcRecruitQ(void)
   {   
   short   iAge, iYear, iFleet;                             
   double  Weight, SumWt, SSq;
   double  TS, Top, Bottom;       
   
   //calculate estimate of q       
   for (iFleet = 0, pTune->IncrFleet(&iFleet); iFleet <= NTuneFleet; pTune->IncrFleet(&iFleet))
      for (iAge = MinTuneAge[iFleet]; iAge <= min(MaxTuneAge[iFleet],Control.LRcrtAge); iAge++)
         {  
         Top = Bottom = 0.0;
         Setqn(iFleet,iAge, 0);
             
         for (iYear = MinTuneYear[iFleet]; iYear <= MaxTuneYear[iFleet]; iYear++)
            if (Theta(iFleet, iAge, iYear))
               {
               Setqn(iFleet,iAge, Getqn(iFleet,iAge)+1);
         
               TS      = Tau(iYear); 
               Bottom += TS;
               Top    += (TS*(log(N[iAge][iYear]) - log(GetCorrectedCPUE(iFleet,iAge,iYear))));
               }
         if (Bottom > 0)
            Setq(iFleet,iAge,(1/exp(Top/Bottom)));
         else
            Setq(iFleet,iAge,0.0);
         }
      
   //calculate estimate of N                                         
   for (iFleet = 0, pTune->IncrFleet(&iFleet); iFleet <= NTuneFleet; pTune->IncrFleet(&iFleet))                                                       
      for (iAge = MinTuneAge[iFleet]; iAge <= min(MaxTuneAge[iFleet],Control.LRcrtAge); iAge++)
         for (iYear = MinTuneYear[iFleet]; iYear <= MaxTuneYear[iFleet]; iYear++)
            if (Getq(iFleet,iAge) > 0.0 && Theta(iFleet, iAge, iYear))
               SetNHat(iFleet,iAge,iYear,  GetCorrectedCPUE(iFleet,iAge,iYear) / Getq(iFleet,iAge, iYear));      
            else
               SetNHat(iFleet,iAge,iYear,  0.0);      
                            
   //calculate estimate of SE                         
   for (iFleet = 0, pTune->IncrFleet(&iFleet); iFleet <= NTuneFleet; pTune->IncrFleet(&iFleet))                                                       
      //Calculate mean NHat
      for (iAge = MinTuneAge[iFleet]; iAge <= min(MaxTuneAge[iFleet],Control.LRcrtAge); iAge++)
         {
         for (iYear = MaxTuneYear[iFleet], SumWt = SSq = 0.0; iYear >= MinTuneYear[iFleet]; iYear--) 
            {                                       
            Weight = Theta(iFleet, iAge, iYear)*Tau(iYear); 
   
            SumWt += Weight;                     
            
            if (Theta(iFleet, iAge, iYear) && N[iAge][iYear] > 0.0)                         
               SSq   += (Weight*SQR( log(N[iAge][iYear]/GetCorrectedCPUE(iFleet,iAge,iYear)) - log(1/Getq(iFleet,iAge)) ) );
            }                      
 
          if (SumWt > 1)
            {
            for (iYear = MinTuneYear[iFleet]; iYear <= MaxTuneYear[iFleet]; iYear++)     
               SetNSE(iFleet,iAge,iYear,(sqrt(SSq/(SumWt-1))*sqrt(1+1/SumWt)));              
            }
         }     
         
   return TRUE;
   }

bool ExtendedSurvivorsAnalysis::XSARecruitsSE(void)
   {
   // VPA User Guide equation (10)
   short iFleet, iAge, iYear;
   double lCPUEBar, TS, SSTS, lSSN, lSSCPUE;
   
   for (iFleet = 0, pTune->IncrFleet(&iFleet); iFleet <= NTuneFleet; pTune->IncrFleet(&iFleet))                                                       
      {
      //Calculate mean NHat
      for (iAge = MinTuneAge[iFleet]; iAge <= max(MinAge,min(Control.LRcrtAge, GetMaxTuneAge(iFleet))); iAge++)
         {
         for (iYear = MinTuneYear[iFleet], SSTS = lCPUEBar = 0.0; iYear <= MaxTuneYear[iFleet]; iYear++)                                                       
            {                                                             
            if (Theta(iFleet, iAge, iYear))
               {
               TS        = Theta(iFleet, iAge, iYear)*Tau(iYear);
               SSTS     += TS;
               lCPUEBar += log(GetCorrectedCPUE(iFleet,iAge,iYear))*TS;
               }
            }     
         lCPUEBar /= SSTS;
                  
         for (iYear = MinTuneYear[iFleet], SSTS = lSSN = lSSCPUE = 0.0; iYear <= MaxTuneYear[iFleet]; iYear++)                                                       
            {
            if (Theta(iFleet, iAge, iYear) && N[iAge][iYear] > 0.0)
               {
               TS       = Theta(iFleet, iAge, iYear)*Tau(iYear);
               SSTS    += TS;
               lSSN    += TS*SQR(log(GetNHat(iFleet,iAge,iYear)) -  log(N[iAge][iYear]));
               lSSCPUE += TS*SQR(log(GetCorrectedCPUE(iFleet,iAge,iYear)) -  lCPUEBar);
               }
            }
         
         for (iYear = MinTuneYear[iFleet]; iYear <= MaxTuneYear[iFleet]; iYear++)
            {
            if (Theta(iFleet, iAge, iYear) && SSTS > 2.0)
               SetNSE(iFleet,iAge,iYear, sqrt(lSSN/(SSTS-2))
                                                *sqrt(1 + 1/SSTS + SQR(log(GetCorrectedCPUE(iFleet,iAge,iYear)) -  lCPUEBar)/lSSCPUE)
                                               );
            }
        }            
     }    
   
   return true;
   };

bool ExtendedSurvivorsAnalysis::iFleetQ(void)
   {
   short iFleet, iAge, iYear, iTargetYear, i;
   double TS, Top, Bottom;             

   //equation (11) VPA user guide by Year
   switch(Control.SmoothQ)
      {
      case 1:

         for (iFleet = 0, pTune->IncrFleet(&iFleet); iFleet <= NTuneFleet; pTune->IncrFleet(&iFleet))
            {
            for (iAge = max(MinTuneAge[iFleet],Control.LRcrtAge+1); iAge <= min(MaxTuneAge[iFleet],Control.ConQAge); iAge++)
               for (iTargetYear = MinTuneYear[iFleet]; iTargetYear <= MaxTuneYear[iFleet]; iTargetYear++)
                  {  
                  Top = Bottom = 0.0;
                  Setqn(iFleet,iAge, iTargetYear, 0);
                  for (iYear = max(MinTuneYear[iFleet], iTargetYear-Control.TSRange); iYear <= min(MaxTuneYear[iFleet], iTargetYear+Control.TSRange); iYear++)
                     if (Theta(iFleet, iAge, iYear) && (TS = Tau(iYear, iTargetYear)) > 0.0)  
                        {
                        Bottom += TS;
                        Top    += TS*(log(N[iAge][iYear]) - log(GetCorrectedCPUE(iFleet,iAge,iYear)));
                           
                        (qn[iFleet][iAge][iTargetYear])++; 
                        }
                     
                  if (Bottom > 0) 
                     Setq(iFleet,iAge, iTargetYear, (1/exp(Top/Bottom)));
                  else
                     Setq(iFleet,iAge, iTargetYear, 0.0); 

                  //not yet calculated
                  SetqSE(iFleet,iAge, iTargetYear, 0.0);
                  }
            }                                           
      break;
         
      case 2:
  
      break;
      
      case 0: default:

         for (iFleet = 0, pTune->IncrFleet(&iFleet); iFleet <= NTuneFleet; pTune->IncrFleet(&iFleet))
            {
            for (iAge = max(MinTuneAge[iFleet],Control.LRcrtAge+1); iAge <= min(MaxTuneAge[iFleet],Control.ConQAge); iAge++)
               {  
               Top = Bottom = 0.0;
               Setqn(iFleet,iAge, 0);
               for (iYear = MinTuneYear[iFleet]; iYear <= MaxTuneYear[iFleet]; iYear++)
                  if (Theta(iFleet, iAge, iYear) && (TS = Tau(iYear)) > 0.0)  
                     {
                     Bottom += TS;
                     Top    += TS*(log(N[iAge][iYear]) - log(GetCorrectedCPUE(iFleet,iAge,iYear)));
                  
                     (qn[iFleet][iAge][MinTuneYear[iFleet]])++; 
                     }
                  
               if (Bottom > 0) 
                  Setq(iFleet,iAge,1/exp(Top/Bottom));
               else
                  Setq(iFleet,iAge, 0.0); 

               //not yet calculated
               SetqSE(iFleet,iAge, 0.0);
               }
            }
                                           
         break;
         }
      
   return true;
   }

bool  ExtendedSurvivorsAnalysis::CalcQ(void)
   {                                        
   //equation (4) VPA user guide
   short iAge, iYear, iFleet;                                        
                
   iFleetQ();
   
   for (iFleet = 0, pTune->IncrFleet(&iFleet); iFleet <= NTuneFleet; pTune->IncrFleet(&iFleet))                                                       
      for (iAge = max(Control.LRcrtAge+1, GetMinTuneAge(iFleet)); iAge <= MaxTuneAge[iFleet]; iAge++)
         for (iYear = MinTuneYear[iFleet]; iYear <= MaxTuneYear[iFleet]; iYear++)
            {
            if (Getq(iFleet,iAge, iYear) > 0)
               SetNHat(iFleet,iAge,iYear,  GetCorrectedCPUE(iFleet,iAge,iYear) / Getq(iFleet,iAge,iYear));      
            else
               SetNHat(iFleet,iAge,iYear,  0.0);      
         }                             

   XSANHatQSE();

   return true;
   }  

bool ExtendedSurvivorsAnalysis::XSANHatQSE(void)
   {
   short iFleet, iAge, iYear, iTargetYear, n;
               
   double Weight, SumWt, SSq;

   //equation (12) VPA user guide
   switch(Control.SmoothQ)
      {
      case 1:  

         for (iFleet = 0, pTune->IncrFleet(&iFleet); iFleet <= NTuneFleet; pTune->IncrFleet(&iFleet))                                                       
            //Calculate mean NHat
            for (iAge = max(MinAge,max(Control.LRcrtAge+1, GetMinTuneAge(iFleet))); iAge <= MaxTuneAge[iFleet]; iAge++)
               {
               for (iYear = MinTuneYear[iFleet], n=0; iYear <= MaxTuneYear[iFleet]; iYear++)                                                       
                  if (Theta(iFleet, iAge, iYear))
                     n++; 

               for (iTargetYear = MaxTuneYear[iFleet]; iTargetYear >= MinTuneYear[iFleet]; iTargetYear--) 
                  {
                  for (iYear = max(MinTuneYear[iFleet], iTargetYear-Control.TSRange), SumWt = SSq = 0.0; iYear <= min(MaxTuneYear[iFleet], iTargetYear+Control.TSRange); iYear++)
                     if (Theta(iFleet, iAge, iYear))
                        {
                        Weight = Tau(iYear, iTargetYear);
                        SumWt += Weight;                     
                        SSq   += (Weight*SQR( log(N[iAge][iYear]/GetCorrectedCPUE(iFleet,iAge,iYear)) - log(1/Getq(iFleet,iAge, iYear)) ) );
                        }
                  
                  if (SumWt > 1)
                     SetNSE(iFleet,iAge,iTargetYear,(sqrt(SSq/(SumWt-1))*sqrt(1+1/SumWt)));                
                  }
               }   
      break;
      
      case 2:

         for (iFleet = 0, pTune->IncrFleet(&iFleet); iFleet <= NTuneFleet; pTune->IncrFleet(&iFleet))                                                       
            //Calculate mean NHat
            for (iAge = max(MinAge,max(Control.LRcrtAge+1, GetMinTuneAge(iFleet))); iAge <= MaxTuneAge[iFleet]; iAge++)
               {
               for (iYear = MinTuneYear[iFleet], n=0; iYear <= MaxTuneYear[iFleet]; iYear++)                                                       
                  if (Theta(iFleet, iAge, iYear))
                     n++; 

               for (iYear = MaxTuneYear[iFleet], SumWt = SSq = 0.0; iYear >= MinTuneYear[iFleet]; iYear--) 
                  if (Theta(iFleet, iAge, iYear))
                     {
                     SumWt++;                     
                     SSq   += (SQR( log(N[iAge][iYear]/GetCorrectedCPUE(iFleet,iAge,iYear)) - log(1/Getq(iFleet,iAge))));
                     }
               
               for (iYear = MinTuneYear[iFleet]; iYear <= MaxTuneYear[iFleet]; iYear++)                                                       
                  if (SumWt > 1)
                     SetNSE(iFleet,iAge,iYear,(sqrt(SSq/(SumWt-1))*sqrt(1+1/SumWt)));                
               }
      break;


      case 0: default:  
      
         for (iFleet = 0, pTune->IncrFleet(&iFleet); iFleet <= NTuneFleet; pTune->IncrFleet(&iFleet))                                                       
            //Calculate mean NHat
            for (iAge = max(MinAge,max(Control.LRcrtAge+1, GetMinTuneAge(iFleet))); iAge <= MaxTuneAge[iFleet]; iAge++)
               {
               for (iYear = MinTuneYear[iFleet], n=0; iYear <= MaxTuneYear[iFleet]; iYear++)                                                       
                  if (Theta(iFleet, iAge, iYear))
                     n++; 

               for (iYear = MaxTuneYear[iFleet], SumWt = SSq = 0.0; iYear >= MinTuneYear[iFleet]; iYear--) 
                  if (Theta(iFleet, iAge, iYear))
                     {
                     Weight = Tau(iYear);
                     SumWt += Weight;                     
                     SSq   += (Weight*SQR( log(N[iAge][iYear]/GetCorrectedCPUE(iFleet,iAge,iYear)) - log(1/Getq(iFleet,iAge))));
                     }
               
               for (iYear = MinTuneYear[iFleet]; iYear <= MaxTuneYear[iFleet]; iYear++)                                                       
                  if (SumWt > 1)
                     SetNSE(iFleet,iAge,iYear,(sqrt(SSq/(SumWt-1))*sqrt(1+1/SumWt)));                
               }
      break;
      }

                  
   return true;
   };      
          
bool ExtendedSurvivorsAnalysis::Initialise(void)
   {      

   if (!TestFlag)
      {
      //Max iAge
      for (short iYear=MinYear; iYear<MaxYear; iYear++)
         N[MaxAge+1][iYear+1] = max(1.0,pCatch->GetCatch(MaxAge, iYear));
   
      //Max iYear   
      for (short iAge =MinAge; iAge<= MaxAge; iAge++)
         N[iAge+1][MaxYear+1] = max(1.0,pCatch->GetCatch(iAge, MaxYear));
      }
   else
      {
      //Max iAge
      for (short iYear=MinYear; iYear<MaxYear; iYear++)
         N[MaxAge+1][iYear+1] = N[MaxAge][iYear]*exp(-F[MaxAge][iYear]-M[MaxAge][iYear]);
   
      //Max iYear   
      for (short iAge =MinAge; iAge<= MaxAge; iAge++)
         N[iAge+1][MaxYear+1] = N[iAge][MaxYear]*exp(-F[iAge][MaxYear]-M[iAge][MaxYear]);
      }

   return true;
   }
   
bool  ExtendedSurvivorsAnalysis::RunVPAFromXSA(void)
   {   
   short iAge;

   if (Control.VPA)
      RunVPA();
   else
      RunCohortAnalysis();
 
   //Calc mean change in Terminal F                       
   for (iAge = MinAge, MSS=0.0; iAge<= MaxAge; iAge++)
      {
      MSS += SQR(PrevF[1][iAge] - F[iAge][MaxYear]); 
      PrevF[1][iAge] = PrevF[2][iAge];
      PrevF[2][iAge] = F[iAge][MaxYear];
      } 
         
   MSS /= (MaxAge - MinAge + 1);

   return true;
   } 
                                    
bool  ExtendedSurvivorsAnalysis::XSATermPop(void)
   {
   if (Control.Shk2N == TRUE )
      Shrink2PopMean();
      
   if (Control.Shk2F == TRUE) 
      Shrink2FMean();
                
   CalcTermPops();
       
   //OutputWts = false;
   //if (OutputWts)
   //   OutputTermPopWts();
 
   return true;                                                                                                     
   }                                                                

bool ExtendedSurvivorsAnalysis::CalcTermPops(void)
   {                
   //open file
   //ofstream outFile("c:\\temp\\xsaNew.csv", ios::out);
   //ofstream outFile2("c:\\temp\\xsaNewFS.csv", ios::out);
   //ofstream outFile3("c:\\temp\\xsaNewU.csv", ios::out);
   //outFile  << "iFleet, Age, Year, YearClass, Type, Sumtop, SumBottom" nl;
   //outFile2 << "iYearClass, Param.Shk2F, GetFShrinkMean, Param.FShkSE" nl;
   //outFile3 << "iFleet, iAge, iYear, iYearClass, pTune->GetMaxAge(iFleet), pTune->GetMinAge(iFleet), pTune->GetMaxYear(iFleet), pTune->GetMinYear(iFleet) ,Theta(iFleet, iAge, iYear), Getqn(iFleet, iAge)" nl;
                 

   short iAge, iFleet, iYear, iYearClass, RecruitAge;
   double SumTop, SumBottom,
          Wt;           

   //equation (13) VPA user guide

   double ECF, ECZ;

   static int Iter = 0;

   Iter++;
   for (iYearClass=MinYear-MaxAge; iYearClass <= MaxYear-MinAge; iYearClass++)
      {
      SumTop         =
      SumBottom      = 0.0;
      ECF            =                
      ECZ            = 1.0;

      for (iAge = min(MaxAge, MaxYear - iYearClass); iAge >= max(MinAge, MinYear - iYearClass); iAge--)
         {
         iYear = iYearClass + iAge; 
         ECF *= exp(F[iAge][iYear]);               
         ECZ *= exp(F[iAge][iYear] + GetM(iAge, iYear));
         for (iFleet=0, pTune->IncrFleet(&iFleet); iFleet<=NTuneFleet; pTune->IncrFleet(&iFleet))
            {
            if (pTune->GetMaxAge(iFleet)  >= iAge  && pTune->GetMinAge(iFleet)  <= iAge &&
                pTune->GetMaxYear(iFleet) >= iYear && pTune->GetMinYear(iFleet) <= iYear  )
                {
                //outFile3 << iFleet cs << iAge cs << iYear cs << iYearClass cs << pTune->GetMaxAge(iFleet) cs << pTune->GetMinAge(iFleet) cs << pTune->GetMaxYear(iFleet)cs << pTune->GetMinYear(iFleet) cs << Theta(iFleet, iAge, iYear) cs << Getqn(iFleet, iAge) nl;
                if (Theta(iFleet, iAge, iYear) && Getqn(iFleet, iAge) > 1)
                  {
                  if (Control.SmoothQ)
                     Wt = Weight(iFleet, iAge, iYear, iYearClass+MaxAge, ECF);                                                                          
                  else
                     Wt = Weight(iFleet, iAge, iYear, ECF);                                                                          

                  SumTop    += Wt*(log(GetNHat(iFleet,iAge,iYear)) - log(ECZ)); 
                  SumBottom += Wt;
                  
                  //outFile << iFleet cs << iAge cs << iYear cs << iYearClass cs << "U" cs << SumTop cs << SumBottom nl;
                  }
               }
            }       
         }

      //outFile2 << iYearClass cs << Control.Shk2F cs << GetFShrinkMean(iYearClass) cs << Control.FShkSE nl;
      if (Control.Shk2F == TRUE && GetFShrinkMean(iYearClass) > 0 && Control.FShkSE > 0.0)
         {
         SumTop    += log(GetFShrinkMean(iYearClass))/(SQR(Control.FShkSE));
         SumBottom += 1/(SQR(Control.FShkSE));
         //outFile << "" cs << "" cs << "" cs <<  iYearClass cs << "Shk2F" cs << SumTop cs << SumBottom nl;
         }

      iAge = min(MaxAge, MaxYear - iYearClass);
      RecruitAge = MaxYear - iYearClass;                                                 

      if (Control.Shk2N == TRUE && RecruitAge <= Control.LRcrtAge && GetPopShrinkMean(RecruitAge) > 0.0 && PopShrinkSE[RecruitAge] > 0.0)             
        {
        SumTop    += log(GetPopShrinkMean(RecruitAge))/(SQR(GetPopShrinkSE(RecruitAge)));                        
        SumBottom += 1/(SQR(GetPopShrinkSE(RecruitAge)));
        //outFile << "" cs << "" cs << ""  cs << iYearClass cs << "Shk2N" cs << SumTop cs << SumBottom nl;
        }
              
      if (SumBottom > 0)
         N[iAge + 1][iYearClass + iAge + 1] = exp(SumTop/SumBottom);
      else
         N[iAge + 1][iYearClass + iAge + 1] = 0.0;
      } 
   
   return true;
   }           
      
bool ExtendedSurvivorsAnalysis::Shrink2PopMean(void)
   //Done for recruiting iAges only
   {                     
   short iYear, iAge, n;
   double tau, Top, Bottom;
     
   //equation (15) VPA User Guide                         
   for (iAge = min(MaxAge,Control.LRcrtAge); iAge >= MinAge; iAge--)
      {
      for (iYear=MinYear, Top = Bottom = 0.0; iYear<= MaxYear; iYear++)
         {
         tau     = Tau(iYear); 
         if (tau > 0.0 && N[iAge+1][iYear] > 0)
            {
            Top    += tau*log(N[iAge+1][iYear]);
            Bottom += tau;  
            }
         }                                  
      if (Bottom > 0)   
         SetPopShrinkMean(iAge,exp(Top/Bottom));  
      else
         SetPopShrinkMean(iAge,0.0);   
      }
                                  
   //equation (16) VPA User Guide                                
   for (iAge = min(MaxAge,Control.LRcrtAge); iAge >= MinAge; iAge--)
      {
      for (iYear=MinYear, Top = Bottom = 0.0, n=0; iYear<= MaxYear; iYear++)
         {    
         n++; 
         tau     = Tau(iYear); 
         Bottom += tau;
         if (N[iAge+1][iYear] > 0 && GetPopShrinkMean(iAge) > 0)
            Top    += tau*SQR(log(N[iAge+1][iYear]) - log(GetPopShrinkMean(iAge)));
         }                        
      if (Bottom-1 > 0)
         SetPopShrinkSE(iAge,sqrt(Top/(Bottom-1)));  
      } 
      
  return TRUE;                                                
  }    

bool ExtendedSurvivorsAnalysis::Shrink2FMean(void)
   //Done for all iAges
   {                     
   short iYearClass, iAge, iYear;

   FBarShrink();

   for (iYearClass = MaxYear-MinAge; iYearClass >= MinYear-MaxAge; iYearClass--)
   //for (iYearClass = MinYear-MaxAge; iYearClass >= MaxYear-MinAge; iYearClass++)
      {
      if (iYearClass >= MaxYear - MaxAge)
         {
         //equation (18) VPA User Guide                         
         iAge   = MaxYear - iYearClass;
        
         if (Control.Shk2F && Control.Shk2FYr > 0 && F[iAge+1][MaxYear+1] > 0.0)
            SetFShrinkMean(iYearClass, pCatch->GetCatch(iAge, MaxYear)
                                       *(F[iAge+1][MaxYear+1] + GetM(iAge, MaxYear))
                                       /(F[iAge+1][MaxYear+1]  
                                       *(exp(F[iAge+1][MaxYear+1] + GetM(iAge, MaxYear)) - 1)));
         else
            SetFShrinkMean(iYearClass, 0.0);
         }
      else                                    
         {
         //equation (20) VPA User Guide                         
         iYear  = MaxAge + iYearClass;  
         
         if (Control.Shk2F && Control.Shk2FAge > 0)
            {
            double eCumZ = 1.0;

            for (short i = 0; i > Control.Shk2FTarget; i--)
               eCumZ *= exp(-F[MaxAge+i][iYear+i] - GetM(MaxAge+i, iYear+i));
            
            SetFShrinkMean(iYearClass, eCumZ
                                       *pCatch->GetCatch(MaxAge+Control.Shk2FTarget, iYear+Control.Shk2FTarget)
                                       *(F[MaxAge+1][iYear+1+Control.Shk2FTarget] + GetM(MaxAge+Control.Shk2FTarget, iYear+Control.Shk2FTarget))
                                       /(F[MaxAge+1][iYear+1+Control.Shk2FTarget]  
                                       *(exp(F[MaxAge+1][iYear+1+Control.Shk2FTarget] + GetM(MaxAge+Control.Shk2FTarget, iYear+Control.Shk2FTarget)) - 1)));
            }
         else
            SetFShrinkMean(iYearClass, 0.0);
         
         }
      }                                
  return TRUE;                  
  }      

bool ExtendedSurvivorsAnalysis::FBarShrink(void)
   {               
   short iYear, iAge, j, k;       

   //equation (17) VPA User Guide                         
   if (Control.Shk2FYr > 0)
      for (iAge = MinAge; iAge < MaxAge; iAge++)
         {
         for (iYear = max(MinYear, MaxYear-Control.Shk2FYr), j=0, F[iAge+1][MaxYear+1] = 0.0; iYear < MaxYear; iYear++)
            {
            if (F[iAge][iYear] > 0.0)
               {
               j++;
               if (GM) 
                  F[iAge+1][MaxYear+1] += log(F[iAge][iYear]);
               else
                  F[iAge+1][MaxYear+1] += F[iAge][iYear];
               }
            }
         if (j > 0)
            {
            F[iAge+1][MaxYear+1] /= j;

            if (GM) 
               F[iAge+1][MaxYear+1]  = exp(F[iAge+1][MaxYear+1]);
            }
         }
                 
   //equation (19) VPA User Guide                         
   //if (Control.Shk2FAge > 0)
      for (iYear = MinYear; iYear <= MaxYear; iYear++)
         {
         for (iAge = max(MinAge, MaxAge-Control.Shk2FAge), k=0, F[MaxAge+1][iYear+1] = 0.0; iAge < MaxAge; iAge++)
            {
            if (F[iAge][iYear] > 0.0)
               {
               k++;                                
               if (GM) 
                  F[MaxAge+1][iYear+1] += log(F[iAge][iYear]);
               else
                  F[MaxAge+1][iYear+1] += F[iAge][iYear];
               }
            }
         if (k > 0)
            {
            F[MaxAge+1][iYear+1] /= k;
            if (GM) 
               F[MaxAge+1][iYear+1]  = exp(F[MaxAge+1][iYear+1]);
            }
         }
                                   
   return true;
   }               
                             
//Weights
void ExtendedSurvivorsAnalysis::SetTau(void)
   {
   double DistanceFromPresent;

   flallocArray(&TSWeight, MinYear, MaxYear);
      
   for (int iYear=MinYear; iYear<=MaxYear; iYear++)
      {
      DistanceFromPresent = MaxYear - iYear;
   
      if (DistanceFromPresent <  0 || DistanceFromPresent >= Control.TSRange)
         TSWeight[iYear] = 0.0;                                      
      else if (DistanceFromPresent == 0 || Control.TSPower == 0)
         TSWeight[iYear] = 1.0;
      else
         TSWeight[iYear] = min(1,max(0,pow(1.0-pow(DistanceFromPresent/Control.TSRange,Control.TSPower),Control.TSPower))); 
      }          
   }

double ExtendedSurvivorsAnalysis::Tau(short iYear)
   {       
   if (iYear > MaxYear || iYear <= MaxYear -  Control.TSRange || iYear < MinYear)
      return 0.0;
   else
      return  TSWeight[iYear]; 
   }          
                             
double ExtendedSurvivorsAnalysis::Tau(short iYear, short TargetYear)
   {
   double DistanceFromPresent;
   
   DistanceFromPresent = TargetYear - iYear;

   if (DistanceFromPresent < 0.0)
      DistanceFromPresent *= -1.0;
   
   if (DistanceFromPresent <  0 || DistanceFromPresent >= Control.TSRange)
      return 0.0;                                      
   
   if (DistanceFromPresent == 0 || Control.TSPower == 0)
      return 1.0;
         
   return  min(1,max(0,pow(1.0-pow(DistanceFromPresent/Control.TSRange,Control.TSPower),Control.TSPower))); 
   }          

bool ExtendedSurvivorsAnalysis::Theta(short iFleet, short iAge, short iYear)
   {          
   if (iYear  < MinTuneYear[iFleet] || iYear > MaxTuneYear[iFleet] ||
       iAge   < MinTuneAge[ iFleet] || iAge  > MaxTuneAge[iFleet]  ||
       iFleet < 0 || iFleet > NTuneFleet && pTune->Tuning(iFleet))
      return FALSE; 

    //if (pTune->GetCPUE(iFleet, iAge, iYear) > 0.0 && (Iters>1 ? !_isnan(GetNHat(iFleet,iAge,iYear)) && _finite(GetNHat(iFleet,iAge,iYear)) : true))
    if (pTune->GetCPUE(iFleet, iAge, iYear) > 0.0) // && (Iters>1 ? !_isnan(GetNHat(iFleet,iAge,iYear)) && _finite(GetNHat(iFleet,iAge,iYear)) : true))
       return TRUE;
    else
       return FALSE;
    }                 
                                    
double  ExtendedSurvivorsAnalysis::Weight(short iFleet, short iAge, short iYear, double ECF)
   { 
   if (Control.NSEWts == NSE_WTS_AUTOMATIC)
      if (Theta(iFleet, iAge, iYear))
         return Tau(iYear) / (SQR(max(GetNSE(iFleet,iAge,iYear), Control.MinNSE))*ECF);
      else
         return 0.0;
   else
      if (Theta(iFleet, iAge, iYear))
         return  Tau(iYear) / (SQR(max(GetManualNSE(iFleet, iAge, iYear), Control.MinNSE))*ECF);
      else
         return 0.0;   
   }

double  ExtendedSurvivorsAnalysis::Weight(short iFleet, short iAge, short iYear, short TargetYear, double ECF)
   { 
   if (Theta(iFleet, iAge, iYear))
      return   Tau(iYear, TargetYear) / (SQR(max(GetNSE(iFleet,iAge,iYear), Control.MinNSE))*ECF);
   else
      return 0.0;                                                              
   }

short  ExtendedSurvivorsAnalysis::GetNTuneFleet(void)
   {
   return NTuneFleet;
   }
   
short  ExtendedSurvivorsAnalysis::GetMinTuneAge(short iFleet)    
   {             
   return MinTuneAge[min(NTuneFleet,max(1,iFleet))]; 
   } 
   
short  ExtendedSurvivorsAnalysis::GetMaxTuneAge(short iFleet)    
   { 
   return MaxTuneAge[min(NTuneFleet,max(1,iFleet))]; 
   }
   
short  ExtendedSurvivorsAnalysis::GetMinTuneYear(short iFleet) 
   { 
   return MinTuneYear[min(NTuneFleet,max(1,iFleet))]; 
   }
   
short  ExtendedSurvivorsAnalysis::GetMaxTuneYear(short iFleet)  
   { 
   return MaxTuneYear[min(NTuneFleet,max(1,iFleet))]; 
   } 

double ExtendedSurvivorsAnalysis::GetNHat(short iFleet, short iAge, short iYear)  
   {                
   if (iFleet < 0 || iFleet > NTuneFleet && pTune->Tuning(iFleet)  ||
       iYear  < MinTuneYear[iFleet] || iYear > MaxTuneYear[iFleet] ||
       iAge   < MinTuneAge[ iFleet] || iAge  > MaxTuneAge[iFleet]    )
      return 0.0; 

   return min(NHat[iFleet][iAge][iYear],exp(30.0)); 
   } 

double ExtendedSurvivorsAnalysis::GetNSE(short iFleet, short iAge, short iYear)  
   {
   if (iYear  < MinTuneYear[iFleet] || iYear > MaxTuneYear[iFleet] ||
       iAge   < MinTuneAge[ iFleet] || iAge  > MaxTuneAge[iFleet]  ||
       iFleet < 0 || iFleet > NTuneFleet && pTune->Tuning(iFleet))
      return 99.99; 

   if (NHat[iFleet][iAge][iYear] >= exp(30.0))
      return 99.99; 
   else
      return NSE[iFleet][iAge][iYear]; 
   } 

double ExtendedSurvivorsAnalysis::GetManualNSE(short iFleet, short iAge, short iYear)  
   {                
   return ManualNSE[max(iFleet,1)][max(0,iAge)][iYear]; // * GetNHat(iFleet, iAge, iYear);                   
   }     

double ExtendedSurvivorsAnalysis::GetTermPopSE(short iYearClass)  
   {                
   short i;
   
   if (iYearClass < MinYear-MaxAge || iYearClass > MaxYear-MinAge)
      i = 1;
      
   return TermPopSE[iYearClass]; 
   } 

double ExtendedSurvivorsAnalysis::GetPopShrinkMean(short iAge)
   {
   return PopShrinkMean[iAge]; 
   }   
   
void ExtendedSurvivorsAnalysis::SetPopShrinkMean(short iAge, double Value)
   {
   if (iAge < MinAge || iAge > Control.LRcrtAge)
      return;
   
   PopShrinkMean[iAge] = Value; 
   }   

void ExtendedSurvivorsAnalysis::SetPopShrinkSE(short iAge, double Value)
   {
   if (iAge < MinAge || iAge > Control.LRcrtAge)
      return;
   
   PopShrinkSE[iAge] = Value; 
   }   
   
double ExtendedSurvivorsAnalysis::GetPopShrinkSE(short iAge)
   {
   return PopShrinkSE[min(iAge,max(MinAge,Control.LRcrtAge))]; 
   }   

double ExtendedSurvivorsAnalysis::GetCorrectedCPUE(short iFleet, short iAge, short iYear)  
   {         
   if (iFleet < 1                   || iFleet > NTuneFleet         || 
       iAge   < MinTuneAge[iFleet]  || iAge   > MaxTuneAge[iFleet] || 
       iYear  < MinTuneYear[iFleet] || iYear  > MaxTuneYear[iFleet]  )
      return 0.0;
   else
      return CorrectedCPUE[iFleet][iAge][iYear]; 
   } 

double ExtendedSurvivorsAnalysis::GetFShrinkMean(short iYearClass)  
   {             
   if (iYearClass < MinYear-MaxAge || iYearClass > MaxYear-MinAge)
      return 0.0;
  
   return FShrinkMean[iYearClass]; 
   } 

void ExtendedSurvivorsAnalysis::SetFShrinkMean(short iYearClass, double Value)  
   {                
   if (iYearClass < MinYear-MaxAge || iYearClass > MaxYear-MinAge)
      return;
      
   FShrinkMean[iYearClass] = Value; 
   } 

void ExtendedSurvivorsAnalysis::SetN(short iAge, short iYear, double Value)
   {
   if (iAge   < MinAge  || iAge   > MaxAge ||             
       iYear  < MinYear || iYear  > MaxYear+1 )
       return;

   N[iAge][iYear] = Value;
   }                                   

void ExtendedSurvivorsAnalysis::SetNHat(short iFleet, short iAge, short iYear, double Value)
   {
   if (iFleet < 1                   || iFleet > NTuneFleet         || 
       iAge   < MinTuneAge[iFleet]  || iAge   > MaxTuneAge[iFleet] ||             
       iYear  < MinTuneYear[iFleet] || iYear  > MaxTuneYear[iFleet] )
       NHat[iFleet][iAge][iYear] = 0.0;

   NHat[iFleet][iAge][iYear] = Value;
   }                                   

void ExtendedSurvivorsAnalysis::SetNSE(short iFleet, short iAge, short iYear, double Value)
   {
   if (iFleet < 1                   || iFleet > NTuneFleet         || 
       iAge   < MinTuneAge[iFleet]  || iAge   > MaxTuneAge[iFleet] ||             
       iYear  < MinTuneYear[iFleet] || iYear  > MaxTuneYear[iFleet] )
       return;

   NSE[iFleet][iAge][iYear] = Value;
   }                                   

void ExtendedSurvivorsAnalysis::Setq(short iFleet, short iAge, double Value)
   {
   if (iFleet < 1                   || iFleet > NTuneFleet        || 
       iAge   < MinTuneAge[iFleet]  || iAge   > MaxTuneAge[iFleet]  )
       return;

   q[iFleet][iAge][MinTuneYear[iFleet]] = Value;
   }                                   
   
void ExtendedSurvivorsAnalysis::Setq_b(short iFleet, short iAge, double Value)
   {
   if (iFleet < 1                   || iFleet > NTuneFleet        || 
       iAge   < MinTuneAge[iFleet]  || iAge   > MaxTuneAge[iFleet]  )
       return;

   q_b[iFleet][iAge][MinTuneYear[iFleet]] = Value;
   }

void ExtendedSurvivorsAnalysis::Setqn(short iFleet, short iAge, short Value)
   {
   if (iFleet < 1                   || iFleet > NTuneFleet        || 
       iAge   < MinTuneAge[iFleet]  || iAge   > MaxTuneAge[iFleet]  )
       return;

   qn[iFleet][iAge][MinTuneYear[iFleet]] = Value;
   }

void ExtendedSurvivorsAnalysis::SetqSE(short iFleet, short iAge, double Value)
   {
   if (iFleet < 1                   || iFleet > NTuneFleet        || 
       iAge   < MinTuneAge[iFleet]  || iAge   > MaxTuneAge[iFleet]  )
       return;

   qSE[iFleet][iAge][MinTuneYear[iFleet]] = Value;
   }

void ExtendedSurvivorsAnalysis::Setq(short iFleet, short iAge, short iYear, double Value)
   {
   if (iFleet < 1                   || iFleet > NTuneFleet         || 
       iAge   < MinTuneAge[iFleet]  || iAge   > MaxTuneAge[iFleet] ||             
       iYear  < MinTuneYear[iFleet] || iYear  > MaxTuneYear[iFleet]  )
       return;

   q[iFleet][iAge][iYear] = Value;
   }                                   
   
void ExtendedSurvivorsAnalysis::Setq_b(short iFleet, short iAge, short iYear, double Value)
   {
   if (iFleet < 1                   || iFleet > NTuneFleet         || 
       iAge   < MinTuneAge[iFleet]  || iAge   > MaxTuneAge[iFleet] ||             
       iYear  < MinTuneYear[iFleet] || iYear  > MaxTuneYear[iFleet]  )
       return;

   q_b[iFleet][iAge][iYear] = Value;
   }

void ExtendedSurvivorsAnalysis::Setqn(short iFleet, short iAge, short iYear, short Value)
   {
      if (iFleet < 1                || iFleet > NTuneFleet         || 
       iAge   < MinTuneAge[iFleet]  || iAge   > MaxTuneAge[iFleet] ||             
       iYear  < MinTuneYear[iFleet] || iYear  > MaxTuneYear[iFleet]  )
       return;

   qn[iFleet][iAge][iYear] = Value;
   }

void ExtendedSurvivorsAnalysis::SetqSE(short iFleet, short iAge, short iYear, double Value)
   {
   if (iFleet < 1                   || iFleet > NTuneFleet         || 
       iAge   < MinTuneAge[iFleet]  || iAge   > MaxTuneAge[iFleet] ||             
       iYear  < MinTuneYear[iFleet] || iYear  > MaxTuneYear[iFleet]  )
       return;

   qSE[iFleet][iAge][iYear] = Value;
   }

void ExtendedSurvivorsAnalysis::CalcMeanF(short AvForNYrs, short MinAgeForOverallF,  short MaxAgeForOverallF)
   {  
   short iAge;
   double SumOverallF, SumNewOverallF;
   
   for (iAge = MinAgeForOverallF, SumOverallF = 0.0; iAge<=MaxAgeForOverallF; iAge++)
      SumOverallF += F[iAge][MaxYear];
   
   for (iAge =MinAge; iAge<=GetMaxCatchAge(); iAge++)
      {
      MeanF[iAge] = 0.0;
      for (short iYear= MaxYear; iYear>MaxYear-AvForNYrs; iYear--)
         MeanF[iAge] +=  F[iAge][iYear];
      MeanF[iAge]/=AvForNYrs;
      }                     

   for (iAge = MinAgeForOverallF, SumNewOverallF = 0.0; iAge<=MaxAgeForOverallF; iAge++)
      SumNewOverallF += F[iAge][MaxYear];
      
   for (iAge =MinAge; iAge<=GetMaxCatchAge(); iAge++)
      MeanF[iAge] /= SumOverallF/SumNewOverallF;
   } 
   
double ExtendedSurvivorsAnalysis::GetMeanF(short iAge)
   {  
   return MeanF[max(min(GetMaxCatchAge(),iAge),MinAge)];       
   }   
   
short ExtendedSurvivorsAnalysis::GetTSRange(void)
   {
   return Control.TSRange;
   }
   
short ExtendedSurvivorsAnalysis::GetTSPower(void)
   {
   return Control.TSPower;
   }
   
short ExtendedSurvivorsAnalysis::GetLRcrtAge(void)
   {
   return Control.LRcrtAge;
   }
   
short ExtendedSurvivorsAnalysis::GetConQAge(void)
   {
   return max(Control.ConQAge,Control.LRcrtAge+2);
   }
   
short ExtendedSurvivorsAnalysis::GetShk2N(void)
   {
   return Control.Shk2N;
   }
   
short ExtendedSurvivorsAnalysis::GetShk2F(void)
   {
   return Control.Shk2F;
   }
   
short ExtendedSurvivorsAnalysis::GetShk2FYr(void)
   {
   return Control.Shk2FYr;
   }
   
short ExtendedSurvivorsAnalysis::GetShk2FAge(void)
   {
   return Control.Shk2FAge;
   }
   
short ExtendedSurvivorsAnalysis::GetMaxIters(void)
   {
   return Control.MaxIters;
   }
   
short ExtendedSurvivorsAnalysis::GetPlusGroup(void)
   {
   return Control.PlusGroup;
   }
          
double ExtendedSurvivorsAnalysis::GetMSSTol(void)
   {
   return Control.MSSTol;
   }
   
double ExtendedSurvivorsAnalysis::GetMinNSE(void)
   {
   return Control.MinNSE;
   }
   
double ExtendedSurvivorsAnalysis::GetFShkSE(void)
   {
   return Control.FShkSE;
   }
   
double ExtendedSurvivorsAnalysis::GetFPreSpwn(void)
   {
   return Control.FPreSpwn;
   }
   
double ExtendedSurvivorsAnalysis::GetMPreSpwn(void)
   {
   return Control.MPreSpwn;   
   }                            
                          
double ExtendedSurvivorsAnalysis::Getq(short iFleet, short iAge)
   {
   if (iFleet < 1 || iFleet > NTuneFleet || min(iAge,Control.ConQAge)  < MinTuneAge[iFleet] )  
      return 0.0;

   return q[iFleet][min(iAge,Control.ConQAge)][MinTuneYear[iFleet]];
   }

double ExtendedSurvivorsAnalysis::Getq_b(short iFleet, short iAge)
   {
   if (iFleet < 1                   || iFleet > NTuneFleet ||
       iAge   < MinTuneAge[iFleet]  || iAge   > Control.LRcrtAge) 
      return 1.0;

   return q_b[iFleet][iAge][MinTuneYear[iFleet]];
   } 
   
double ExtendedSurvivorsAnalysis::GetqSE(short iFleet, short iAge)
   {
   if (iFleet < 1                   || iFleet > NTuneFleet ||
       iAge   < MinTuneAge[iFleet]  || iAge   > MaxTuneAge[iFleet]) 
      return 0;

   return qSE[iFleet][min(iAge,Control.ConQAge)][MinTuneYear[iFleet]];
   } 
   
short ExtendedSurvivorsAnalysis::Getqn(short iFleet, short iAge)
   {
   if (iFleet < 1                   || iFleet > NTuneFleet ||
       iAge   < MinTuneAge[iFleet]  || iAge   > MaxTuneAge[iFleet]) 
      return 0;

   return qn[iFleet][min(iAge,Control.ConQAge)][MinTuneYear[iFleet]];
   }   

double ExtendedSurvivorsAnalysis::Getq(short iFleet, short iAge, short iYear)
   {
   if (Control.SmoothQ)
      return q[min(max(iFleet,1),NTuneFleet)][min(max(iAge,MinTuneAge[iFleet]),Control.ConQAge)][min(max(iYear,MinTuneYear[iFleet]),MaxTuneYear[iFleet])];
   else 
      return Getq(iFleet, iAge);
   }

double ExtendedSurvivorsAnalysis::Getq_b(short iFleet, short iAge, short iYear)
   {
   if (Control.SmoothQ)
      return q_b[min(max(iFleet,1),NTuneFleet)][min(max(iAge,MinTuneAge[iFleet]),Control.ConQAge)][min(max(iYear,MinTuneYear[iFleet]),MaxTuneYear[iFleet])];
   else 
      return Getq_b(iFleet, iAge);
   } 
   
double ExtendedSurvivorsAnalysis::GetqSE(short iFleet, short iAge, short iYear)
   {
   if (Control.SmoothQ)
      return qSE[min(max(iFleet,1),NTuneFleet)][min(max(iAge,MinTuneAge[iFleet]),Control.ConQAge)][min(max(iYear,MinTuneYear[iFleet]),MaxTuneYear[iFleet])];
   else 
      return GetqSE(iFleet, iAge);
   } 
   
short ExtendedSurvivorsAnalysis::Getqn(short iFleet, short iAge, short iYear)
   {
   if (Control.SmoothQ)
      return qn[min(max(iFleet,1),NTuneFleet)][min(max(iAge,MinTuneAge[iFleet]),Control.ConQAge)][min(max(iYear,MinTuneYear[iFleet]),MaxTuneYear[iFleet])];
   else 
      return Getqn(iFleet, iAge);
   }   
             
LPDOUBLE ExtendedSurvivorsAnalysis::GetRecruits()
   {  
   return &N[MinAge][MinYear];       
   } 

void ExtendedSurvivorsAnalysis::CalcPlusGroupNF(void)
   {
   short iYear;
   
   double PlusGroupCatch;
   
   if (Control.PlusGroup != MaxAge+1)
      return;
   
   for (iYear=MinYear; iYear<= MaxYear; iYear++)
      {              
      PlusGroupCatch = 0.0;
            
      if (F[Control.PlusGroup-1][iYear] > 0.0)
         F[Control.PlusGroup][iYear] = F[Control.PlusGroup-1][iYear];                           
      else
         F[Control.PlusGroup][iYear] = F[max(MinAge,Control.PlusGroup-2)][iYear];                           
      
      for (short iAge = Control.PlusGroup; iAge <= pCatch->GetMaxAge(); iAge++)
         PlusGroupCatch += pCatch->GetCatch(iAge, iYear);
      
      if (F[Control.PlusGroup][iYear] > 0.0)
         N[Control.PlusGroup][iYear] =  PlusGroupCatch*(F[Control.PlusGroup][iYear]+GetM(Control.PlusGroup, iYear))
                                     /(F[Control.PlusGroup][iYear]*(1-exp(-F[Control.PlusGroup][iYear]-GetM(Control.PlusGroup, iYear))));  
      }
   
   N[Control.PlusGroup][MaxYear+1] += N[Control.PlusGroup][MaxYear]*exp(-F[Control.PlusGroup][iYear]-GetM(Control.PlusGroup-1, iYear));  
   }

double ExtendedSurvivorsAnalysis::GetCPUEHat(short iFleet, short iYear, short iAge)
   {
   double ReturnValue;

   if (pTune->GetEffort(iFleet, iYear) > 0.0)
      ReturnValue = pTune->GetEffort(iFleet, iYear) * pow(Getq(iFleet, iAge, iYear) * GetN(iAge, iYear), Getq_b(iFleet, iAge))
                     * AvCatch(iAge, iYear, pTune->GetStartFishing(iFleet), pTune->GetEndFishing(iFleet));
   else
      ReturnValue = 0.0;
   
   return ReturnValue;
   }
   
short ExtendedSurvivorsAnalysis::GetMaxCatchAge()
   {
   if (GetPlusGroup() == MaxAge+1)     
      return MaxAge+1;
   else
      return MaxAge;
   }                               

void ExtendedSurvivorsAnalysis::CalcWeightedLogMeanN(short iYearClass, double *pMean, double *pSE)
   {
   short  iAge, iYear;
  
   double Top, Bottom;

   *pMean = 
   *pSE   = 0.0;

   iAge = MaxYear - iYearClass;

   //equation (15) VPA User Guide                         
   for (iYear = max(MinYear,MaxYear-Control.TSRange), Top = Bottom = 0.0; iYear <= MaxYear; iYear++)
      if (Tau(iYear) > 0.0 && N[iAge+1][iYear] > 0)
         {
         Top    += Tau(iYear)*log(N[iAge][iYear]);
         Bottom += Tau(iYear);  
         }
 
   if (Bottom > 0)   
      *pMean = exp(Top/Bottom);  
   else
      *pMean = 0.0;   
                                  
   //equation (16) VPA User Guide                                
   for (iYear=MinYear, Top = Bottom = 0.0; iYear<= MaxYear; iYear++)
      if (N[iAge][iYear] > 0.0 && *pMean > 0.0)
         Top += Tau(iYear)*SQR(log(N[iAge][iYear]) - log(*pMean));
      
   if (Bottom-1 > 0)
      *pSE = sqrt(Top/(Bottom-1));
   else
      *pSE = 0.0;   
   }

void ExtendedSurvivorsAnalysis::RunRetro(short NYr, LP2DOUBLE Fecundity, LP2SHORT FBarRanges, short NumFBar, LP2SHORT NBarRanges, short NumNBar, LP2DOUBLE RetroSSB,LP2DOUBLE RetroR,LP3DOUBLE RetroNBar,LP3DOUBLE RetroFBar)
   {   
   short i, iAge, iRetro, _MaxYear = MaxYear, iYear;

   for (iRetro=1; iRetro <= NYr && MaxYear > MinYear + 10; MaxYear--, iRetro++)
      {
      Run();

      for (iYear=MinYear; iYear<=MaxYear; iYear++)
         {
         RetroR[iRetro][iYear] = N[MinAge][iYear];
         
         for (iAge = MinAge, RetroSSB[iRetro][iYear] = 0.0; iAge <= min(GetPlusGroup(),MaxAge+1); iAge++)
            RetroSSB[iRetro][iYear] += N[iAge][iYear]*Fecundity[iAge][iYear];

         for (i = 1; i <= NumNBar; i++)
            {
            for (iAge = max(GetPlusGroup(), NBarRanges[i][1]), RetroNBar[iRetro][iYear][i] = 0.0; iAge <= min(MaxAge, min(GetPlusGroup(), NBarRanges[i][2])); iAge++)
               RetroNBar[iRetro][iYear][i] += N[iAge][iYear];
            RetroNBar[iRetro][iYear][i] /= (NBarRanges[i][2]-NBarRanges[i][1]+1);
            }

         for (i = 1; i <= NumFBar; i++)
            {
            for (iAge = max(MinAge, FBarRanges[i][1]), RetroFBar[iRetro][iYear][i] = 0.0; iAge <= min(MaxAge, min(GetPlusGroup(), FBarRanges[i][2])); iAge++)
               RetroFBar[iRetro][iYear][i] += F[iAge][iYear];
            RetroFBar[iRetro][iYear][i] /= (FBarRanges[i][2]-FBarRanges[i][1]+1);
            }
         }
      }

   MaxYear = _MaxYear;
   }    

double ExtendedSurvivorsAnalysis::GetLogLikelihood(void)
   {
   short iFleet, iYear, iAge;

   DegreesOfFreedom = NObs = NParam = 0;
   LogLikelihood = 0.0;

      
   for (iFleet=0, pTune->IncrFleet(&iFleet); iFleet <= NTuneFleet; pTune->IncrFleet(&iFleet))
      {
      short  n = -1;
  
      double Temp = 0.0;
      for (iAge = GetMinTuneAge(iFleet); iAge < min(Control.ConQAge,GetMaxTuneAge(iFleet)); iAge++)
         {
         for (short iYear = MinTuneYear[iFleet]; iYear <= MaxTuneYear[iFleet]; iYear++)
            if (Theta(iFleet, iAge, iYear))
               {
               n++;
               Temp += SQR(log(GetCorrectedCPUE(iFleet, iAge, iYear)/N[iAge][iYear]) - log(Getq(iFleet,iAge)));
               }
      
         if (n>1)
            SetqSE(iFleet, iAge, sqrt(Temp/n));
         }
      
      for (iAge = GetMinTuneAge(iFleet); iAge <= min(Control.ConQAge,GetMaxTuneAge(iFleet)); iAge++)
         if (iAge <= Control.LRcrtAge)
            NParam += 2;
         else if (iAge < Control.ConQAge)
            NParam++;

      for (iYear = GetMinTuneYear(iFleet); iYear <= GetMaxTuneYear(iFleet); iYear++)
         for (iAge = GetMinTuneAge(iFleet); iAge < min(Control.ConQAge,GetMaxTuneAge(iFleet)); iAge++)
            if (Theta(iFleet, iAge, iYear))
               {
      			LogLikelihood -= SQR(GetResidual(iFleet,iAge,iYear))/(2.0*SQR(GetqSE(iFleet,iAge)))+log(GetqSE(iFleet,iAge));
               
               NObs++;
               }
      }

   DegreesOfFreedom = NObs - NParam;

   AIC = -2.0*LogLikelihood + 2.0 * NParam;

   return LogLikelihood;
   }

double ExtendedSurvivorsAnalysis::GetRSquared(void)
   {
   short iFleet, iYear, iAge;

   double   Top    = 0.0,
            Bottom = 0.0;

   for (iFleet=0, pTune->IncrFleet(&iFleet); iFleet<=NTuneFleet; pTune->IncrFleet(&iFleet))
      {
		double MnObs  = 0.0;

      for (iAge = GetMinTuneAge(iFleet); iAge <= GetMaxTuneAge(iFleet); iAge++)
         {
         for (iYear = GetMinTuneYear(iFleet); iYear <= GetMaxTuneYear(iFleet); iYear++)
            if (Theta(iFleet, iAge, iYear))
               MnObs += GetCorrectedCPUE(iFleet,iAge,iYear)*Tau(iYear);

         MnObs /= Tau(iYear);
         for (iYear = GetMinTuneYear(iFleet); iYear <= GetMaxTuneYear(iFleet); iYear++)
            if (Theta(iFleet, iAge, iYear))
               {
               Top    += SQR(GetResidual(iFleet,iAge,iYear));
               Bottom += SQR(GetCorrectedCPUE(iFleet,iAge,iYear)-MnObs);
               }
         }
      }

   return 1.0-Top/Bottom;
   }

bool ExtendedSurvivorsAnalysis::AllocArrays(void)
   {
   short iFleet;

   NTuneFleet = pTune->GetNFleet();   

   flallocArray(&MinTuneAge,  1, NTuneFleet);
   flallocArray(&MaxTuneAge,  1, NTuneFleet);
   flallocArray(&MinTuneYear, 1, NTuneFleet);
   flallocArray(&MaxTuneYear, 1, NTuneFleet);   
   flallocArray(&InternalSE,  0, NTuneFleet, MinYear-MaxAge, MaxYear-MinAge);   
   flallocArray(&ExternalSE,  0, NTuneFleet, MinYear-MaxAge, MaxYear-MinAge);   
   
   flallocArray(&M,  MinAge, MaxAge+1, MinYear, MaxYear);
   flallocArray(&N,  MinAge, MaxAge+1, MinYear, MaxYear+1);
   flallocArray(&F,  MinAge, MaxAge+1, MinYear, MaxYear+1);

   flallocArray(&PrevF, 1, 2, MinAge, MaxAge);
   for (short iAge = MinAge; iAge<= MaxAge; iAge++)
      {
      PrevF[1][iAge] = 1.0;
      PrevF[2][iAge] = 0.0;
      }

   SetTau();
   
   flallocArray(&TermPopSE,    MinYear-MaxAge, MaxYear-MinAge);

   flallocArray(&FShrinkMean,  MinYear-MaxAge, MaxYear-MinAge);
   
   if (Control.LRcrtAge >= MinAge)
      {
      flallocArray(&PopShrinkMean,MinAge, Control.LRcrtAge);
      flallocArray(&PopShrinkSE,  MinAge, Control.LRcrtAge);
      } 
   
   flallocArray(&CorrectedCPUE, 1, NTuneFleet);
   flallocArray(&NHat,          1, NTuneFleet);
   flallocArray(&NSE,           1, NTuneFleet);
   flallocArray(&ManualNSE,     1, NTuneFleet);
   flallocArray(&q,             1, NTuneFleet);
   flallocArray(&q_b,           1, NTuneFleet);
   flallocArray(&qSE,           1, NTuneFleet);
   flallocArray(&qn,            1 ,NTuneFleet);
                                       
   for (iFleet=0, pTune->IncrFleet(&iFleet); iFleet <= NTuneFleet; pTune->IncrFleet(&iFleet))
      {
      MinTuneAge[iFleet]  = max(MinAge, pTune->GetMinAge(iFleet));
      MaxTuneAge[iFleet]  = min(MaxAge, pTune->GetMaxAge(iFleet));
      MinTuneYear[iFleet] = max(MinYear,pTune->GetMinYear(iFleet));
      MaxTuneYear[iFleet] = min(MaxYear,pTune->GetMaxYear(iFleet));

      flallocArray(&NHat[iFleet],          MinTuneAge[iFleet], MaxTuneAge[iFleet], MinTuneYear[iFleet], MaxTuneYear[iFleet]);
      flallocArray(&NSE[iFleet],           MinTuneAge[iFleet], MaxTuneAge[iFleet], MinTuneYear[iFleet], MaxTuneYear[iFleet]);
      flallocArray(&ManualNSE[iFleet],  MinTuneAge[iFleet], MaxTuneAge[iFleet], MinTuneYear[iFleet], MaxTuneYear[iFleet]);
      flallocArray(&CorrectedCPUE[iFleet], MinTuneAge[iFleet], MaxTuneAge[iFleet], MinTuneYear[iFleet], MaxTuneYear[iFleet]);
      flallocArray(&q[iFleet],             MinTuneAge[iFleet], MaxTuneAge[iFleet], MinTuneYear[iFleet], MaxTuneYear[iFleet]);
      flallocArray(&q_b[iFleet],           MinTuneAge[iFleet], MaxTuneAge[iFleet], MinTuneYear[iFleet], MaxTuneYear[iFleet]);
      flallocArray(&qSE[iFleet],           MinTuneAge[iFleet], MaxTuneAge[iFleet], MinTuneYear[iFleet], MaxTuneYear[iFleet]);
      flallocArray(&qn[iFleet],            MinTuneAge[iFleet], MaxTuneAge[iFleet], MinTuneYear[iFleet], MaxTuneYear[iFleet]);
      }                                              
      
   flallocArray(&MeanF, MinAge, GetMaxCatchAge());      
   
   flallocArray(&ExtraVariable, MinAge, MaxAge, MinYear, MaxYear);
   
   for (int iYrCls = MinYear-MaxAge; iYrCls <= MaxYear-MinAge; iYrCls++)
      {
      TermPopSE[iYrCls]   =  0.0;
      SetFShrinkMean(iYrCls, 0.0); 
   
      for (iFleet=0; iFleet <= NTuneFleet; iFleet++)
         InternalSE[iFleet][iYrCls]  =
         ExternalSE[iFleet][iYrCls]  = 0.0;
      }
   
   flallocArray(&LogLikelihoodRecord,    1, Control.MaxIters);
   
   ArraysSet = TRUE;

   return TRUE;
   }                              

void ExtendedSurvivorsAnalysis::FreeAllocArrays(void)
   {
   if(!ArraysSet)
      return;
    
   short iFleet;
  
   for (iFleet=0, pTune->IncrFleet(&iFleet); iFleet <= NTuneFleet; pTune->IncrFleet(&iFleet))
      {
      free_flallocArray(NHat[iFleet],          MinTuneAge[iFleet], MaxTuneAge[iFleet], MinTuneYear[iFleet], MaxTuneYear[iFleet]);
      free_flallocArray(NSE[iFleet],           MinTuneAge[iFleet], MaxTuneAge[iFleet], MinTuneYear[iFleet], MaxTuneYear[iFleet]);
      free_flallocArray(ManualNSE[iFleet],     MinTuneAge[iFleet], MaxTuneAge[iFleet], MinTuneYear[iFleet], MaxTuneYear[iFleet]);
      free_flallocArray(CorrectedCPUE[iFleet], MinTuneAge[iFleet], MaxTuneAge[iFleet], MinTuneYear[iFleet], MaxTuneYear[iFleet]);
      free_flallocArray(q[iFleet],             MinTuneAge[iFleet], MaxTuneAge[iFleet], MinTuneYear[iFleet], MaxTuneYear[iFleet]);
      free_flallocArray(q_b[iFleet],           MinTuneAge[iFleet], MaxTuneAge[iFleet], MinTuneYear[iFleet], MaxTuneYear[iFleet]);
      free_flallocArray(qSE[iFleet],           MinTuneAge[iFleet], MaxTuneAge[iFleet], MinTuneYear[iFleet], MaxTuneYear[iFleet]);
      free_flallocArray(qn[iFleet],            MinTuneAge[iFleet], MaxTuneAge[iFleet], MinTuneYear[iFleet], MaxTuneYear[iFleet]);
      }

   free_flallocArray(NHat,          1, NTuneFleet);
   free_flallocArray(NSE,           1, NTuneFleet);   
   free_flallocArray(ManualNSE,     1, NTuneFleet);   
   free_flallocArray(CorrectedCPUE, 1, NTuneFleet);
   free_flallocArray(q,             1, NTuneFleet);
   free_flallocArray(q_b,           1, NTuneFleet);
   free_flallocArray(qSE,           1, NTuneFleet);
   free_flallocArray(qn,            1, NTuneFleet);

   free_flallocArray(PrevF, 1, 2, MinAge, MaxAge);
   
   free_flallocArray(TSWeight, MinYear, MaxYear);

   free_flallocArray(TermPopSE,    MinYear-MaxAge, MaxYear-MinAge);
   free_flallocArray(FShrinkMean,  MinYear-MaxAge, MaxYear-MinAge);
   if (Control.LRcrtAge >= MinAge)
      {
      free_flallocArray(PopShrinkMean,MinAge, Control.LRcrtAge);
      free_flallocArray(PopShrinkSE,  MinAge, Control.LRcrtAge);
      }
   
   free_flallocArray(MeanF,  MinAge, GetMaxCatchAge());
                     
   free_flallocArray(InternalSE,  0, NTuneFleet, MinYear-MaxAge, MaxYear-MinAge);   
   free_flallocArray(ExternalSE,  0, NTuneFleet, MinYear-MaxAge, MaxYear-MinAge);   
   free_flallocArray(MinTuneAge,  1, NTuneFleet);
   free_flallocArray(MaxTuneAge,  1, NTuneFleet);
   free_flallocArray(MinTuneYear, 1, NTuneFleet);
   free_flallocArray(MaxTuneYear, 1, NTuneFleet);   

   free_flallocArray(ExtraVariable, MinAge, MaxAge, MinYear, MaxYear);
   
   free_flallocArray(LogLikelihoodRecord,  1, Control.MaxIters);

   free_flallocArray(N,  MinAge, MaxAge+1, MinYear, MaxYear+1);
   free_flallocArray(F,  MinAge, MaxAge+1, MinYear, MaxYear+1);
   free_flallocArray(M,  MinAge, MaxAge+1, MinYear, MaxYear);  
   
   ArraysSet = FALSE;
   }   
   
ExtendedSurvivorsAnalysis::~ExtendedSurvivorsAnalysis()
   { 
   FreeAllocArrays();
   }

void ExtendedSurvivorsAnalysis::Project(double FMult, double TransitionRate, bool FlagFBelowBPA, double MinF, double BPA, double BLim, short NYrSRR, short MinFBarAge, short MaxFBarAge, LP2DOUBLE SWt, LP2DOUBLE Mat, LP2DOUBLE CWt, LPDOUBLE pFBar, LPDOUBLE pABCFBar, LPDOUBLE pSSB, LPDOUBLE pBiomass, LPDOUBLE pABC, LPDOUBLE pRecruits, bool RCT3Flag)
   {
	short iAge,
			iYear; 

   double  Z; 

   LPDOUBLE  Sel, _CWt, _SWt; 
   LP2DOUBLE  _N;
  
   short _MaxAge = MaxAge + (Control.PlusGroup == MaxAge+1 ? 1 : 0);
 
   flallocArray(&Sel,   MinAge, _MaxAge);
   flallocArray(&_SWt,  MinAge, _MaxAge);
   flallocArray(&_CWt,  MinAge, _MaxAge);
   flallocArray(&_N,    MinAge, _MaxAge, 1, 2);

   *pFBar     =  
   *pSSB      = 
   *pBiomass  =
   *pRecruits =
   *pABC      = 0.0; 
  
   for (iAge = MinAge; iAge <= _MaxAge; iAge++)
      {
      if (iAge >= MinFBarAge && iAge <= MaxFBarAge)
         *pFBar += F[iAge][MaxYear];

   	*pSSB     += N[iAge][MaxYear]*SWt[iAge][MaxYear]*Mat[iAge][MaxYear];
      *pBiomass += N[iAge][MaxYear]*SWt[iAge][MaxYear];

      Sel[iAge] = _SWt[iAge] = _CWt[iAge] = _N[iAge][1] = _N[iAge][2] = 0.0;

      for (iYear = MaxYear-2; iYear <= MaxYear; iYear++)
         {
         Sel[iAge]  += F[  iAge][iYear];
         _SWt[iAge] += SWt[iAge][iYear];
         _CWt[iAge] += CWt[iAge][iYear];
         }

      Sel[iAge]  /= 3;
      _SWt[iAge] /= 3;
      _CWt[iAge] /= 3;
      }
   *pFBar /= (MaxFBarAge-MinFBarAge+1);

   double Scaling = 0.0;
   for (iAge = MinFBarAge; iAge <= MaxFBarAge; iAge++)
      Scaling += Sel[iAge];

   Scaling /= (MaxFBarAge-MinFBarAge+1);
   for (iAge = MinAge; iAge <= _MaxAge; iAge++)
      Sel[iAge] /= Scaling;
      
    //Geo mean recruits
   *pRecruits = 0.0;
   int NObs = 0;
   for (iYear = min(max(MaxYear-NYrSRR+1,MinYear),MaxYear); iYear <= MaxYear; iYear++)
   if (N[MinAge][iYear]>0)
      {
      *pRecruits += log(N[MinAge][iYear]);
      NObs++;
      }

   *pRecruits = exp(*pRecruits/NObs);
   _N[MinAge][2] = *pRecruits;
     
   if (RCT3Flag)
      *pRecruits    = _N[MinAge][1] = N[MinAge][MaxYear+1];
   else
      _N[MinAge][1] = *pRecruits;
   
   //N at start of current year
   for (iAge=MinAge+1; iAge<=_MaxAge; iAge++)
	   _N[iAge][1]  = N[iAge-1][MaxYear]*exp(-F[iAge-1][MaxYear]-M[iAge-1][MaxYear]);

   if (_MaxAge == MaxAge+1)
      _N[_MaxAge][1] += N[_MaxAge][MaxYear]*exp(-F[_MaxAge][MaxYear]-M[_MaxAge][MaxYear]);
      
   //N at start of TAC year
   for (iAge=MinAge+1; iAge<=_MaxAge; iAge++)
	   _N[iAge][2]  = _N[iAge-1][1]*exp(-Sel[iAge-1]*(*pFBar)-M[iAge-1][MaxYear]);

   if (_MaxAge == MaxAge+1)
      _N[_MaxAge][2] += _N[_MaxAge][1]*exp(-Sel[_MaxAge]*(*pFBar)-M[_MaxAge][MaxYear]);
   

double t;
t = max(t,0.0);

   //ABC
   if (TransitionRate > 0.0)
      (*pABCFBar) =  max(FMult,(1.0-TransitionRate)*(*pFBar));
   else
      (*pABCFBar) = FMult;

   if (*pSSB < BPA && FlagFBelowBPA==1) 
      (*pABCFBar) *= (*pSSB - BLim)/(BPA - BLim);
   else if (*pSSB < BPA && FlagFBelowBPA < 0 && (*pABCFBar > *pFBar))
      (*pABCFBar) *= (*pSSB - BLim)/(BPA - BLim);

   (*pABCFBar) = max(MinF,(*pABCFBar));

	for (iAge=MinAge; iAge<=_MaxAge; iAge++)
		{
		Z      = Sel[iAge]*(*pABCFBar)+M[iAge][MaxYear];  
     
		*pABC += _N[iAge][2]*Sel[iAge]*(*pABCFBar)/Z*(1.0-exp(-Z))*_CWt[iAge];
      }

   free_flallocArray(Sel,   MinAge, _MaxAge);
	free_flallocArray(_SWt,  MinAge, _MaxAge);
	free_flallocArray(_CWt,  MinAge, _MaxAge);
	free_flallocArray(_N,    MinAge, _MaxAge, 1, 2);
   }

bool ExtendedSurvivorsAnalysis::CheckRange(short _MinAge, short _MaxAge, short _MinYear, short _MaxYear)
   {
   return (_MinAge != MinAge || _MaxAge != MaxAge || _MinYear != MinYear || _MaxYear != MaxYear);
   }
