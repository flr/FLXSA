#include "fl__vpa.hpp"
#include <math.h>

VirtualPopAnalysis::VirtualPopAnalysis(void)     
   {
   Error     = !FiFiErrNull;
   ArraysSet = FALSE;  

   FRatio     = 1.0;
   FRatioFlag = FALSE;
   }

VirtualPopAnalysis::VirtualPopAnalysis(LPSHORT pAge, LPSHORT pYear, LPDOUBLE pM, LPDOUBLE pF, LPSHORT pn, LPSHORT pPlusGroup, LPDOUBLE pFRatio, CatchData *pCatchData)
   {
   Error     = !FiFiErrNull;
   ArraysSet = FALSE;  

   pCatch = pCatchData;

   if (*pFRatio <= 0)
      FRatioFlag = FALSE;
   else
      {
      FRatio = *pFRatio; 
      FRatioFlag = TRUE;
      } 
   
   SetRange(pAge,pYear,*pn);

   //Check array range against catch array range
   CheckRangeWithCatch(*pPlusGroup);
   
   if (!AllocArrays())
      return;           

   //Initialise Arrays in case missing values in input vectors
   InitialiseMatrix(F);
   InitialiseMatrix(M);
   
   //Input Data
   for (short i = 0; i < *pn; i++)
      {
      SetF((short)pAge[i],(short)pYear[i],pF[i]);
      SetM((short)pAge[i],(short)pYear[i],pM[i]);
      }            

   SetTerminalNs();   

   Error = FiFiErrNull;
   }

void VirtualPopAnalysis::SetRange(LPSHORT Age,LPSHORT Year,short n)
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

void VirtualPopAnalysis::InitialiseMatrix(LP2DOUBLE Value)
   {
   for (short iAge=MinAge; iAge<=MaxAge; iAge++)
      for (short iYear=MinYear; iYear<=MaxYear; iYear++)
//          if (!IsBadReadPtr(&Value[iAge][iYear], sizeof(short)))
              Value[iAge][iYear] = 0.0;
   }   

void VirtualPopAnalysis::CheckRangeWithCatch(short PlusGroup)
   {   
   MinAge  = max(MinAge, pCatch->GetMinAge());
   if (PlusGroup > MinAge+2)
      MaxAge  = min(pCatch->GetMaxAge(),PlusGroup-1);
   MinYear = max(MinYear, pCatch->GetMinYear());
   MaxYear = min(MaxYear, pCatch->GetMaxYear());
   }
 
void VirtualPopAnalysis::SetM(short iAge, double Value)  
   {                                                 
   for (short iYear = MinYear; iYear <= MaxYear; iYear++)
      M[min(MaxAge,max(MinAge,iAge))][iYear] = Value;
   }

void VirtualPopAnalysis::SetM(short iAge, short iYear, double Value)  
   {                                                 
   M[min(MaxAge,max(MinAge,iAge))][min(MaxYear,max(MinYear,iYear))] = Value;
   }

void VirtualPopAnalysis::SetF(short iAge, short iYear, double Value)  
   {                                                 
   F[min(MaxAge,max(MinAge,iAge))][min(MaxYear,max(MinYear,iYear))] = Value;
   }

void VirtualPopAnalysis::SetTerminalNs(void)  
   {                                                                                     
   double Z;
   
   //set Ns at end of last year
   for (short iAge=MinAge; iAge<=MaxAge; iAge++)
      {
      Z = F[iAge][MaxYear] + M[iAge][MaxYear];
      
      N[iAge+1][MaxYear+1] = pCatch->GetCatch(iAge, MaxYear)*Z/(F[iAge][MaxYear]*(exp(Z)-1));
      }
   
   //set Ns at end of oldest age
   if (!FRatioFlag)
      for (short iYear=MinYear; iYear<MaxYear; iYear++)
         {
         Z = F[MaxAge][iYear] + M[MaxAge][iYear];
      
         N[MaxAge+1][iYear+1] = pCatch->GetCatch(MaxAge, iYear)*Z/(F[MaxAge][iYear]*(exp(Z)-1));
         }
   }
                                 
short VirtualPopAnalysis::GetMinAge(void)    
   {             
   return MinAge; 
   } 
   
short VirtualPopAnalysis::GetMaxAge(void)    
   { 
   return MaxAge; 
   }
   
short VirtualPopAnalysis::GetMinYear(void) 
   { 
   return MinYear; 
   }
   
short VirtualPopAnalysis::GetMaxYear(void)  
   { 
   return MaxYear; 
   } 

double VirtualPopAnalysis::GetN(short iAge, short iYear)
   {
   if (iAge < MinAge || iAge > MaxAge+1 || iYear < MinYear ||  iYear > MaxYear)
      return 0.0;
      
   return max(0.0,N[iAge][iYear]);
   }
   
double VirtualPopAnalysis::GetF(short iAge, short iYear)
   {
   return max(0.0,F[min(max(iAge,MinAge),MaxAge)][min(max(iYear,MinYear),MaxYear)]);
   }

double VirtualPopAnalysis::GetM(short iAge)
   {
   return max(0.0,M[min(max(iAge,MinAge),MaxAge)][MinYear]);
   }
   
double VirtualPopAnalysis::GetM(short iAge, short iYear)
   {
   return max(0.0,M[min(max(iAge,MinAge),MaxAge)][min(max(iYear,MinYear),MaxYear)]);
   }

bool  VirtualPopAnalysis::RunCohortAnalysis(void)
   {
   short iAge, iYearClass;

   for (iYearClass = MaxYear - MinAge; iYearClass >= MinYear - MaxAge; iYearClass--)
      for (iAge  = min(MaxAge, MaxYear - iYearClass); iAge >= max(MinAge, MinYear - iYearClass); iAge--)
         if (FRatioFlag && iAge == MaxAge && iYearClass < MaxYear - MaxAge)
            {                                   
            F[iAge][iYearClass+iAge] = max(0.0,FRatio*F[iAge-1][iYearClass+iAge]);
              
            N[iAge][iYearClass+iAge] = max(0.0,pCatch->GetCatch(iAge,iYearClass+iAge)
                                   /(F[iAge][iYearClass+iAge]/(F[iAge][iYearClass+iAge] + GetM(iAge, iYearClass+iAge))
                                   *(1.0 - exp(-F[iAge][iYearClass+iAge] - GetM(iAge, iYearClass+iAge)))));   
        
            }
         else
            {
            double t;
            t = N[iAge+1][iYearClass+iAge+1]*exp(GetM(iAge, iYearClass+iAge));
            N[iAge][iYearClass+iAge] = max(0.0, t + pCatch->GetCatch(iAge, iYearClass+iAge) * exp(0.5*GetM(iAge, iYearClass+iAge))); 
            
            if (N[iAge+1][iYearClass+iAge+1] > 0.0)
               F[iAge][iYearClass+iAge] = max(0.0,log(N[iAge][iYearClass+iAge]/N[iAge+1][iYearClass+iAge+1]) - GetM(iAge, iYearClass+iAge));
            else
               F[iAge][iYearClass+iAge] = 0.0;
		      } 
         
   CalcPlusGroupFN();

   return (TRUE);
   } 

bool VirtualPopAnalysis::RunVPA(void)
   {
   short iAge, Iters; 
   short  iYearClass;  
   double f, dfdx;         

   for (iYearClass = MaxYear - MinAge; iYearClass >= MinYear - MaxAge; iYearClass--)
      for (iAge = min(MaxAge, MaxYear - iYearClass); iAge >= max(MinAge, MinYear - iYearClass); iAge--)
         if (FRatioFlag && iAge == MaxAge && iYearClass < MaxYear - MaxAge)
            {                                   
            F[iAge][iYearClass+iAge] = max(0.0,FRatio*F[iAge-1][iYearClass+iAge]);
//            N[iAge][iYearClass+iAge] = max(0.0,pCatch->GetCatch(iAge,iYearClass+iAge)
//                                       /(F[iAge][iYearClass+iAge]/(F[iAge][iYearClass+iAge] + GetM(iAge, iYearClass+iAge))
//                                       *(1.0 - exp(-F[iAge][iYearClass+iAge] - GetM(iAge, iYearClass+iAge)))));   
            N[iAge][iYearClass+iAge] = max(0.0,  pCatch->GetCatch(iAge,iYearClass+iAge)*(F[iAge][iYearClass+iAge] + GetM(iAge, iYearClass+iAge))
                                               /(F[iAge][iYearClass+iAge]*(1.0 - exp(-F[iAge][iYearClass+iAge] - GetM(iAge, iYearClass+iAge)))));   
            }
         else
            {
            if (pCatch->GetCatch(iAge, iYearClass+iAge) <= 0.0)
               {
               N[iAge][iYearClass+iAge] = max(0.0,N[iAge+1][iYearClass+iAge+1]*exp(GetM(iAge, iYearClass+iAge)));
               F[iAge][iYearClass+iAge] = 0.0;
               }
            else
               {
               Iters = 0;  
               
               double t;
               
               t = N[iAge+1][iYearClass+iAge+1]*exp(GetM(iAge, iYearClass+iAge));
               N[iAge][iYearClass+iAge] = max(0.0, t + pCatch->GetCatch(iAge, iYearClass+iAge) * exp(0.5*GetM(iAge, iYearClass+iAge)));
               
               do
                  {
                  Iters++;                                           
                  //do Newton Raphson to estimate N
                  f    = Calcf(GetM(iAge, iYearClass+iAge), pCatch->GetCatch(iAge, iYearClass+iAge), 
                               N[iAge][iYearClass+iAge], N[iAge+1][iYearClass+iAge+1]);
               
                  dfdx = Calcdfdx(GetM(iAge, iYearClass+iAge),
                               N[iAge][iYearClass+iAge], N[iAge+1][iYearClass+iAge+1]);
                  //calc N        
                  N[iAge][iYearClass+iAge] = NewtonRhapson(N[iAge][iYearClass+iAge], f, dfdx);
                  }
               while (fabs(f) >= VPA_TOL && Iters <= VPA_ITS);    
               
               //calc F at iAge
               F[iAge][iYearClass+iAge] =  max(0.0,-log(N[iAge+1][iYearClass+iAge+1]/N[iAge][iYearClass+iAge]) - GetM(iAge, iYearClass+iAge));
               }
            }
            
   CalcPlusGroupFN();

   return TRUE;   
   }   
                                       
double VirtualPopAnalysis::Calcf(double M, double Catch, double N, double N1)
   {                                            
   return (Catch - (1 - M/(log(N)-log(N1)))*(N-N1)); 
   }
   
double VirtualPopAnalysis::Calcdfdx(double M, double N, double N1)
   {       
   return  (-1 -((log(N)-log(N1))*M -M*(N-N1)/N)/log(2*N/N1));
   } 
   
inline double VirtualPopAnalysis::CalcY(double M, double F, double Catch, double N)
   {                                            
   return (exp(-M-F)*Catch/N - F/(F+M)*(1-exp(-F-M))); 
   }
   
inline double VirtualPopAnalysis::CalcdYdF(double M, double F, double N, double N1)
   {       
   return (-exp(-F)*exp(-1.5f*M)*N/N1);   
   } 


double VirtualPopAnalysis::NewtonRhapson(double x, double f, double dfdx)
   {
   return (x - f/dfdx);
   }  
                
void VirtualPopAnalysis::CalcPlusGroupFN(void)
   {
   for (short iYear=MinYear; iYear<= MaxYear; iYear++)
      {
      F[MaxAge+1][iYear] = F[MaxAge][iYear];
      N[MaxAge+1][iYear] = 0.0f;
      for (short iAge = MaxAge+1; iAge<=pCatch->GetMaxAge(); iAge++)
          {
          double Z = F[MaxAge][iYear]+M[MaxAge][iYear];
          N[MaxAge+1][iYear] += pCatch->GetCatch(iAge, iYear)*Z/(F[MaxAge][iYear]*(1-exp(-Z)));
          }
      }
   }
      
bool VirtualPopAnalysis::AllocArrays(void)
   {                           
   if (ArraysSet)
      FreeAllocArrays();  
   
   if (!flallocArray(&N, MinAge, MaxAge+1, MinYear, MaxYear+1))
      return FALSE;
   if (!flallocArray(&F, MinAge, MaxAge+1, MinYear, MaxYear+1))
      {
      free_flallocArray(N, MinAge, MaxAge+1, MinYear, MaxYear+1);
      return FALSE;
      }
   if (!flallocArray(&M, MinAge, MaxAge+1, MinYear, MaxYear))
      {
      free_flallocArray(N, MinAge, MaxAge+1, MinYear, MaxYear+1);
      free_flallocArray(F, MinAge, MaxAge+1, MinYear, MaxYear+1);
      return FALSE;
      } 
      
   for (short iAge = MinAge; iAge <= MaxAge+1; iAge++)
      for (short iYear = MinYear; iYear <= MaxYear+1; iYear++)
         N[iAge][iYear] = 
         F[iAge][iYear] = 0.0;   

   ArraysSet = TRUE;
   
   return TRUE;
   }   

void VirtualPopAnalysis::FreeAllocArrays(void)
   {                                    
   if (!ArraysSet)
      return;

   free_flallocArray(N,  MinAge, MaxAge+1, MinYear, MaxYear+1);
   free_flallocArray(F,  MinAge, MaxAge+1, MinYear, MaxYear+1);
   free_flallocArray(M,  MinAge, MaxAge+1, MinYear, MaxYear);  

   ArraysSet = FALSE;
   }                 
                
VirtualPopAnalysis::~VirtualPopAnalysis(void)
   {                    
   FreeAllocArrays();
   }   
   

/*
Option Explicit

Declare Function VBSetAY Lib "flvb32.DLL" _
   (Ref As Range, AY() As Double) As Integer

Declare Function VBSPR Lib "flvb32.DLL" _
   (OverallF As Double, F() As Double, M() As Double, WtAtAge() As Double, Mat() As Double, PFSpawn As Double, PMSpawn As Double) As Double
   
Declare Function VBVPA Lib "flvb32.DLL" _
   (Catch() As Double, NTermYear() As Double, NTermAge() As Double, M() As Double, N() As Double, F() As Double) As Integer

Declare Function VBSepVPA Lib "flvb32.DLL" _
    (Catch() As Double, FinalF As Double, OldestSel As Double, RefAge As Integer, _
     RefSel As Double, M() As Double, N() As Double, F() As Double, CatchHat() As Double) As Integer

Function EricVPA(Catch As Object, M As Object, EricYear As Integer, FinalF As Double, OldestSel As Double, RefAge As Integer, RefSel As Double, Optional LastTrueAge As Variant) As Variant
   Dim iAge As Integer
   Dim iYear As Integer
   Dim MinAge As Integer
   Dim MaxAge As Integer
   Dim MinYear As Integer
   Dim MaxYear As Integer
   
   Dim CatchMatrix() As Double
   Dim CatchHat() As Double
   Dim MMatrix() As Double
   Dim N() As Double
   Dim F() As Double
   
   Dim EricCatch() As Double
   Dim EricM() As Double
   Dim EricN() As Double
   Dim EricF() As Double
   
   Dim VPACatch() As Double
   Dim VPAM() As Double
   Dim VPAN() As Double
   Dim VPAF() As Double
      
   Dim TermYear() As Double
   Dim TermAgeFRatio() As Double
      
   VBSetAY Catch, CatchMatrix
   VBSetAY M, MMatrix
      
   MinAge = LBound(CatchMatrix, 1)
   MaxAge = UBound(CatchMatrix, 1)
   MinYear = LBound(CatchMatrix, 2)
   MaxYear = UBound(CatchMatrix, 2)
   
   If IsMissing(LastTrueAge) Then
      LastTrueAge = MaxAge
   End If
   
   ReDim EricCatch(MinAge To LastTrueAge, EricYear To MaxYear) As Double
   ReDim EricM(MinAge To LastTrueAge, EricYear To MaxYear) As Double
   ReDim EricN(MinAge To LastTrueAge, EricYear To MaxYear) As Double
   ReDim EricF(MinAge To LastTrueAge, EricYear To MaxYear) As Double
   ReDim VPACatch(MinAge To LastTrueAge, MinYear To EricYear - 1) As Double
   ReDim VPAM(MinAge To LastTrueAge, MinYear To EricYear - 1) As Double
   ReDim VPAN(MinAge To LastTrueAge, MinYear To EricYear - 1) As Double
   ReDim VPAF(MinAge To LastTrueAge, MinYear To EricYear - 1) As Double
   ReDim CatchHat(MinAge To LastTrueAge, EricYear To MaxYear) As Double
    
   For iAge = MinAge To LastTrueAge
      For iYear = MinYear To EricYear - 1
         VPACatch(iAge, iYear) = CatchMatrix(iAge, iYear)
         VPAM(iAge, iYear) = MMatrix(iAge, iYear)
      Next iYear
      For iYear = EricYear To MaxYear
         EricCatch(iAge, iYear) = CatchMatrix(iAge, iYear)
         EricM(iAge, iYear) = MMatrix(iAge, iYear)
      Next iYear
   Next iAge
   
   VBSepVPA EricCatch, FinalF, OldestSel, RefAge, RefSel, EricM, EricN, EricF, CatchHat

   ReDim TermYear(MinAge + 1 To Application.Min(MaxAge, LastTrueAge) + 1) As Double
   ReDim TermAgeFRatio(1 To 1) As Double
   For iAge = MinAge + 1 To LastTrueAge
      TermYear(iAge) = EricN(iAge, EricYear)
   Next
   
   TermYear(LastTrueAge + 1) = CatchMatrix(LastTrueAge, EricYear - 1)
   TermAgeFRatio(1) = 1#
   
   'do VPA to get F in oldest age
   VBVPA VPACatch, TermYear, TermAgeFRatio, VPAM, VPAN, VPAF
   
   Dim TempF As Double
   Dim Z As Double
   TempF = TermAgeFRatio(1) * VPAF(LastTrueAge - 1, EricYear - 1)
   Z = TempF + MMatrix(LastTrueAge, EricYear - 1)
   TermYear(LastTrueAge + 1) = CatchMatrix(LastTrueAge, EricYear - 1) * Exp(-Z) * Z / (TempF * (1 - Exp(-Z)))
   
   VBVPA VPACatch, TermYear, TermAgeFRatio, VPAM, VPAN, VPAF
   
   EricVPA = ReturnEric(EricN, EricF, CatchHat, VPAN, VPAF, CatchMatrix, MMatrix)

End Function

'Creates Excel output
Function ReturnEric(EricN() As Double, EricF() As Double, EricCatchHat() As Double, VPAN() As Double, VPAF() As Double, Catch() As Double, M() As Double) As Variant
    Dim ReturnArray() As Variant
    Dim PlusGroup As Boolean
    
    'Output InputArray
    Dim MinAge As Integer
    Dim MaxAge As Integer
    Dim MinYear As Integer
    Dim MaxYear As Integer
    Dim NAges As Integer
    Dim NYears As Integer
    Dim iAge As Integer
    Dim iYear As Integer
    Dim PlusGroupCatch As Double
    Dim Z As Double
    
    MinAge = LBound(VPAN, 1)
    If UBound(Catch, 1) > UBound(EricN, 1) Then
      MaxAge = UBound(EricN, 1) + 1
      PlusGroup = True
    Else
      MaxAge = UBound(EricN, 1)
      PlusGroup = False
    End If
    
    MinYear = LBound(VPAN, 2)
    MaxYear = UBound(EricN, 2)
    NAges = MaxAge - MinAge + 1
    
    ReDim ReturnArray(0 To NAges * 3 + 5, 0 To MaxYear - MinYear + 1) As Variant
    
    'Age Headings
    For iAge = MinAge To MaxAge
        ReturnArray(iAge, 0) = iAge
        ReturnArray(iAge - MinAge + 3 + NAges, 0) = iAge
        ReturnArray(iAge - MinAge + NAges * 2 + 5, LBound(EricN, 2) - MinYear) = iAge
    Next
    
    'Year Headings
    For iYear = MinYear To MaxYear
        ReturnArray(0, iYear - MinYear + 1) = iYear
        ReturnArray(NAges + 2, iYear - MinYear + 1) = iYear
        If iYear >= LBound(EricN, 2) Then _
           ReturnArray(NAges * 2 + 4, iYear - MinYear + 1) = iYear
    Next
    
    'Titles
    ReturnArray(0, 0) = "N"
    ReturnArray(NAges + 2, 0) = "F"
    ReturnArray(NAges * 2 + 4, LBound(EricN, 2) - MinYear) = "CatchHat"
    
    'Blanks
    For iYear = MinYear To MaxYear + 1
        ReturnArray(NAges + 1, iYear - MinYear) = ""
        ReturnArray(NAges * 2 + 3, iYear - MinYear) = ""
    Next
    For iAge = MinAge - 1 To MaxAge
       For iYear = MinYear To UBound(VPAN, 2)
         ReturnArray(iAge - MinAge + NAges * 2 + 5, iYear - MinYear) = ""
       Next iYear
    Next

    'Data
    For iAge = MinAge To UBound(VPAN, 1)
       For iYear = MinYear To UBound(VPAN, 2)
         ReturnArray(iAge - MinAge + 1, iYear - MinYear + 1) = VPAN(iAge, iYear)
         ReturnArray(iAge - MinAge + NAges + 3, iYear - MinYear + 1) = VPAF(iAge, iYear)
       Next iYear
       For iYear = LBound(EricN, 2) To MaxYear
         ReturnArray(iAge - MinAge + 1, iYear - MinYear + 1) = EricN(iAge, iYear)
         ReturnArray(iAge - MinAge + NAges + 3, iYear - MinYear + 1) = EricF(iAge, iYear)
         ReturnArray(iAge - MinAge + NAges * 2 + 5, iYear - MinYear + 1) = EricCatchHat(iAge, iYear)
       Next iYear
    Next
    
    'Plus group
    If PlusGroup Then
      For iYear = MinYear To UBound(VPAN, 2)
         ReturnArray(MaxAge - MinAge + NAges + 3, iYear - MinYear + 1) = VPAF(MaxAge - 1, iYear)
         ReturnArray(MaxAge - MinAge + 1, iYear - MinYear + 1) = VPAF(MaxAge - 1, iYear)
         PlusGroupCatch = 0
         For iAge = UBound(VPAN, 1) + 1 To UBound(Catch, 1)
            PlusGroupCatch = PlusGroupCatch + Catch(iAge, iYear)
         Next
         Z = VPAF(MaxAge - 1, iYear) + M(MaxAge - 1, iYear)
         ReturnArray(MaxAge - MinAge + 1, iYear - MinYear + 1) = Z / (VPAF(MaxAge - 1, iYear) * (1 - Exp(-Z))) * PlusGroupCatch
      Next iYear
      For iYear = LBound(EricN, 2) To MaxYear
         ReturnArray(MaxAge - MinAge + 1, iYear - MinYear + 1) = 0#
         ReturnArray(MaxAge - MinAge + NAges + 3, iYear - MinYear + 1) = EricF(MaxAge - 1, iYear)
         PlusGroupCatch = 0
         For iAge = UBound(VPAN, 1) + 1 To UBound(Catch, 1)
            PlusGroupCatch = PlusGroupCatch + Catch(iAge, iYear)
         Next
         Z = EricF(MaxAge - 1, iYear) + M(MaxAge - 1, iYear)
         ReturnArray(MaxAge - MinAge + 1, iYear - MinYear + 1) = Z / (EricF(MaxAge - 1, iYear) * (1 - Exp(-Z))) * PlusGroupCatch
      Next iYear
   End If
    
    ReturnEric = ReturnArray
    
End Function

Sub RunVPA(N As Double, F As Double, Nt1 As Double, Catch As Double, M As Double)
   Dim Iters As Integer
   Dim g As Double
   Dim dgdx As Double
   
   Iters = 0
   Do
      Iters = Iters + 1
      'do Newton Raphson to estimate N
      g = Calcf(M, Catch, N, Nt1)
               
      dfdx = Calcdfdx(M, N, Nt1)
      
      'calc N
       N = NewtonRhapson(N, g, dgdx)
                 
   While fabs(g) >= 0.0000000001 And Iters <= 25
               
   'calc F at iAge
   F = Application.Max(0#, -Log(Nt1 / N) - M)
 

End Sub
Function ForwardVPABevHolt(Catch As Object, M As Object, SWt As Object, Mat As Object, SelVector As Object, MVector As Object, SWtVector As Object, MatVector As Object, Steepness As Double, VirginSSB As Double, InitialFMult As Double) As Variant
   On Error Resume Next
   
   Dim i As Integer
   Dim iAge As Integer
   Dim iYear As Integer
   Dim MinAge As Integer
   Dim MaxAge As Integer
   Dim MinYear As Integer
   Dim MaxYear As Integer
   
   Dim CatchMatrix() As Double
   Dim MMatrix() As Double
   Dim SWtMatrix() As Double
   Dim MatMatrix() As Double
   Dim NMatrix() As Double
   Dim FMatrix() As Double
   
   Dim MSteadyState() As Double
   Dim SWtSteadyState() As Double
   Dim MatSteadyState() As Double
   Dim SelSteadyState() As Double
   
   Dim Z As Double
   Dim SSB As Double
   Dim SPRF0 As Double
   Dim DeltaF As Double
      
   VBSetAY Catch, CatchMatrix
   VBSetAY M, MMatrix
   VBSetAY SWt, SWtMatrix
   VBSetAY Mat, MatMatrix
      
   SetVector SelVector, SelSteadyState
   SetVector MVector, MSteadyState
   SetVector SWtVector, SWtSteadyState
   SetVector MatVector, MatSteadyState
   
   SPRF0 = VBSPR(0#, SelSteadyState, MSteadyState, SWtSteadyState, MatSteadyState, 0#, 0#)
   
   MinAge = LBound(CatchMatrix, 1)
   MaxAge = UBound(CatchMatrix, 1)
   MinYear = LBound(CatchMatrix, 2)
   MaxYear = UBound(CatchMatrix, 2)
   
   ReDim NMatrix(MinAge To MaxAge, MinYear To MaxYear) As Double
   ReDim FMatrix(MinAge To MaxAge, MinYear To MaxYear) As Double
   
   NMatrix(MinAge, MinYear) = BevHoltSPR(VBSPR(InitialFMult, SelSteadyState, MSteadyState, SWtSteadyState, MatSteadyState, 0#, 0#), Steepness, VirginSSB, SPRF0)
   For iAge = MinAge + 1 To MaxAge
      NMatrix(iAge, MinYear) = NMatrix(iAge - 1, MinYear) * Exp(-SelSteadyState(iAge - 1) * InitialFMult - MSteadyState(iAge - 1))
   Next
   NMatrix(MaxAge, MinYear) = NMatrix(MaxAge, MinYear) + NMatrix(MaxAge, MinYear) * -Exp(-SelSteadyState(MaxAge) * InitialFMult - MSteadyState(MaxAge)) / (Exp(-SelSteadyState(MaxAge) * InitialFMult - MSteadyState(MaxAge)) - 1#)
   
   For iYear = MinYear + 1 To MinYear + MinAge
      NMatrix(MinAge, iYear) = NMatrix(iAge, MinYear)
   Next
   
   SSB = 0#
   For iAge = MinAge To MaxAge
      SSB = SSB + NMatrix(iAge, MinYear) * SWtMatrix(iAge, MinYear) * MatMatrix(iAge, MinYear)
   Next
   NMatrix(MinAge, MinYear + MinAge) = BevHolt(SSB, Steepness, VirginSSB, SPRF0)
   
   For iYear = MinYear + MinAge - 1 To MaxYear - 1
      SSB = 0#
      For iAge = MinAge To MaxAge - 1
        NMatrix(iAge + 1, iYear + 1) = NMatrix(iAge, iYear) - CatchMatrix(iAge, iYear) * Exp(0.5 * MMatrix(iAge, iYear))
        i = 0
        Do
            DeltaF = g(NMatrix(iAge + 1, iYear + 1), NMatrix(iAge, iYear), MMatrix(iAge, iYear), CatchMatrix(iAge, iYear))
            NMatrix(iAge + 1, iYear + 1) = NMatrix(iAge + 1, iYear + 1) - DeltaF
            i = i + 1
        Loop While i < 100 And Abs(DeltaF) > 0.0000000001
        NMatrix(iAge + 1, iYear + 1) = Application.Max(0#, NMatrix(iAge + 1, iYear + 1))
        FMatrix(iAge, iYear) = -Log(NMatrix(iAge + 1, iYear + 1) / NMatrix(iAge, iYear)) - MMatrix(iAge, iYear)
      
        SSB = SSB + NMatrix(MaxAge, iYear + MinAge) * SWtMatrix(MaxAge, iYear) * MatMatrix(MaxAge, iYear)
      Next
     
      NMatrix(MinAge, iYear + MinAge + 1) = BevHolt(SSB, Steepness, VirginSSB, SPRF0)
    
   Next
  
  ForwardVPABevHolt = ReturnNandF(NMatrix, FMatrix)
   
End Function
Function g(Nt1 As Double, Nt As Double, M As Double, Catch As Double) As Double
   If Nt > 0 Then
      g = Catch - (1# - M / (Log(Nt) - Log(Nt1))) * (Nt - Nt1)
   End If
   g = g
End Function

Function dGdF(Nt1 As Double, Nt As Double, M As Double) As Double
   Dim Expr3 As Double
   Dim Expr5 As Double
   Dim Expr6 As Double
   Dim Expr10 As Double
   
   If Nt And Nt1 > 0 Then
      Expr3 = (Log(Nt)) - (Log(Nt1))
      Expr5 = 1 - (M / Expr3)
      Expr6 = Nt - Nt1
      dGdF = (((M * (1 / Nt)) / (Expr3 ^ 2)) * Expr6) + Expr5
   End If
   
End Function


Sub SetVector(InputRef As Range, OutputArray() As Double)
    Dim i As Integer
    Dim N As Integer
    
    N = InputRef.Cells.Count
    
    ReDim OutputArray(1 To N)
    For i = 1 To N
       OutputArray(i) = InputRef.Cells(i).Value
    Next
End Sub

'Creates Excel output
Function ReturnNandF(N() As Double, F() As Double) As Variant
    Dim NandF() As Variant

    'Output InputArray
    Dim MinAge As Integer
    Dim MaxAge As Integer
    Dim MinYear As Integer
    Dim MaxYear As Integer
    Dim NAges As Integer
    Dim NYears As Integer
    Dim iAge As Integer
    Dim iYear As Integer
    Dim iRow As Integer
    Dim iCol As Integer
    
    MinAge = LBound(N, 1)
    MaxAge = UBound(N, 1)
    MinYear = LBound(N, 2)
    MaxYear = UBound(N, 2)
    NAges = MaxAge - MinAge + 1
    NYears = MaxYear - MinYear + 1
    
    ReDim NandF(0 To NAges * 2 + 2, 0 To NYears) As Variant
    
    'Titles
    NandF(0, 0) = "N"
    NandF(NAges + 1, 0) = ""
    NandF(NAges + 2, 0) = "F"
   'Age Headings
    For iAge = 1 To NAges
        NandF(iAge, 0) = iAge - 1 + MinAge
        NandF(iAge + NAges + 2, 0) = iAge - 1 + MinAge
    Next
    'Year Headings
    For iYear = 1 To NYears
        NandF(0, iYear) = iYear - 1 + MinYear
        NandF(NAges + 1, iYear) = ""
        NandF(NAges + 2, iYear) = iYear - 1 + MinYear
    Next
    'Data
    For iAge = 1 To NAges
        For iYear = 1 To NYears
            NandF(iAge, iYear) = N(MinAge + iAge - 1, MinYear + iYear - 1)
            NandF(iAge + 2 + NAges, iYear) = F(MinAge + iAge - 1, MinYear + iYear - 1)
        Next
    Next

    ReturnNandF = NandF
    
End Function

Function BevHolt(SSB As Double, Steepness As Double, VirginSSB As Double, SPRF0 As Double) As Double
   Dim alpha As Double
   Dim beta As Double
   Dim VirginRecruits As Double
   
   VirginRecruits = VirginSSB / SPRF0
   
   alpha = 4# * Steepness * VirginRecruits / (5# * Steepness - 1#)
   beta = VirginSSB * (1# - Steepness) / (5# * Steepness - 1#)
   
   BevHolt = alpha * SSB / (beta + SSB)

End Function

Function BevHoltSPR(SPR As Double, Steepness As Double, VirginSSB As Double, SPRF0 As Double) As Double
   Dim alpha As Double
   Dim beta As Double
   Dim VirginRecruits As Double
   
   VirginRecruits = VirginSSB / SPRF0
   
   alpha = 4# * Steepness * VirginRecruits / (5# * Steepness - 1#)
   beta = VirginSSB * (1# - Steepness) / (5# * Steepness - 1#)
   
   BevHoltSPR = alpha * (alpha * SPR - beta) / (alpha * SPR)

End Function
Function Calcg(M As Double, Catch As Double, N As Double, Nt1 As Double) As Double
   Calcg = Catch - (1 - M / (Log(N) - Log(Nt1))) * (N - Nt1)
End Function
Function Calcdgdx(M As Double, N As Double, Nt1 As Double) As Double
   Calcdgdx = -1 - ((Log(N) - Log(N1)) * M - M * (N - N1) / N) / Log(2 * N / N1)
End Function
Function NewtonRhapson(x As Double, F As Double, dgdx As Double)
   NewtonRhapson = x - F / dgdx
End Function



*/

