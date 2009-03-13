//
// flrxsa.hpp
//
// Last Change: 10 Ago 2005 11:44
// $Id: flrxsa.hpp,v 1.10 2007/10/22 09:29:13 ltkell Exp $

#ifndef _INC_flrxsa
#define _INC_flrxsa

#include "math.h"
#include <stdio.h>

#include "fl__r.hpp"
#include "fl__xsa.hpp"

#define NXSAPARAMS  15
 
class CatchDataR : public CatchData
{
public:
   
   CatchDataR(SEXP);  
};                  


class TuningFleetsR: public TuningFleets
{
public:
                                                  
TuningFleetsR(SEXP, CatchData *, short TuningWindow=-1);
};                  


class ExtendedSurvivorsAnalysisR : public ExtendedSurvivorsAnalysis
{
public:                                             
   ExtendedSurvivorsAnalysisR(SEXP, SEXP, CatchData *, TuningFleets *);
   
   SEXP ReturnQ(short,short);
   SEXP ReturnQRes(short);
   SEXP ReturnIndexHat(short);
   SEXP ReturnIndexVar(short);
   SEXP ReturnCorrectedCPUE(short iFleet, short _MaxYear = -1);
   SEXP ReturnControl(void);

    
   SEXP Return(int _MaxYear=-1);
   SEXP ReturnSimple(int _MaxYear=-1);
   SEXP ReturnWts(bool ScaleFlag = true, bool LogTransformFlag = true);
   SEXP ReturnInternalSE(void);
   SEXP ReturnExternalSE(void);

   SEXP dataframe(void);
   bool OutputTermPopWts(LP2DOUBLE *, int *, int *);
   void InputNF(SEXP);
};                  

#endif /* _INC_flrxsa */
                                                  
