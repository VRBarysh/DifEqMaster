//---------------------------------------------------------------------------

#ifndef Equ25DEl05H
#define Equ25DEl05H

#include "Dump25DElModel.h"
#include "Dump25DElForm.h"
#include "Equ25DEl02.h"
#include "Equ25DEl03.h"
#include "SimplePSO.h"
//---------------------------------------------------------------------------

class cEqu25DEl_PSOOptimizer : public cDifEqu_DataBankIO
{
public:
  cEqu25DEl_PSOOptimizer();
  ~cEqu25DEl_PSOOptimizer();

  void InitReport(cDifEquReport *hReport);
  void Report(cDifEquReport *hReport);
  void FinalReport(cDifEquReport *hReport);
  void StepRoutine();

  void ReportMode2(cDifEquReport *hReport);

  void LoadInitData(cDifEqu_InitDataRec *hInitData);

  void SetHEquation(cEqu25DElModel *phEquation)
    { hEquation=phEquation; }
protected:
  void NormalizeKappaData(long double *hData);
  void NormalizeDeltaZData(long double *hData);
  long double FuckFunction(long double *hData);

cDifEqu_InitDataRec_25DElModel    InitData;
CPSOSimple                        PSOOptimizer;
CPSOInitComp                      *hPointData;
long double                       *hVariable,BestOutValue,BestTotalValue,*hDeltaZData;
long double                       BestOutTime,BestTotalTime,*hKappaData,*hMaxKappaData;
cEqu25DElModel                    *hEquation;
long double                       MaxKappa[1000],MaxOut;
int                               tStep,vStep,nTSteps;
};

#endif
