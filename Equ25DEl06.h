//---------------------------------------------------------------------------

#ifndef Equ25DEl06H
#define Equ25DEl06H
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

#ifndef Equ25DEl04H
#define Equ25DEl04H
//---------------------------------------------------------------------------
#include "Equ25DEl02.h"
//---------------------------------------------------------------------------
class cEquSimpleLBVSpeedDistr : public cEqu25DElModel
{
public:
  cEquVerticalASpeedDistr();
  ~cEquVerticalASpeedDistr();

  virtual void InitReport(cDifEquReport *hReport);
  virtual void Report(cDifEquReport *hReport);
  virtual void FinalReport(cDifEquReport *hReport);

  virtual void StepRoutine();
  void StepRoutineRunge4();
  void StepLOV1();
  void StepRoutineRunge2_SimplePhaseModel();

  void FillSplineArray();
  void FillSplineCoefArray(long double *hData,
          long double *hSplineCoefB,long double *hSplineCoefC,long double *hSplineCoefD);


  virtual void LoadInitData(cDifEqu_InitDataRec *hInitData);

  long double SimplePhaseDeSync(int nZStep);
  long double SimplePhaseKappa(int nZStep);
//  virtual long double kappa(int zstep, int tstep) {return (1); }
  inline long double Right0();
  inline long double Right1();

  long double EADesync(int nZStep);
  void InitEls3(int rMul);

  TCplxLong Calc_I(long double * hEls);

  void SavePData(long double *hPdata);
  void SaveDataArray(long double x1, long double x2, long double *hData, int n, AnsiString filename);
  void LoadEAData();
  void MakeNormalP();

  long double CalcEnergy();

  void MakeSpeedDistibution();

protected:
  int       NMovingEls,MaxAPos;
  long double *hAvgP,maxW,maxdeltaTS,*hdeltaTS,intA,intW,*hNormalP;
  long double *hEAData,*EATempIData,EACoef;
  long double *hElk3,*hElk4,*hEldata;
  long double *hElSpeedFracs,*hElSpeed;
  long double BaseRandom,lastRandom,lastRandom1,lastRandom2,*hRandomData;
  TCplxLong *hAk1,*hAk2,*hABack,*hIBack,*hAdata,*hIdata,*hAdata1,*hIdata1;
  TCplxLong *hInterpolated;
  long double *hSplineCoefBRe,*hSplineCoefCRe,*hSplineCoefDRe,*htmpSplineData;
  long double *hSplineCoefBIm,*hSplineCoefCIm,*hSplineCoefDIm;
  long double workZ;
  int workElNum,workDeltaNum;
  TCplxLong workA;
  int nMovingEls,SplineRequest;
};

#endif

#endif
