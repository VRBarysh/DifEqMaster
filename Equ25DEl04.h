//---------------------------------------------------------------------------

#ifndef Equ25DEl04H
#define Equ25DEl04H
//---------------------------------------------------------------------------
#include "Equ25DEl02.h"
//---------------------------------------------------------------------------
class cEquVerticalASpeedDistr : public cEqu25DElModel
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
  void StepLBVMoving();
  void StepLBVMovingNoEls();
  void StepLBVMovingElWave();
  void StepLOVMoving();
  void StepRoutineRunge2_SimplePhaseModel();
  void StepIMode();
  void StepLBVMoving_EA();

  void FillSplineArray();
  void FillSplineCoefArray(long double *hData,
          long double *hSplineCoefB,long double *hSplineCoefC,long double *hSplineCoefD);


  virtual void LoadInitData(cDifEqu_InitDataRec *hInitData);

  long double SimplePhaseDeSync(int nZStep);
  long double SimplePhaseKappa(int nZStep);
//  virtual long double kappa(int zstep, int tstep) {return (1); }
  inline long double Right0();
  inline long double Right1();
  void CheckIfElStopped(long double *hEl,int nEl);

  double ElDroppedFractionCurrent() {
    return( 1.0-((double)nMovingEls)/InitData.nEls/InitData.nSpeedFracs );
  }

  double ElDroppedFractionTotal() {
    return( ((double)nTotalDroppedEls)/InitData.nEls/InitData.nSpeedFracs/nStepNumber );
  }

  long double EADesync(int nZStep);
  void InitEls3(int rMul);

  TCplxLong Calc_I(long double * hEls);

  void SavePData(long double *hPdata);
  void SaveDataArray(long double x1, long double x2, long double *hData, int n, AnsiString filename);
  void LoadEAData();
  void MakeNormalP();

  void InitTSliceSave();
  void TSliceSave();
  void TSliceZSaveA();
  void TSliceZSave(int num);
  void StopTSliceSave();

  long double GetOutA()
//    { return( InitData.fUseLOVModel ? abs(hAdata[0]) :
//                                      abs(hAdata[InitData.nSteps]) ); }
    {
      return  abs( hOutA[nStepNumber-1] );
    }             
  virtual long double GetMaxA() {
    long double Res=abs(hAdata[0]),tmp;
    for(int i=1;i<=InitData.nSteps;i++) {
      tmp=abs(hAdata[i]);
      Res = Res > tmp ? Res : tmp;
    }
    return( Res );
  }

  int EACalcMaxIPos(long double delta);
  long double EACalcMaxI(long double delta);
  long double EACalcMaxA(long double Coef);

  void MakeSpeedDistibution();
  long double CalcEnergy() { return CalcEnergyMaxI(); }
  long double CalcEnergyMaxI();

protected:
  int       NMovingEls;
  long double *hAvgP,maxW,maxdeltaTS,*hdeltaTS,intA,intW,*hNormalP;
  long double *hEAData,*EATempIData,EACoef;
  long double *hElk3,*hElk4,*hEldata;
  long double *hElSpeedFracs,*hElSpeed;
  long double BaseRandom,lastRandom,lastRandom1,lastRandom2,*hRandomData;
  TCplxLong *hAk1,*hAk2,*hABack,*hIBack,*hAdata,*hIdata,*hAdata1,*hIdata1,*hIk1,*hIk2,*hIk3,*hIk4;
  TCplxLong *hInterpolated;
  long double *hSplineCoefBRe,*hSplineCoefCRe,*hSplineCoefDRe,*htmpSplineData;
  long double *hSplineCoefBIm,*hSplineCoefCIm,*hSplineCoefDIm;
  long double workZ,EATime,TSliceSaveZ;
  int workElNum,workDeltaNum;
  TCplxLong workA;
  int nMovingEls,SplineRequest,nTotalDroppedEls;
  int TSliceSaveFlag;
  ofstream *hOutStream;
};

#endif
