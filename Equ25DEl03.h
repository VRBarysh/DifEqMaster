//---------------------------------------------------------------------------

#ifndef Equ25DEl03H
#define Equ25DEl03H

#include "Equ25DEl02.h"
//---------------------------------------------------------------------------
class cEqu2PhaseElectron {
public:
cEqu2PhaseElectron operator+(cEqu2PhaseElectron El);
cEqu2PhaseElectron operator*(long double mul);

long double tI,tS,W;
TCplxLong P,P1,P2;
int Present;
};

class cEqu2PhaseElModel : public cEqu25DElModel
{
public:
  cEqu2PhaseElModel();
  ~cEqu2PhaseElModel();

  virtual void InitReport(cDifEquReport *hReport);
  virtual void Report(cDifEquReport *hReport);
  virtual void FinalReport(cDifEquReport *hReport);

  virtual void StepRoutine();
  void StepRoutineRunge2();
  void StepRoutineRunge4();
  void StepRoutineRunge2_SimplePhaseModel();

  virtual void InitEls2(int rMul);
  virtual void LoadInitData(cDifEqu_InitDataRec *hInitData);

  cEqu2PhaseElectron RightPartEl2(long double El_TI,long double El_TS,
                             TCplxLong El_P,long double El_W,TCplxLong A,int nZStep);
  long double SimplePhaseDeSync(int nZStep);
  long double SimplePhaseKappa(int nZStep);

  long double EADesync(int nZStep);

  void SavePData(long double *hPdata);
  void SaveDataArray(long double x1, long double x2, long double *hData, int n, AnsiString filename);
  void LoadEAData();
  void MakeNormalP();

  void SaveEnergyAdaptCalibrationData();
  int EACalcMaxIPos(long double delta);
  long double CalcEnergy();

protected:
  cEqu2PhaseElectron *hEl2,*hElNew2,*hElk12,*hElk22;
  int       NMovingEls;
  long double *hAvgP,maxW,maxdeltaTS,*hdeltaTS,intA,intW,*hNormalP;
  long double *hEAData;
  long double BaseRandom,lastRandom,lastRandom1,lastRandom2,*hRandomData;
  TCplxLong *hAk1,*hAk2,*hABack,*hIBack;
};

#endif
