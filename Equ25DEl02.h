//---------------------------------------------------------------------------

#ifndef Equ25DEl02H
#define Equ25DEl02H

#include <vcl.h>
#include "Dump25DElModel.h"
#include "Dump25DElForm.h"

#define KAPPAMODE_SIMPLE      0
#define KAPPAMODE_TOMSK       1
#define KAPPAMODE_LINEAR      2
#define KAPPAMODE_LINEARFRACS 3
#define KAPPAMODE_CONSTFRACS  4
#define KAPPAMODE_SPLINE      5
#define KAPPAMODE_FOURIER     6
#define KAPPAMODE_SIN2        7
#define KAPPAMODE_MOVEA       8
#define KAPPAMODE_MOVEASIN2    9

#define LDHALF ((long double)0.5)

//---------------------------------------------------------------------------

class cEqu25DElModel : public cDifEqu_DataBankIO
{
public:
  cEqu25DElModel();
  virtual ~cEqu25DElModel();

  virtual void InitReport(cDifEquReport *hReport);
  virtual void Report(cDifEquReport *hReport);
  virtual void FinalReport(cDifEquReport *hReport);

  TCplxLong GetBoundI() { return( t<InitData.rTime ? BoundIValue : TCplxLong(0,0)  ); }
  long double r() { return( t<InitData.rTime ? InitData.r : 0  ); }
  TCplxLong GetAInput() { return(TCplxLong(0,0)); }

  virtual long double KappaEXT(long double z,long double t=0,cDifEqu_InitDataRec_25DElModel *hInitData=NULL);
  virtual long double Kappa(long double z, long double t=0);
  virtual long double AvgKappa(long double t=0,cDifEqu_InitDataRec_25DElModel *hInitData=NULL);
  virtual void InitKappa();
  virtual int LoadKappaData(cDifEqu_InitDataRec_25DElModel *hData=NULL);
  virtual int LoadDeltaZData(cDifEqu_InitDataRec_25DElModel *hData=NULL);

  virtual long double Delta();
  virtual long double DeltaZ();
  long double DeltaZEXT(long double z,long double t,cDifEqu_InitDataRec_25DElModel *hInitData);
  virtual long double AInput();
  
  virtual void CalcI();
  virtual void StepRoutine();
  virtual void StepDataMove();
  virtual void FillAData();

  virtual long double GetOutA() { return( abs(hA[InitData.nSteps]) ); }
  virtual long double GetMaxA();


  virtual void InitEls(int nElFrac,int rMul);
  virtual void LoadInitData(cDifEqu_InitDataRec *hInitData);
  void SetHDump(cDump25DEl *phDump) { hDump=phDump; }

  long double RightPartEl0(long double El0,long double El1,TCplxLong A);
  long double RightPartEl1(long double El0,long double El1,TCplxLong A);
  long double *GetHElTime() { return(&InitData.ElTime); }

  virtual int EACalcMaxIPos(long double delta);
  virtual void EAMakeDelta();
  virtual long double CalcEnergy();

  int MaxAchieved() { return Got1stMax; }

protected:
  cDifEqu_InitDataRec_25DElModel InitData;
  long double *hEl,*hElNew;
  int       *hElPresent,*hNMovingEls;
  int       nZStepNumber,MaxAPos;
  long double *hElk1,*hElk2,*EATempIData,EACoef;
  long double KappaC,KappaZInv,dzKappa,z,zDeltaZBack,zBack;
  TCplxLong *hA,*hI,*hOutA;
  TCplxLong BoundIValue;
  cDump25DEl *hDump;
  int Got1stMax;
  long double tmpRealArray[10];
};

#endif
 