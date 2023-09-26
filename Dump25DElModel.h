//---------------------------------------------------------------------------

#ifndef Dump25DElModelH
#define Dump25DElModelH

#include <vcl.h>
#include "EquBaseMaster.h"
#include "EquBaseDump.h"
#include "Gr3DPlane.h"
#include "SimplePSO.h"
//---------------------------------------------------------------------------

class DrawDumpInitData {
public:
int nXSteps;
long double t,dx;
bool DoubleSided,DrawLines;
int Skip;
};

class cDifEqu_InitDataRec_25DElModel : public cDifEqu_InitDataRec
{
public:
  cDifEqu_InitDataRec_25DElModel() { fDeleteKappaData=0; fDeleteDeltaData=0;
                                     fDeleteDeltaZData=0;
                                     hKappaData=NULL; hDeltaData=NULL; }
int nSteps,nEls,ElDirection,ADirection,ElsOnStart,nEASteps,nSpeedFracs,KappaMode;
long double r,rTime,a0,a0Time,L,ElTime,ICoef,deltaI,deltaS,rndCoef,Ai_kappa,nKappaData;
long double KappaV1,KappaV2,KappaV3,KappaV4,KappaV5,KappaV6,KappaV7,KappaV8,KappaV9,KappaV10,KappaV11;
TCplxLong AI;
long double EAStartTime,EAnSteps,EADelta,EAEndTime,EADeltaZ0,EAkz,SpeedDistribFactor;
long double TSliceSaveTime,TSliceSaveZ;
bool UseSimplePhaseModel,NoRandom,UseEnergyAdapt,MoveAData,OutMaxA,fLoadKappaData,fDeleteKappaData;
bool fLoadDeltaZData,fDeleteDeltaZData;
bool fUseLOVModel,fNeed1stPeak,fDeleteDeltaData;
long double *hKappaData,Nu;
long double deltaV1,deltaV2,AInputV1,AInputV2,AInputV3;
int PSOMode,nKappaSteps,ElMode,AInputMode;
int DeltaMode,DeltaZMode,nDeltaZSteps;
long double *hDeltaData,*hDeltaZData,DeltaZMax;
CPSOInitData PSOData;
AnsiString KappaDataFileName;
};

class cDump25DEl : public cEquBaseDump
{
public:
  cDump25DEl();
  ~cDump25DEl() {}

  void LoadInitData(cDifEqu_InitDataRec_25DElModel *hInit);
  void Dump(TCplxLong *hAData,TCplxLong *hIData,long double *hElData,int *hElPresDAta);
  void SaveCURR(char *hFileName);
  void SaveXTSlice(DrawDumpInitData DrawInit);
  void SaveXZSlice(DrawDumpInitData DrawInit);
  void SaveZTSlice(DrawDumpInitData DrawInit);

  void Show() {}
  void MakeXProfileCoef(long double pdx, int pnXSteps);
  void DrawXProfileCoef(Graphics::TBitmap *hBitmap,cGr3DInitData Gr3DInit,DrawDumpInitData DrawInit);
  void DrawZSlice(Graphics::TBitmap *hBitmap,cGr3DInitData Gr3DInit,DrawDumpInitData DrawInit);
  void DrawZTSlice(Graphics::TBitmap *hBitmap,cGr3DInitData Gr3DInit,DrawDumpInitData DrawInit);
  void DrawXTSlice(Graphics::TBitmap *hBitmap,cGr3DInitData Gr3DInit,DrawDumpInitData DrawInit);
  void DrawXZSlice(Graphics::TBitmap *hBitmap,cGr3DInitData Gr3DInit,DrawDumpInitData DrawInit);
  TCplxLong GetAVal(int ZStep,int XStep,int TStep);
  double GetMaxTime() {return(Step*Init.dt);}
  double GetMaxZ() {return(Init.nSteps*Init.dt);}


cDifEqu_InitDataRec_25DElModel Init;
int DumpEls,nTSteps,Step,nXSteps;
TCplxLong *hA,*hI,*hXProfileA,*hXProfileB;
long double *hEl,dx;
int *hElPresent;

int DumpChanged;
protected:
};


/*
class cDumpStorage
{
public:
  cDumpStorage(){ hADump=hIDump=NULL; hElDump=hMaxADump=NULL;
                  nEls=nSteps=nStepsT=nOptiSteps1=nOptiSteps2=nOptiSteps3=0;
                  Dumping=DUMPING_NONE; DumpReady=DumpRequired=DumpElectrons=false; }
  ~cDumpStorage(){ delete [] hADump; delete [] hIDump;
                   delete [] hElDump; delete [] hMaxADump; }

  void Init_Diffur(int pnSteps, int pnEls, REAL_TYPE pL, REAL_TYPE ptMax);
  void Init_Opti(int Type,
            int pnOptiSteps1,int pnOptiSteps2=1,int pnOptiSteps3=1);

  void DumpOpti(REAL_TYPE MaxA);
  void DumpDiffur( complex<REAL_TYPE> *hA,
                   complex<REAL_TYPE> *hI, REAL_TYPE *hEls);
  bool Ready() { return(DumpReady); }

  void Save_AUDump(char *hFileName,float time);
  void Save_AIDump(char *hFileName,float time);
  void Save_AIUDump(char *hFileName,float time);

  void SaveSurfaceZT(char *hFileName);
  void SaveCURR(char *hFileName);

  void Draw_AUDump(TImage *hImage,float time);
  void Draw_AIUDump(TImage *hImage,float time);

  complex<float> GetAValue(REAL_TYPE z,REAL_TYPE t);
  complex<float> GetIValue(REAL_TYPE z,REAL_TYPE t);

  float GetEndTime() {return(EndTime);}
  float GetL() { return(L); }
  double GetMaxA() { return(MaxA); }
  double GetMaxI() { return(MaxI); }

  void SetRequiredState(bool pDumpRequired,bool pDumpElectrons)
    { DumpRequired=pDumpRequired; DumpElectrons=pDumpElectrons; }

protected:
complex<float>        *hADump,*hIDump;
float                 *hElDump;
float                 *hMaxADump;
REAL_TYPE             EndTime,L;
int                   nSteps,nStepsT,nEls;
int                   nOptiSteps1,nOptiSteps2,nOptiSteps3;
int                   Step,StepT,OptiStep1,OptiStep2,OptiStep3;
int                   Dumping;
bool                  DumpReady,DumpRequired,DumpElectrons;
double                MaxI,MaxA,MinU,MaxU;
};
*/


#endif
