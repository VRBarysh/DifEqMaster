//---------------------------------------------------------------------------

#ifndef Dump2WaveElModelH
#define Dump2WaveElModelH
//---------------------------------------------------------------------------
#include <vcl.h>
#include "EquBaseMaster.h"
#include "EquBaseDump.h"
#include "Gr3DPlane.h"
#include "SimplePSO.h"
//---------------------------------------------------------------------------

#define REAL_TYPE double

class cDifEqu_InitDataRec_2WaveElModel : public cDifEqu_InitDataRec
{
public:
  cDifEqu_InitDataRec_2WaveElModel() { fDeleteKappaData=0; fDeleteDeltaData=0;
                                     fDeleteDeltaZData=0;
                                     hKappaData=NULL; hDeltaData=NULL; }
int nSteps,nEls,ElDirection,ADirection,ElsOnStart,nEASteps,nSpeedFracs,KappaMode,nAmSteps,nStepsT;
REAL_TYPE r,rTime,a0,a0Time,L,ElTime,ICoef,deltaI,deltaS,rndCoef,Ai_kappa,nKappaData,ACoef;
REAL_TYPE KappaV1,KappaV2,KappaV3,KappaV4,KappaV5,KappaV6,KappaV7,KappaV8,KappaV9,KappaV10,KappaV11;
complex<REAL_TYPE> AI;
REAL_TYPE EAStartTime,EAnSteps,EADelta,EAEndTime,EADeltaZ0,EAkz,SpeedDistribFactor;
REAL_TYPE TSliceSaveTime,TSliceSaveZ;
bool UseSimplePhaseModel,NoRandom,UseEnergyAdapt,MoveAData,OutMaxA,fLoadKappaData,fDeleteKappaData;
bool fLoadDeltaZData,fDeleteDeltaZData;
bool fUseLOVModel,fNeed1stPeak,fDeleteDeltaData;
REAL_TYPE *hKappaData,Nu;
REAL_TYPE deltaV1,deltaV2,AInputV1,AInputV2,AInputV3;
int PSOMode,nKappaSteps,ElMode,AInputMode,nThreads;
int DeltaMode,DeltaZMode,nDeltaZSteps;
REAL_TYPE *hDeltaData,*hDeltaZData,DeltaZMax;
CPSOInitData PSOData;
AnsiString KappaDataFileName;
int RandSeed;
};

#endif
