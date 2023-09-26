//---------------------------------------------------------------------------

#ifndef EquOptic2D_X_01H
#define EquOptic2D_X_01H
//---------------------------------------------------------------------------
#include "EquTaskThread.h"
#include "EquOpticTask.h"
#include "EquOpticTask_2BraggsSolidMedia.h"
#include "EquOpticTask_2D1D2DSolidMedia.h"
#include "Dump2DOpticModel.h"
#include "FFTUnit.h"

class EnergyDataType {
public:
OPTIC_TYPE A,R,Out,OutXp,OutXm,OutZp,OutZm,OutGenZ,OutGenX,E,PRel,OutZ2,OutXp2,OutXm2;
};

class cEquOptic2D_X_base : public cDifEqu_DataBankIO {
public:
  cEquOptic2D_X_base() {
    InitData.nThreads=nEqus=0; hTaskEvent=NULL; hTaskCritSection=NULL;
    hAXp=hAXm=hAZp=hAZm=NULL; hM=NULL; hReportBuf=NULL; t=0; hFFTSrc=hFFTDest=NULL;
    hAbsStore=NULL; hEnergyData=NULL;
    hAXpOutData=hAZpOutData=NULL;
    InitData.nStepsX=InitData.nStepsZ=InitData.nStepsT=0;
    fDebug=0;
  }

  ~cEquOptic2D_X_base() {
    for(int i=0;i<nEqus;i++) delete hEqus[i];
    for(int i=0;i<InitData.nThreads;i++) delete hThreads[i];
    delete hTaskEvent; delete hTaskCritSection;  delete [] hReportBuf;
    delete [] hAXp; delete [] hAXm; delete [] hAZp; delete [] hAZm; delete [] hM;
    delete [] hFFTSrc; delete [] hFFTDest; delete [] hAbsStore; delete [] hEnergyData;
    delete [] hAXpOutData; delete [] hAZpOutData; 
    hFFTSrc=hFFTDest=NULL; hAbsStore=NULL; hEnergyData=NULL;
    InitData.nThreads=nEqus=0; hTaskEvent=NULL; hTaskCritSection=NULL;
    hAXpOutData=hAZpOutData=NULL;
    hAXp=hAXm=hAZp=hAZm=NULL;  hReportBuf=NULL; hM=NULL; t=0;
  }

  virtual void InitReport(cDifEquReport *hReport);
  virtual void Report(cDifEquReport *hReport);
  virtual void FinalReport(cDifEquReport *hReport);

  virtual void SaveADataMaxE();
  virtual void FinalSave();

  virtual void StepRoutine();
  virtual void DoTasks();
  virtual void StepRoutineEulerNoMedia();
  virtual void StepRoutineRunge2NoMedia();
  
  virtual void FillBorders(int iLayer);
  virtual void PrepareMedia();

  virtual void AddNoise();
  virtual void PurgeNoise();
  virtual void PurgeRoutine1();
  virtual void PurgeRoutine2();
  virtual void ForceLinear();
  inline OPTIC_TYPE MicrowaveModePulse();

  virtual void CalcEnergy();
  virtual void CalcOutEnergy();

  virtual void SaveBitmaps();
  int TimeToSaveBitmaps();
  void PrepareBitmapFolders();

  virtual int Get2DColor(int Number2D, double x, double y);
  virtual void Prepare2DColor(int Number2D);
  void PrepareColorData();

  virtual void LoadInitData(cDifEqu_InitDataRec *hInitData);

  inline int IndexXZ(int ix,int iz) { return((InitData.nStepsZ+2)*ix+iz); }
  inline int Index_ZX(int iz,int ix) { return((InitData.nStepsZ+2)*ix+iz); }

protected:
cDifEqu_InitDataRec_Optic2D_T<OPTIC_TYPE> InitData;
complex<OPTIC_TYPE> *hAXp,*hAXm,*hAZp,*hAZm;  // +0 - base
cOptic2D_Media_X    *hM;                      // +Size - k1
                                              // +2*Size - k2
OPTIC_TYPE EnergyXP,EnergyXM,EnergyZP,EnergyZM,EnergyZPlast,ESumLast,EDelayBuf, Frequency, lastFrequency;
OPTIC_TYPE EnergyXPOut,EnergyXMOut,EnergyZPOut,EnergyZMOut,EnergyR,EnergyQ,EnergyA0,EnergyR0,EnergyUPOut;
OPTIC_TYPE          *hReportBuf;
EnergyDataType      *hEnergyData;
OPTIC_TYPE          maxAXp,maxAXm,maxAZp,maxAZm,maxPXp,maxPXm,maxPZp,maxPZm;
OPTIC_TYPE          maxmaxAXp,maxmaxAXm,maxmaxAZp,maxmaxAZm;
OPTIC_TYPE          normAXp,normAXm,normAZp,normAZm,ComplexBitmapNorm;
OPTIC_TYPE          LastReportTime,LastMaxOut,LastBitmapSaveTime;
OPTIC_TYPE          decA,decZP;
complex<OPTIC_TYPE> LastA;
int PurgeNoiseTime, fDebug;

complex<OPTIC_TYPE> *hFFTSrc,*hFFTDest;
CFFTrans<OPTIC_TYPE> FFT;
int FFTSizeFull,nBitmap,MakeGlobalNormBitmap,MakeComplexBitmap;

cOptic2D_AbsStore *hAbsStore;

complex<OPTIC_TYPE> *hAXpOutData, *hAZpOutData;

int ThreadCounter,nEqus,MaxESaved;
cEquOpticBaseTask *hEqus[MAX_THREADS+1];
TTaskThread       *hThreads[MAX_THREADS];
TEvent            *hTaskEvent;
TCriticalSection  *hTaskCritSection;
int               ColorData[256],ColorDataGray[256];
int               iDebug1,iDebug2,iDebug3,iDebug4,iDebug5;
int               iSpectrumX,iSpectrumZ;  
};

#endif
