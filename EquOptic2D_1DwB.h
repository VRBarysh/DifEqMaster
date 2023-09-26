//---------------------------------------------------------------------------

#ifndef EquOptic2D_1DwBH
#define EquOptic2D_1DwBH
//---------------------------------------------------------------------------

#include "Dump2DOpticModel.h"
#include "EquOptic2D_X_01.h"
#include "FFTUnit.h"

class cOptic2D_Media_2Way_B {
public:
  /*
  inline cOptic2D_Media_X operator+(cOptic2D_Media_X m);
  inline cOptic2D_Media_X operator+=(cOptic2D_Media_X m);
  inline cOptic2D_Media_X operator*(OPTIC_TYPE v);
  inline cOptic2D_Media_X operator*=(OPTIC_TYPE v);
  inline void Add(cOptic2D_Media_X &m);
  inline void Mul(OPTIC_TYPE v);
  */
  inline void M1M2v(cOptic2D_Media_2Way_B &m1,cOptic2D_Media_2Way_B &m2,OPTIC_TYPE v) {
    PZp=m1.PZp+m2.PZp*v; PZm=m1.PZm+m2.PZm*v; PB=m1.PB+m2.PB*v;
    R0=m1.R0+m2.R0*v; RZ=m1.RZ+m2.RZ*v;
    R0Generated=m1.R0Generated+m2.R0Generated*v;
  }

  complex<OPTIC_TYPE> PZp,PZm,PB,RZ; // r2z, r2x, rz+x, rz-x
  OPTIC_TYPE R0,R0Generated;   // r0
};


class cEquOptic2D_1DwB : public cDifEqu_DataBankIO {
public:
  cEquOptic2D_1DwB() {
    InitData.nThreads=nEqus=0; hTaskEvent=NULL; hTaskCritSection=NULL;
    hB=hAZp=hAZm=NULL; hM=NULL; hReportBuf=NULL; t=0; hFFTSrc=hFFTDest=NULL;
    hAbsStore=NULL; hEnergyData=NULL;
    hAZpOutData=NULL;
    InitData.nStepsX=InitData.nStepsZ=InitData.nStepsT=0;
    hAZpD=hAZmD=hBD=NULL;
    hMD=NULL;
  }

  ~cEquOptic2D_1DwB() {
    delete [] hReportBuf;
    delete [] hAZp; delete [] hAZm; delete [] hB; delete [] hM;
    delete [] hFFTSrc; delete [] hFFTDest; delete [] hAbsStore;
    delete [] hEnergyData; delete [] hAZpOutData;
    delete [] hAZpD; delete [] hAZmD; delete [] hBD; delete [] hMD;
    hFFTSrc=hFFTDest=NULL; hAbsStore=NULL; hEnergyData=NULL;
    InitData.nThreads=nEqus=0; hAZpOutData=NULL;
    hB=hAZp=hAZm=NULL;  hReportBuf=NULL; hM=NULL; t=0;
    hAZpD=hAZmD=hBD=NULL;
    hMD=NULL;
  }

  void InitReport(cDifEquReport *hReport);
  void Report(cDifEquReport *hReport);
  void FinalReport(cDifEquReport *hReport);

  virtual void FinalSave();

  virtual void StepRoutine();
  virtual void StepRoutineEuler();
  virtual void StepRoutineRunge2();
  virtual void StepRoutineRunge2Matrix();

  virtual void StepDebug();

  virtual void PrepareMedia();
  virtual void FillBorders();

  virtual void CalcEnergy();
  virtual void CalcOutEnergy();

  virtual void LoadInitData(cDifEqu_InitDataRec *hInitData);

protected:
cDifEqu_InitDataRec_Optic2D_T<OPTIC_TYPE> InitData;
complex<OPTIC_TYPE> *hAZp,*hAZm,*hB;          // +0 - base
cOptic2D_Media_2Way_B    *hM;                 // +Size - k1
                                              // +2*Size - k2

complex<OPTIC_TYPE> *hAZpD,*hAZmD,*hBD; 
cOptic2D_Media_2Way_B    *hMD;

// InitData.AGenX - коэффициент связи среды с волной B !
OPTIC_TYPE EnergyZP,EnergyZM,EnergyZPlast,ESumLast,EDelayBuf;
OPTIC_TYPE EnergyZPOut,EnergyZMOut,EnergyR,EnergyQ,EnergyA0,EnergyR0,EnergyUPOut,EnergyB,EnergyGenZ,EnergyGenX,EnergyQBack;
OPTIC_TYPE          PGenR, PGenA, POutA, PRel, PGenZ, PGenX;
OPTIC_TYPE          *hReportBuf;
EnergyDataType      *hEnergyData;
OPTIC_TYPE          maxAXp,maxAXm,maxAZp,maxAZm,maxPXp,maxPXm,maxPZp,maxPZm;
OPTIC_TYPE          maxmaxAXp,maxmaxAXm,maxmaxAZp,maxmaxAZm;
OPTIC_TYPE          normAXp,normAXm,normAZp,normAZm,ComplexBitmapNorm;
OPTIC_TYPE          LastReportTime,LastMaxOut,LastBitmapSaveTime;
int PurgeNoiseTime;

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
};

#endif
