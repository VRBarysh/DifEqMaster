//---------------------------------------------------------------------------

#ifndef EquOptic2D_1DwB_Mod1H
#define EquOptic2D_1DwB_Mod1H
//---------------------------------------------------------------------------

#include "EquOptic2D_1DwB.h"

class cOptic2D_Media_2Way_B_Mod1 {
public:
  /*
  inline cOptic2D_Media_X operator+(cOptic2D_Media_X m);
  inline cOptic2D_Media_X operator+=(cOptic2D_Media_X m);
  inline cOptic2D_Media_X operator*(OPTIC_TYPE v);
  inline cOptic2D_Media_X operator*=(OPTIC_TYPE v);
  inline void Add(cOptic2D_Media_X &m);
  inline void Mul(OPTIC_TYPE v);
  */
  inline void M1M2v(cOptic2D_Media_2Way_B_Mod1 &m1,cOptic2D_Media_2Way_B_Mod1 &m2,OPTIC_TYPE v) {
    R0=m1.R0+m2.R0*v;
    R0Generated=m1.R0Generated+m2.R0Generated*v;
  }

  OPTIC_TYPE R0,R0Generated;   // r0
};


// ------------------------ Много продольных волн с разными поперечными числами 
class cEquOptic2D_1DwB_Mod1 : public cDifEqu_DataBankIO {
public:
  cEquOptic2D_1DwB_Mod1() {
    InitData.nThreads=nEqus=0; hTaskEvent=NULL; hTaskCritSection=NULL;
    hB=hAZp=hAZm=NULL; hM=NULL; hReportBuf=NULL; t=0; hFFTSrc=hFFTDest=NULL;
    hAbsStore=NULL; hEnergyData=NULL;
    hAZpOutData=NULL;
    InitData.nStepsX=InitData.nStepsZ=InitData.nStepsT=0;
    hAZpD=hAZmD=hBD=NULL; hMatrix=NULL; htmpM=NULL;
    hMD=NULL; hDelta=NULL; hEAHarmP=NULL; hDeltaMatrixCoef=NULL;
  }

  ~cEquOptic2D_1DwB_Mod1() {
    delete [] hReportBuf;
    delete [] hAZp; delete [] hAZm; delete [] hB; delete [] hM;
    delete [] hFFTSrc; delete [] hFFTDest; delete [] hAbsStore;
    delete [] hEnergyData; delete [] hAZpOutData;
    delete [] hAZpD; delete [] hAZmD; delete [] hBD; delete [] hMD;
    delete [] hDelta; delete [] hEAHarmP; delete [] htmpM;
    delete [] hDeltaMatrixCoef; delete [] hMatrix;
    hFFTSrc=hFFTDest=NULL; hAbsStore=NULL; hEnergyData=NULL;
    InitData.nThreads=nEqus=0; hAZpOutData=NULL;
    hB=hAZp=hAZm=NULL;  hReportBuf=NULL; hM=NULL; t=0;
    hAZpD=hAZmD=hBD=NULL; hMatrix=NULL; htmpM=NULL;
    hMD=NULL; hDelta=NULL; hEAHarmP=NULL; hDeltaMatrixCoef=NULL;
  }

  void InitReport(cDifEquReport *hReport);
  void Report(cDifEquReport *hReport);
  void FinalReport(cDifEquReport *hReport);

  virtual void FinalSave();

  virtual void StepRoutine();
  virtual void StepRoutineEuler();
  virtual void StepRoutineRunge2();
  virtual void StepRoutineRunge2Matrix();
  virtual void StepRoutineRunge2MatrixExp();
  virtual void StepRoutineRunge2MatrixExpRx();

  virtual void PrepareMedia();
  virtual void FillBorders();

  virtual void CalcEnergy();
  virtual void CalcOutEnergy();

  complex<OPTIC_TYPE> Ap(int iz, double x); 

  virtual void LoadInitData(cDifEqu_InitDataRec *hInitData);
  inline int I(int i,int iz,int ix)
    { return i*InitData.nStepsX*InitData.SizeZ+InitData.SizeZ*ix+iz; }
  inline int IM(int i, int iz, int ix)
    { return i*InitData.nStepsX2*InitData.SizeZ+InitData.SizeZ*ix+iz; }

  inline double MediaCoef(int ix, int ix2) {
    return cos(M_PI*ix*ix2/((InitData.nStepsX)*(InitData.nStepsX)));
  }
  inline complex<OPTIC_TYPE> K1p(int ix)
    {  return hDeltaMatrixCoef[4*ix]; }
  inline complex<OPTIC_TYPE> K2p(int ix)
    {  return hDeltaMatrixCoef[4*ix+1]; }
  inline complex<OPTIC_TYPE> K1m(int ix)
    {  return hDeltaMatrixCoef[4*ix+2]; }
  inline complex<OPTIC_TYPE> K2m(int ix)
    {  return hDeltaMatrixCoef[4*ix+3]; }


protected:
cDifEqu_InitDataRec_Optic2D_T<OPTIC_TYPE> InitData;
complex<OPTIC_TYPE> *hDelta,*hDeltaMatrixCoef;
complex<OPTIC_TYPE> *hAZp,*hAZm,*hB;          // +0 - base
cOptic2D_Media_2Way_B_Mod1    *hM,*htmpM;     // +Size - k1
                                              // +2*Size - k2
                                              // +3*Size - new

complex<OPTIC_TYPE> *hAZpD,*hAZmD,*hBD,*hMatrix;
cOptic2D_Media_2Way_B_Mod1    *hMD;

// InitData.AGenX - коэффициент связи среды с волной B !
OPTIC_TYPE EnergyZP,EnergyZM,EnergyZPlast,ESumLast,EDelayBuf;
OPTIC_TYPE EnergyZPOut,EnergyZMOut,EnergyR,EnergyQ,EnergyA0,EnergyR0,EnergyUPOut,EnergyB,EnergyGenZ,EnergyGenX,EnergyQBack;
OPTIC_TYPE          PGenR, PGenA, POutA, PRel, PGenZ, PGenX;
OPTIC_TYPE          *hEAHarmP;
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
