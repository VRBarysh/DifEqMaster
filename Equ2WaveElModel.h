//---------------------------------------------------------------------------

#ifndef Equ2WaveElModelH
#define Equ2WaveElModelH
//---------------------------------------------------------------------------

#include "EquTaskThread.h"
#include "Equ2WaveTask.h"
#include "FFTUnit.h"
#include "SplineUnit.h"

/*
class EnergyDataType {
public:
OPTIC_TYPE A,R,Out,OutXp,OutXm,OutZp,OutZm,OutGenZ,OutGenX,E;
};
*/

class cEnergyData_2WaveEl {
public:
REAL_TYPE E,Ap,Am,Out;
};

class cEqu2WaveElModel : public cDifEqu_DataBankIO {
public:
  cEqu2WaveElModel() {
    InitData.nThreads=nEqus=0; hTaskEvent=NULL; hTaskCritSection=NULL;
    hAp=hAm=hAhalf=hI=NULL; hEl=NULL; t=0; nStepNumber=0; hReportBuf=NULL;
    InitData.nSteps=InitData.nEls=InitData.nStepsT=InitData.nAmSteps=0;
    htmpAp=htmpAm=NULL; hZ=NULL; hApInter=hAmInter=NULL;
    for(int i=0;i<6;i++) hSCoef[i]=NULL; hOutA=NULL;
    hAmDump=hApDump=NULL; hEnergyData=NULL;
  }

  ~cEqu2WaveElModel() {
    for(int i=0;i<nEqus;i++) delete hEqus[i];
    for(int i=0;i<InitData.nThreads;i++) delete hThreads[i];
    delete hTaskEvent; delete hTaskCritSection; delete [] hReportBuf;
    delete [] hAp; delete [] hAm; delete [] hAhalf; delete [] hEl; delete [] hI;
    delete [] htmpAp; delete [] htmpAm; delete [] hZ; delete [] hOutA;
    for(int i=0;i<6;i++) delete [] hSCoef[i]; delete [] hApInter; delete [] hAmInter;
    delete [] hAmDump; delete [] hApDump; delete [] hEnergyData;
    hEnergyData=NULL;
    InitData.nThreads=nEqus=0; hTaskEvent=NULL; hTaskCritSection=NULL;
    hAp=hAm=hAhalf=hI=NULL; hEl=NULL; t=0; nStepNumber=0; hReportBuf=NULL;
    InitData.nSteps=InitData.nEls=InitData.nStepsT=InitData.nAmSteps=0;
    htmpAp=htmpAm=NULL; hZ=NULL; hApInter=hAmInter=NULL; hAmDump=hApDump=NULL;
    for(int i=0;i<6;i++) hSCoef[i]=NULL; hOutA=NULL;
  }

  void InitReport(cDifEquReport *hReport);
  void Report(cDifEquReport *hReport);
  void FinalReport(cDifEquReport *hReport);

  virtual void FinalSave();

  virtual void StepRoutine();
  virtual void DoTasks();

  virtual void LoadInitData(cDifEqu_InitDataRec *hInitData);

  void PrepareMedia();
  void CalcEnergy();


  // inline REAL_TYPE Fz(int iz) { return 1; }
  inline REAL_TYPE Fz(REAL_TYPE z) { return 1; }

protected:
cDifEqu_InitDataRec_2WaveElModel InitData;
complex<REAL_TYPE> *hAp,*hAm,*hAhalf;      // +0 - base
REAL_TYPE          *hEl,*hZ;               // +Size - k1
                                           // +2*Size - k2
complex<REAL_TYPE> *hI,*htmpAp,*htmpAm,*hApInter,*hAmInter,*hOutA;

complex<float>              *hAmDump,*hApDump;  
REAL_TYPE          *hSCoef[6];

REAL_TYPE          *hReportBuf;
CSpline<REAL_TYPE>    SplineAp,SplineAm;

cEnergyData_2WaveEl *hEnergyData;

REAL_TYPE         EAp,EAm,EOut,EE;

int ThreadCounter,nEqus,MaxESaved;
cEqu2WaveBaseTask *hEqus[MAX_THREADS+1];
TTaskThread       *hThreads[MAX_THREADS];
TEvent            *hTaskEvent;
TCriticalSection  *hTaskCritSection;
};

#endif
