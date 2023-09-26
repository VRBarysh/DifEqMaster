//---------------------------------------------------------------------------

#ifndef EquOptic2D_X_03H
#define EquOptic2D_X_03H

#include "EquOpticTask_2Way.h"
#include "EquOptic2D_X_01.h"
//---------------------------------------------------------------------------

//           mirror - media - bragg 
class cEquOptic2D_X_Mod1 : public cEquOptic2D_X_base {
public:
  cEquOptic2D_X_Mod1() {
    InitData.nThreads=nEqus=0; hTaskEvent=NULL; hTaskCritSection=NULL;
    hAXp=hAXm=hAZp=hAZm=NULL; hM=NULL; hReportBuf=NULL; t=0; hFFTSrc=hFFTDest=NULL;
    hAbsStore=NULL; hEnergyData=NULL; hAZpMedia=hAZmMedia=NULL;
    InitData.nStepsX=InitData.nStepsZ=InitData.nStepsT=0; hM2=NULL;
  }

  ~cEquOptic2D_X_Mod1() {
    for(int i=0;i<nEqus;i++) delete hEqus[i];
    for(int i=0;i<InitData.nThreads;i++) delete hThreads[i];
    delete hTaskEvent; delete hTaskCritSection;  delete [] hReportBuf;
    delete [] hAXp; delete [] hAXm; delete [] hAZp; delete [] hAZm; delete [] hM;
    delete [] hFFTSrc; delete [] hFFTDest; delete [] hAbsStore; delete [] hEnergyData;
    delete [] hAZpMedia; delete [] hAZmMedia; delete [] hM2;
    hFFTSrc=hFFTDest=NULL; hAbsStore=NULL; hEnergyData=NULL; hAZpMedia=hAZmMedia=NULL;
    InitData.nThreads=nEqus=0; hTaskEvent=NULL; hTaskCritSection=NULL;
    hAXp=hAXm=hAZp=hAZm=NULL;  hReportBuf=NULL; hM=NULL; t=0; hM2=NULL;
  }

  void InitReport(cDifEquReport *hReport);
  void Report(cDifEquReport *hReport);
  void FinalReport(cDifEquReport *hReport);

  //void SaveADataMaxE();
  void FinalSave();

  void StepRoutine();
  void DoTasks();
  void StepRoutineEulerNoMedia();
  void StepRoutineRunge2NoMedia();
  
  void FillBorders(int iLayer);
  void PrepareMedia();

  void CalcEnergy();
  void CalcOutEnergy();

  virtual void LoadInitData(cDifEqu_InitDataRec *hInitData);

  inline int IndexXZ_M(int ix,int iz) { return((InitData.nStepsZB+2)*ix+iz); }

protected:
complex<OPTIC_TYPE> *hAZpMedia,*hAZmMedia;     // +2*Size - k2
cOptic2D_Media_X_2Way *hM2;

OPTIC_TYPE EnergyXP,EnergyZPMedia,EnergyZMMedia;
};

#endif
