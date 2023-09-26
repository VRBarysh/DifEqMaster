//---------------------------------------------------------------------------

#ifndef EquOptic_X_04H
#define EquOptic_X_04H
//---------------------------------------------------------------------------

#include "EquOpticTask_BMB_holes.h"
#include "EquOptic2D_X_03.h"

class cEquOptic2D_X_Mod2 : public cEquOptic2D_X_Mod1 {
public:
  cEquOptic2D_X_Mod2() {
    InitData.nThreads=nEqus=0; hTaskEvent=NULL; hTaskCritSection=NULL;
    hAXp=hAXm=hAZp=hAZm=NULL; hM=NULL; hReportBuf=NULL; t=0; hFFTSrc=hFFTDest=NULL;
    hAbsStore=NULL; hEnergyData=NULL; hAZpMedia=hAZmMedia=NULL;
    InitData.nStepsX=InitData.nStepsZ=InitData.nStepsT=0; hM2=NULL;
    hAXp1=NULL; hAXm1=NULL; hAZp1=NULL; hAZm1=NULL;
    hAXp2=NULL; hAXm2=NULL; hAZp2=NULL; hAZm2=NULL;
    fDebug=0;
  }

  ~cEquOptic2D_X_Mod2() {
    for(int i=0;i<nEqus;i++) delete hEqus[i];
    for(int i=0;i<InitData.nThreads;i++) delete hThreads[i];
    delete hTaskEvent; delete hTaskCritSection;  delete [] hReportBuf;
    delete [] hM;
    delete [] hFFTSrc; delete [] hFFTDest; delete [] hAbsStore; delete [] hEnergyData;
    delete [] hAZpMedia; delete [] hAZmMedia; delete [] hM2;
    delete [] hAXp1; delete [] hAXm1; delete [] hAZp1; delete [] hAZm1;
    delete [] hAXp2; delete [] hAXm2; delete [] hAZp2; delete [] hAZm2;
    hFFTSrc=hFFTDest=NULL; hAbsStore=NULL; hEnergyData=NULL; hAZpMedia=hAZmMedia=NULL;
    InitData.nThreads=nEqus=0; hTaskEvent=NULL; hTaskCritSection=NULL;
    hAXp=hAXm=hAZp=hAZm=NULL;  hReportBuf=NULL; hM=NULL; t=0; hM2=NULL;
    hAXp1=NULL; hAXm1=NULL; hAZp1=NULL; hAZm1=NULL;
    hAXp2=NULL; hAXm2=NULL; hAZp2=NULL; hAZm2=NULL;
    iDebug1=0;
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
  inline int IndexXZ_L(int ix,int iz) {
    int s1=(InitData.nStepsZ+2)*(InitData.nStepsX+2);
    int s2=(InitData.nStepsZ+2)*(InitData.nStepsX+2)+(InitData.nStepsZB+2)*(InitData.nStepsX+2);
    int nz1=InitData.nStepsZ+2;
    int nz2=InitData.nStepsZ+InitData.nStepsZB+4;
    if(iz<nz1) return IndexXZ(ix,iz);
    if(iz<nz2) return IndexXZ_M(ix,iz-nz1)+s1;
    return IndexXZ(ix,iz-nz2)+s2;
  }

  int isBragg(int ix);
  OPTIC_TYPE ActivePartSize();

  virtual int Get2DColor(int Number2D, double x, double y);
  virtual void Prepare2DColor(int Number2D);

protected:

complex<OPTIC_TYPE> *hAXp1,*hAXm1,*hAZp1,*hAZm1;
complex<OPTIC_TYPE> *hAXp2,*hAXm2,*hAZp2,*hAZm2;

complex<OPTIC_TYPE> lastAforFreq[128];
OPTIC_TYPE DeltaArray[128];
};

#endif
