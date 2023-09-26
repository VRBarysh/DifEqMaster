//---------------------------------------------------------------------------

#ifndef EquOptic2D_X_MultiMediaH
#define EquOptic2D_X_MultiMediaH
//---------------------------------------------------------------------------
#include "EquOptic2D_X_01.h"

class cEquOptic2D_X_MultiMedia : public cEquOptic2D_X_base {
public:
  cEquOptic2D_X_MultiMedia() {
    InitData.nThreads=nEqus=0; hTaskEvent=NULL; hTaskCritSection=NULL;
    hAXp=hAXm=hAZp=hAZm=NULL; hM=NULL; hReportBuf=NULL; t=0; hFFTSrc=hFFTDest=NULL;
    hAbsStore=NULL; hEnergyData=NULL;
    hAXpOutData=hAZpOutData=NULL;
    InitData.nStepsX=InitData.nStepsZ=InitData.nStepsT=0;
    fDebug=0;
  }

  ~cEquOptic2D_X_MultiMedia() {
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

  virtual void Report(cDifEquReport *hReport);
  virtual void FinalReport(cDifEquReport *hReport);

  virtual void FinalSave();

  virtual void FillBorders(int iLayer);
  virtual void PrepareMedia();
  virtual void ForceLinear();

  virtual void CalcEnergy();
  virtual void CalcOutEnergy();

  virtual void LoadInitData(cDifEqu_InitDataRec *hInitData);

  inline int IndexXZ(int ix,int iz) { return((InitData.nStepsZ+2)*ix+iz); }
  inline int Index_ZX(int iz,int ix) { return((InitData.nStepsZ+2)*ix+iz); }
  inline int IndexXZM(int ix,int iz,int im) { return( ((InitData.nStepsZ+2)*ix+iz)*InitData.nMedia+im ); }

  OPTIC_TYPE MediaDistrNoNorm(int im);
  OPTIC_TYPE MediaDistrNormCoef();
  void TransformMediaCoef();

  OPTIC_TYPE MicrowaveModePulse() { return 0; }

protected:

};

#endif
