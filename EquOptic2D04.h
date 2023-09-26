//---------------------------------------------------------------------------

#ifndef EquOptic2D04H
#define EquOptic2D04H

#include "EquBaseMaster.h"
#include "FFTUnit.h"
#include "Dump2DOpticModel.h"
#include "EquOptic2D01.h"
#include "EquOptic2D02.h"
#include "EquOptic2D03.h"
//---------------------------------------------------------------------------

class cEquOptic2D_01 : public cDifEqu_DataBankIO {
public:
  cEquOptic2D_01() : cDifEqu_DataBankIO() { hMedia=hMediaNew=NULL;
                  hMediaBragg=hMediaBraggNew=NULL;
                  hFFTSrc=hFFTDest=hApData=hAmData=NULL;
                  InitData.MediaSize=InitData.nStepsT=InitData.nStepsZ=InitData.nStepsX=0;
                  InitData.MediaSizeB=InitData.nStepsZB=InitData.nStepsDelay=0;  }
  ~cEquOptic2D_01() { delete [] hMedia; delete [] hMediaNew; delete [] hMediaBragg;
                   delete [] hMediaBraggNew; delete [] hFFTSrc; delete [] hFFTDest;
                   delete [] hApData; delete [] hAmData;
                   hMedia=hMediaNew=NULL;
                   hMediaBragg=hMediaBraggNew=NULL;
                   hFFTSrc=hFFTDest=hApData=hAmData=NULL;
                   InitData.MediaSize=InitData.nStepsT=InitData.nStepsZ=InitData.nStepsX=0;
                   InitData.MediaSizeB=InitData.nStepsZB=InitData.nStepsDelay=0; }
  void StepRoutine();
  void FillBorders();
  void PrepareMedia();

  TCplxLong BorderAZm(int i) {
      return( hMedia[Index(1,i)].AZp*polar(InitData.Z0RefCoef,InitData.Z0RefPhase) ); }
  TCplxLong BorderAZp(int i) {
      return( hApData[Index((InitData.DelayPos+1)%InitData.nStepsDelay,i)] ); }
  TCplxLong BorderAZmBragg(int i) {
      return( hAmData[Index((InitData.DelayPos+1)%InitData.nStepsDelay,i)] ); }
  TCplxLong BorderAZpBragg(int i) {
      return( 0 ); }
  TCplxLong BorderAXmBragg(int i) { if (InitData.LockX)
      return( hMediaBragg[Index(i,InitData.nStepsX)].AXm ); else return(0); }
  TCplxLong BorderAXpBragg(int i) { if (InitData.LockX)
      return( hMediaBragg[Index(i,1)].AXp ); else return(0); }
  void CalcEnergy();
  void CalcOutEnergy();

  void InitReport(cDifEquReport *hReport);
  void Report(cDifEquReport *hReport);
  void FinalReport(cDifEquReport *hReport);

  void SaveAData();

  virtual void LoadInitData(cDifEqu_InitDataRec *hInitData);

  void SetHDump(cDump2DOpticModel *phDump) { hDump=phDump;}
  inline int Index(int iz,int ix)
    {return(InitData.SizeX*iz+ix);}

//  virtual void FillAData( TCplxLong *hP );

  TCplxLong RightPartP(TCplxLong P,long double R, TCplxLong A);
  long double RightPartR(TCplxLong P,long double R, TCplxLong A);

protected:
long double EnergyZP,EnergyZM,EnergyZPlast,ESumLast;
long double EnergyXPB,EnergyXMB,EnergyZPB,EnergyZMB,EnergyZPT,EnergyZMT;
long double decZP,LastReportTime,LastDecTime;
long double EnergyR0;
long double EnergyXPOut,EnergyXMOut,EnergyZPOut,EnergyZMOut;
cOptic2D_Media_2Waves *hMedia,*hMediaNew;
cOptic2D_Media_BraggOnly *hMediaBragg,*hMediaBraggNew;
TCplxLong *hApData,*hAmData;
cDifEqu_InitDataRec_Optic2D InitData;
TCplxLong *hFFTSrc,*hFFTDest;
CFFTrans<long double> FFT;
int FFTSizeFull;
cDump2DOpticModel *hDump;
};


#endif
