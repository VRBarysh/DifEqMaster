//---------------------------------------------------------------------------

#ifndef EquOptic2DH
#define EquOptic2DH

#include "EquBaseMaster.h"
#include "FFTUnit.h"
#include "Dump2DOpticModel.h"
//---------------------------------------------------------------------------

class cEquOptic2D : public cDifEqu_DataBankIO
{
public:
  cEquOptic2D() : cDifEqu_DataBankIO() { hMedia=hMediaNew=hMediaK=NULL; hFFTSrc=hFFTDest=NULL;
                  InitData.MediaSize=InitData.nStepsT=InitData.nStepsZ=InitData.nStepsX=0;
                  hAbsStore=NULL; hMirrorDelayBuf=NULL; hEnergyData=NULL;
                  nReportNumber=nStepNumber=0; }
  ~cEquOptic2D() { delete [] hMedia; delete [] hMediaNew; delete [] hMediaK;
                   delete [] hFFTSrc; delete [] hFFTDest; delete [] hAbsStore;
                   delete [] hMirrorDelayBuf; delete [] hEnergyData;
                   hMedia=hMediaNew=hMediaK=NULL; hFFTSrc=hFFTDest=NULL; hAbsStore=NULL;
                   hMirrorDelayBuf=NULL; hEnergyData=NULL;
                   InitData.MediaSize=InitData.nStepsT=InitData.nStepsZ=InitData.nStepsX=0; }
  void StepRoutine();
  void FillBorders();
  void PrepareMedia();
  void CalcMirrorDelay();

//  TCplxLong BorderAZm(int i) { if (0)
//      return( hMedia[Index(InitData.nStepsZ,i)].AZm );
//      else return(InitData.Z0RefCoef*hMedia[Index(1,i)].AZp); }
//      { return(polar(1,(long double)i*M_PI*0.25/InitData.nStepsX)); }
  TCplxLong BorderAZm(int i) {
      if(InitData.UseMirrorDelay)
        return( polar(InitData.Z0RefCoef,InitData.Z0RefPhase)*hMirrorDelayBuf[Index((InitData.DelayPos+1)%InitData.nStepsDelay,i)] );
      else return(polar(InitData.Z0RefCoef,InitData.Z0RefPhase)*hMedia[Index(1,i)].AZp); }
  TCplxLong BorderAZp(int i) { if (0)
      return( hMedia[Index(1,i)].AZp );
      else return(InitData.Z1RefCoef*hMedia[Index(InitData.nStepsZ,i)].AZm); }
  TCplxLong BorderAXm(int i) { if (InitData.LockX)
      return( hMedia[Index(i,InitData.nStepsX)].AXm ); else return(0); }
  TCplxLong BorderAXp(int i) { if (InitData.LockX)
      return( hMedia[Index(i,1)].AXp ); else return(0); }

  int Get2DColor(int Number2D, double x, double y);
  void Prepare2DColor(int Number2D);

  void CalcEnergy();
  void CalcOutEnergy();

  void InitReport(cDifEquReport *hReport);
  void Report(cDifEquReport *hReport);
  void FinalReport(cDifEquReport *hReport);

  virtual void SaveAData();

  virtual void LoadInitData(cDifEqu_InitDataRec *hInitData);

  void SetHDump(cDump2DOpticModel *phDump) { hDump=phDump;}
  inline int Index(int iz,int ix)
    {return(InitData.SizeX*iz+ix);}

//  virtual void FillAData( TCplxLong *hP );

  TCplxLong RightPartP(TCplxLong P,long double R, TCplxLong A);
  long double RightPartR(TCplxLong P,long double R, TCplxLong A);

protected:
TCplxLong *A,*P,*Pk1,*Pk2,*Pnew;
long double *R,*Rk1,*Rk2,*Rnew;
long double EnergyXP,EnergyXM,EnergyZP,EnergyZM,EnergyZPlast,ESumLast,EDelayBuf;
long double *hEnergyData;
long double decZP,LastReportTime,LastDecTime,DiffEnergy,EDiffOut,EWrong;
long double EnergyR0,EnergyRZ,EnergyRX,EnergyRp,EnergyRm;
long double EnergyXPOut,EnergyXMOut,EnergyZPOut,EnergyZMOut;
long double maxAZm,maxAZp,maxAXp,maxAXm;
cOptic2D_Media *hMedia,*hMediaNew,*hMediaK;
cDifEqu_InitDataRec_Optic2D InitData;
TCplxLong *hFFTSrc,*hFFTDest,*hMirrorDelayBuf;
CFFTrans<long double> FFT;
int FFTSizeFull;
cDump2DOpticModel *hDump;
cOptic2D_AbsStore *hAbsStore;
};

#endif
