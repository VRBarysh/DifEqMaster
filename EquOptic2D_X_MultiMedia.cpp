//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquOptic2D_X_MultiMedia.h"
#include "EquOpticTask_MultiMedia.h"
#include <fstream.h>

//---------------------------------------------------------------------------

#pragma package(smart_init)

void cEquOptic2D_X_MultiMedia::Report(cDifEquReport *hReport) {
  nReportNumber++;

  for(int iz=1; iz<=InitData.nStepsZ; iz++)
    hReportBuf[iz-1]=abs(hAZp[IndexXZ(InitData.SizeX/2,iz)]);
  hReport->Graph[0].GetData(0,InitData.nStepsZ,hReportBuf);
  for(int ix=1; ix<=InitData.nStepsX; ix++)
    hReportBuf[ix-1]=abs(hAZp[IndexXZ(ix,InitData.SizeZ/2)]);
  hReport->Graph[1].GetData(0,InitData.nStepsX,hReportBuf);
  for(int iz=1; iz<=InitData.nStepsZ; iz++)
    hReportBuf[iz-1]=abs(hAXm[IndexXZ(InitData.SizeX/2,iz)]);
  hReport->Graph[2].GetData(0,InitData.nStepsZ,hReportBuf);
/*  for(int ix=1; ix<=InitData.nStepsX; ix++)
    hReportBuf[ix-1]=abs(hAXm[IndexXZ(ix,InitData.SizeZ/2)]);
  hReport->Graph[3].GetData(0,InitData.nStepsX,hReportBuf);*/

  for(int im=0; im<InitData.nMedia; im++)
    hReportBuf[im]= MediaDistrNoNorm(im);
  hReport->Graph[3].GetData(0,InitData.nMedia,hReportBuf);

  hReport->nVals=1;
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
  CalcEnergy();
  hReport->nVals=1;
  OPTIC_TYPE EnergyA=EnergyXP+EnergyXM+EnergyZP+EnergyZM;
  OPTIC_TYPE dtReport=1,decA,decZP;
  if(t-LastReportTime>0.0000001) dtReport=t-LastReportTime;
  if((ESumLast>       0.00000000000000000000001)&&
     (EnergyZPlast>   0.00000000000000000000001)&&
     (EnergyA>        0.00000000000000000000001)&&
     (EnergyZP>       0.00000000000000000000001)) {
    decA=0.5*log(EnergyA/ESumLast)/dtReport;
    decZP=0.5*log(EnergyZP/EnergyZPlast)/dtReport;
  } else decA=decZP=0;

  hReport->ValText[hReport->nVals]="Average R ";
  hReport->Val[hReport->nVals++]=EnergyR/InitData.nStepsX/InitData.nStepsZ/InitData.dt/InitData.dt;

  hReport->ValText[hReport->nVals]="Dec A "; hReport->Val[hReport->nVals++]=decA;
  hReport->ValText[hReport->nVals]="Dec ZP "; hReport->Val[hReport->nVals++]=decZP;

  hReport->ValText[hReport->nVals]="Freq ";
  hReport->Val[hReport->nVals++]=Frequency;

  hReport->ValText[hReport->nVals]="Media Norm ";
  hReport->Val[hReport->nVals++]=MediaDistrNormCoef();
  hReport->ValText[hReport->nVals]="Media Max ";
  hReport->Val[hReport->nVals++]=MediaDistrNoNorm(InitData.nMedia/2)*MediaDistrNormCoef();
  lastFrequency = Frequency;

  if(InitData.Coaxial) {
    hReport->ValText[hReport->nVals]="EoutZ/EGen ";
    hReport->Val[hReport->nVals++]=(hEnergyData[nStepNumber-1].OutZp+hEnergyData[nStepNumber-1].OutZm)/EnergyA/InitData.AGenCoef;
  }

  if(nStepNumber) {
    hReport->ValText[hReport->nVals]="OutZp ";
    hReport->Val[hReport->nVals++]=hEnergyData[nStepNumber-1].OutZp/InitData.dt;
    hReport->ValText[hReport->nVals]="OutZ2 ";
    hReport->Val[hReport->nVals++]=hEnergyData[nStepNumber-1].OutZ2/InitData.dt;
    hReport->ValText[hReport->nVals]="OutXp ";
    hReport->Val[hReport->nVals++]=hEnergyData[nStepNumber-1].OutXp/InitData.dt;
    hReport->ValText[hReport->nVals]="OutXp2 ";
    hReport->Val[hReport->nVals++]=hEnergyData[nStepNumber-1].OutXp2/InitData.dt;
  }

  FFT.TransForward(hFFTDest,hFFTSrc);
  complex<OPTIC_TYPE> tmp;
  for(int i=0;i<(FFTSizeFull >> 1);i++) {
    tmp=hFFTDest[i];
    hFFTDest[i]=hFFTDest[i+(FFTSizeFull >> 1)];
    hFFTDest[i+(FFTSizeFull >> 1)]=tmp;
  }
  hReport->Graph[4].GetDataAbs(0,(FFTSizeFull >> InitData.FFTScale),
             hFFTDest+(FFTSizeFull>>1)-(FFTSizeFull>> (1+InitData.FFTScale) ) );
  hReport->Graph[5].GetDataAbs(0,FFTSizeFull,hFFTSrc);

  hReport->Graph[6].Number2D=0;
  hReport->Graph[6].Size=hReport->Graph[6].MemSize;
  hReport->Graph[7].Number2D=1;
  hReport->Graph[7].Size=hReport->Graph[7].MemSize;
  hReport->Graph[8].Number2D=2;
  hReport->Graph[8].Size=hReport->Graph[8].MemSize;
  hReport->Graph[9].Number2D=3;
  hReport->Graph[9].Size=hReport->Graph[9].MemSize;
  hReport->Graph[10].Number2D=4;
  hReport->Graph[10].Size=hReport->Graph[10].MemSize;
  hReport->Graph[11].Number2D=5;
  hReport->Graph[11].Size=hReport->Graph[11].MemSize;
  hReport->Graph[12].Number2D=6;
  hReport->Graph[12].Size=hReport->Graph[12].MemSize;
  hReport->Graph[13].Number2D=7;
  hReport->Graph[13].Size=hReport->Graph[13].MemSize;


  EnergyZPlast=EnergyZP; ESumLast=EnergyA; LastReportTime=t;

}

void cEquOptic2D_X_MultiMedia::FinalReport(cDifEquReport *hReport) {
  hReport->nVals=0;
  hReport->ValText[hReport->nVals]="Average R ";
  hReport->Val[hReport->nVals++]=EnergyR/InitData.nStepsX/InitData.nStepsZB/InitData.dt/InitData.dt;
//  hReport->ValText[hReport->nVals]="Energy wrong ";
//  OPTIC_TYPE EOut=EnergyXPOut+EnergyXMOut+EnergyZPOut+EnergyZMOut+EnergyUPOut;
//  hReport->Val[hReport->nVals++]=(EnergyR0-EnergyR)+EnergyQ+(EnergyA0-EnergyA)-EOut;
  if(nStepNumber) {
    hReport->ValText[hReport->nVals]="OutZp ";
    hReport->Val[hReport->nVals++]=hEnergyData[nStepNumber-1].OutZp/InitData.dt;
    hReport->ValText[hReport->nVals]="OutZ2 ";
    hReport->Val[hReport->nVals++]=hEnergyData[nStepNumber-1].OutZ2/InitData.dt;
    hReport->ValText[hReport->nVals]="OutXp ";
    hReport->Val[hReport->nVals++]=hEnergyData[nStepNumber-1].OutXp/InitData.dt;
    hReport->ValText[hReport->nVals]="OutXp2 ";
    hReport->Val[hReport->nVals++]=hEnergyData[nStepNumber-1].OutXp2/InitData.dt;
  }

  hReport->nGraphs=2;
  for(int i=0;i<nStepNumber;i++) hReportBuf[i]=hEnergyData[i].Out;
  hReport->Graph[0].GetData(0,nStepNumber,hReportBuf);
  for(int i=0;i<nStepNumber;i++) hReportBuf[i]=hEnergyData[i].R;
  hReport->Graph[1].GetData(0,nStepNumber,hReportBuf);

  hReport->ValText[hReport->nVals]="Max P "; hReport->Val[hReport->nVals++]=hReport->Graph[0].Max;
  if(InitData.Coaxial) {
    hReport->ValText[hReport->nVals]="EoutZ/EGen ";
    hReport->Val[hReport->nVals++]=(hEnergyData[nStepNumber-1].OutZp+hEnergyData[nStepNumber-1].OutZm)/(EnergyZP+EnergyZM)/InitData.AGenCoef;
  }
  FinalSave();
}

void cEquOptic2D_X_MultiMedia::FinalSave() {
  long double tmp1;
  ofstream FileOutAZp = ofstream("AZp.dat",ios::out);
  //ofstream FileOutAZp_ph = ofstream("AZP_ph.dat",ios::out);
  //ofstream FileOutAZm_ph = ofstream("AZM_ph.dat",ios::out);
  ofstream FileOutAXp = ofstream("AXp.dat",ios::out);
  ofstream FileOutAZm = ofstream("AZm.dat",ios::out);
  ofstream FileOutAXm = ofstream("AXm.dat",ios::out);
  ofstream FileOutAXMicro = ofstream("AXMicro.dat",ios::out);
  ofstream FileOutAZMicro = ofstream("AZMicro.dat",ios::out);
  ofstream FileOutAXArg = ofstream("AXArg.dat",ios::out);
  ofstream FileOutAZArg = ofstream("AZArg.dat",ios::out);
  ofstream FileOutAXpt = ofstream("AXpA(t).dat",ios::out);
  ofstream FileOutAZpt = ofstream("AZpA(t).dat",ios::out);
  ofstream FileOutAXpArgt = ofstream("AXpArg(t).dat",ios::out);
  ofstream FileOutAZpArgt = ofstream("AZpArg(t).dat",ios::out);
  ofstream FileOutR0 = ofstream("R0.dat",ios::out);
  ofstream FileOutRZ = ofstream("RZ.dat",ios::out);
  ofstream FileOutRX = ofstream("RX.dat",ios::out);
  ofstream FileOutRp = ofstream("Rp.dat",ios::out);
  ofstream FileOutRm = ofstream("Rm.dat",ios::out);
  ofstream FileOutRpACoef = ofstream("RpACoef.dat",ios::out);
  ofstream FileOutRmACoef = ofstream("RmACoef.dat",ios::out);
  //ofstream FileOutEt = ofstream("Et.dat",ios::out);
  //ofstream FileOutEt2= ofstream("Et2.dat",ios::out);
  ofstream FileOutEtOut = ofstream("EtOut.dat",ios::out);
  ofstream FileOutEtOut2= ofstream("EtOut2.dat",ios::out);
  ofstream FileOutRt = ofstream("Rt.dat",ios::out);
  ofstream FileOutRtNorm = ofstream("RtNorm.dat",ios::out);
  ofstream FileOutAt = ofstream("At.dat",ios::out);
  ofstream FileOutEXp = ofstream("EOutXp.dat",ios::out);
  ofstream FileOutEXm = ofstream("EOutXm.dat",ios::out);
  ofstream FileOutEZp = ofstream("EOutZp.dat",ios::out);
  ofstream FileOutEZm = ofstream("EOutZm.dat",ios::out);
  ofstream FileOutSpectrum = ofstream("Spectrum.dat",ios::out);
  complex<OPTIC_TYPE> Ap,Am;
  for(int i=0;i<InitData.nStepsZ;i++)
    for(int j=0;j<InitData.nStepsX;j++) {
      FileOutAZm << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hAZm[Index_ZX(i+1,j+1)]) << endl;
      FileOutAXm << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hAXm[Index_ZX(i+1,j+1)]) << endl;
      FileOutAZp << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hAZp[Index_ZX(i+1,j+1)]) << endl;
      FileOutAXp << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hAXp[Index_ZX(i+1,j+1)]) << endl;

      Ap=hAXp[Index_ZX(i+1,j+1)]; Am=hAXm[Index_ZX(i+1,j+1)];

      FileOutAXMicro << InitData.dt*i << " " << InitData.dt*j
      /* C */    << " " << abs(Am)+abs(Ap)
      /* D */    << " " << abs(Am)*abs(Am)+abs(Ap)*abs(Ap)
      /* E */    << " " << sqrt(abs(Am)*abs(Am)+abs(Ap)*abs(Ap))
      /* F */    << " " << sqrt(Am.real()*Am.real()+Am.imag()*Am.imag()+
                                Ap.real()*Ap.real()+Ap.imag()*Ap.imag())
      /* G */    << " " << sqrt(Am.real()*Am.real()+Am.imag()*Am.imag()+
                                Ap.real()*Ap.real()+Ap.imag()*Ap.imag()-(Am*Ap).real())
                 << endl;
/*
      FileOutAZMicro << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hAZm[Index_ZX(i+1,j+1)])+
                           abs(hAZp[Index_ZX(i+1,j+1)])
                 << endl;

      FileOutRX << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hM[Index_ZX(i+1,j+1)].RX) << endl;
      FileOutRZ << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hM[Index_ZX(i+1,j+1)].RZ) << endl;
      FileOutRp << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hM[Index_ZX(i+1,j+1)].Rp) << endl;
      FileOutRm << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hM[Index_ZX(i+1,j+1)].Rm) << endl;
      FileOutR0 << InitData.dt*i << " " << InitData.dt*j
                 << " " << hM[Index_ZX(i+1,j+1)].R0 << endl;

      FileOutRpACoef << InitData.dt*i << " " << InitData.dt*j << " "
                 << abs(complex<OPTIC_TYPE>(InitData.ACoef)+
                 InitData.beta*InitData.PCoef*hM[Index_ZX(i+1,j+1)].Rp) << endl;
      FileOutRpACoef << InitData.dt*i << " " << InitData.dt*j << " "
                 << abs(complex<OPTIC_TYPE>(InitData.ACoef)+
                 InitData.beta*InitData.PCoef*hM[Index_ZX(i+1,j+1)].Rm) << endl;*/
    }

//  for(int i=0;i<nReportNumber;i++) {
//    FileOutEt << i*t/nReportNumber << " " << hEnergyData[i].A << endl;
//  }
  OPTIC_TYPE tmpR;
  for(int i=0;i<nStepNumber;i++) {
    FileOutEXp << i*t/nStepNumber << " " << hEnergyData[i].OutXp << endl;
    FileOutEXm << i*t/nStepNumber << " " << hEnergyData[i].OutXm << endl;
    FileOutEZp << i*t/nStepNumber << " " << hEnergyData[i].OutZp << endl;
    FileOutEZm << i*t/nStepNumber << " " << hEnergyData[i].OutZm << endl;
    //FileOutEt  << i*t/nStepNumber << " " << hEnergyData[i].E << endl;
    tmpR=hEnergyData[i].OutXp+hEnergyData[i].OutXm+
         hEnergyData[i].OutZp+hEnergyData[i].OutZm;
    FileOutEtOut  << i*t/nStepNumber << " " << tmpR << endl;
    FileOutRt  << i*t/nStepNumber << " " << hEnergyData[i].R << endl;
    FileOutRtNorm  << i*t/nStepNumber << " " <<
      hEnergyData[i].R/(InitData.dt*InitData.dt*InitData.nStepsX*InitData.nStepsZ) << endl;
    if(InitData.RCoef>0.00000000001)
      FileOutEtOut2 << i*t/nStepNumber << " " << tmpR/InitData.RCoef << endl;

    for(int ix=0; ix<InitData.nStepsX; ix++) {
      FileOutAZpt << InitData.dt*i << " "<< InitData.dt*ix << " " <<
                abs(hAZpOutData[i*InitData.nStepsX+ix]) << endl;
      if(abs(hAZpOutData[i*InitData.nStepsX+ix])>0.000001) {
        FileOutAZpArgt << InitData.dt*i << " "<< InitData.dt*ix << " " <<
                  arg(hAZpOutData[i*InitData.nStepsX+ix]) << endl;
      } else {
        FileOutAZpArgt << InitData.dt*i << " "<< InitData.dt*ix << " "
                  << 0 <<endl;
      }
    }
    for(int iz=0; iz<InitData.nStepsZ; iz++) {
      FileOutAXpt << InitData.dt*i << " "<< InitData.dt*iz << " " <<
                abs(hAXpOutData[i*InitData.nStepsZ+iz]) << endl;
      if(abs(hAXpOutData[i*InitData.nStepsZ+iz])>0.000001) {
        FileOutAXpArgt << InitData.dt*i << " "<< InitData.dt*iz << " " <<
                  arg(hAXpOutData[i*InitData.nStepsZ+iz]) << endl;
      } else {
        FileOutAXpArgt << InitData.dt*i << " "<< InitData.dt*iz << " "
                  << 0 <<endl;
      }
    }
  }

  FFT.TransForward(hFFTDest,hFFTSrc);
  complex<OPTIC_TYPE> tmp;
  for(int i=0;i<(FFTSizeFull >> 1);i++) {
    tmp=hFFTDest[i];
    hFFTDest[i]=hFFTDest[i+(FFTSizeFull >> 1)];
    hFFTDest[i+(FFTSizeFull >> 1)]=tmp;
  }
  for(int i=0;i<FFTSizeFull;i++) {
    tmp1=(i-(FFTSizeFull>>1))*2.0*M_PI/FFTSizeFull/InitData.dt;
    FileOutSpectrum << tmp1 << " " << abs(hFFTDest[i]) << endl;
  }
  for(int i=0;i<FFTSizeFull;i++) {
    tmp1=t-(FFTSizeFull-i+1)*InitData.dt;
    FileOutAt << tmp1 << " " << abs(hFFTSrc[i]) << endl;
  }
}

void cEquOptic2D_X_MultiMedia::LoadInitData(cDifEqu_InitDataRec *hInitData) {
  cDifEqu_InitDataRec_Optic2D_T<OPTIC_TYPE> *hInit=
             (cDifEqu_InitDataRec_Optic2D_T<OPTIC_TYPE> *)hInitData;
  int newSize=(hInit->nStepsX+2)*(hInit->nStepsZ+2);
  int MaxSize = hInit->nStepsX > hInit->nStepsZ ? hInit->nStepsX+2 : hInit->nStepsZ+2;
  MaxSize = MaxSize > (hInit->tMax/hInit->dt+1) ? MaxSize : (hInit->tMax/hInit->dt+1);
  int nStepsT=hInit->tMax/hInit->dt+1;
  if(!hTaskCritSection) hTaskCritSection = new TCriticalSection;
  if(!hTaskEvent) hTaskEvent = new TEvent(NULL,true,false,"");
  delete [] hEnergyData;
  hEnergyData = new EnergyDataType[hInit->tMax/hInit->dt+1];
  delete [] hAXpOutData; delete [] hAZpOutData;
  hAXpOutData = new complex<OPTIC_TYPE>[hInit->nStepsZ*nStepsT];
  hAZpOutData = new complex<OPTIC_TYPE>[hInit->nStepsX*nStepsT];
  if( (newSize!=(InitData.nStepsX+2)*(InitData.nStepsZ+2)) || (hInit->nMedia!=InitData.nMedia) ) {
    delete [] hReportBuf;
    delete [] hAXp; delete [] hAXm; delete [] hAZp; delete [] hAZm;
    delete [] hM; delete [] hAbsStore;
    hAXp = new complex<OPTIC_TYPE>[newSize*3];
    hAXm = new complex<OPTIC_TYPE>[newSize*3];
    hAZp = new complex<OPTIC_TYPE>[newSize*3];
    hAZm = new complex<OPTIC_TYPE>[newSize*3];
    hM   = new cOptic2D_Media_X[newSize*3*hInit->nMedia];
    hReportBuf = new OPTIC_TYPE[MaxSize];
    hAbsStore = new cOptic2D_AbsStore[ newSize ];
  }
  for(int i=0;i<hInit->nThreads+1-nEqus;i++) {
    hEqus[i] = new cEquOpticTask_MultiMedia;
  }
  for(int i=0;i<hInit->nThreads-InitData.nThreads;i++) {
    hThreads[i] = new TTaskThread(true);
    hThreads[i]->Priority = tpLower;
    /*
        tpIdle	The thread executes only when the system is idle. The system will not interrupt other threads to execute a thread with tpIdle priority.
        tpLowest	The thread's priority is two points below normal.
        tpLower	The thread's priority is one point below normal.
        tpNormal	The thread has normal priority.
        tpHigher	The thread's priority is one point above normal.
        tpHighest	The thread's priority is two points above normal.
        tpTimeCritical	The thread gets highest priority.
    */
  }
  InitData=*hInit; nEqus=InitData.nThreads+1;
  InitData.SizeX=InitData.nStepsX+2;  InitData.SizeZ=InitData.nStepsZ+2;
  TransformMediaCoef();
  int xPart=InitData.nStepsX/nEqus;
  int xRest=InitData.nStepsX - nEqus*xPart;
  for(int i=0;i<nEqus;i++) {
    hEqus[i]->InitData=InitData;
    hEqus[i]->dtT=InitData.dt;
    hEqus[i]->hAXp=hAXp; hEqus[i]->hAXm=hAXm;
    hEqus[i]->hAZp=hAZp; hEqus[i]->hAZm=hAZm;
    hEqus[i]->hM=hM;
    hEqus[i]->ix1st=1+i*xPart + (i<xRest?i:xRest);
    hEqus[i]->ixLast=(i+1)*xPart + (i<xRest?(i+1):xRest);
    hEqus[i]->FillMatrixes();
  }
  for(int i=0;i<InitData.nThreads;i++) {
    hThreads[i]->Init(hEqus[i+1],hTaskEvent,hTaskCritSection,&ThreadCounter);
  }
  t=0; dt=InitData.dt; tMax=InitData.tMax;

  delete [] hFFTSrc; delete [] hFFTDest; FFTSizeFull = 1 << hInit->FFTSize;
  hFFTSrc = new complex<OPTIC_TYPE>[FFTSizeFull];
  hFFTDest = new complex<OPTIC_TYPE>[FFTSizeFull];
  for (int i=0;i<FFTSizeFull;i++) hFFTSrc[i]=hFFTDest[i]=0;
  FFT.MakeRevTable(InitData.FFTSize);

  PrepareMedia();
  FillBorders(0);
  EnergyXPOut=EnergyXMOut=EnergyZPOut=EnergyZMOut=EnergyQ=0;
  CalcEnergy(); EnergyA0=EnergyXP+EnergyXM+EnergyZP+EnergyZM; EnergyR0=EnergyR;
  ESumLast=EnergyZPlast=0; LastReportTime=0; LastMaxOut=0; EnergyUPOut=0;
  LastA=0; Frequency=100;
  nStepNumber=0; t=0; MaxESaved=0; LastBitmapSaveTime=-1; nBitmap=0;
  maxmaxAXp=maxmaxAXm=maxmaxAZp=maxmaxAZm=0;
  normAXp=normAXm=normAZp=normAZm=100;
  MakeGlobalNormBitmap=MakeComplexBitmap=0;
  PurgeNoiseTime=0;
  ComplexBitmapNorm=100;
  if(InitData.BitmapNorm) {
    MakeGlobalNormBitmap=1;
    normAXp=normAXm=normAZp=normAZm=InitData.BitmapNorm;
  }
  if(InitData.SaveBitmapInterval) PrepareBitmapFolders();
  PrepareColorData();
  }

void cEquOptic2D_X_MultiMedia::PrepareMedia() {
  OPTIC_TYPE r1,r2,r3,r4;
  for(int ix=1;ix<=InitData.nStepsX;ix++)
    for(int iz=1;iz<=InitData.nStepsZ;iz++) {
      hAXp[IndexXZ(ix,iz)]=0;
      hAXm[IndexXZ(ix,iz)]=0;
      if(!InitData.MicrowaveModeLength) {
        hAZp[IndexXZ(ix,iz)]=sin(M_PI*(iz-1)/(InitData.nStepsZ-1))*
          (InitData.A0+
           InitData.A1*sin(M_PI*(ix-1)/(InitData.nStepsX-1))+
           InitData.A2*sin(2.0*M_PI*(ix-1)/(InitData.nStepsX-1)) );
      } else hAZp[IndexXZ(ix,iz)]=0;
      hAZm[IndexXZ(ix,iz)]=0;
      if(InitData.rPolNoise) {
        if((iz%2)) {
          r1=2.0*M_PI*rand()/RAND_MAX;
          r2=2.0*M_PI*rand()/RAND_MAX;
          r3=2.0*M_PI*rand()/RAND_MAX;
          r4=2.0*M_PI*rand()/RAND_MAX;
        }
        hAXp[IndexXZ(ix,iz)]+=polar(InitData.rPolNoise,r1);
        hAXm[IndexXZ(ix,iz)]+=polar(InitData.rPolNoise,r2);
        hAZp[IndexXZ(ix,iz)]+=polar(InitData.rPolNoise,r3);
        hAZm[IndexXZ(ix,iz)]+=polar(InitData.rPolNoise,r4);
      }
      for(int im=0; im<InitData.nMedia; im++) {
        hM[IndexXZM(ix,iz,im)].PXp=0;//polar(InitData.rPolNoise,2.0*M_PI*rand()/RAND_MAX);
        hM[IndexXZM(ix,iz,im)].PXm=0;//polar(InitData.rPolNoise,2.0*M_PI*rand()/RAND_MAX);
        hM[IndexXZM(ix,iz,im)].PZp=0;//polar(InitData.rPolNoise,2.0*M_PI*rand()/RAND_MAX);
        hM[IndexXZM(ix,iz,im)].PZm=0;//polar(InitData.rPolNoise,2.0*M_PI*rand()/RAND_MAX);

        hM[IndexXZM(ix,iz,im)].RX=hM[IndexXZM(ix,iz,im)].RZ=0;
        hM[IndexXZM(ix,iz,im)].Rp=hM[IndexXZM(ix,iz,im)].Rm=0;
        hM[IndexXZM(ix,iz,im)].R0=InitData.StartingInv;
      }
    }
  complex<OPTIC_TYPE> tmpC;
  if(InitData.rPolNoise<0) {
    for(int ix=1;ix<=InitData.nStepsX;ix++)
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        tmpC=hAXp[IndexXZ(ix,iz)]*OPTIC_TYPE(4)
            +hAXp[IndexXZ(ix-1,iz)]
            +hAXp[IndexXZ(ix,iz-1)]
            +hAXp[IndexXZ(ix+1,iz)]
            +hAXp[IndexXZ(ix,iz+1)];
        hM[IndexXZM(ix,iz,InitData.nMedia/2)].PXp=tmpC*0.125;
        tmpC=hAXm[IndexXZ(ix,iz)]*OPTIC_TYPE(4)
            +hAXm[IndexXZ(ix-1,iz)]
            +hAXm[IndexXZ(ix,iz-1)]
            +hAXm[IndexXZ(ix+1,iz)]
            +hAXm[IndexXZ(ix,iz+1)];
        hM[IndexXZM(ix,iz,InitData.nMedia/2)].PXm=tmpC*0.125;

        tmpC=hAZp[IndexXZ(ix,iz)]*OPTIC_TYPE(4)
            +hAZp[IndexXZ(ix-1,iz)]
            +hAZp[IndexXZ(ix,iz-1)]
            +hAZp[IndexXZ(ix+1,iz)]
            +hAZp[IndexXZ(ix,iz+1)];
        hM[IndexXZM(ix,iz,InitData.nMedia/2)].PZp=tmpC*0.125;
        tmpC=hAZm[IndexXZ(ix,iz)]*OPTIC_TYPE(4)
            +hAZm[IndexXZ(ix-1,iz)]
            +hAZm[IndexXZ(ix,iz-1)]
            +hAZm[IndexXZ(ix+1,iz)]
            +hAZm[IndexXZ(ix,iz+1)];
        hM[IndexXZM(ix,iz,InitData.nMedia/2)].PZm=tmpC*0.125;
      }
    for(int ix=1;ix<=InitData.nStepsX;ix++)
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        hAXp[IndexXZ(ix,iz)]=hM[IndexXZM(ix,iz,InitData.nMedia/2)].PXp;
        hAXm[IndexXZ(ix,iz)]=hM[IndexXZM(ix,iz,InitData.nMedia/2)].PXm;
        hAZp[IndexXZ(ix,iz)]=hM[IndexXZM(ix,iz,InitData.nMedia/2)].PZp;
        hAZm[IndexXZ(ix,iz)]=hM[IndexXZM(ix,iz,InitData.nMedia/2)].PZm;
        hM[IndexXZM(ix,iz,InitData.nMedia/2)].PXp=hM[IndexXZM(ix,iz,InitData.nMedia/2)].PXm=0;
        hM[IndexXZM(ix,iz,InitData.nMedia/2)].PZp=hM[IndexXZM(ix,iz,InitData.nMedia/2)].PZm=0;
      }
  }
}

void cEquOptic2D_X_MultiMedia::ForceLinear() {
  if(InitData.fForceLinear)
    for(int iz=0; iz<=InitData.nStepsZ+1; iz++)
      for(int ix=0; ix<=InitData.nStepsX+1; ix++) {
        for(int im=0; im<InitData.nMedia; im++) {
          hM[IndexXZM(ix,iz,im)].R0=InitData.Q;
          hM[IndexXZM(ix,iz,im)].Rm=InitData.QGrid;
          hM[IndexXZM(ix,iz,im)].Rp=InitData.QGrid;
          hM[IndexXZM(ix,iz,im)].RX=hM[IndexXZM(ix,iz,im)].RZ=0;
        }
      }
}

void cEquOptic2D_X_MultiMedia::FillBorders(int iLayer) {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  complex<OPTIC_TYPE> phase=polar(1,InitData.Z0RefPhase*M_PI);
  if(!InitData.LockX) {
                                        // --------------- Open Conditions + reflections from borders
    for(int ix=0;ix<=InitData.nStepsX+1;ix++) {
      hAXp[iLayer*Size+IndexXZ(ix,0)]=hAXp[iLayer*Size+IndexXZ(ix,1)];
      hAXm[iLayer*Size+IndexXZ(ix,0)]=hAXm[iLayer*Size+IndexXZ(ix,1)];
      hAZp[iLayer*Size+IndexXZ(ix,0)]=InitData.Z0RefCoef*hAZm[iLayer*Size+IndexXZ(ix,1)];

      hAZm[iLayer*Size+IndexXZ(ix,0)]=0;
      hAXp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=hAXp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];
      hAXm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=hAXm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];
      hAZp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;
      hAZm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=InitData.Z0RefCoef*hAZp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];

      if(!InitData.fNoMedia) {
        for(int im=0; im<InitData.nMedia; im++) {
          hM[IndexXZM(ix,0,im)].PXp=hM[IndexXZM(ix,1,im)].PXp;
          hM[IndexXZM(ix,0,im)].PXm=hM[IndexXZM(ix,1,im)].PXm;
          hM[IndexXZM(ix,0,im)].PZp=0;
          hM[IndexXZM(ix,0,im)].PZm=0;
          hM[IndexXZM(ix,InitData.nStepsZ+1,im)].PXp=hM[IndexXZM(ix,InitData.nStepsZ,im)].PXp;
          hM[IndexXZM(ix,InitData.nStepsZ+1,im)].PXm=hM[IndexXZM(ix,InitData.nStepsZ,im)].PXm;
          hM[IndexXZM(ix,InitData.nStepsZ+1,im)].PZp=0;
          hM[IndexXZM(ix,InitData.nStepsZ+1,im)].PZm=0;
        }
      }

    }
    for(int iz=0;iz<=InitData.nStepsZ+1;iz++) {
      hAXp[iLayer*Size+IndexXZ(0,iz)]=InitData.Z0RefCoef*hAXp[iLayer*Size+IndexXZ(1,iz)];
      hAXm[iLayer*Size+IndexXZ(0,iz)]=0;
      hAZp[iLayer*Size+IndexXZ(0,iz)]=hAZp[iLayer*Size+IndexXZ(1,iz)];
      hAZm[iLayer*Size+IndexXZ(0,iz)]=hAZm[iLayer*Size+IndexXZ(1,iz)];
      hAXp[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=0;
      hAXm[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=InitData.Z0RefCoef*hAXp[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];;
      hAZp[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=hAZp[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
      hAZm[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=hAZm[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];

      if(!InitData.fNoMedia) {
        for(int im=0; im<InitData.nMedia; im++) {
          hM[IndexXZM(0,iz,im)].PXp=0;
          hM[IndexXZM(0,iz,im)].PXm=0;
          hM[IndexXZM(0,iz,im)].PZp=hM[IndexXZM(1,iz,im)].PZp;
          hM[IndexXZM(0,iz,im)].PZm=hM[IndexXZM(1,iz,im)].PZm;
          hM[IndexXZM(InitData.nStepsX+1,iz,im)].PXp=0;
          hM[IndexXZM(InitData.nStepsX+1,iz,im)].PXm=0;
          hM[IndexXZM(InitData.nStepsX+1,iz,im)].PZp=hM[IndexXZM(InitData.nStepsX,iz,im)].PZp;
          hM[IndexXZM(InitData.nStepsX+1,iz,im)].PZm=hM[IndexXZM(InitData.nStepsX,iz,im)].PZm;
        }
      }

    }
  } else {
    if(InitData.Coaxial) {
                                            // --------------- Coaxial Conditions
      for(int ix=0;ix<=InitData.nStepsX+1;ix++) {
        hAXp[iLayer*Size+IndexXZ(ix,0)]=0;//hAXp[iLayer*Size+IndexXZ(ix,1)];
        hAXm[iLayer*Size+IndexXZ(ix,0)]=0;//hAXm[iLayer*Size+IndexXZ(ix,1)];
        hAZp[iLayer*Size+IndexXZ(ix,0)]=0;
        hAZm[iLayer*Size+IndexXZ(ix,0)]=0;
        hAXp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=hAXp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];
        hAXm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=hAXm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];
        hAZp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;
        hAZm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;
        if(!InitData.fNoMedia) {
          for(int im=0; im<InitData.nMedia; im++) {
            hM[IndexXZM(ix,0,im)].PXp=0;//hM[iLayer*Size+IndexXZ(ix,1)].PXp;
            hM[IndexXZM(ix,0,im)].PXm=0;//hM[iLayer*Size+IndexXZ(ix,1)].PXm;
            hM[IndexXZM(ix,0,im)].PZp=0;
            hM[IndexXZM(ix,0,im)].PZm=0;
            hM[IndexXZM(ix,InitData.nStepsZ+1,im)].PXp=0;//hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)].PXp;
            hM[IndexXZM(ix,InitData.nStepsZ+1,im)].PXm=0;//hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)].PXm;
            hM[IndexXZM(ix,InitData.nStepsZ+1,im)].PZp=0;
            hM[IndexXZM(ix,InitData.nStepsZ+1,im)].PZm=0;
          }
        }
      }
      for(int iz=0;iz<=InitData.nStepsZ+1;iz++) {
        hAXp[iLayer*Size+IndexXZ(0,iz)]=hAXp[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
        hAXm[iLayer*Size+IndexXZ(0,iz)]=hAXm[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
        hAZp[iLayer*Size+IndexXZ(0,iz)]=hAZp[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
        hAZm[iLayer*Size+IndexXZ(0,iz)]=hAZm[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
        hAXp[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=hAXp[iLayer*Size+IndexXZ(1,iz)];
        hAXm[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=hAXm[iLayer*Size+IndexXZ(1,iz)];
        hAZp[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=hAZp[iLayer*Size+IndexXZ(1,iz)];
        hAZm[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=hAZm[iLayer*Size+IndexXZ(1,iz)];
        if(!InitData.fNoMedia) {
          for(int im=0; im<InitData.nMedia; im++) {
            hM[IndexXZM(0,iz,im)].PXp=hM[IndexXZM(InitData.nStepsX,iz,im)].PXp;
            hM[IndexXZM(0,iz,im)].PXm=hM[IndexXZM(InitData.nStepsX,iz,im)].PXm;
            hM[IndexXZM(0,iz,im)].PZp=hM[IndexXZM(InitData.nStepsX,iz,im)].PZp;
            hM[IndexXZM(0,iz,im)].PZm=hM[IndexXZM(InitData.nStepsX,iz,im)].PZm;
            hM[IndexXZM(InitData.nStepsX+1,iz,im)].PXp=hM[IndexXZM(1,iz,im)].PXp;
            hM[IndexXZM(InitData.nStepsX+1,iz,im)].PXm=hM[IndexXZM(1,iz,im)].PXm;
            hM[IndexXZM(InitData.nStepsX+1,iz,im)].PZp=hM[IndexXZM(1,iz,im)].PZp;
            hM[IndexXZM(InitData.nStepsX+1,iz,im)].PZm=hM[IndexXZM(1,iz,im)].PZm;
          }
        }
      }
    } else {
                                                // --------------- Mirror Conditions
                        // --------------- Z0RefCoef, Z1RefCoef along AXpm borders
      for(int ix=0;ix<=InitData.nStepsX+1;ix++) {
        hAXp[iLayer*Size+IndexXZ(ix,0)]=hAXp[iLayer*Size+IndexXZ(ix,1)];
        hAXm[iLayer*Size+IndexXZ(ix,0)]=hAXm[iLayer*Size+IndexXZ(ix,1)];
        hAZp[iLayer*Size+IndexXZ(ix,0)]=0;
        hAZm[iLayer*Size+IndexXZ(ix,0)]=0;
        hAXp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=hAXp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];
        hAXm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=hAXm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];
        hAZp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;
        hAZm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;
        if(!InitData.fNoMedia) {
          for(int im=0; im<InitData.nMedia; im++) {
            hM[IndexXZM(ix,0,im)].PXp=hM[IndexXZM(ix,1,im)].PXp;
            hM[IndexXZM(ix,0,im)].PXm=hM[IndexXZM(ix,1,im)].PXm;
            hM[IndexXZM(ix,0,im)].PZp=0;
            hM[IndexXZM(ix,0,im)].PZm=0;
            hM[IndexXZM(ix,InitData.nStepsZ+1,im)].PXp=hM[IndexXZM(ix,InitData.nStepsZ,im)].PXp;
            hM[IndexXZM(ix,InitData.nStepsZ+1,im)].PXm=hM[IndexXZM(ix,InitData.nStepsZ,im)].PXm;
            hM[IndexXZM(ix,InitData.nStepsZ+1,im)].PZp=0;
            hM[IndexXZM(ix,InitData.nStepsZ+1,im)].PZm=0;
          }
        }
      }
      for(int iz=0;iz<=InitData.nStepsZ+1;iz++) {
        hAXp[iLayer*Size+IndexXZ(0,iz)]=InitData.Z0RefCoef*phase*hAXm[iLayer*Size+IndexXZ(1,iz)];
        hAXm[iLayer*Size+IndexXZ(0,iz)]=InitData.Z0RefCoef*phase*hAXp[iLayer*Size+IndexXZ(1,iz)];
        hAZp[iLayer*Size+IndexXZ(0,iz)]=InitData.Z0RefCoef*phase*hAZm[iLayer*Size+IndexXZ(1,iz)];
        hAZm[iLayer*Size+IndexXZ(0,iz)]=InitData.Z0RefCoef*phase*hAZp[iLayer*Size+IndexXZ(1,iz)];
        hAXp[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=InitData.Z1RefCoef*phase*hAXm[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
        hAXm[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=InitData.Z1RefCoef*phase*hAXp[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
        hAZp[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=InitData.Z1RefCoef*phase*hAZm[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
        hAZm[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=InitData.Z1RefCoef*phase*hAZp[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
        if(!InitData.fNoMedia) {
          for(int im=0; im<InitData.nMedia; im++) {
            hM[IndexXZM(0,iz,im)].PXp=InitData.Z0RefCoef*phase*hM[IndexXZM(1,iz,im)].PXm;
            hM[IndexXZM(0,iz,im)].PXm=InitData.Z0RefCoef*phase*hM[IndexXZM(1,iz,im)].PXp;
            hM[IndexXZM(0,iz,im)].PZp=InitData.Z0RefCoef*phase*hM[IndexXZM(1,iz,im)].PZm;
            hM[IndexXZM(0,iz,im)].PZm=InitData.Z0RefCoef*phase*hM[IndexXZM(1,iz,im)].PZp;
            hM[IndexXZM(InitData.nStepsX+1,iz,im)].PXp=InitData.Z1RefCoef*phase*hM[IndexXZM(InitData.nStepsX,iz,im)].PXm;
            hM[IndexXZM(InitData.nStepsX+1,iz,im)].PXm=InitData.Z1RefCoef*phase*hM[IndexXZM(InitData.nStepsX,iz,im)].PXp;
            hM[IndexXZM(InitData.nStepsX+1,iz,im)].PZp=InitData.Z1RefCoef*phase*hM[IndexXZM(InitData.nStepsX,iz,im)].PZm;
            hM[IndexXZM(InitData.nStepsX+1,iz,im)].PZm=InitData.Z1RefCoef*phase*hM[IndexXZM(InitData.nStepsX,iz,im)].PZp;
          }
        }
      }
    }
  }
  if(InitData.MicrowaveModeLength) {
    for(int ix=0;ix<=InitData.nStepsX+1;ix++) {
      hAZp[iLayer*Size+IndexXZ(ix,0)]=MicrowaveModePulse();
    }
  }
}

OPTIC_TYPE cEquOptic2D_X_MultiMedia::MediaDistrNoNorm(int im) {
  int tm=InitData.nMedia/2;
  OPTIC_TYPE tmpC=InitData.DeltaCoef/InitData.MediaGaussCoef;
  return exp(-tmpC*tmpC*(im-tm)*(im-tm));
}

OPTIC_TYPE cEquOptic2D_X_MultiMedia::MediaDistrNormCoef() {
  OPTIC_TYPE tmp=0, tmpC=InitData.DeltaCoef/InitData.MediaGaussCoef;
  int tm=InitData.nMedia/2;
  for(int im=0; im<InitData.nMedia; im++) tmp+=exp(-tmpC*tmpC*(im-tm)*(im-tm));
  return 1.0/tmp;
}

void cEquOptic2D_X_MultiMedia::TransformMediaCoef() {
  // Gauss exp(-delta*delta/MediaGaussCoef/MediaGaussCoef) ends at exp(-4)
  int nm = InitData.nMedia/2;
  if(InitData.nMedia>=3) {
    InitData.DeltaCoef = InitData.MediaGaussCoef/nm;
  } else InitData.DeltaCoef = 0;
}

void cEquOptic2D_X_MultiMedia::CalcEnergy() {
  int iI;
  EnergyXP=EnergyXM=EnergyZP=EnergyZM=EnergyR=EnergyQ=0;
  return;                                                                               // !!!!!!!!!!!!!!!!!!!!!!!!!!! return!!!
  for(int ix=1;ix<=InitData.nStepsX;ix++)
    for(int iz=1;iz<=InitData.nStepsZ;iz++) {

      iI=IndexXZ(ix,iz);
      EnergyXP+=hAXp[iI].real()*hAXp[iI].real()+hAXp[iI].imag()*hAXp[iI].imag();
      EnergyXM+=hAXm[iI].real()*hAXm[iI].real()+hAXm[iI].imag()*hAXm[iI].imag();
      EnergyZP+=hAZp[iI].real()*hAZp[iI].real()+hAZp[iI].imag()*hAZp[iI].imag();
      EnergyZM+=hAZm[iI].real()*hAZm[iI].real()+hAZm[iI].imag()*hAZm[iI].imag();
      for(int im=0; im<InitData.nMedia; im++) {
        EnergyR+=hM[IndexXZM(ix,iz,im)].R0*MediaDistrNoNorm(im);
        EnergyQ+=hM[IndexXZM(ix,iz,im)].R0Generated*MediaDistrNoNorm(im);
      }
    }
  EnergyXP*=0.5*InitData.dt*InitData.dt;
  EnergyXM*=0.5*InitData.dt*InitData.dt;
  EnergyZP*=0.5*InitData.dt*InitData.dt;
  EnergyZM*=0.5*InitData.dt*InitData.dt;
  EnergyR*=InitData.dt*InitData.dt*MediaDistrNormCoef();
  EnergyQ*=InitData.dt*InitData.dt*MediaDistrNormCoef();
}

void cEquOptic2D_X_MultiMedia::CalcOutEnergy() {
  int iI;
  OPTIC_TYPE tmpXP=0,tmpXM=0,tmpZP=0,tmpZM=0,tmpGenZ=0,tmpGenX=0,tmpR=0;
  for(int ix=1;ix<=InitData.nStepsX;ix++) {
    tmpZP+= hAZp[IndexXZ(ix,InitData.nStepsZ)].real()
           *hAZp[IndexXZ(ix,InitData.nStepsZ)].real()
           +hAZp[IndexXZ(ix,InitData.nStepsZ)].imag()
           *hAZp[IndexXZ(ix,InitData.nStepsZ)].imag();
    tmpZM+=  hAZm[IndexXZ(ix,1)].real()
            *hAZm[IndexXZ(ix,1)].real()
            +hAZm[IndexXZ(ix,1)].imag()
            *hAZm[IndexXZ(ix,1)].imag();
  }
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    tmpXP+= hAXp[IndexXZ(InitData.nStepsX,iz)].real()
           *hAXp[IndexXZ(InitData.nStepsX,iz)].real()
           +hAXp[IndexXZ(InitData.nStepsX,iz)].imag()
           *hAXp[IndexXZ(InitData.nStepsX,iz)].imag();
    tmpXM+= hAXm[IndexXZ(1,iz)].real()
           *hAXm[IndexXZ(1,iz)].real()
           +hAXm[IndexXZ(1,iz)].imag()
           *hAXm[IndexXZ(1,iz)].imag();
  }
  if((InitData.AGenZ)||(InitData.AGenX)) {
    for(int ix=1;ix<=InitData.nStepsX;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        tmpGenZ+=hAZp[IndexXZ(ix,iz)].real()
                *hAZp[IndexXZ(ix,iz)].real()
                +hAZp[IndexXZ(ix,iz)].imag()
                *hAZp[IndexXZ(ix,iz)].imag()
                +hAZm[IndexXZ(ix,iz)].real()
                *hAZm[IndexXZ(ix,iz)].real()
                +hAZm[IndexXZ(ix,iz)].imag()
                *hAZm[IndexXZ(ix,iz)].imag();
        tmpGenX+=hAXp[IndexXZ(ix,iz)].real()
                *hAXp[IndexXZ(ix,iz)].real()
                +hAXp[IndexXZ(ix,iz)].imag()
                *hAXp[IndexXZ(ix,iz)].imag()
                +hAXm[IndexXZ(ix,iz)].real()
                *hAXm[IndexXZ(ix,iz)].real()
                +hAXm[IndexXZ(ix,iz)].imag()
                *hAXm[IndexXZ(ix,iz)].imag();
      }
    }
  }
  for(int ix=1;ix<=InitData.nStepsX;ix++)
    for(int iz=1;iz<=InitData.nStepsZ;iz++)
      for(int im=0; im<InitData.nMedia; im++) {
        iI=IndexXZM(ix,iz,im);
        tmpR+=hM[iI].R0*MediaDistrNoNorm(im);
      }

  hEnergyData[nStepNumber].R=tmpR*InitData.dt*InitData.dt*MediaDistrNormCoef();
  tmpXP*=0.5*InitData.dt*InitData.dt;  EnergyXPOut+=tmpXP;
  tmpXM*=0.5*InitData.dt*InitData.dt;  EnergyXMOut+=tmpXM;
  tmpZP*=0.5*InitData.dt*InitData.dt;  EnergyZPOut+=tmpZP;
  tmpZM*=0.5*InitData.dt*InitData.dt;  EnergyZMOut+=tmpZM;
  tmpGenZ*=-InitData.dt*InitData.dt*InitData.dt*InitData.AGenZ;
  tmpGenX*=-InitData.dt*InitData.dt*InitData.dt*InitData.AGenX;
  EnergyUPOut+=tmpGenZ+tmpGenX;
  hEnergyData[nStepNumber].OutXp=tmpXP/InitData.dt;
  hEnergyData[nStepNumber].OutXm=tmpXM/InitData.dt;
  hEnergyData[nStepNumber].OutZp=tmpZP/InitData.dt;
  hEnergyData[nStepNumber].OutZm=tmpZM/InitData.dt;
  hEnergyData[nStepNumber].OutGenZ=tmpGenZ/InitData.dt;
  hEnergyData[nStepNumber].Out=(tmpXP+tmpXM+tmpZP+tmpZM)/InitData.dt;
  if( (hEnergyData[nStepNumber].Out>InitData.MaxV[0])&&
      (LastMaxOut>hEnergyData[nStepNumber].Out)&&
      (!MaxESaved))
              SaveADataMaxE();
  LastMaxOut=hEnergyData[nStepNumber].Out;

  complex<OPTIC_TYPE> C1=hAZp[IndexXZ(InitData.nStepsX/3,InitData.nStepsZ/3)];
  complex<OPTIC_TYPE> C2;
  Frequency = 100;
  if( abs(LastA)>0.000000001 ) {
    C2 = C1/LastA;
    if( abs(C2)>0.000000001 ) Frequency = arg(C2);
  }
  LastA = C1;
}




