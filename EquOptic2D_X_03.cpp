//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquOptic2D_X_03.h"
#include <fstream.h>

//---------------------------------------------------------------------------

#pragma package(smart_init)

void cEquOptic2D_X_Mod1::InitReport(cDifEquReport *hReport) {
  int Lines[EQU_MAX_GRAPHS];
  for(int i=0;i<EQU_MAX_GRAPHS;i++) Lines[i]=1;
  Lines[0]=1; Lines[1]=1; Lines[2]=1; Lines[3]=1; Lines[4]=1;
  hReport->LoadText(hReportText);
  hReport->SetGraphSize(EQU_GRAPH_SIZE,Lines);
  hReport->SetMaxTo1();
//  CalcEnergy();
//  nReportNumber++;
//  hEnergyData[nReportNumber-1]=EnergyXP+EnergyXM+EnergyZP+EnergyZM+EDelayBuf;
}
void cEquOptic2D_X_Mod1::Report(cDifEquReport *hReport) {
//  nReportNumber++;
  for(int iz=1; iz<=InitData.nStepsZ; iz++)
    hReportBuf[iz-1]=abs(hAZp[IndexXZ(InitData.SizeX/2,iz)]);
  hReport->Graph[0].GetData(0,InitData.nStepsZ,hReportBuf);
  for(int ix=1; ix<=InitData.nStepsX; ix++)
    hReportBuf[ix-1]=abs(hAZp[IndexXZ(ix,InitData.SizeZ/2)]);
  hReport->Graph[1].GetData(0,InitData.nStepsX,hReportBuf);
  for(int iz=1; iz<=InitData.nStepsZ; iz++)
    hReportBuf[iz-1]=abs(hAXm[IndexXZ(InitData.SizeX/2,iz)]);
  hReport->Graph[2].GetData(0,InitData.nStepsZ,hReportBuf);
  for(int ix=1; ix<=InitData.nStepsX; ix++)
    hReportBuf[ix-1]=abs(hAXm[IndexXZ(ix,InitData.SizeZ/2)]);
  hReport->Graph[3].GetData(0,InitData.nStepsX,hReportBuf);
  hReport->nVals=1;
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
  CalcEnergy();
  hReport->nVals=3;
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
  hReport->ValText[1]="Energy A "; hReport->Val[1]=EnergyA;
  hReport->ValText[2]="Energy R "; hReport->Val[2]=EnergyR;
  hReport->ValText[hReport->nVals]="Energy wrong ";
  OPTIC_TYPE EOut=EnergyXPOut+EnergyXMOut+EnergyZPOut+EnergyZMOut+EnergyUPOut;
  hReport->Val[hReport->nVals++]=(EnergyR0-EnergyR)+EnergyQ+(EnergyA0-EnergyA)-EOut;
  /*
  hReport->ValText[hReport->nVals]="Energy out X ";
  hReport->Val[hReport->nVals++]=EnergyXPOut+EnergyXMOut;
  hReport->ValText[hReport->nVals]="Energy out Z ";
  hReport->Val[hReport->nVals++]=EnergyZPOut+EnergyZMOut;
  hReport->ValText[hReport->nVals]="Dec A "; hReport->Val[hReport->nVals++]=decA;
  hReport->ValText[hReport->nVals]="Dec ZP "; hReport->Val[hReport->nVals++]=decZP;
  hReport->ValText[hReport->nVals]="Energy out Up ";
  hReport->Val[hReport->nVals++]=EnergyUPOut;
  hReport->ValText[hReport->nVals]="Max Max AXp ";
  hReport->Val[hReport->nVals++]=maxmaxAXp;
  hReport->ValText[hReport->nVals]="Max Max AZp ";
  hReport->Val[hReport->nVals++]=maxmaxAZp;
  hReport->ValText[hReport->nVals]="Max Max AXm ";
  hReport->Val[hReport->nVals++]=maxmaxAXm;
  hReport->ValText[hReport->nVals]="Max Max AZm ";
  hReport->Val[hReport->nVals++]=maxmaxAZm;
*/
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

//  hReport->Graph[10].GetData(0,nReportNumber,hEnergyData);
//  hReport->SetMaxTo1();
  /*
  hReport->nVals=3;
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
  hReport->ValText[1]="Max I "; hReport->Val[1]=hReport->Graph[1].Max;
  hReport->ValText[2]="Max Gamma "; hReport->Val[2]=hReport->Graph[2].Max;
  hReport->Graph[2].Max=2; */
//  LastReportTime=t;
  EnergyZPlast=EnergyZP; ESumLast=EnergyA; LastReportTime=t;
}

void cEquOptic2D_X_Mod1::FinalReport(cDifEquReport *hReport) {
  hReport->nVals=0; hReport->nGraphs=1;
  /*
  for(int i=0;i<nStepNumber;i++) hReportBuf[i]=hEnergyData[i].Out;
  hReport->Graph[0].GetData(0,nStepNumber,hReportBuf);
  */
/*  hReport->Graph[0].GetDataAbs(0,nStepNumber,hOutA);
  hReport->GraphText[0]="Out A";
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;*/
  FinalSave();
}

void cEquOptic2D_X_Mod1::FinalSave() {
  long double tmp1;
  ofstream FileOutAZp = ofstream("AZp.dat",ios::out);
  //ofstream FileOutAZp_ph = ofstream("AZP_ph.dat",ios::out);
  //ofstream FileOutAZm_ph = ofstream("AZM_ph.dat",ios::out);
  ofstream FileOutAXp = ofstream("AXp.dat",ios::out);
  ofstream FileOutAZm = ofstream("AZm.dat",ios::out);
  ofstream FileOutAXm = ofstream("AXm.dat",ios::out);
  ofstream FileOutR0 = ofstream("R0.dat",ios::out);
  ofstream FileOutRZ = ofstream("RZ.dat",ios::out);
  ofstream FileOutRX = ofstream("RX.dat",ios::out);
  ofstream FileOutRp = ofstream("Rp.dat",ios::out);
  ofstream FileOutRm = ofstream("Rm.dat",ios::out);
  //ofstream FileOutEt = ofstream("Et.dat",ios::out);
  //ofstream FileOutEt2= ofstream("Et2.dat",ios::out);
  ofstream FileOutEtOut = ofstream("EtOut.dat",ios::out);
  ofstream FileOutEtOut2= ofstream("EtOut2.dat",ios::out);
  ofstream FileOutAt = ofstream("At.dat",ios::out);
  ofstream FileOutEXp = ofstream("EOutXp.dat",ios::out);
  ofstream FileOutEXm = ofstream("EOutXm.dat",ios::out);
  ofstream FileOutEZp = ofstream("EOutZp.dat",ios::out);
  ofstream FileOutEZm = ofstream("EOutZm.dat",ios::out);
  ofstream FileOutSpectrum = ofstream("Spectrum.dat",ios::out);
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
/*
      FileOutRX << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hM[Index_ZX(i+1,j+1)].RX) << endl;
      FileOutRZ << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hM[Index_ZX(i+1,j+1)].RZ) << endl;
      FileOutRp << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hM[Index_ZX(i+1,j+1)].Rp) << endl;
      FileOutRm << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hM[Index_ZX(i+1,j+1)].Rm) << endl;
      FileOutR0 << InitData.dt*i << " " << InitData.dt*j
                 << " " << hM[Index_ZX(i+1,j+1)].R0 << endl;*/
    }

//  for(int i=0;i<nReportNumber;i++) {
//    FileOutEt << i*t/nReportNumber << " " << hEnergyData[i].A << endl;
//  }
  OPTIC_TYPE tmpR;
  /*
  for(int i=0;i<nStepNumber;i++) {
    FileOutEXp << i*t/nStepNumber << " " << hEnergyData[i].OutXp << endl;
    FileOutEXm << i*t/nStepNumber << " " << hEnergyData[i].OutXm << endl;
    FileOutEZp << i*t/nStepNumber << " " << hEnergyData[i].OutZp << endl;
    FileOutEZm << i*t/nStepNumber << " " << hEnergyData[i].OutZm << endl;
    //FileOutEt  << i*t/nStepNumber << " " << hEnergyData[i].E << endl;
    tmpR=hEnergyData[i].OutXp+hEnergyData[i].OutXm+
         hEnergyData[i].OutZp+hEnergyData[i].OutZm;
    FileOutEtOut  << i*t/nStepNumber << " " << tmpR << endl;
    if(InitData.RCoef>0.00000000001)
      FileOutEtOut2 << i*t/nStepNumber << " " << tmpR/InitData.RCoef << endl;
  }
  */
  
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

void cEquOptic2D_X_Mod1::LoadInitData(cDifEqu_InitDataRec *hInitData) {
  cDifEqu_InitDataRec_Optic2D_T<OPTIC_TYPE> *hInit=
             (cDifEqu_InitDataRec_Optic2D_T<OPTIC_TYPE> *)hInitData;
  int newSize=(hInit->nStepsX+2)*(hInit->nStepsZ+2);
  int newSizeMedia=(hInit->nStepsX+2)*(hInit->nStepsZB+2);
  if(!hTaskCritSection) hTaskCritSection = new TCriticalSection;
  if(!hTaskEvent) hTaskEvent = new TEvent(NULL,true,false,"");
  delete [] hEnergyData;
  hEnergyData = new EnergyDataType[hInit->tMax/hInit->dt+1];
  if(newSize!=(InitData.nStepsX+2)*(InitData.nStepsZ+2)) {
    delete [] hReportBuf;
    delete [] hAXp; delete [] hAXm; delete [] hAZp; delete [] hAZm;
    delete [] hAbsStore;
    hAXp = new complex<OPTIC_TYPE>[newSize*3];
    hAXm = new complex<OPTIC_TYPE>[newSize*3];
    hAZp = new complex<OPTIC_TYPE>[newSize*3];
    hAZm = new complex<OPTIC_TYPE>[newSize*3];
    hReportBuf = new OPTIC_TYPE[newSize*3];
    hAbsStore = new cOptic2D_AbsStore[ newSize ];
  }
  if(newSizeMedia!=(InitData.nStepsX+2)*(InitData.nStepsZB+2)) {
    delete [] hM2;
    delete [] hAZpMedia; delete [] hAZmMedia;
    hM2   = new cOptic2D_Media_X_2Way[newSizeMedia*3];
    hAZpMedia = new complex<OPTIC_TYPE>[newSizeMedia*3];
    hAZmMedia = new complex<OPTIC_TYPE>[newSizeMedia*3];
  }
  for(int i=0;i<hInit->nThreads+1-nEqus;i++) {
    hEqus[i] = new cEquOpticBaseTask_2Way;
  }
  for(int i=0;i<hInit->nThreads-InitData.nThreads;i++) {
    hThreads[i] = new TTaskThread(true);
    hThreads[i]->Priority = tpLowest;
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
  int xPart=InitData.nStepsX/nEqus;
  int xRest=InitData.nStepsX - nEqus*xPart;
  for(int i=0;i<nEqus;i++) {
    hEqus[i]->InitData=InitData;
    hEqus[i]->dtT=InitData.dt;
    hEqus[i]->hAXp=hAXp; hEqus[i]->hAXm=hAXm;
    hEqus[i]->hAZp=hAZp; hEqus[i]->hAZm=hAZm;
    ((cEquOpticBaseTask_2Way *)hEqus[i])->hAZpMedia=hAZpMedia;
    ((cEquOpticBaseTask_2Way *)hEqus[i])->hAZmMedia=hAZmMedia;
    ((cEquOpticBaseTask_2Way *)hEqus[i])->hM2=hM2;
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
  nStepNumber=0; t=0; MaxESaved=0; LastBitmapSaveTime=-1; nBitmap=0;
  maxmaxAXp=maxmaxAXm=maxmaxAZp=maxmaxAZm=0;
  normAXp=normAXm=normAZp=normAZm=100;
  MakeGlobalNormBitmap=MakeComplexBitmap=0;
  ComplexBitmapNorm=100;
  if(InitData.BitmapNorm) {
    MakeGlobalNormBitmap=1;
    normAXp=normAXm=normAZp=normAZm=InitData.BitmapNorm;
  }
  if(InitData.SaveBitmapInterval) PrepareBitmapFolders();
  PrepareColorData();
}

void cEquOptic2D_X_Mod1::PrepareMedia() {
  for(int ix=1;ix<=InitData.nStepsX;ix++) {
    for(int iz=1;iz<=InitData.nStepsZ;iz++) {
      hAXp[IndexXZ(ix,iz)]=0;
      hAXm[IndexXZ(ix,iz)]=0;
      hAZp[IndexXZ(ix,iz)]=0;
      /*
      hAZp[IndexXZ(ix,iz)]=sin(M_PI*(iz-1)/(InitData.nStepsZ-1))*
        (InitData.A0+
         InitData.A1*sin(M_PI*(ix-1)/(InitData.nStepsX-1))+
         InitData.A2*sin(2.0*M_PI*(ix-1)/(InitData.nStepsX-1)) );
      */
      hAZm[IndexXZ(ix,iz)]=hAZp[IndexXZ(ix,iz)];
    }
    
    for(int iz=1;iz<=InitData.nStepsZB;iz++) {
      hAZpMedia[IndexXZ_M(ix,iz)]=sin(M_PI*(iz-1)/(InitData.nStepsZB-1))*
        (InitData.A0+
         InitData.A1*sin(M_PI*(ix-1)/(InitData.nStepsX-1))+
         InitData.A2*sin(2.0*M_PI*(ix-1)/(InitData.nStepsX-1)) );
      hAZmMedia[IndexXZ_M(ix,iz)]=0;
      hM2[IndexXZ_M(ix,iz)].PZp=polar(InitData.rPolNoise,2.0*M_PI*rand()/RAND_MAX);
      hM2[IndexXZ_M(ix,iz)].PZm=polar(InitData.rPolNoise,2.0*M_PI*rand()/RAND_MAX);
      hM2[IndexXZ_M(ix,iz)].RZ=0;
      hM2[IndexXZ_M(ix,iz)].R0=InitData.StartingInv;
    }
  }
}

void cEquOptic2D_X_Mod1::FillBorders(int iLayer) {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  int SizeMedia=(InitData.nStepsX+2)*(InitData.nStepsZB+2);
  if(!InitData.LockX) {
    for(int ix=0;ix<=InitData.nStepsX+1;ix++) {
      hAXp[iLayer*Size+IndexXZ(ix,0)]=hAXp[iLayer*Size+IndexXZ(ix,1)];
      hAXm[iLayer*Size+IndexXZ(ix,0)]=hAXm[iLayer*Size+IndexXZ(ix,1)];
      hAZp[iLayer*Size+IndexXZ(ix,0)]=hAZpMedia[iLayer*SizeMedia+IndexXZ_M(ix,InitData.nStepsZB)];
      hAZm[iLayer*Size+IndexXZ(ix,0)]=0;
      hAXp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=hAXp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];
      hAXm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=hAXm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];
      hAZp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;
      hAZm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;
/*
      if(!InitData.fNoMedia) {
        hM[iLayer*Size+IndexXZ(ix,0)].PXp=hM[iLayer*Size+IndexXZ(ix,1)].PXp;
        hM[iLayer*Size+IndexXZ(ix,0)].PXm=hM[iLayer*Size+IndexXZ(ix,1)].PXm;
        hM[iLayer*Size+IndexXZ(ix,0)].PZp=0;
        hM[iLayer*Size+IndexXZ(ix,0)].PZm=0;
        hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PXp=hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)].PXp;
        hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PXm=hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)].PXm;
        hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PZp=0;
        hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PZm=0;
      }
/*
    }
    for(int iz=0;iz<=InitData.nStepsZ+1;iz++) {
      hAXp[iLayer*Size+IndexXZ(0,iz)]=0;
      hAXm[iLayer*Size+IndexXZ(0,iz)]=0;
      hAZp[iLayer*Size+IndexXZ(0,iz)]=hAZp[iLayer*Size+IndexXZ(1,iz)];
      hAZm[iLayer*Size+IndexXZ(0,iz)]=hAZm[iLayer*Size+IndexXZ(1,iz)];
      hAXp[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=0;
      hAXm[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=0;
      hAZp[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=hAZp[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
      hAZm[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=hAZm[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
/*
      if(!InitData.fNoMedia) {
        hM[iLayer*Size+IndexXZ(0,iz)].PXp=0;
        hM[iLayer*Size+IndexXZ(0,iz)].PXm=0;
        hM[iLayer*Size+IndexXZ(0,iz)].PZp=hM[iLayer*Size+IndexXZ(1,iz)].PZp;
        hM[iLayer*Size+IndexXZ(0,iz)].PZm=hM[iLayer*Size+IndexXZ(1,iz)].PZm;
        hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PXp=0;
        hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PXm=0;
        hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PZp=hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PZp;
        hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PZm=hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PZm;
      }
*/
    }
  } else {
    if(InitData.Coaxial) {
      for(int ix=0;ix<=InitData.nStepsX+1;ix++) {
        hAXp[iLayer*Size+IndexXZ(ix,0)]=hAXp[iLayer*Size+IndexXZ(ix,1)];
        hAXm[iLayer*Size+IndexXZ(ix,0)]=hAXm[iLayer*Size+IndexXZ(ix,1)];
        hAZp[iLayer*Size+IndexXZ(ix,0)]=hAZpMedia[iLayer*SizeMedia+IndexXZ_M(ix,InitData.nStepsZB)];
        hAZm[iLayer*Size+IndexXZ(ix,0)]=0;
        hAXp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=hAXp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];
        hAXm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=hAXm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];
        hAZp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;
        hAZm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;
/*
        if(!InitData.fNoMedia) {
          hM[iLayer*Size+IndexXZ(ix,0)].PXp=hM[iLayer*Size+IndexXZ(ix,1)].PXp;
          hM[iLayer*Size+IndexXZ(ix,0)].PXm=hM[iLayer*Size+IndexXZ(ix,1)].PXm;
          hM[iLayer*Size+IndexXZ(ix,0)].PZp=0;
          hM[iLayer*Size+IndexXZ(ix,0)].PZm=0;
          hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PXp=hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)].PXp;
          hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PXm=hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)].PXm;
          hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PZp=0;
          hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PZm=0;
        }
*/
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
/*
        if(!InitData.fNoMedia) {
          hM[iLayer*Size+IndexXZ(0,iz)].PXp=hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PXp;
          hM[iLayer*Size+IndexXZ(0,iz)].PXm=hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PXm;
          hM[iLayer*Size+IndexXZ(0,iz)].PZp=hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PZp;
          hM[iLayer*Size+IndexXZ(0,iz)].PZm=hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PZm;
          hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PXp=hM[iLayer*Size+IndexXZ(1,iz)].PXp;
          hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PXm=hM[iLayer*Size+IndexXZ(1,iz)].PXm;
          hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PZp=hM[iLayer*Size+IndexXZ(1,iz)].PZp;
          hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PZm=hM[iLayer*Size+IndexXZ(1,iz)].PZm;
        }
*/
      }
    } else {
      for(int ix=0;ix<=InitData.nStepsX+1;ix++) {
        hAXp[iLayer*Size+IndexXZ(ix,0)]=hAXp[iLayer*Size+IndexXZ(ix,1)];
        hAXm[iLayer*Size+IndexXZ(ix,0)]=hAXm[iLayer*Size+IndexXZ(ix,1)];
        hAZp[iLayer*Size+IndexXZ(ix,0)]=hAZpMedia[iLayer*SizeMedia+IndexXZ_M(ix,InitData.nStepsZB)];
        hAZm[iLayer*Size+IndexXZ(ix,0)]=0;
        hAXp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=hAXp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];
        hAXm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=hAXm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];
        hAZp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;
        hAZm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;
/*
        if(!InitData.fNoMedia) {
          hM[iLayer*Size+IndexXZ(ix,0)].PXp=hM[iLayer*Size+IndexXZ(ix,1)].PXp;
          hM[iLayer*Size+IndexXZ(ix,0)].PXm=hM[iLayer*Size+IndexXZ(ix,1)].PXm;
          hM[iLayer*Size+IndexXZ(ix,0)].PZp=0;
          hM[iLayer*Size+IndexXZ(ix,0)].PZm=0;
          hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PXp=hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)].PXp;
          hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PXm=hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)].PXm;
          hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PZp=0;
          hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PZm=0;
        }
*/
      }
      for(int iz=0;iz<=InitData.nStepsZ+1;iz++) {
        hAXp[iLayer*Size+IndexXZ(0,iz)]=InitData.Z0RefCoef*hAXm[iLayer*Size+IndexXZ(1,iz)];
        hAXm[iLayer*Size+IndexXZ(0,iz)]=InitData.Z0RefCoef*hAXp[iLayer*Size+IndexXZ(1,iz)];
        hAZp[iLayer*Size+IndexXZ(0,iz)]=InitData.Z0RefCoef*hAZm[iLayer*Size+IndexXZ(1,iz)];
        hAZm[iLayer*Size+IndexXZ(0,iz)]=InitData.Z0RefCoef*hAZp[iLayer*Size+IndexXZ(1,iz)];
        hAXp[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=InitData.Z1RefCoef*hAXm[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
        hAXm[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=InitData.Z1RefCoef*hAXp[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
        hAZp[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=InitData.Z1RefCoef*hAZm[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
        hAZm[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=InitData.Z1RefCoef*hAZp[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
/*
        if(!InitData.fNoMedia) {
          hM[iLayer*Size+IndexXZ(0,iz)].PXp=InitData.Z0RefCoef*hM[iLayer*Size+IndexXZ(1,iz)].PXm;
          hM[iLayer*Size+IndexXZ(0,iz)].PXm=InitData.Z0RefCoef*hM[iLayer*Size+IndexXZ(1,iz)].PXp;
          hM[iLayer*Size+IndexXZ(0,iz)].PZp=InitData.Z0RefCoef*hM[iLayer*Size+IndexXZ(1,iz)].PZm;
          hM[iLayer*Size+IndexXZ(0,iz)].PZm=InitData.Z0RefCoef*hM[iLayer*Size+IndexXZ(1,iz)].PZp;
          hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PXp=InitData.Z1RefCoef*hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PXm;
          hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PXm=InitData.Z1RefCoef*hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PXp;
          hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PZp=InitData.Z1RefCoef*hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PZm;
          hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PZm=InitData.Z1RefCoef*hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PZp;
        } */
      }
    }
  }
  for(int ix=0;ix<=InitData.nStepsX+1;ix++) {
    hAZpMedia[iLayer*SizeMedia+IndexXZ_M(ix,0)]=hAZmMedia[iLayer*SizeMedia+IndexXZ_M(ix,1)];
    hAZmMedia[iLayer*SizeMedia+IndexXZ_M(ix,0)]=0;
    hAZpMedia[iLayer*SizeMedia+IndexXZ_M(ix,InitData.nStepsZB+1)]=0;
    hAZmMedia[iLayer*SizeMedia+IndexXZ_M(ix,InitData.nStepsZB+1)]=hAZm[iLayer*Size+IndexXZ(ix,1)];
    if(!InitData.fNoMedia) {
      hM2[iLayer*SizeMedia+IndexXZ_M(ix,0)].PZp=hM2[iLayer*SizeMedia+IndexXZ_M(ix,1)].PZm;
      hM2[iLayer*SizeMedia+IndexXZ_M(ix,0)].PZm=0;
      hM2[iLayer*SizeMedia+IndexXZ_M(ix,InitData.nStepsZB+1)].PZp=0;
      hM2[iLayer*SizeMedia+IndexXZ_M(ix,InitData.nStepsZB+1)].PZm=0;
    }
  }
}

void cEquOptic2D_X_Mod1::StepRoutine() {
/*
  if(InitData.UseFastP)
    StepRoutineEulerNoMedia();
  else
    StepRoutineRunge2NoMedia();
*/
  StepRoutineRunge2NoMedia();
}

void cEquOptic2D_X_Mod1::DoTasks() {
  hTaskCritSection->Acquire();
  ThreadCounter=InitData.nThreads;
  hTaskEvent->ResetEvent();
  hTaskCritSection->Release();
  for(int i=0;i<InitData.nThreads;i++) {
    hThreads[i]->Resume();
  }
  hEqus[0]->Step();
  if(InitData.nThreads) {
    if (hTaskEvent->WaitFor(100000) != wrSignaled) {
      Application->MessageBoxA("Thread wait time is off", "OMG",MB_OK);
    }
  }
}

void cEquOptic2D_X_Mod1::StepRoutineEulerNoMedia() {
  FillBorders(0);
  DoTasks();
  DoTasks();
}

void cEquOptic2D_X_Mod1::StepRoutineRunge2NoMedia() {
  PurgeNoise();
  FillBorders(0);
  DoTasks();
  DoTasks();
  DoTasks();
  //CalcOutEnergy();
  //CalcEnergy();
  //hEnergyData[nStepNumber].E=EnergyXP+EnergyXM+EnergyZP+EnergyZM;
  //SaveBitmaps();
  for(int i=0;i<FFTSizeFull-1;i++) hFFTSrc[i]=hFFTSrc[i+1];
  hFFTSrc[FFTSizeFull-1]=hAZp[InitData.SizeZ*InitData.SizeX/2+InitData.SizeX/2];
//  hFFTSrc[FFTSizeFull-1]=hMedia[(InitData.SizeZ-1)*InitData.SizeX+InitData.SizeX/2].AZm;
}

void cEquOptic2D_X_Mod1::CalcEnergy() {
  int iI;
  EnergyXP=EnergyXM=EnergyZP=EnergyZM=EnergyR=EnergyQ=0;
  for(int ix=1;ix<=InitData.nStepsX;ix++)
    for(int iz=1;iz<=InitData.nStepsZ;iz++) {
      iI=IndexXZ(ix,iz);
      EnergyXP+=hAXp[iI].real()*hAXp[iI].real()+hAXp[iI].imag()*hAXp[iI].imag();
      EnergyXM+=hAXm[iI].real()*hAXm[iI].real()+hAXm[iI].imag()*hAXm[iI].imag();
      EnergyZP+=hAZp[iI].real()*hAZp[iI].real()+hAZp[iI].imag()*hAZp[iI].imag();
      EnergyZM+=hAZm[iI].real()*hAZm[iI].real()+hAZm[iI].imag()*hAZm[iI].imag();
//      EnergyR+=hM2[iI].R0;
//      EnergyQ+=hM2[iI].R0Generated;
    }
  for(int ix=1;ix<=InitData.nStepsX;ix++)
    for(int iz=1;iz<=InitData.nStepsZB;iz++) {
      iI=IndexXZ_M(ix,iz);
      EnergyZP+=hAZpMedia[iI].real()*hAZpMedia[iI].real()+
                hAZpMedia[iI].imag()*hAZpMedia[iI].imag();
      EnergyZM+=hAZmMedia[iI].real()*hAZmMedia[iI].real()+
                hAZmMedia[iI].imag()*hAZmMedia[iI].imag();
      EnergyR+=hM2[iI].R0;
      EnergyQ+=hM2[iI].R0Generated;
    }
  EnergyXP*=0.5*InitData.dt*InitData.dt;
  EnergyXM*=0.5*InitData.dt*InitData.dt;
  EnergyZP*=0.5*InitData.dt*InitData.dt;
  EnergyZM*=0.5*InitData.dt*InitData.dt;
  EnergyR*=InitData.dt*InitData.dt;
  EnergyQ*=InitData.dt*InitData.dt;
}

void cEquOptic2D_X_Mod1::CalcOutEnergy() {
  int iI;
  OPTIC_TYPE tmpXP=0,tmpXM=0,tmpZP=0,tmpZM=0,tmpGenZ=0,tmpGenX=0;
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
  tmpXP*=0.5*InitData.dt*InitData.dt;  EnergyXPOut+=tmpXP;
  tmpXM*=0.5*InitData.dt*InitData.dt;  EnergyXMOut+=tmpXM;
  tmpZP*=0.5*InitData.dt*InitData.dt;  EnergyZPOut+=tmpZP;
  tmpZM*=0.5*InitData.dt*InitData.dt;  EnergyZMOut+=tmpZM;
  tmpGenZ*=-InitData.dt*InitData.dt*InitData.dt*InitData.AGenZ;
  tmpGenX*=-InitData.dt*InitData.dt*InitData.dt*InitData.AGenX;
  EnergyUPOut+=tmpGenZ+tmpGenX;
  hEnergyData[nStepNumber].OutXp=tmpXP;
  hEnergyData[nStepNumber].OutXm=tmpXM;
  hEnergyData[nStepNumber].OutZp=tmpZP;
  hEnergyData[nStepNumber].OutZm=tmpZM;
  hEnergyData[nStepNumber].OutGenZ=tmpGenZ;
  hEnergyData[nStepNumber].Out=tmpXP+tmpXM+tmpZP+tmpZM;
  if( (hEnergyData[nStepNumber].Out>InitData.MaxV[0])&&
      (LastMaxOut>hEnergyData[nStepNumber].Out)&&
      (!MaxESaved))
              SaveADataMaxE();
  LastMaxOut=hEnergyData[nStepNumber].Out;
}

/*
void cEquOptic2D_X_Mod1::SaveADataMaxE() {
  ofstream FileOutData = ofstream("DataMaxE.txt",ios::out);
  FileOutData << "t = " << t << endl;
  FileOutData << "SizeZ = " << InitData.dt*InitData.nStepsZ << endl;
  FileOutData << "SizeX = " << InitData.dt*InitData.nStepsX << endl;
  FileOutData << "dt = " << InitData.nStepsZ << endl;
  FileOutData << "a = " << InitData.ACoef << endl;
  FileOutData << "beta = " << InitData.beta << endl;
  FileOutData << "Q = " << InitData.Q << endl;
  FileOutData << "PCoef = " << InitData.PCoef << endl;
  FileOutData << "RCoef = " << InitData.RCoef << endl;
  ofstream FileOutAZp = ofstream("AZpMaxE.dat",ios::out);
  ofstream FileOutAXp = ofstream("AXpMaxE.dat",ios::out);
  ofstream FileOutAZm = ofstream("AZmMaxE.dat",ios::out);
  ofstream FileOutAXm = ofstream("AXmMaxE.dat",ios::out);

  ofstream FileOutAZpF = ofstream("AZpFMaxE.dat",ios::out);
  ofstream FileOutAXpF = ofstream("AXpFMaxE.dat",ios::out);
  ofstream FileOutAZmF = ofstream("AZmFMaxE.dat",ios::out);
  ofstream FileOutAXmF = ofstream("AXmFMaxE.dat",ios::out);
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
    }
  for(int iz=1;iz<InitData.nStepsZ;iz++) {
    FileOutAXpF << InitData.dt*iz << " " << arg(hAXp[IndexXZ(InitData.nStepsX,iz)]) << endl;
    FileOutAXmF << InitData.dt*iz << " " <<arg(hAXm[IndexXZ(1,iz)]) << endl;
  }
  for(int ix=1;ix<InitData.nStepsZ;ix++) {
    FileOutAZpF << InitData.dt*ix << " " <<arg(hAZp[IndexXZ(ix,InitData.nStepsZ)]) << endl;
    FileOutAZmF << InitData.dt*ix << " " <<arg(hAZm[IndexXZ(ix,1)]) << endl;
  }
  MaxESaved=1;
}

void cEquOptic2D_X_Mod1::SaveBitmaps() {
  int tmpInt;
  if(!InitData.SaveBitmapInterval) return;
  if( (LastBitmapSaveTime!=-1)&&
    (t-LastBitmapSaveTime < InitData.SaveBitmapInterval) ) return;

  AnsiString FileName=AnsiString(nBitmap);
  while(FileName.Length()<8) FileName="0"+FileName;
  FileName+=".bmp";

  Graphics::TBitmap *hBitmap = new Graphics::TBitmap;
  hBitmap->Width=InitData.nStepsZ; hBitmap->Height=InitData.nStepsX;

  OPTIC_TYPE maxAXp=0.0000000000000001,tmpV;
  OPTIC_TYPE maxAZp=maxAXp,maxAZm=maxAXp,maxAXm=maxAXp,maxR0=maxAXp;
  for(int ix=0;ix<InitData.nStepsX;ix++) {
    for(int iz=0;iz<InitData.nStepsZ;iz++) {
      tmpV=hAXp[IndexXZ(ix+1,iz+1)].real()*hAXp[IndexXZ(ix+1,iz+1)].real()+
           hAXp[IndexXZ(ix+1,iz+1)].imag()*hAXp[IndexXZ(ix+1,iz+1)].imag();
      maxAXp = maxAXp < tmpV ? tmpV : maxAXp;
      tmpV=hAXm[IndexXZ(ix+1,iz+1)].real()*hAXm[IndexXZ(ix+1,iz+1)].real()+
           hAXm[IndexXZ(ix+1,iz+1)].imag()*hAXm[IndexXZ(ix+1,iz+1)].imag();
      maxAXm = maxAXm < tmpV ? tmpV : maxAXm;
      tmpV=hAZp[IndexXZ(ix+1,iz+1)].real()*hAZp[IndexXZ(ix+1,iz+1)].real()+
           hAZp[IndexXZ(ix+1,iz+1)].imag()*hAZp[IndexXZ(ix+1,iz+1)].imag();
      maxAZp = maxAZp < tmpV ? tmpV : maxAZp;
      tmpV=hAZm[IndexXZ(ix+1,iz+1)].real()*hAZm[IndexXZ(ix+1,iz+1)].real()+
           hAZm[IndexXZ(ix+1,iz+1)].imag()*hAZm[IndexXZ(ix+1,iz+1)].imag();
      maxAZm = maxAZm < tmpV ? tmpV : maxAZm;
      tmpV=hM[IndexXZ(ix+1,iz+1)].R0;
      maxR0 = maxR0 < tmpV ? tmpV : maxR0;
    }
  }
  maxAXp=sqrt(maxAXp); maxAXm=sqrt(maxAXm);
  maxAZp=sqrt(maxAZp); maxAZm=sqrt(maxAZm);

  maxmaxAXp = maxmaxAXp < maxAXp ? maxAXp : maxmaxAXp;
  maxmaxAXm = maxmaxAXm < maxAXm ? maxAXm : maxmaxAXm;
  maxmaxAZp = maxmaxAZp < maxAZp ? maxAZp : maxmaxAZp;
  maxmaxAZm = maxmaxAZm < maxAZm ? maxAZm : maxmaxAZm;

  for(int ix=0;ix<InitData.nStepsX;ix++) {
    for(int iz=0;iz<InitData.nStepsZ;iz++) {
      hBitmap->Canvas->Pixels[iz][ix]=ColorData[int(255*abs(hAXp[IndexXZ(ix+1,iz+1)])/maxAXp)];
    }
  }
  hBitmap->SaveToFile("Bitmaps\\AXp\\"+FileName);

  for(int ix=0;ix<InitData.nStepsX;ix++) {
    for(int iz=0;iz<InitData.nStepsZ;iz++) {
      hBitmap->Canvas->Pixels[iz][ix]=ColorData[int(255*abs(hAXm[IndexXZ(ix+1,iz+1)])/maxAXm)];
    }
  }
  hBitmap->SaveToFile("Bitmaps\\AXm\\"+FileName);

  for(int ix=0;ix<InitData.nStepsX;ix++) {
    for(int iz=0;iz<InitData.nStepsZ;iz++) {
      hBitmap->Canvas->Pixels[iz][ix]=ColorData[int(255*abs(hAZp[IndexXZ(ix+1,iz+1)])/maxAZp)];
    }
  }
  hBitmap->SaveToFile("Bitmaps\\AZp\\"+FileName);

  for(int ix=0;ix<InitData.nStepsX;ix++) {
    for(int iz=0;iz<InitData.nStepsZ;iz++) {
      hBitmap->Canvas->Pixels[iz][ix]=ColorData[int(255*abs(hAZm[IndexXZ(ix+1,iz+1)])/maxAZm)];
    }
  }
  hBitmap->SaveToFile("Bitmaps\\AZm\\"+FileName);

  for(int ix=0;ix<InitData.nStepsX;ix++) {
    for(int iz=0;iz<InitData.nStepsZ;iz++) {
      tmpInt=int(255*abs(hM[IndexXZ(ix+1,iz+1)].R0)/maxR0);
      tmpInt = tmpInt > 255 ? 255 : tmpInt;
      tmpInt = tmpInt < 0 ? 0 : tmpInt;
      hBitmap->Canvas->Pixels[iz][ix]=ColorData[int(255*abs(hM[IndexXZ(ix+1,iz+1)].R0)/maxR0)];    }
  }
  hBitmap->SaveToFile("Bitmaps\\R0\\"+FileName);

  if(MakeGlobalNormBitmap) {
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      for(int iz=0;iz<InitData.nStepsZ;iz++) {
        tmpInt=int(255*abs(hAXp[IndexXZ(ix+1,iz+1)])/normAXp);
        tmpInt = tmpInt > 255 ? 255 : tmpInt;
        hBitmap->Canvas->Pixels[iz][ix]=ColorData[tmpInt];
      }
    }
    hBitmap->SaveToFile("Bitmaps\\AXpNorm\\"+FileName);

    for(int ix=0;ix<InitData.nStepsX;ix++) {
      for(int iz=0;iz<InitData.nStepsZ;iz++) {
        tmpInt=int(255*abs(hAXm[IndexXZ(ix+1,iz+1)])/normAXm);
        tmpInt = tmpInt > 255 ? 255 : tmpInt;
        hBitmap->Canvas->Pixels[iz][ix]=ColorData[tmpInt];
      }
    }
    hBitmap->SaveToFile("Bitmaps\\AXmNorm\\"+FileName);

    for(int ix=0;ix<InitData.nStepsX;ix++) {
      for(int iz=0;iz<InitData.nStepsZ;iz++) {
        tmpInt=int(255*abs(hAZp[IndexXZ(ix+1,iz+1)])/normAZp);
        tmpInt = tmpInt > 255 ? 255 : tmpInt;
        hBitmap->Canvas->Pixels[iz][ix]=ColorData[tmpInt];
      }
    }
    hBitmap->SaveToFile("Bitmaps\\AZpNorm\\"+FileName);

    for(int ix=0;ix<InitData.nStepsX;ix++) {
      for(int iz=0;iz<InitData.nStepsZ;iz++) {
        tmpInt=int(255*abs(hAZm[IndexXZ(ix+1,iz+1)])/normAZm);
        tmpInt = tmpInt > 255 ? 255 : tmpInt;
        hBitmap->Canvas->Pixels[iz][ix]=ColorData[tmpInt];
      }
    }
    hBitmap->SaveToFile("Bitmaps\\AZmNorm\\"+FileName);
  }

  int ComplexBaseLine=InitData.nStepsX*3+19;
  if (MakeComplexBitmap) {
    hBitmap->Width=InitData.nStepsZ; hBitmap->Height=InitData.nStepsX*3+20;
    */
    /*
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      for(int iz=0;iz<InitData.nStepsZ;iz++) {
        hBitmap->Canvas->Pixels[iz][ix]=ColorData[int(255*abs(hAXp[IndexXZ(ix+1,iz+1)])/maxAXp)];
      }
    }
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      for(int iz=0;iz<InitData.nStepsZ;iz++) {
        tmpInt=int(255*abs(hAXp[IndexXZ(ix+1,iz+1)])/normAXp);
        tmpInt = tmpInt > 255 ? 255 : tmpInt;
        hBitmap->Canvas->Pixels[iz][InitData.nStepsX+10+ix]=ColorData[tmpInt];
      }
    }
    hBitmap->Canvas->Pen->Color=clBlack;
    hBitmap->Canvas->MoveTo(0,ComplexBaseLine-InitData.nStepsX);
    hBitmap->Canvas->LineTo(0,ComplexBaseLine);
    hBitmap->Canvas->LineTo(InitData.nStepsZ,ComplexBaseLine);
    int nSt=InitData.tMax/InitData.dt;
    for(int i=0;i<nStepNumber;i++) {
      tmpInt=(hEnergyData[i].OutXp+hEnergyData[i].OutXm)*0.5*InitData.nStepsX/ComplexBitmapNorm;
      hBitmap->Canvas->Pixels[i*InitData.nStepsZ/nSt][ComplexBaseLine-tmpInt]=clGreen;
    }
    */
    /*
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      for(int iz=0;iz<InitData.nStepsZ;iz++) {
        hBitmap->Canvas->Pixels[iz][ix]=ColorData[int(255*abs(hAZp[IndexXZ(ix+1,iz+1)])/maxAZp)];
      }
    }
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      for(int iz=0;iz<InitData.nStepsZ;iz++) {
        tmpInt=int(255*abs(hAZp[IndexXZ(ix+1,iz+1)])/normAZp);
        tmpInt = tmpInt > 255 ? 255 : tmpInt;
        hBitmap->Canvas->Pixels[iz][InitData.nStepsX+10+ix]=ColorData[tmpInt];
      }
    }
    hBitmap->Canvas->Pen->Color=clBlack;
    hBitmap->Canvas->MoveTo(0,ComplexBaseLine-InitData.nStepsX);
    hBitmap->Canvas->LineTo(0,ComplexBaseLine);
    hBitmap->Canvas->LineTo(InitData.nStepsZ,ComplexBaseLine);
    int nSt=InitData.tMax/InitData.dt;
    for(int i=0;i<nStepNumber;i++) {
      tmpInt=(hEnergyData[i].OutZp+hEnergyData[i].OutZm)*0.5*InitData.nStepsX/ComplexBitmapNorm;
      hBitmap->Canvas->Pixels[i*InitData.nStepsZ/nSt][ComplexBaseLine-tmpInt]=clGreen;
    }
    hBitmap->Canvas->Pen->Color=clRed;
    hBitmap->Canvas->MoveTo((nStepNumber-1)*InitData.nStepsZ/nSt,ComplexBaseLine-InitData.nStepsX);
    hBitmap->Canvas->LineTo((nStepNumber-1)*InitData.nStepsZ/nSt,ComplexBaseLine);
    hBitmap->SaveToFile("Bitmaps\\Complex\\"+FileName);
  }

  LastBitmapSaveTime=t;  nBitmap++;
  delete hBitmap;
}     */

/*
void cEquOptic2D_X_Mod1::PrepareBitmapFolders() {
  TSearchRec SearchRec;
  int findV;
  CreateDir("Bitmaps");
  CreateDir("Bitmaps\\AXp");
  CreateDir("Bitmaps\\AXm");
  CreateDir("Bitmaps\\AZp");
  CreateDir("Bitmaps\\AZm");
  CreateDir("Bitmaps\\R0");
  CreateDir("Bitmaps\\AXpNorm");
  CreateDir("Bitmaps\\AXmNorm");
  CreateDir("Bitmaps\\AZpNorm");
  CreateDir("Bitmaps\\AZmNorm");
  CreateDir("Bitmaps\\Complex");
  findV = FindFirst("Bitmaps\\AXp\\*.*",0,SearchRec);
  while(!findV) {
    DeleteFile("Bitmaps\\AXp\\"+SearchRec.Name);
    findV=FindNext(SearchRec);
  }
  FindClose(SearchRec);
  findV = FindFirst("Bitmaps\\AXm\\*.*",0,SearchRec);
  while(!findV) {
    DeleteFile("Bitmaps\\AXm\\"+SearchRec.Name);
    findV=FindNext(SearchRec);
  }
  FindClose(SearchRec);
  findV = FindFirst("Bitmaps\\AZp\\*.*",0,SearchRec);
  while(!findV) {
    DeleteFile("Bitmaps\\AZp\\"+SearchRec.Name);
    findV=FindNext(SearchRec);
  }
  FindClose(SearchRec);
  findV = FindFirst("Bitmaps\\AZm\\*.*",0,SearchRec);
  while(!findV) {
    DeleteFile("Bitmaps\\AZm\\"+SearchRec.Name);
    findV=FindNext(SearchRec);
  }
  FindClose(SearchRec);
  findV = FindFirst("Bitmaps\\R0\\*.*",0,SearchRec);
  while(!findV) {
    DeleteFile("Bitmaps\\R0\\"+SearchRec.Name);
    findV=FindNext(SearchRec);
  }
  FindClose(SearchRec);
  findV = FindFirst("Bitmaps\\AXpNorm\\*.*",0,SearchRec);
  while(!findV) {
    DeleteFile("Bitmaps\\AXpNorm\\"+SearchRec.Name);
    findV=FindNext(SearchRec);
  }
  FindClose(SearchRec);
  findV = FindFirst("Bitmaps\\AXmNorm\\*.*",0,SearchRec);
  while(!findV) {
    DeleteFile("Bitmaps\\AXmNorm\\"+SearchRec.Name);
    findV=FindNext(SearchRec);
  }
  FindClose(SearchRec);
  findV = FindFirst("Bitmaps\\AZpNorm\\*.*",0,SearchRec);
  while(!findV) {
    DeleteFile("Bitmaps\\AZpNorm\\"+SearchRec.Name);
    findV=FindNext(SearchRec);
  }
  FindClose(SearchRec);
  findV = FindFirst("Bitmaps\\AZmNorm\\*.*",0,SearchRec);
  while(!findV) {
    DeleteFile("Bitmaps\\AZmNorm\\"+SearchRec.Name);
    findV=FindNext(SearchRec);
  }
  FindClose(SearchRec);
  findV = FindFirst("Bitmaps\\Complex\\*.*",0,SearchRec);
  while(!findV) {
    DeleteFile("Bitmaps\\Complex\\"+SearchRec.Name);
    findV=FindNext(SearchRec);
  }
  FindClose(SearchRec);
}

int cEquOptic2D_X_Mod1::TimeToSaveBitmaps() {
  if(!InitData.SaveBitmapInterval) return 0;
  if( (LastBitmapSaveTime!=-1)&&
    (t-LastBitmapSaveTime < InitData.SaveBitmapInterval) ) return 0;
  return 1;
}

*/

