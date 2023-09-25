//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquOptic2D_1DwB.h"
#include "fstream.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)


void cEquOptic2D_1DwB::InitReport(cDifEquReport *hReport) {
  int Lines[EQU_MAX_GRAPHS];
  for(int i=0;i<EQU_MAX_GRAPHS;i++) Lines[i]=1;
  Lines[0]=1; Lines[1]=1; Lines[2]=1; Lines[3]=1; Lines[4]=1;
  hReport->LoadText(hReportText);
  hReport->SetGraphSize(EQU_GRAPH_SIZE,Lines);
  hReport->SetMaxTo1();
}

void cEquOptic2D_1DwB::Report(cDifEquReport *hReport) {
  nReportNumber++;

  for(int iz=1; iz<=InitData.nStepsZ; iz++)
    hReportBuf[iz-1]=abs(hAZp[iz]);
  hReport->Graph[0].GetData(0,InitData.nStepsZ,hReportBuf);
  for(int iz=1; iz<=InitData.nStepsZ; iz++)
    hReportBuf[iz-1]=abs(hAZm[iz]);
  hReport->Graph[1].GetData(0,InitData.nStepsZ,hReportBuf);
  for(int iz=1; iz<=InitData.nStepsZ; iz++)
    hReportBuf[iz-1]=abs(hB[iz]);
  hReport->Graph[2].GetData(0,InitData.nStepsZ,hReportBuf);
  for(int iz=1; iz<=InitData.nStepsZ; iz++)
    hReportBuf[iz-1]=abs(hM[iz].R0);
  hReport->Graph[3].GetData(0,InitData.nStepsZ,hReportBuf);

  hReport->nVals=1;
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
  hReport->ValText[hReport->nVals]="Max B ";
  hReport->Val[hReport->nVals++]=hReport->Graph[2].Max;
  CalcEnergy();
  OPTIC_TYPE EnergyA=EnergyZP+EnergyZM+EnergyB;
  double gen = (EnergyGenX+EnergyGenZ)*InitData.dt;
  double ea = EnergyZP+EnergyZM+EnergyB;
  double out = EnergyZPOut+EnergyZMOut;
  double deltaEA=ea-EnergyA0;
  double deltaER=EnergyR-EnergyR0;
  double eq = EnergyQ-EnergyQBack;
  double balance=deltaEA-gen+out-eq+deltaER;
  hReport->ValText[hReport->nVals]="EA ";
  hReport->Val[hReport->nVals++]=EnergyZP+EnergyZM;
  hReport->ValText[hReport->nVals]="EB ";
  hReport->Val[hReport->nVals++]=EnergyB;
  hReport->ValText[hReport->nVals]="ER ";
  hReport->Val[hReport->nVals++]=EnergyR;
  hReport->ValText[hReport->nVals]="delta EA ";
  hReport->Val[hReport->nVals++]=EnergyA-EnergyA0;
  hReport->ValText[hReport->nVals]="A balance ";
  hReport->Val[hReport->nVals++]=deltaEA-gen+out;
  hReport->ValText[hReport->nVals]="E balance ";
  hReport->Val[hReport->nVals++]=balance;
  hReport->ValText[hReport->nVals]="R0 balance ";
  hReport->Val[hReport->nVals++]=deltaER-eq;
//  hReport->ValText[hReport->nVals]="P Wrong ";
//  hReport->Val[hReport->nVals++]=PGenR+PGenA-POutA;
//  hReport->ValText[hReport->nVals]="Tmp Wrong ";
//  hReport->Val[hReport->nVals++]=EnergyA+EnergyR-EnergyA0-EnergyQ-EnergyR0;
  
  hReport->ValText[hReport->nVals]="Inc ZP ";
  OPTIC_TYPE tmpInc=0;
  if(EnergyZPlast > 0.000000001)
    tmpInc=2.0*log(EnergyZP/EnergyZPlast)/InitData.dt;
  hReport->Val[hReport->nVals++]=tmpInc;


  /*

  hReport->nVals=3;

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

//  hReport->Graph[10].GetData(0,nReportNumber,hEnergyData);
//  hReport->SetMaxTo1();
  /*
  hReport->nVals=3;
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
  hReport->ValText[1]="Max I "; hReport->Val[1]=hReport->Graph[1].Max;
  hReport->ValText[2]="Max Gamma "; hReport->Val[2]=hReport->Graph[2].Max;
  hReport->Graph[2].Max=2; */
  EnergyZPlast=EnergyZP; ESumLast=EnergyA; LastReportTime=t;
}

void cEquOptic2D_1DwB::FinalReport(cDifEquReport *hReport) {
  hReport->nVals=0; hReport->nGraphs=2;
  for(int i=0;i<nStepNumber;i++) hReportBuf[i]=hEnergyData[i].Out;
  hReport->Graph[0].GetData(0,nStepNumber,hReportBuf);
  for(int i=0;i<nStepNumber;i++) hReportBuf[i]=hEnergyData[i].R;
  hReport->Graph[1].GetData(0,nStepNumber,hReportBuf);
/*  hReport->Graph[0].GetDataAbs(0,nStepNumber,hOutA);
  hReport->GraphText[0]="Out A";
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;*/

  hReport->nVals=1;
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
  hReport->ValText[hReport->nVals]="Max B ";
  hReport->Val[hReport->nVals++]=hReport->Graph[2].Max;
  CalcEnergy();
  OPTIC_TYPE EnergyA=EnergyZP+EnergyZM+EnergyB;
  double gen = (EnergyGenX+EnergyGenZ)*InitData.dt;
  double ea = EnergyZP+EnergyZM+EnergyB;
  double out = EnergyZPOut+EnergyZMOut;
  double deltaEA=ea-EnergyA0;
  double deltaER=EnergyR-EnergyR0;
  double eq = EnergyQ-EnergyQBack;
  double balance=deltaEA-gen+out-eq+deltaER;
  hReport->ValText[hReport->nVals]="EA ";
  hReport->Val[hReport->nVals++]=EnergyZP+EnergyZM;
  hReport->ValText[hReport->nVals]="EB ";
  hReport->Val[hReport->nVals++]=EnergyB;
  hReport->ValText[hReport->nVals]="ER ";
  hReport->Val[hReport->nVals++]=EnergyR;
  hReport->ValText[hReport->nVals]="delta EA ";
  hReport->Val[hReport->nVals++]=EnergyA-EnergyA0; 
  hReport->ValText[hReport->nVals]="A balance ";
  hReport->Val[hReport->nVals++]=deltaEA-gen+out;
  hReport->ValText[hReport->nVals]="E balance ";
  hReport->Val[hReport->nVals++]=balance;
  hReport->ValText[hReport->nVals]="R0 balance ";
  hReport->Val[hReport->nVals++]=deltaER-eq;
  hReport->ValText[hReport->nVals]="P Out B ";
  hReport->Val[hReport->nVals++]=PGenA;
  hReport->ValText[hReport->nVals]="P Out A ";
  hReport->Val[hReport->nVals++]=POutA;
  hReport->ValText[hReport->nVals]="P Rel ";
  hReport->Val[hReport->nVals++]=PRel;

  FinalSave();
}

void cEquOptic2D_1DwB::FinalSave() {
  long double tmp1;

  ofstream FileOutAZp = ofstream("AZp.dat",ios::out);
  ofstream FileOutAZm = ofstream("AZm.dat",ios::out);
  ofstream FileOutB = ofstream("B.dat",ios::out);
  ofstream FileOutR0 = ofstream("R0.dat",ios::out);
  ofstream FileOutRt = ofstream("Rt.dat",ios::out);
  ofstream FileOutGenXt = ofstream("GenXt.dat",ios::out);
  ofstream FileOutRelt = ofstream("PRelt.dat",ios::out);
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    FileOutAZp << iz*InitData.dt << "  " << abs(hAZp[iz]) << endl;
    FileOutAZm << iz*InitData.dt << "  " << abs(hAZm[iz]) << endl;
    FileOutB << iz*InitData.dt << "  " << abs(hB[iz]) << endl;
    FileOutR0 << iz*InitData.dt << "  " << hM[iz].R0 << endl;
  }
  for(int i=0;i<nStepNumber;i++) {
    FileOutRt << i*InitData.dt << " " << hEnergyData[i].R << endl;
    FileOutGenXt << i*InitData.dt << " " << hEnergyData[i].OutGenX << endl;
    FileOutRelt << i*InitData.dt << " "  << hEnergyData[i].PRel << endl;
  }
}

void cEquOptic2D_1DwB::LoadInitData(cDifEqu_InitDataRec *hInitData) {
  cDifEqu_InitDataRec_Optic2D_T<OPTIC_TYPE> *hInit=
             (cDifEqu_InitDataRec_Optic2D_T<OPTIC_TYPE> *)hInitData;
  int nStepsT=hInit->tMax/hInit->dt+1;
  OPTIC_TYPE MaxSize = hInit->nStepsZ+2 > nStepsT ? hInit->nStepsZ+2 : nStepsT;
  delete [] hEnergyData;
  hEnergyData = new EnergyDataType[hInit->tMax/hInit->dt+1];
  delete [] hAZpOutData;
  hAZpOutData = new complex<OPTIC_TYPE>[nStepsT];
  if(InitData.nStepsZ!=hInit->nStepsZ) {
    delete [] hReportBuf;
    delete [] hB; delete [] hAZp; delete [] hAZm;
    delete [] hM; delete [] hAbsStore;
    delete [] hAZpD; delete [] hAZmD; delete [] hBD; delete [] hMD;
    hB = new  complex<OPTIC_TYPE>[(hInit->nStepsZ+2)*3];
    hAZp = new complex<OPTIC_TYPE>[(hInit->nStepsZ+2)*3];
    hAZm = new complex<OPTIC_TYPE>[(hInit->nStepsZ+2)*3];
    hM   = new cOptic2D_Media_2Way_B[(hInit->nStepsZ+2)*3];

    hBD = new  complex<OPTIC_TYPE>[(hInit->nStepsZ+2)*4];
    hAZpD = new complex<OPTIC_TYPE>[(hInit->nStepsZ+2)*4];
    hAZmD = new complex<OPTIC_TYPE>[(hInit->nStepsZ+2)*4];
    hMD   = new cOptic2D_Media_2Way_B[(hInit->nStepsZ+2)*4];
    hReportBuf = new OPTIC_TYPE[MaxSize];
  }

  InitData=*hInit;
  InitData.SizeZ=hInit->nStepsZ+2;
  t=0; dt=InitData.dt; tMax=InitData.tMax;

  delete [] hFFTSrc; delete [] hFFTDest; FFTSizeFull = 1 << hInit->FFTSize;
  hFFTSrc = new complex<OPTIC_TYPE>[FFTSizeFull];
  hFFTDest = new complex<OPTIC_TYPE>[FFTSizeFull];
  for (int i=0;i<FFTSizeFull;i++) hFFTSrc[i]=hFFTDest[i]=0;
  FFT.MakeRevTable(InitData.FFTSize);

  PrepareMedia();
  FillBorders();
  EnergyZPOut=EnergyZMOut=EnergyQ=0;
  CalcEnergy(); EnergyA0=EnergyZP+EnergyZM+EnergyB; EnergyR0=EnergyR;
  ESumLast=EnergyZPlast=0; LastReportTime=0; LastMaxOut=0; EnergyUPOut=0;
  EnergyGenZ=EnergyGenX=0; EnergyQBack=0;
  nStepNumber=0; t=0; MaxESaved=0; LastBitmapSaveTime=-1; nBitmap=0;
  maxmaxAZp=maxmaxAZm=0;
  normAZp=normAZm=100;
  MakeGlobalNormBitmap=MakeComplexBitmap=0;
  PurgeNoiseTime=0;
}

void cEquOptic2D_1DwB::PrepareMedia() {
  OPTIC_TYPE r1,r2,r3,r4;
  for(int iz=0;iz<InitData.nStepsZ+2;iz++) {
    hAZp[iz]=0; hAZm[iz]=0; hB[iz]=0;
    hM[iz].PZp=0;
    hM[iz].PZm=0;
    hM[iz].PB=0;
    hM[iz].RZ=0;
    hM[iz].R0=InitData.StartingInv;
    hM[iz].R0Generated=0;
    hAZp[iz]+=InitData.A0*sin(M_PI*iz/2/InitData.nStepsZ)+
              InitData.A1*sin(M_PI*iz/InitData.nStepsZ)+
              InitData.A2*sin(4*M_PI*iz/InitData.nStepsZ);
    hAZm[iz]+=InitData.A0*sin(M_PI*iz/InitData.nStepsZ)+
              InitData.A1*sin(2*M_PI*iz/InitData.nStepsZ)+
              InitData.A2*sin(3*M_PI*iz/InitData.nStepsZ);
    if(InitData.rPolNoise) {
      r1=2.0*M_PI*rand()/RAND_MAX;
      r2=2.0*M_PI*rand()/RAND_MAX;
      hAZp[iz]+=polar(InitData.rPolNoise,r1);
      hAZm[iz]+=polar(InitData.rPolNoise,r2);
    }
  }
}

void cEquOptic2D_1DwB::FillBorders() {
  int sz=InitData.SizeZ;
  hAZp[0]=0;
  hAZm[InitData.nStepsZ+1]=0;
  hAZp[sz]=0;
  hAZm[InitData.nStepsZ+1+sz]=0;
  hM[0].PZp=hM[0].PZm=0;
  hM[InitData.nStepsZ+1].PZp=hM[InitData.nStepsZ+1].PZm=0;
}

void cEquOptic2D_1DwB::StepRoutine() {

  if (InitData.fNoMedia)
    StepRoutineRunge2Matrix();
  else
    StepRoutineRunge2();

//  StepRoutineRunge2();
//  StepRoutineRunge2Matrix();
}

void cEquOptic2D_1DwB::StepDebug() {
/*
  int SizeM=InitData.SizeZ*sizeof(cOptic2D_Media_2Way_B);
  int SizeA=InitData.SizeZ*sizeof(complex<OPTIC_TYPE>);
  int sz=InitData.SizeZ;
  memcpy(hMD,hM,SizeM);
  memcpy(hAZpD,hAZp,SizeA);
  memcpy(hAZmD,hAZm,SizeA);
  memcpy(hBD,hBp,SizeA);
  StepRoutineEuler();

  memcpy(hMD+sz,hM,SizeM);
  memcpy(hAZpD+sz,hAZp,SizeA);
  memcpy(hAZmD+sz,hAZm,SizeA);
  memcpy(hBD+sz,hBp,SizeA);

  memcpy(hMD,hMD,SizeM);
  memcpy(hAZp,hAZpD,SizeA);
  memcpy(hAZm,hAZmD,SizeA);
  memcpy(hB,hBpD,SizeA);
  StepRoutineRunge2();

  memcpy(hMD+2*sz,hM,SizeM);
  memcpy(hAZpD+2*sz,hAZp,SizeA);
  memcpy(hAZmD+2*sz,hAZm,SizeA);
  memcpy(hBD+2*sz,hBp,SizeA);

  memcpy(hMD,hMD,SizeM);
  memcpy(hAZp,hAZpD,SizeA);
  memcpy(hAZm,hAZmD,SizeA);
  memcpy(hB,hBpD,SizeA);
  StepRoutineRunge2Matrix();

  memcpy(hMD+3*sz,hM,SizeM);
  memcpy(hAZpD+3*sz,hAZp,SizeA);
  memcpy(hAZmD+3*sz,hAZm,SizeA);
  memcpy(hBD+3*sz,hBp,SizeA);

  OPTIC_TYPE EZP0=0, EZP1=0, EZP2=0, EZP3=0, EZM0=0, EZM1=0, EZM2=0, EZM3=0,
  OPTIC_TYPE EB0=0, EB1=0, EB2=0, EB3=0;
  PGenR=PGenA=0;
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    EB+=hB[iz].real()*hB[iz].real()+hB[iz].imag()*hB[iz].imag();
    EZP+=hAZp[iz].real()*hAZp[iz].real()+hAZp[iz].imag()*hAZp[iz].imag();
    EZM+=hAZm[iz].real()*hAZm[iz].real()+hAZm[iz].imag()*hAZm[iz].imag();
    EnergyR+=hM[iz].R0;
    EnergyQ+=hM[iz].R0Generated;
    PGenA+=InitData.AGenX*(hB[iz].real()*hB[iz].real()+hB[iz].imag()*hB[iz].imag());
    PGenA+=InitData.AGenZ*(hAZp[iz].real()*hAZp[iz].real()+hAZp[iz].imag()*hAZp[iz].imag());
    PGenA+=InitData.AGenZ*(hAZm[iz].real()*hAZm[iz].real()+hAZm[iz].imag()*hAZm[iz].imag());
    PGenR+=InitData.Q-hM[iz].R0*InitData.RCoef;
  }
  PGenA*=InitData.dt;
  PGenR*=InitData.dt;
  POutA=hAZp[InitData.nStepsZ].real()*hAZp[InitData.nStepsZ].real()+
                hAZp[InitData.nStepsZ].imag()*hAZp[InitData.nStepsZ].imag();
  POutA+=hAZm[1].real()*hAZm[1].real()+hAZm[1].imag()*hAZm[1].imag();
  POutA*=0.5;
  PRel=-PGenR+InitData.Q*InitData.dt*InitData.nStepsZ;

  EnergyB*=0.5*InitData.dt;
  EnergyZP*=0.5*InitData.dt;
  EnergyZM*=0.5*InitData.dt;
  EnergyR*=InitData.dt;
  EnergyQ*=InitData.dt;
  */
}

void cEquOptic2D_1DwB::StepRoutineEuler() {
  complex<OPTIC_TYPE> tmpPZp,tmpPZm,tmpPB;
  OPTIC_TYPE dt2=InitData.dt;
  complex<OPTIC_TYPE> ii=complex<OPTIC_TYPE>(0,1),tmpAZp,tmpAZm,tmpB;
  cOptic2D_Media_2Way_B tmpM;
  OPTIC_TYPE two=2;
  int sz=InitData.SizeZ;
                        // InitData.AGenCoef - коэффициент связи среды с волной B !
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    if(InitData.UseFastP) {
      hM[iz].PZp=tmpPZp=(hAZp[iz]*hM[iz].R0*two+hAZm[iz]*hM[iz].RZ)*InitData.beta;
      hM[iz].PZm=tmpPZm=(hAZm[iz]*hM[iz].R0*two+hAZp[iz]*conj(hM[iz].RZ))*InitData.beta;
      hM[iz].PB =tmpPB=(hB[iz]*hM[iz].R0*two*InitData.AGenCoef)*InitData.beta;
    }
    hM[iz+sz].PZp=(hAZp[iz]*hM[iz].R0*two+hAZm[iz]*hM[iz].RZ)*InitData.beta*dt2;
    hM[iz+sz].PZm=(hAZm[iz]*hM[iz].R0*two+hAZp[iz]*conj(hM[iz].RZ))*InitData.beta*dt2;
    hM[iz+sz].PB=(hB[iz]*hM[iz].R0*two*InitData.AGenCoef)*InitData.beta*dt2;
    hM[iz+sz].R0 =(InitData.Q-InitData.RCoef*hM[iz].R0
                       -(hAZp[iz]*conj(hM[iz].PZp)
                        +hAZm[iz]*conj(hM[iz].PZm)
                        +hB[iz]*conj(hM[iz].PB)
                       ).real() )*dt2;
    hM[iz+sz].RZ=(-hM[iz].RZ*InitData.RCoef
                        -(hAZp[iz]*conj(hM[iz].PZm)
                        +hM[iz].PZm*conj(hAZp[iz])) )*dt2;
    hM[iz+sz].R0Generated=(InitData.Q-InitData.RCoef*hM[iz].R0)*dt2;
    hAZp[iz+sz]=dt2*hM[iz].PZp+dt2*InitData.AGenZ*hAZp[iz]+ii*InitData.ACoef*hB[iz]*dt2;
    hAZm[iz+sz]=dt2*hM[iz].PZm+dt2*InitData.AGenZ*hAZm[iz]+ii*InitData.ACoef*hB[iz]*dt2;
    hB[iz+sz]=(ii*InitData.ACoef*(hAZp[iz]+hAZm[iz])+hM[iz].PB+InitData.AGenX*hB[iz])*dt2;
  }

double dAP7=abs(hAZp[7]+hAZp[7+sz])*abs(hAZp[7]+hAZp[7+sz])-abs(hAZp[7])*abs(hAZp[7]);
double dAM7=abs(hAZm[7]+hAZm[7+sz])*abs(hAZm[7]+hAZm[7+sz])-abs(hAZm[7])*abs(hAZm[7]);
double dR7=hM[7+sz].R0;
double dd=abs(dAP7+dAM7+dR7);
/*
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    hM[iz].M1M2v(hM[iz],hM[iz+sz],1);
    hAZp[iz]=hAZp[iz]+hAZp[iz+sz];
    hAZm[iz]=hAZm[iz]+hAZm[iz+sz];
    hB[iz]=hB[iz]+hB[iz+sz];
  }
*/

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    hM[iz].M1M2v(hM[iz],hM[iz+sz],1);
    hAZp[InitData.nStepsZ+1-iz]=hAZp[InitData.nStepsZ-iz]+hAZp[InitData.nStepsZ-iz+sz];
    hAZm[iz]=hAZm[iz+1]+hAZm[iz+1+sz];
    hB[iz]=hB[iz]+hB[iz+sz];
  }

  FillBorders();
  CalcOutEnergy();
  hAZpOutData[nStepNumber]=hAZp[InitData.nStepsZ-1];
  for(int i=0;i<FFTSizeFull-1;i++) hFFTSrc[i]=hFFTSrc[i+1];
  hFFTSrc[FFTSizeFull-1]=hAZp[InitData.nStepsZ/2];
}




void cEquOptic2D_1DwB::StepRoutineRunge2Matrix() {
  complex<OPTIC_TYPE> tmpPZp,tmpPZm,tmpPB;
  OPTIC_TYPE dt2=0.5*InitData.dt;
  complex<OPTIC_TYPE> ii=complex<OPTIC_TYPE>(0,1),tmpAZp,tmpAZm,tmpB,tAZpM,tAZmM,tBM;
  complex<OPTIC_TYPE> BCoef=complex<OPTIC_TYPE>(InitData.AGenX,InitData.Delta);
  complex<OPTIC_TYPE> q=ii*InitData.ACoef*dt2;
  OPTIC_TYPE two=2;
  complex<OPTIC_TYPE> q1=complex<OPTIC_TYPE>(1,0)/(OPTIC_TYPE(1)-q*q*two);
  cOptic2D_Media_2Way_B tmpM;
  int sz=InitData.SizeZ;

                        // InitData.AGenCoef - коэффициент связи среды с волной B !

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    if(InitData.UseFastP) {
      hM[iz].PZp=tmpPZp=(hAZp[iz]*hM[iz].R0*two+hAZm[iz]*hM[iz].RZ)*InitData.beta;
      hM[iz].PZm=tmpPZm=(hAZm[iz]*hM[iz].R0*two+hAZp[iz]*conj(hM[iz].RZ))*InitData.beta;
      hM[iz].PB =tmpPB=hB[iz]*hM[iz].R0*two*InitData.AGenCoef;
    }
    hM[iz+sz].PZp=(hAZp[iz]*hM[iz].R0*two+hAZm[iz]*hM[iz].RZ)*InitData.beta*dt2;
    hM[iz+sz].PZm=(hAZm[iz]*hM[iz].R0*two+hAZp[iz]*conj(hM[iz].RZ))*InitData.beta*dt2;
    hM[iz+sz].PB=(hB[iz]*hM[iz].R0*two*InitData.AGenCoef)*dt2;
    hM[iz+sz].R0 =(InitData.Q-InitData.RCoef*hM[iz].R0
                       -(hAZp[iz]*conj(hM[iz].PZp)
                        +hAZm[iz]*conj(hM[iz].PZm)
                        +hB[iz]*conj(hM[iz].PB)
                       ).real() )*dt2;
    hM[iz+sz].RZ=(-hM[iz].RZ*InitData.RCoef
                        -(hAZp[iz]*conj(hM[iz].PZm)
                        +hM[iz].PZp*conj(hAZm[iz])) )*dt2;
    hM[iz+sz].R0Generated=(InitData.Q-InitData.RCoef*hM[iz].R0)*dt2;
    hAZp[iz+sz]=dt2*hM[iz].PZp+dt2*InitData.AGenZ*hAZp[iz]+ii*InitData.ACoef*hB[iz]*dt2;
    hAZm[iz+sz]=dt2*hM[iz].PZm+dt2*InitData.AGenZ*hAZm[iz]+ii*InitData.ACoef*hB[iz]*dt2;
    hB[iz+sz]=(ii*InitData.ACoef*(hAZp[iz]+hAZm[iz])+hM[iz].PB+hB[iz]*BCoef)*dt2;
  }

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    tmpAZp = hAZp[iz-1]+hAZp[iz+sz-1];
    tmpAZm = hAZm[iz+1]+hAZm[iz+sz+1];
    tmpB   = hB[iz]+hB[iz+sz];
    hB[iz+2*sz]=q1*(tmpB+q*(tmpAZp+tmpAZm));
    hAZp[iz+2*sz]=tmpAZp+q*hB[iz+2*sz];
    hAZm[iz+2*sz]=tmpAZm+q*hB[iz+2*sz];
  }

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    tmpM.M1M2v(hM[iz],hM[iz+sz],2);
    tmpAZp = hAZp[iz-1]+hAZp[iz+sz-1]*two;
    tmpAZm = hAZm[iz+1]+hAZm[iz+sz+1]*two;
    tmpB   = hB[iz]+hB[iz+sz]*two;
    if(InitData.UseFastP) {
      tmpM.PZp=(tmpAZp*tmpM.R0*two
             +tmpAZm*tmpM.RZ)*InitData.beta;
      tmpM.PZm=(tmpAZm*tmpM.R0*two
             +tmpAZp*conj(tmpM.RZ))*InitData.beta;
      tmpM.PB =tmpB*tmpM.R0*two*InitData.AGenCoef;
    }
    hM[iz+2*sz].PZp=(tmpAZp*tmpM.R0*two+tmpAZm*tmpM.RZ)*InitData.beta*dt2;
    hM[iz+2*sz].PZm=(tmpAZm*tmpM.R0*two+tmpAZp*conj(tmpM.RZ))*InitData.beta*dt2;
    hM[iz+2*sz].PB=(tmpB*tmpM.R0*two*InitData.AGenCoef)*dt2;
    hM[iz+2*sz].R0 =(InitData.Q-InitData.RCoef*tmpM.R0
                       -(tmpAZp*conj(tmpM.PZp)+tmpAZm*conj(tmpM.PZm)).real() )*dt2;
    hM[iz+2*sz].RZ=(-tmpM.RZ*InitData.RCoef
                        -(tmpAZp*conj(tmpM.PZm)
                        +tmpM.PZp*conj(tmpAZm)) )*dt2;
    hM[iz+2*sz].R0Generated=(InitData.Q-InitData.RCoef*tmpM.R0)*dt2;

    hAZp[iz+2*sz]+=dt2*tmpM.PZp+dt2*InitData.AGenZ*tmpAZp;
    hAZm[iz+2*sz]+=dt2*tmpM.PZm+dt2*InitData.AGenZ*tmpAZm;
    hB[iz+2*sz]+=(tmpM.PB+tmpB*BCoef)*dt2;
  }

/*
for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    hM[iz].PZp=tmpPZp=(hAZp[iz]*hM[iz].R0*two)*InitData.beta;
    hM[iz].PZm=tmpPZm=(hAZm[iz]*hM[iz].R0*two)*InitData.beta;
    hM[iz].PB =tmpPB=hB[iz]*hM[iz].R0*two*InitData.AGenCoef;
    hM[iz+sz].R0 =(InitData.Q-InitData.RCoef*hM[iz].R0
                       -(hAZp[iz]*conj(hM[iz].PZp)
                        +hAZm[iz]*conj(hM[iz].PZm)
                        +hB[iz]*conj(hM[iz].PB)
                       ).real() )*dt2;
    hM[iz+sz].RZ=(-hM[iz].RZ*InitData.RCoef
                        -(hAZp[iz]*conj(hM[iz].PZm)
                        +hM[iz].PZp*conj(hAZm[iz])) )*dt2;
    hM[iz+sz].R0Generated=(InitData.Q-InitData.RCoef*hM[iz].R0)*dt2;
    hAZp[iz+sz]=dt2*hM[iz].PZp+dt2*InitData.AGenZ*hAZp[iz]+ii*InitData.ACoef*hB[iz]*dt2;
    hAZm[iz+sz]=dt2*hM[iz].PZm+dt2*InitData.AGenZ*hAZm[iz]+ii*InitData.ACoef*hB[iz]*dt2;
    hB[iz+sz]=(ii*InitData.ACoef*(hAZp[iz]+hAZm[iz])+hM[iz].PB+hB[iz]*BCoef)*dt2;

    hAZp[iz+sz]=dt2*InitData.AGenZ*hAZp[iz]+ii*InitData.ACoef*hB[iz]*dt2;
    hAZm[iz+sz]=dt2*InitData.AGenZ*hAZm[iz]+ii*InitData.ACoef*hB[iz]*dt2;
    hB[iz+sz]=(ii*InitData.ACoef*(hAZp[iz]+hAZm[iz])+hB[iz]*BCoef)*dt2;
  }

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    tmpAZp = hAZp[iz-1]+hAZp[iz+sz-1];
    tmpAZm = hAZm[iz+1]+hAZm[iz+sz+1];
    tmpB   = hB[iz]+hB[iz+sz];
    hB[iz+2*sz]=q1*(tmpB+q*(tmpAZp+tmpAZm));
    hAZp[iz+2*sz]=tmpAZp+q*hB[iz+2*sz];
    hAZm[iz+2*sz]=tmpAZm+q*hB[iz+2*sz];
  }

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    tmpM.M1M2v(hM[iz],hM[iz+sz],2);
    tmpAZp = hAZp[iz-1]+hAZp[iz+sz-1]*two;
    tmpAZm = hAZm[iz+1]+hAZm[iz+sz+1]*two;
    tmpB   = hB[iz]+hB[iz+sz]*two;

    tmpM.PZp=(tmpAZp*tmpM.R0*two)*InitData.beta;
    tmpM.PZm=(tmpAZm*tmpM.R0*two)*InitData.beta;
    tmpM.PB =tmpB*tmpM.R0*two*InitData.AGenCoef;
    hM[iz+2*sz].R0 =(InitData.Q-InitData.RCoef*tmpM.R0
                       -(tmpAZp*conj(tmpM.PZp)+tmpAZm*conj(tmpM.PZm)).real() )*dt2;
    hM[iz+2*sz].RZ=(-tmpM.RZ*InitData.RCoef
                        -(tmpAZp*conj(tmpM.PZm)
                        +tmpM.PZp*conj(tmpAZm)) )*dt2;
    hM[iz+2*sz].R0Generated=(InitData.Q-InitData.RCoef*tmpM.R0)*dt2;

    hAZp[iz+2*sz]+=dt2*InitData.AGenZ*tmpAZp;
    hAZm[iz+2*sz]+=dt2*InitData.AGenZ*tmpAZm;
    hB[iz+2*sz]+=(tmpB*BCoef)*dt2;
  }
*/
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    tmpM.M1M2v(hM[iz],hM[iz+sz],1);
    hM[iz].M1M2v(tmpM,hM[iz+2*sz],1);
    hAZp[iz]=hAZp[iz+2*sz];
    hAZm[iz]=hAZm[iz+2*sz];
    hB[iz]=hB[iz+2*sz];
    /*
    hAZp[InitData.nStepsZ+1-iz]=hAZp[InitData.nStepsZ-iz]
          +hAZp[InitData.nStepsZ-iz+sz]+hAZp[InitData.nStepsZ+1-iz+2*sz];
    hAZm[iz]=hAZm[iz+1]+hAZm[iz+1+sz]+hAZm[iz+2*sz];
    hB[iz]=hB[iz]+hB[iz+sz]+hB[iz+2*sz];
    */
  }

  CalcEnergy();
  hEnergyData[nStepNumber].OutGenX=-PGenX;
  hEnergyData[nStepNumber].OutGenZ=-PGenZ;
  hEnergyData[nStepNumber].PRel=PRel;

  EnergyGenX+=PGenX;
  EnergyGenZ+=PGenZ;

  FillBorders();
  CalcOutEnergy();

  double gen = (EnergyGenX+EnergyGenZ)*InitData.dt;
  double ea = EnergyZP+EnergyZM+EnergyB;
  double out = EnergyZPOut+EnergyZMOut;
  double deltaEA=ea-EnergyA0;
  double deltaER=EnergyR-EnergyR0;
  double eq = EnergyQ-EnergyQBack;
  double balance=deltaEA-gen+out-eq+deltaER;

  int fDebugReset=0;
  if(fDebugReset) {
    EnergyA0=EnergyZP+EnergyZM+EnergyB;
    EnergyR0=EnergyR;
    EnergyQBack=EnergyQ;
    EnergyQ=EnergyZPOut=EnergyZMOut=EnergyGenZ=EnergyGenX=0;
  }


  hAZpOutData[nStepNumber]=hAZp[InitData.nStepsZ-1];
  for(int i=0;i<FFTSizeFull-1;i++) hFFTSrc[i]=hFFTSrc[i+1];
  hFFTSrc[FFTSizeFull-1]=hAZp[InitData.nStepsZ/2];
}





void cEquOptic2D_1DwB::StepRoutineRunge2() {
  complex<OPTIC_TYPE> tmpPZp,tmpPZm,tmpPB;
  OPTIC_TYPE dt2=0.5*InitData.dt;
  complex<OPTIC_TYPE> ii=complex<OPTIC_TYPE>(0,1),tmpAZp,tmpAZm,tmpB;
  complex<OPTIC_TYPE> BCoef=complex<OPTIC_TYPE>(InitData.AGenX,InitData.Delta);
  cOptic2D_Media_2Way_B tmpM;
  OPTIC_TYPE two=2;
  int sz=InitData.SizeZ;
  //CalcEnergy();
  OPTIC_TYPE ER=EnergyR,EQ=EnergyQ,EB=EnergyB,EZp=EnergyZP,EZm=EnergyZM;
  OPTIC_TYPE dEo=abs(hAZp[InitData.nStepsZ])*abs(hAZp[InitData.nStepsZ])+abs(hAZm[1])*abs(hAZm[1]);


                        // InitData.AGenCoef - коэффициент связи среды с волной B !
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    if(InitData.UseFastP) {
      hM[iz].PZp=tmpPZp=(hAZp[iz]*hM[iz].R0*two+hAZm[iz]*hM[iz].RZ)*InitData.beta;
      hM[iz].PZm=tmpPZm=(hAZm[iz]*hM[iz].R0*two+hAZp[iz]*conj(hM[iz].RZ))*InitData.beta;
      hM[iz].PB =tmpPB=hB[iz]*hM[iz].R0*two*InitData.AGenCoef;
    }
    hM[iz+sz].PZp=(hAZp[iz]*hM[iz].R0*two+hAZm[iz]*hM[iz].RZ)*InitData.beta*dt2;
    hM[iz+sz].PZm=(hAZm[iz]*hM[iz].R0*two+hAZp[iz]*conj(hM[iz].RZ))*InitData.beta*dt2;
    hM[iz+sz].PB=(hB[iz]*hM[iz].R0*two*InitData.AGenCoef)*dt2;
    hM[iz+sz].R0 =(InitData.Q-InitData.RCoef*hM[iz].R0
                       -(hAZp[iz]*conj(hM[iz].PZp)
                        +hAZm[iz]*conj(hM[iz].PZm)
                        +hB[iz]*conj(hM[iz].PB)
                       ).real() )*dt2;
    hM[iz+sz].RZ=(-hM[iz].RZ*InitData.RCoef
                        -(hAZp[iz]*conj(hM[iz].PZm)
                        +hM[iz].PZp*conj(hAZm[iz])) )*dt2;
    hM[iz+sz].R0Generated=(InitData.Q-InitData.RCoef*hM[iz].R0)*dt2;
    hAZp[iz+sz]=dt2*hM[iz].PZp+dt2*InitData.AGenZ*hAZp[iz]+ii*InitData.ACoef*hB[iz]*dt2;
    hAZm[iz+sz]=dt2*hM[iz].PZm+dt2*InitData.AGenZ*hAZm[iz]+ii*InitData.ACoef*hB[iz]*dt2;
    hB[iz+sz]=(ii*InitData.ACoef*(hAZp[iz]+hAZm[iz])+hM[iz].PB+hB[iz]*BCoef)*dt2;

    /*
    hAZp[iz+sz]=ii*InitData.ACoef*hB[iz]*dt2;
    hAZm[iz+sz]=ii*InitData.ACoef*hB[iz]*dt2;
    hB[iz+sz]=(ii*InitData.ACoef*(hAZp[iz]+hAZm[iz]))*dt2;
    */
  }

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    /*
    tmpAZp = hAZp[iz]+hAZp[iz+sz]*two;
    tmpAZm = hAZm[iz]+hAZm[iz+sz]*two;
    tmpB   = hB[iz]+hB[iz+sz]*two;
    */
    tmpM.M1M2v(hM[iz],hM[iz+sz],2);

    tmpAZp = hAZp[iz-1]+hAZp[iz+sz-1]*two;
    tmpAZm = hAZm[iz+1]+hAZm[iz+sz+1]*two;
    tmpB   = hB[iz]+hB[iz+sz]*two;
    

    if(InitData.UseFastP) {
      tmpM.PZp=(tmpAZp*tmpM.R0*two
             +tmpAZm*tmpM.RZ)*InitData.beta;
      tmpM.PZm=(tmpAZm*tmpM.R0*two
             +tmpAZp*conj(tmpM.RZ))*InitData.beta;
      tmpM.PB =tmpB*tmpM.R0*two*InitData.AGenCoef;
    }
    hM[iz+2*sz].PZp=(tmpAZp*tmpM.R0*two+tmpAZm*tmpM.RZ)*InitData.beta*dt2;
    hM[iz+2*sz].PZm=(tmpAZm*tmpM.R0*two+tmpAZp*conj(tmpM.RZ))*InitData.beta*dt2;
    hM[iz+2*sz].PB=(tmpB*tmpM.R0*two*InitData.AGenCoef)*dt2;
    hM[iz+2*sz].R0 =(InitData.Q-InitData.RCoef*tmpM.R0
                       -(tmpAZp*conj(tmpM.PZp)+tmpAZm*conj(tmpM.PZm)).real() )*dt2;
    hM[iz+2*sz].RZ=(-tmpM.RZ*InitData.RCoef
                        -(tmpAZp*conj(tmpM.PZm)
                        +tmpM.PZp*conj(tmpAZm)) )*dt2;
    hM[iz+2*sz].R0Generated=(InitData.Q-InitData.RCoef*tmpM.R0)*dt2;
    hAZp[iz+2*sz]=dt2*tmpM.PZp+dt2*InitData.AGenZ*tmpAZp+ii*InitData.ACoef*tmpB*dt2;
    hAZm[iz+2*sz]=dt2*tmpM.PZm+dt2*InitData.AGenZ*tmpAZm+ii*InitData.ACoef*tmpB*dt2;
    hB[iz+2*sz]=(ii*InitData.ACoef*(tmpAZp+tmpAZm)+tmpM.PB+tmpB*BCoef)*dt2;
    /*
    hAZp[iz+2*sz]=ii*InitData.ACoef*tmpB*dt2;
    hAZm[iz+2*sz]=ii*InitData.ACoef*tmpB*dt2;
    hB[iz+2*sz]=(ii*InitData.ACoef*(tmpAZp+tmpAZm))*dt2;
    */
  }
  /*
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    tmpM.M1M2v(hM[iz],hM[iz+sz],1);
    hM[iz].M1M2v(tmpM,hM[iz+2*sz],1);
    hAZp[iz]=hAZp[iz]+hAZp[iz+sz]+hAZp[iz+2*sz];
    hAZm[iz]=hAZm[iz]+hAZm[iz+sz]+hAZm[iz+2*sz];
    hB[iz]=hB[iz]+hB[iz+sz]+hB[iz+2*sz];
  }
  */

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    tmpM.M1M2v(hM[iz],hM[iz+sz],1);
    hM[iz].M1M2v(tmpM,hM[iz+2*sz],1);
    hAZp[InitData.nStepsZ+1-iz]=hAZp[InitData.nStepsZ-iz]
          +hAZp[InitData.nStepsZ-iz+sz]+hAZp[InitData.nStepsZ+1-iz+2*sz];
    hAZm[iz]=hAZm[iz+1]+hAZm[iz+1+sz]+hAZm[iz+2*sz];
    hB[iz]=hB[iz]+hB[iz+sz]+hB[iz+2*sz];
  }

  CalcEnergy();
  hEnergyData[nStepNumber].OutGenX=-PGenX;
  hEnergyData[nStepNumber].OutGenZ=-PGenZ;
  hEnergyData[nStepNumber].PRel=PRel;

  dEo+=abs(hAZp[InitData.nStepsZ])*abs(hAZp[InitData.nStepsZ])+abs(hAZm[1])*abs(hAZm[1]);
  dEo*=0.25*InitData.dt;
  OPTIC_TYPE E1=EnergyZP+EnergyZM+EnergyB+EnergyQ;
  OPTIC_TYPE E0=EZp+EZm+EB+EQ;
  OPTIC_TYPE dE=E1+dEo-E0;
  OPTIC_TYPE dEr=dE/E0;
  OPTIC_TYPE dEr1=dE/E0/InitData.dt;
  OPTIC_TYPE dEr2=dE/dEo;
    
  FillBorders();
  CalcOutEnergy();
  hAZpOutData[nStepNumber]=hAZp[InitData.nStepsZ-1];
  for(int i=0;i<FFTSizeFull-1;i++) hFFTSrc[i]=hFFTSrc[i+1];
  hFFTSrc[FFTSizeFull-1]=hAZp[InitData.nStepsZ/2];
}

void cEquOptic2D_1DwB::CalcEnergy() {
  EnergyB=EnergyZP=EnergyZM=EnergyR=EnergyQ=0;
  PGenR=PGenA=PGenZ=PGenX=0;
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    EnergyB+=hB[iz].real()*hB[iz].real()+hB[iz].imag()*hB[iz].imag();
    EnergyZP+=hAZp[iz].real()*hAZp[iz].real()+hAZp[iz].imag()*hAZp[iz].imag();
    EnergyZM+=hAZm[iz].real()*hAZm[iz].real()+hAZm[iz].imag()*hAZm[iz].imag();
    EnergyR+=hM[iz].R0;
    EnergyQ+=hM[iz].R0Generated;
    PGenX+=InitData.AGenX*(hB[iz].real()*hB[iz].real()+hB[iz].imag()*hB[iz].imag());
    PGenZ+=InitData.AGenZ*(hAZp[iz].real()*hAZp[iz].real()+hAZp[iz].imag()*hAZp[iz].imag());
    PGenZ+=InitData.AGenZ*(hAZm[iz].real()*hAZm[iz].real()+hAZm[iz].imag()*hAZm[iz].imag());
    PGenR+=InitData.Q-hM[iz].R0*InitData.RCoef;
  }
  PGenX*=InitData.dt;
  PGenZ*=InitData.dt;
  PGenA=PGenX+PGenZ;
  PGenR*=InitData.dt;
  POutA=hAZp[InitData.nStepsZ].real()*hAZp[InitData.nStepsZ].real()+
                hAZp[InitData.nStepsZ].imag()*hAZp[InitData.nStepsZ].imag();
  POutA+=hAZm[1].real()*hAZm[1].real()+hAZm[1].imag()*hAZm[1].imag();
  POutA*=0.5;
  PRel=-PGenR+InitData.Q*InitData.dt*InitData.nStepsZ;

  EnergyB*=0.5*InitData.dt;
  EnergyZP*=0.5*InitData.dt;
  EnergyZM*=0.5*InitData.dt;
  EnergyR*=InitData.dt;
  EnergyQ*=InitData.dt;
}

void cEquOptic2D_1DwB::CalcOutEnergy() {
  OPTIC_TYPE tmpXP=0,tmpXM=0,tmpZP=0,tmpZM=0,tmpGenZ=0,tmpGenX=0,tmpR=0;
  tmpZP+= hAZp[InitData.nStepsZ].real()*hAZp[InitData.nStepsZ].real()
         +hAZp[InitData.nStepsZ].imag()*hAZp[InitData.nStepsZ].imag();
  tmpZM+=  hAZm[1].real()*hAZm[1].real()
          +hAZm[1].imag()*hAZm[1].imag();
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    tmpR+=hM[iz].R0;
  }
  hEnergyData[nStepNumber].R=tmpR*InitData.dt;
//  tmpZP*=0.5*InitData.dt;  EnergyZPOut+=tmpZP;
//  tmpZM*=0.5*InitData.dt;  EnergyZMOut+=tmpZM;
  tmpZP*=0.5*InitData.dt;  EnergyZPOut+=tmpZP;
  tmpZM*=0.5*InitData.dt;  EnergyZMOut+=tmpZM;
  hEnergyData[nStepNumber].OutZp=tmpZP/InitData.dt;
  hEnergyData[nStepNumber].OutZm=tmpZM/InitData.dt;
  hEnergyData[nStepNumber].OutGenZ=tmpGenZ/InitData.dt;
  hEnergyData[nStepNumber].Out=(tmpZP+tmpZM)/InitData.dt;
  /*
  if( (hEnergyData[nStepNumber].Out>InitData.MaxV[0])&&
      (LastMaxOut>hEnergyData[nStepNumber].Out)&&
      (!MaxESaved))
              SaveADataMaxE();
  LastMaxOut=hEnergyData[nStepNumber].Out;
  */
}

