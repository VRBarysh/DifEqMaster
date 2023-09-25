//---------------------------------------------------------------------------


#pragma hdrstop

#include <fstream>
#include "EquOptic2D04.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

void cEquOptic2D_01::InitReport(cDifEquReport *hReport) {
  int Lines[EQU_MAX_GRAPHS];
  for(int i=0;i<EQU_MAX_GRAPHS;i++) Lines[i]=1;
  Lines[0]=1; Lines[1]=1; Lines[2]=1; Lines[3]=1; Lines[4]=1;
  hReportText->nGraphs=10;
  hReportText->GraphText[0]="AZm(Z)";
  hReportText->GraphText[1]="AZp(Z)";
  hReportText->GraphText[2]="AZmBragg(Z)";
  hReportText->GraphText[3]="AZpBragg(Z)";
  hReportText->GraphText[6]="AZm(X)";
  hReportText->GraphText[7]="AZp(X)";
  hReportText->GraphText[8]="AZmBragg(X)";
  hReportText->GraphText[9]="AZpBragg(X)";
  hReportText->GraphText[4]="Spectrum AZp";
  hReportText->GraphText[5]="AZp(t)";
  hReport->LoadText(hReportText);
  hReport->SetGraphSize(EQU_GRAPH_SIZE,Lines);
  hReport->SetMaxTo1();
}
void cEquOptic2D_01::Report(cDifEquReport *hReport) {
  int BufSize=InitData.nStepsX > InitData.nStepsZ ? InitData.nStepsX : InitData.nStepsZ;
  BufSize=BufSize > InitData.nStepsZB ? BufSize : InitData.nStepsZB;
  BufSize=BufSize > InitData.nStepsDelay ? BufSize : InitData.nStepsDelay;
  long double *hData = new long double [BufSize]; 
  for(int i=1; i<=InitData.nStepsZ; i++)
    hData[i-1]=abs(hMedia[Index(i,InitData.SizeX/2)].AZm);
  hReport->Graph[0].GetData(0,InitData.nStepsZ,hData);
  for(int i=1; i<=InitData.nStepsZ; i++)
    hData[i-1]=abs(hMedia[Index(i,InitData.SizeX/2)].AZp);
  hReport->Graph[1].GetData(0,InitData.nStepsZ,hData);
  for(int i=1; i<=InitData.nStepsZB; i++)
    hData[i-1]=abs(hMediaBragg[Index(i,InitData.SizeX/2)].AZm);
  hReport->Graph[2].GetData(0,InitData.nStepsZB,hData);
  for(int i=1; i<=InitData.nStepsZB; i++)
    hData[i-1]=abs(hMediaBragg[Index(i,InitData.SizeX/2)].AZp);
  hReport->Graph[3].GetData(0,InitData.nStepsZB,hData);

  for(int i=1; i<=InitData.nStepsX; i++)
    hData[i-1]=abs(hMedia[Index(InitData.SizeZ/2,i)].AZm);
  hReport->Graph[6].GetData(0,InitData.nStepsX,hData);
  for(int i=1; i<=InitData.nStepsX; i++)
    hData[i-1]=abs(hMedia[Index(InitData.SizeZ/2,i)].AZp);
  hReport->Graph[7].GetData(0,InitData.nStepsX,hData);
  for(int i=1; i<=InitData.nStepsX; i++)
    hData[i-1]=abs(hMediaBragg[Index(InitData.SizeZB/2,i)].AZm);
  hReport->Graph[8].GetData(0,InitData.nStepsX,hData);
  for(int i=1; i<=InitData.nStepsX; i++)
    hData[i-1]=abs(hMediaBragg[Index(InitData.SizeZB/2,i)].AZp);
  hReport->Graph[9].GetData(0,InitData.nStepsX,hData);

  hReport->nVals=9;
  CalcEnergy();
  hReport->ValText[0]="Energy sum "; hReport->Val[0]=EnergyZP+EnergyZM+
        EnergyXPB+EnergyXMB+EnergyZPB+EnergyZMB+EnergyZPT+EnergyZMT;
  hReport->ValText[1]="Energy media "; hReport->Val[1]=EnergyZP+EnergyZM;
  hReport->ValText[2]="Energy bragg "; hReport->Val[2]=EnergyXPB+EnergyXMB+EnergyZPB+EnergyZMB;
  hReport->ValText[3]="Energy ZP "; hReport->Val[3]=EnergyZP;
  hReport->ValText[4]="Energy ZM "; hReport->Val[4]=EnergyZM;
  hReport->ValText[5]="Energy sum+out "; hReport->Val[5]=EnergyZP+EnergyZM+
        EnergyXPB+EnergyXMB+EnergyZPB+EnergyZMB+EnergyZPT+EnergyZMT+
        EnergyXPOut+EnergyXMOut+EnergyZPOut+EnergyZMOut;
  hReport->ValText[6]="Energy R "; hReport->Val[6]=EnergyR0;
  hReport->ValText[7]="Energy sum+out+R "; hReport->Val[7]=
        EnergyR0+
        EnergyZP+EnergyZM+EnergyXPB+EnergyXMB+
        EnergyZPB+EnergyZMB+EnergyZPT+EnergyZMT+
        EnergyXPOut+EnergyXMOut+EnergyZPOut+EnergyZMOut;
  hReport->ValText[8]="ZP inc "; hReport->Val[8]=decZP;

  FFT.TransForward(hFFTDest,hFFTSrc);
  TCplxLong tmp;
  for(int i=0;i<(FFTSizeFull >> 1);i++) {
    tmp=hFFTDest[i];
    hFFTDest[i]=hFFTDest[i+(FFTSizeFull >> 1)];
    hFFTDest[i+(FFTSizeFull >> 1)]=tmp;
  }
  /*
  hReport->Graph[4].GetDataAbs(0,FFTSizeFull,hFFTDest);
  hReport->Graph[5].GetDataAbs(0,FFTSizeFull,hFFTSrc);
  */
//  hReport->Graph[2].GetDataRe(0,FFTSizeFull,hFFTDest);
//  hReport->Graph[2].Max = hReport->Graph[2].Max > -hReport->Graph[2].Min ?
//                          hReport->Graph[2].Max : -hReport->Graph[2].Min;
//  hReport->Graph[3].GetDataAbs(0,FFTSizeFull,hFFTDest);
//  hReport->Graph[3].Max = hReport->Graph[3].Max > -hReport->Graph[3].Min ?
//                          hReport->Graph[3].Max : -hReport->Graph[3].Min;
  hReport->Val[0]=hReport->Graph[2].Max;
  hReport->Val[1]=hReport->Graph[3].Max;
  hReport->Graph[4].GetDataAbs(0,(FFTSizeFull >> InitData.FFTScale),
             hFFTDest+(FFTSizeFull>>1)-(FFTSizeFull>> (1+InitData.FFTScale) ) );
  hReport->Graph[5].GetDataAbs(0,FFTSizeFull,hFFTSrc);
//  hReport->SetMaxTo1();
  /*
  hReport->nVals=3;
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
  hReport->ValText[1]="Max I "; hReport->Val[1]=hReport->Graph[1].Max;
  hReport->ValText[2]="Max Gamma "; hReport->Val[2]=hReport->Graph[2].Max;
  hReport->Graph[2].Max=2; */
  LastReportTime=t;
  delete [] hData;
}
void cEquOptic2D_01::FinalReport(cDifEquReport *hReport) {
  hReport->nVals=0; hReport->nGraphs=0;
/*  hReport->Graph[0].GetDataAbs(0,nStepNumber,hOutA);
  hReport->GraphText[0]="Out A";
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;*/
  SaveAData();
}

void cEquOptic2D_01::LoadInitData(cDifEqu_InitDataRec *hInitData)
{
  cDifEqu_InitDataRec_Optic2D *hInit=(cDifEqu_InitDataRec_Optic2D *)hInitData;
  int nTStepsNew=hInit->tMax/hInit->dt;
  hInit->MediaSize=(hInit->nStepsX+2)*(hInit->nStepsZ+2);
  hInit->MediaSizeB=(hInit->nStepsX+2)*(hInit->nStepsZB+2);
  if ( InitData.nStepsT!=nTStepsNew ) {
//    delete [] hOutA; hOutA=new TCplxLong[nTStepsNew+1];
  }
  if( InitData.MediaSize!=hInit->MediaSize )  {
    delete [] hMedia; delete [] hMediaNew;
    hMedia = new cOptic2D_Media_2Waves[ hInit->MediaSize ];
    hMediaNew = new cOptic2D_Media_2Waves[ hInit->MediaSize ];
  }

  if( InitData.MediaSizeB!=hInit->MediaSizeB )  {
    delete [] hMediaBragg; delete [] hMediaBraggNew;
    hMediaBragg = new cOptic2D_Media_BraggOnly[ hInit->MediaSizeB ];
    hMediaBraggNew = new cOptic2D_Media_BraggOnly[ hInit->MediaSizeB ];
  }

  if ((InitData.nStepsDelay*(InitData.nStepsX+2))!=(hInit->nStepsDelay*(hInit->nStepsX+2))) {
    delete [] hApData; delete [] hAmData;
    hApData = new TCplxLong[hInit->nStepsDelay*(hInit->nStepsX+2)];
    hAmData = new TCplxLong[hInit->nStepsDelay*(hInit->nStepsX+2)];
  }

  delete [] hFFTSrc; delete [] hFFTDest; FFTSizeFull = 1 << hInit->FFTSize;
  hFFTSrc = new TCplxLong[FFTSizeFull];
  hFFTDest = new TCplxLong[FFTSizeFull];
  for (int i=0;i<FFTSizeFull;i++) hFFTSrc[i]=hFFTDest[i]=0;

  InitData=*hInit; nStepNumber=0; InitData.nStepsT=nTStepsNew;
  t=0; tMax=InitData.tMax; dt=InitData.dt; InitData.DelayPos=0;
  if(InitData.UseFastP) InitData.PCoefFast=1.0/InitData.PCoef;
  InitData.SizeX=InitData.nStepsX+2; InitData.SizeZ=InitData.nStepsZ+2;
  InitData.SizeZB=InitData.nStepsZB+2;
  FFT.MakeRevTable(InitData.FFTSize);
  PrepareMedia(); EnergyZPlast=0; ESumLast=0; LastReportTime=LastDecTime=0;
  EnergyXPOut=EnergyXMOut=EnergyZPOut=EnergyZMOut=0;
}

void cEquOptic2D_01::PrepareMedia() {
  for(int i=0;i<=InitData.nStepsZ+1;i++)
    for(int j=0;j<=InitData.nStepsX+1;j++) {
      hMedia[Index(i,j)].AZp=InitData.A0*sin(2*M_PI*i/InitData.nStepsZ)+
         InitData.A1*sin(M_PI*i/InitData.nStepsZ)*sin(M_PI*j/InitData.nStepsX);
      hMedia[Index(i,j)].AZm=InitData.A0*sin(3*M_PI*i/InitData.nStepsZ)+
         InitData.A1*sin(M_PI*i/InitData.nStepsZ)*sin(M_PI*j/InitData.nStepsX);
   /*
      hMedia[i*InitData.SizeX+j].AZp=InitData.A0;
      hMedia[i*InitData.SizeX+j].AZm=InitData.A0;
   */
      hMedia[Index(i,j)].PZp=0;
      hMedia[Index(i,j)].PZm=0;
      hMedia[Index(i,j)].R0=1;
      hMedia[Index(i,j)].RZ=0;
      hMedia[Index(i,j)].hInitData=&InitData;
    }
  for(int i=0;i<=InitData.nStepsZB+1;i++)
    for(int j=0;j<=InitData.nStepsX+1;j++) {
      hMediaBragg[Index(i,j)].AZp=InitData.A0*sin(M_PI*i/InitData.nStepsZ)+
         InitData.A1*sin(M_PI*i/InitData.nStepsZ)*sin(M_PI*j/InitData.nStepsX);
      hMediaBragg[Index(i,j)].AZm=InitData.A0*sin(4*M_PI*i/InitData.nStepsZ)+
         InitData.A1*sin(M_PI*i/InitData.nStepsZ)*sin(M_PI*j/InitData.nStepsX);
      hMediaBragg[Index(i,j)].AXp=InitData.A0*sin(M_PI*i/InitData.nStepsZ)+
         InitData.A1*sin(M_PI*i/InitData.nStepsZ)*sin(M_PI*j/InitData.nStepsX);
      hMediaBragg[Index(i,j)].AXm=InitData.A0*sin(M_PI*i/InitData.nStepsZ)+
         InitData.A1*sin(M_PI*i/InitData.nStepsZ)*sin(M_PI*j/InitData.nStepsX);
      hMediaBragg[Index(i,j)].hInitData=&InitData;
    }
  for(int i=0;i<InitData.nStepsDelay;i++)
    for(int j=0;j<=InitData.nStepsX+1;j++) {
      hApData[Index(i,j)]=0;
      hAmData[Index(i,j)]=0;
    }
  FillBorders();
}

void cEquOptic2D_01::FillBorders() {
  for(int i=0;i<InitData.SizeX;i++) {
    hAmData[Index(InitData.DelayPos,i)]=hMedia[Index(InitData.SizeZ-2,i)].AZm;
    hApData[Index(InitData.DelayPos,i)]=hMediaBragg[Index(1,i)].AZp;
  }
  InitData.DelayPos=(InitData.DelayPos+1)%InitData.nStepsDelay;
  for(int i=1;i<InitData.SizeX-1;i++) {
    hMedia[Index(0,i)]=hMedia[Index(1,i)];
    hMedia[Index(0,i)].AZm=BorderAZm(i);
    hMedia[Index(InitData.SizeZ-1,i)]=hMedia[Index(InitData.SizeZ-2,i)];
    hMedia[Index(InitData.SizeZ-1,i)].AZp=BorderAZp(i);
  }
  for(int i=1;i<InitData.SizeZB-1;i++) {
    hMediaBragg[Index(i,0)]=hMediaBragg[Index(i,1)];
    hMediaBragg[Index(i,0)].AXm=BorderAXmBragg(i);
    hMediaBragg[Index(i,InitData.SizeX-1)]=hMediaBragg[Index(i,InitData.SizeX-2)];
    hMediaBragg[Index(i,InitData.SizeX-1)].AXp=BorderAXpBragg(i);
  }
  for(int i=1;i<InitData.SizeX-1;i++) {
    hMediaBragg[Index(0,i)]=hMediaBragg[Index(1,i)];
    hMediaBragg[Index(0,i)].AZm=BorderAZmBragg(i);
    hMediaBragg[Index(InitData.SizeZB-1,i)]=hMediaBragg[Index(InitData.SizeZB-2,i)];
    hMediaBragg[Index(InitData.SizeZB-1,i)].AZp=BorderAZpBragg(i);
  }
}

void cEquOptic2D_01::StepRoutine() {
  for(int i=1;i<=InitData.nStepsZ;i++)
    for(int j=1;j<=InitData.nStepsX;j++)
      hMedia[Index(i,j)].MakeNewSeveralSteps(hMediaNew+Index(i,j),i,j);
  for(int i=1;i<=InitData.nStepsZ;i++)
    for(int j=1;j<=InitData.nStepsX;j++)
      hMedia[Index(i,j)]=hMediaNew[Index(i,j)];
  for(int i=1;i<=InitData.nStepsZB;i++)
    for(int j=1;j<=InitData.nStepsX;j++)
      hMediaBragg[Index(i,j)].MakeNewSeveralSteps(hMediaBraggNew+Index(i,j),i,j);
  for(int i=1;i<=InitData.nStepsZB;i++)
    for(int j=1;j<=InitData.nStepsX;j++)
      hMediaBragg[Index(i,j)]=hMediaBraggNew[Index(i,j)];
  CalcOutEnergy();
  FillBorders();
  for(int i=0;i<=InitData.nStepsZB;i++)
    for(int j=0;j<=InitData.nStepsX;j++) {
      hMediaBragg[Index(i,j)].AXp=hMediaBragg[Index(i,j+1)].AXp;
      hMediaBragg[Index(i,j)].AZp=hMediaBragg[Index(i+1,j)].AZp;
    }
  for(int i=InitData.nStepsZB+1;i>=1;i--)
    for(int j=InitData.nStepsX+1;j>=1;j--) {
      hMediaBragg[Index(i,j)].AXm=hMediaBragg[Index(i,j-1)].AXm;
      hMediaBragg[Index(i,j)].AZm=hMediaBragg[Index(i-1,j)].AZm;
    }
  for(int i=0;i<=InitData.nStepsZ;i++)
    for(int j=0;j<=InitData.nStepsX;j++) {
      hMedia[Index(i,j)].AZp=hMedia[Index(i+1,j)].AZp;
    }
  for(int i=InitData.nStepsZ+1;i>=1;i--)
    for(int j=InitData.nStepsX+1;j>=1;j--) {
      hMedia[Index(i,j)].AZm=hMedia[Index(i-1,j)].AZm;
    }
  for(int i=0;i<FFTSizeFull-1;i++) hFFTSrc[i]=hFFTSrc[i+1];
  hFFTSrc[FFTSizeFull-1]=hMedia[(InitData.SizeZ-2)*InitData.SizeX+InitData.SizeX/2].AZp;
}

void cEquOptic2D_01::CalcOutEnergy(){
  long double tmpEnergyXP=0,tmpEnergyXM=0,tmpEnergyZP=0,tmpEnergyZM=0;
  for(int i=InitData.nStepsZB;i>=1;i--) {
    tmpEnergyXM+=abs(hMediaBragg[Index(i,InitData.nStepsX)].AXm)*abs(hMediaBragg[Index(i,InitData.nStepsX)].AXm);
    tmpEnergyXP+=abs(hMediaBragg[Index(i,1)].AXp)*abs(hMediaBragg[Index(i,1)].AXp);
  }
  for(int j=InitData.nStepsX;j>=1;j--) {
    tmpEnergyZM+=abs(hMediaBragg[Index(InitData.nStepsZB,j)].AZm)*abs(hMediaBragg[Index(InitData.nStepsZB,j)].AZm);
    tmpEnergyZP+=abs(hMedia[Index(1,j)].AZp)*abs(hMedia[Index(1,j)].AZp);
  }
  EnergyXPOut+=InitData.dt*InitData.dt*tmpEnergyXP;
  EnergyXMOut+=InitData.dt*InitData.dt*tmpEnergyXM;
  EnergyZPOut+=InitData.dt*InitData.dt*tmpEnergyZP*(1.0-InitData.Z0RefCoef*InitData.Z0RefCoef);
  EnergyZMOut+=InitData.dt*InitData.dt*tmpEnergyZM;
}

void cEquOptic2D_01::CalcEnergy(){
  EnergyXPB=EnergyXMB=EnergyZPB=EnergyZMB=EnergyR0=EnergyZP=EnergyZM=EnergyZPT=EnergyZMT=0;
  for(int i=InitData.nStepsZ;i>=1;i--)
    for(int j=InitData.nStepsX;j>=1;j--) {
      EnergyZM+=abs(hMedia[i*InitData.SizeX+j].AZm)*abs(hMedia[i*InitData.SizeX+j].AZm);
      EnergyZP+=abs(hMedia[i*InitData.SizeX+j].AZp)*abs(hMedia[i*InitData.SizeX+j].AZp);
      EnergyR0+=hMedia[i*InitData.SizeX+j].R0;
    }
  for(int i=InitData.nStepsZB;i>=1;i--)
    for(int j=InitData.nStepsX;j>=1;j--) {
      EnergyZMB+=abs(hMediaBragg[i*InitData.SizeX+j].AZm)*abs(hMediaBragg[i*InitData.SizeX+j].AZm);
      EnergyZPB+=abs(hMediaBragg[i*InitData.SizeX+j].AZp)*abs(hMediaBragg[i*InitData.SizeX+j].AZp);
      EnergyXMB+=abs(hMediaBragg[i*InitData.SizeX+j].AXm)*abs(hMediaBragg[i*InitData.SizeX+j].AXm);
      EnergyXPB+=abs(hMediaBragg[i*InitData.SizeX+j].AXp)*abs(hMediaBragg[i*InitData.SizeX+j].AXp);
    }
  for(int i=0;i<InitData.nStepsDelay;i++)
    for(int j=InitData.nStepsX;j>=1;j--) {
      EnergyZMT+=abs(hAmData[i*InitData.SizeX+j])*abs(hAmData[i*InitData.SizeX+j]);
      EnergyZPT+=abs(hApData[i*InitData.SizeX+j])*abs(hApData[i*InitData.SizeX+j]);
    }
  EnergyZP*=InitData.dt*InitData.dt;
  EnergyZM*=InitData.dt*InitData.dt;
  EnergyR0*=InitData.dt*InitData.dt*2;
  EnergyXPB*=InitData.dt*InitData.dt;
  EnergyXMB*=InitData.dt*InitData.dt;
  EnergyZPB*=InitData.dt*InitData.dt;
  EnergyZMB*=InitData.dt*InitData.dt;
  EnergyZPT*=InitData.dt*InitData.dt;
  EnergyZMT*=InitData.dt*InitData.dt;
  long double ESum=EnergyZP+EnergyZM+EnergyXPB+EnergyXMB+EnergyZPB+EnergyZMB+EnergyZPT+EnergyZMT;
  if (t-LastDecTime>0.00000000000000001) {
    if( (ESumLast>0.000000000000000001)&&(ESum>0.000000000000000001))
      decZP=0.5*logl((ESum)/ESumLast)/((t-LastDecTime));
    else decZP=0;
    LastDecTime=t;
    EnergyZPlast=EnergyZP;
    ESumLast=ESum;
  }
}

void cEquOptic2D_01::SaveAData(){
  ofstream FileOutAZp = ofstream("AZpM.dat",ios::out);
  ofstream FileOutAZpB = ofstream("AZpB.dat",ios::out);
  ofstream FileOutAXpB = ofstream("AXpB.dat",ios::out);
  ofstream FileOutAZm = ofstream("AZmM.dat",ios::out);
  ofstream FileOutAZmB = ofstream("AZmB.dat",ios::out);
  ofstream FileOutAXmB = ofstream("AXmB.dat",ios::out);
  ofstream FileOutR0 = ofstream("R0.dat",ios::out);
  ofstream FileOutRZ = ofstream("RZ.dat",ios::out);
  for(int i=0;i<InitData.nStepsZ;i++)
    for(int j=0;j<InitData.nStepsX;j++) {
      FileOutAZm << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hMedia[Index(i+1,j+1)].AZm) << endl;
      FileOutAZp << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hMedia[Index(i+1,j+1)].AZp) << endl;
      FileOutRZ << InitData.dt*i << " " << InitData.dt*j
                 << " " << hMedia[Index(i+1,j+1)].RZ << endl;
      FileOutR0 << InitData.dt*i << " " << InitData.dt*j
                 << " " << hMedia[Index(i+1,j+1)].R0 << endl;
    }
  for(int i=0;i<InitData.nStepsZB;i++)
    for(int j=0;j<InitData.nStepsX;j++) {
      FileOutAZmB << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hMediaBragg[Index(i+1,j+1)].AZm) << endl;
      FileOutAXmB << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hMediaBragg[Index(i+1,j+1)].AXm) << endl;
      FileOutAZpB << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hMediaBragg[Index(i+1,j+1)].AZp) << endl;
      FileOutAXpB << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hMediaBragg[Index(i+1,j+1)].AXp) << endl;
    }
}
