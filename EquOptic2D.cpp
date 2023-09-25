//---------------------------------------------------------------------------


#pragma hdrstop

#include <fstream>
#include "EquOptic2D.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

void cEquOptic2D::Prepare2DColor(int Number2D) {
  switch(Number2D) {
    case 0:
      maxAZp=0.000000000000001;
      for(int i=0;i<=InitData.nStepsZ+1;i++)
        for(int j=0;j<=InitData.nStepsX+1;j++) {
          hAbsStore[Index(i,j)].AZp=abs(hMedia[Index(i,j)].AZp);
          if(hAbsStore[Index(i,j)].AZp>maxAZp) maxAZp=hAbsStore[Index(i,j)].AZp;
        }
    break;
    case 1:
      maxAZm=0.000000000000001;
      for(int i=0;i<=InitData.nStepsZ+1;i++)
        for(int j=0;j<=InitData.nStepsX+1;j++) {
          hAbsStore[Index(i,j)].AZm=abs(hMedia[Index(i,j)].AZm);
          if(hAbsStore[Index(i,j)].AZm>maxAZm) maxAZm=hAbsStore[Index(i,j)].AZm;
        }
    break;
    case 2:
      maxAXp=0.000000000000001;
      for(int i=0;i<=InitData.nStepsZ+1;i++)
        for(int j=0;j<=InitData.nStepsX+1;j++) {
          hAbsStore[Index(i,j)].AXp=abs(hMedia[Index(i,j)].AXp);
          if(hAbsStore[Index(i,j)].AXp>maxAXp) maxAXp=hAbsStore[Index(i,j)].AXp;
        }
    break;
    case 3:
      maxAXm=0.000000000000001;
      for(int i=0;i<=InitData.nStepsZ+1;i++)
        for(int j=0;j<=InitData.nStepsX+1;j++) {
          hAbsStore[Index(i,j)].AXm=abs(hMedia[Index(i,j)].AXm);
          if(hAbsStore[Index(i,j)].AXm>maxAXm) maxAXm=hAbsStore[Index(i,j)].AXm;
        }
    break;
  }
}

int cEquOptic2D::Get2DColor(int Number2D, double x, double y) {
  int i0=x*(InitData.nStepsZ-1)+1;
  int i1=i0+1;
  int j0=y*(InitData.nStepsX-1)+1;
  int j1=j0+1;
  double c0=x*(InitData.nStepsZ-1)+1-i0;
  double c1=y*(InitData.nStepsX-1)+1-j0;
  double c00=c0*c1;
  double c01=c0*(1.0-c1);
  double c10=(1.0-c0)*c1;
  double c11=(1.0-c0)*(1.0-c1);
  double tmpinv;
  switch(Number2D) {
    case 0:
      tmpinv=1.0/maxAZp;
//      return(tmpinv*255*(hAbsStore[Index(i0,j0)].AZp*c00+hAbsStore[Index(i0,j0+1)].AZp*c01+
//                  hAbsStore[Index(i0+1,j0)].AZp*c10+hAbsStore[Index(i0+1,j0+1)].AZp*c11));
      return(tmpinv*255*(hAbsStore[Index(i0,j0)].AZp));
    break;
    case 1:
      tmpinv=1.0/maxAZm;
      return(tmpinv*255*(hAbsStore[Index(i0,j0)].AZm*c00+hAbsStore[Index(i0,j0+1)].AZm*c01+
                  hAbsStore[Index(i0+1,j0)].AZm*c10+hAbsStore[Index(i0+1,j0+1)].AZm*c11));
    break;
    case 2:
      tmpinv=1.0/maxAXp;
      return(tmpinv*255*(hAbsStore[Index(i0,j0)].AXp*c00+hAbsStore[Index(i0,j0+1)].AXp*c01+
                  hAbsStore[Index(i0+1,j0)].AXp*c10+hAbsStore[Index(i0+1,j0+1)].AXp*c11));
    break;
    case 3:
      tmpinv=1.0/maxAXm;
      return(tmpinv*255*(hAbsStore[Index(i0,j0)].AXm*c00+hAbsStore[Index(i0,j0+1)].AXm*c01+
                  hAbsStore[Index(i0+1,j0)].AXm*c10+hAbsStore[Index(i0+1,j0+1)].AXm*c11));
    break;
    default:
      return(255*x*y);
  }
}

void cEquOptic2D::InitReport(cDifEquReport *hReport) {
  int Lines[EQU_MAX_GRAPHS];
  for(int i=0;i<EQU_MAX_GRAPHS;i++) Lines[i]=1;
  Lines[0]=1; Lines[1]=1; Lines[2]=1; Lines[3]=1; Lines[4]=1;
  hReport->LoadText(hReportText);
  hReport->SetGraphSize(EQU_GRAPH_SIZE,Lines);
  hReport->SetMaxTo1();
  CalcEnergy();
  nReportNumber++;
  hEnergyData[nReportNumber-1]=EnergyXP+EnergyXM+EnergyZP+EnergyZM+EDelayBuf;
}
void cEquOptic2D::Report(cDifEquReport *hReport) {
  nReportNumber++;
  long double *hData = new long double
     [ InitData.nStepsX > InitData.nStepsZ ? InitData.nStepsX : InitData.nStepsZ];
  for(int i=1; i<=InitData.nStepsZ; i++)
    hData[i-1]=abs(hMedia[Index(i,InitData.SizeX/2)].AZm);
  hReport->Graph[0].GetData(0,InitData.nStepsZ,hData);
  for(int i=1; i<=InitData.nStepsX; i++)
    hData[i-1]=abs(hMedia[Index(InitData.SizeZ/2,i)].AZm);
  hReport->Graph[1].GetData(0,InitData.nStepsX,hData);
  for(int i=1; i<=InitData.nStepsZ; i++)
    hData[i-1]=abs(hMedia[Index(i,InitData.SizeX/2)].AXm);
  hReport->Graph[2].GetData(0,InitData.nStepsZ,hData);
  for(int i=1; i<=InitData.nStepsX; i++)
    hData[i-1]=abs(hMedia[Index(InitData.SizeZ/2,i)].AXm);
  hReport->Graph[3].GetData(0,InitData.nStepsX,hData);
  hReport->nVals=13;
  CalcEnergy();
  hEnergyData[nReportNumber-1]=EnergyXP+EnergyXM+EnergyZP+EnergyZM+EDelayBuf;
  hReport->ValText[0]="Energy sum "; hReport->Val[0]=EnergyXP+EnergyXM+EnergyZP+EnergyZM+EDelayBuf;
  hReport->ValText[1]="Energy XP "; hReport->Val[1]=EnergyXP;
  hReport->ValText[2]="Energy XM "; hReport->Val[2]=EnergyXM;
  hReport->ValText[3]="Energy ZP "; hReport->Val[3]=EnergyZP;
  hReport->ValText[4]="Energy ZM "; hReport->Val[4]=EnergyZM;
  hReport->ValText[5]="Energy sum+out "; hReport->Val[5]=
        EnergyXP+EnergyXM+EnergyZP+EnergyZM+EDelayBuf+
        EnergyXPOut+EnergyXMOut+EnergyZPOut+EnergyZMOut+EDiffOut;
  hReport->ValText[6]="Energy R "; hReport->Val[6]=EnergyR0;
  hReport->ValText[7]="Energy sum+out+R "; hReport->Val[7]=
        EnergyR0+
        EnergyXP+EnergyXM+EnergyZP+EnergyZM+
        EnergyXPOut+EnergyXMOut+EnergyZPOut+EnergyZMOut;
  hReport->ValText[8]="ZP inc "; hReport->Val[8]=decZP;
  hReport->ValText[9]="Dif R Coef "; hReport->Val[9]=DiffEnergy;
  hReport->ValText[10]="Out bragg ";
  hReport->Val[10]=EnergyXPOut+EnergyXMOut+EnergyZPOut+EnergyZMOut;
  hReport->ValText[11]="Out diff "; hReport->Val[11]=EDiffOut;
  hReport->ValText[12]="EWrong "; hReport->Val[12]=InitData.EWrong;

  FFT.TransForward(hFFTDest,hFFTSrc);
  TCplxLong tmp;
  for(int i=0;i<(FFTSizeFull >> 1);i++) {
    tmp=hFFTDest[i];
    hFFTDest[i]=hFFTDest[i+(FFTSizeFull >> 1)];
    hFFTDest[i+(FFTSizeFull >> 1)]=tmp;
  }
//  hReport->Graph[4].GetDataAbs(0,FFTSizeFull,hFFTDest);
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
  hReport->Graph[10].GetData(0,nReportNumber,hEnergyData);
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
void cEquOptic2D::FinalReport(cDifEquReport *hReport) {
  hReport->nVals=0; hReport->nGraphs=0;
/*  hReport->Graph[0].GetDataAbs(0,nStepNumber,hOutA);
  hReport->GraphText[0]="Out A";
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;*/
  SaveAData();
}

void cEquOptic2D::LoadInitData(cDifEqu_InitDataRec *hInitData)
{
  cDifEqu_InitDataRec_Optic2D *hInit=(cDifEqu_InitDataRec_Optic2D *)hInitData;
  int nTStepsNew=hInit->tMax/hInit->dt;
  hInit->MediaSize=(hInit->nStepsX+2)*(hInit->nStepsZ+2);
  if ( InitData.nStepsT!=nTStepsNew ) {
//    delete [] hOutA; hOutA=new TCplxLong[nTStepsNew+1];
      delete [] hEnergyData; hEnergyData = new long double[nTStepsNew+1];
  }
  if( InitData.MediaSize!=hInit->MediaSize )  {
    delete [] hMedia; delete [] hMediaK; delete [] hMediaNew; delete [] hAbsStore;
    hMedia = new cOptic2D_Media[ hInit->MediaSize ];
    hMediaK = new cOptic2D_Media[ hInit->MediaSize*4 ];
    hMediaNew = new cOptic2D_Media[ hInit->MediaSize ];
    hAbsStore = new cOptic2D_AbsStore[ hInit->MediaSize ];
  }
  
  if ((InitData.nStepsDelay*(InitData.nStepsX+2))!=(hInit->nStepsDelay*(hInit->nStepsX+2))) {
    delete [] hMirrorDelayBuf;
    hMirrorDelayBuf = new TCplxLong[hInit->nStepsDelay*(hInit->nStepsX+2)];
  }

  delete [] hFFTSrc; delete [] hFFTDest; FFTSizeFull = 1 << hInit->FFTSize;
  hFFTSrc = new TCplxLong[FFTSizeFull];
  hFFTDest = new TCplxLong[FFTSizeFull];
  for (int i=0;i<FFTSizeFull;i++) hFFTSrc[i]=hFFTDest[i]=0;

  InitData=*hInit; nStepNumber=nReportNumber=0; InitData.nStepsT=nTStepsNew;
  t=0; tMax=InitData.tMax; dt=InitData.dt;
  if(InitData.UseFastP) InitData.PCoefFast=1.0/InitData.PCoef;
  InitData.SizeX=InitData.nStepsX+2; InitData.SizeZ=InitData.nStepsZ+2;
  FFT.MakeRevTable(InitData.FFTSize); InitData.DelayPos=0;
  PrepareMedia(); EnergyZPlast=0; ESumLast=0; LastReportTime=LastDecTime=0;
  EnergyXPOut=EnergyXMOut=EnergyZPOut=EnergyZMOut=EDiffOut=InitData.EWrong=0;
}

void cEquOptic2D::PrepareMedia() {
  for(int i=0;i<=InitData.nStepsZ+1;i++)
    for(int j=0;j<=InitData.nStepsX+1;j++) {
      hMedia[Index(i,j)].AZp=InitData.A0*sin(M_PI*i/InitData.nStepsZ)+
         InitData.A1*sin(M_PI*i/InitData.nStepsZ)*sin(M_PI*j/InitData.nStepsX)+
         InitData.A2*sin(M_PI*i/InitData.nStepsZ)*sin(2*M_PI*j/InitData.nStepsX);
      hMedia[Index(i,j)].AZm=InitData.A0*sin(M_PI*i/InitData.nStepsZ)+
         InitData.A1*sin(M_PI*i/InitData.nStepsZ)*sin(M_PI*j/InitData.nStepsX)+
         InitData.A2*sin(M_PI*i/InitData.nStepsZ)*sin(2*M_PI*j/InitData.nStepsX);
      hMedia[Index(i,j)].AXp=InitData.A0*sin(M_PI*i/InitData.nStepsZ)+
         InitData.A1*sin(M_PI*i/InitData.nStepsZ)*sin(M_PI*j/InitData.nStepsX)+
         InitData.A2*sin(M_PI*i/InitData.nStepsZ)*sin(2*M_PI*j/InitData.nStepsX);
      hMedia[Index(i,j)].AXm=InitData.A0*sin(M_PI*i/InitData.nStepsZ)+
         InitData.A1*sin(M_PI*i/InitData.nStepsZ)*sin(M_PI*j/InitData.nStepsX)+
         InitData.A2*sin(M_PI*i/InitData.nStepsZ)*sin(2*M_PI*j/InitData.nStepsX);
      /*
      hMedia[Index(i,j)].AZp=InitData.A0+
         InitData.A1*sin(M_PI*(j-1)/(InitData.nStepsX-1))+
         InitData.A2*sin(2*M_PI*(j-1)/(InitData.nStepsX-1));
      hMedia[Index(i,j)].AZm=InitData.A0+
         InitData.A1*sin(M_PI*(j-1)/(InitData.nStepsX-1))+
         InitData.A2*sin(2*M_PI*(j-1)/(InitData.nStepsX-1));
      hMedia[Index(i,j)].AXp=InitData.A0+
         InitData.A1*sin(M_PI*(j-1)/(InitData.nStepsX-1))+
         InitData.A2*sin(2*M_PI*(j-1)/(InitData.nStepsX-1));
      hMedia[Index(i,j)].AXm=InitData.A0+
         InitData.A1*sin(M_PI*(j-1)/(InitData.nStepsX-1))+
         InitData.A2*sin(2*M_PI*(j-1)/(InitData.nStepsX-1));*/
   /*
      hMedia[i*InitData.SizeX+j].AZp=InitData.A0;
      hMedia[i*InitData.SizeX+j].AZm=InitData.A0;
      hMedia[i*InitData.SizeX+j].AXp=0;
      hMedia[i*InitData.SizeX+j].AXm=0; */

      hMedia[Index(i,j)].PZp=InitData.rPolNoise*sin(M_PI*i/InitData.nStepsZ)+
         InitData.rPolNoise*sin(M_PI*i/InitData.nStepsZ)*sin(M_PI*j/InitData.nStepsX);
      hMedia[Index(i,j)].PZm=InitData.rPolNoise*sin(M_PI*i/InitData.nStepsZ)+
         InitData.rPolNoise*sin(M_PI*i/InitData.nStepsZ)*sin(M_PI*j/InitData.nStepsX);;
      hMedia[Index(i,j)].PXp=InitData.rPolNoise*sin(M_PI*i/InitData.nStepsZ)+
         InitData.rPolNoise*sin(M_PI*i/InitData.nStepsZ)*sin(M_PI*j/InitData.nStepsX);;
      hMedia[Index(i,j)].PXm=InitData.rPolNoise*sin(M_PI*i/InitData.nStepsZ)+
         InitData.rPolNoise*sin(M_PI*i/InitData.nStepsZ)*sin(M_PI*j/InitData.nStepsX);;
      hMedia[Index(i,j)].R0=1;
      hMedia[Index(i,j)].RZ=0;
      hMedia[Index(i,j)].RX=0;
      hMedia[Index(i,j)].Rp=0;
      hMedia[Index(i,j)].Rm=0;
      hMedia[Index(i,j)].hInitData=&InitData;
    }
  for(int i=0;i<InitData.nStepsDelay;i++)
    for(int j=0;j<=InitData.nStepsX+1;j++) {
      hMirrorDelayBuf[Index(i,j)]=0;
         /*InitData.A0+
         InitData.A1*sin(M_PI*(j-1)/(InitData.nStepsX-1))+
         InitData.A2*sin(2*M_PI*(j-1)/(InitData.nStepsX-1));;*/
    }
  for(int i=0;i<InitData.nStepsT+1;i++)
    hEnergyData[i]=0;
  FillBorders();
}

void cEquOptic2D::CalcMirrorDelay() {
  TCplxLong Buf[1000];
  TCplxLong Buf1[1000];

  long double Energy0=0,Energy1=0,Energy2=0;
  for(int i=1;i<=InitData.nStepsX;i++)
    Energy0+=hMirrorDelayBuf[Index((InitData.DelayPos+1)%InitData.nStepsDelay,i)].real()*
              hMirrorDelayBuf[Index((InitData.DelayPos+1)%InitData.nStepsDelay,i)].real()+
              hMirrorDelayBuf[Index((InitData.DelayPos+1)%InitData.nStepsDelay,i)].imag()*
              hMirrorDelayBuf[Index((InitData.DelayPos+1)%InitData.nStepsDelay,i)].imag();
  for(int i=1;i<=InitData.nStepsX;i++) {
    Buf[i]=0;
    for(int j=1;j<=InitData.nStepsX;j++) {
      Buf[i]+=hMirrorDelayBuf[Index((InitData.DelayPos+1)%InitData.nStepsDelay,j)]*
             polar(1,InitData.DiffractionCoef*InitData.dt*InitData.dt*(i-j)*(i-j));
    }
//    Buf[i]*=InitData.DiffractionCoef/((long double)(M_PI*InitData.nStepsX));
  }
  for(int i=1;i<=InitData.nStepsX;i++) {
    hMirrorDelayBuf[Index((InitData.DelayPos+1)%InitData.nStepsDelay,i)]=0;
    for(int j=1;j<=InitData.nStepsX;j++) {
      hMirrorDelayBuf[Index((InitData.DelayPos+1)%InitData.nStepsDelay,i)]+=
                    Buf[j]*polar(1,InitData.DiffractionCoef*InitData.dt*InitData.dt*(i-j)*(i-j));
    }
    hMirrorDelayBuf[Index((InitData.DelayPos+1)%InitData.nStepsDelay,i)]*=
              InitData.dt*InitData.dt*InitData.DiffractionCoef/((long double)(M_PI));
    Energy2+=hMirrorDelayBuf[Index((InitData.DelayPos+1)%InitData.nStepsDelay,i)].real()*
            hMirrorDelayBuf[Index((InitData.DelayPos+1)%InitData.nStepsDelay,i)].real()+
            hMirrorDelayBuf[Index((InitData.DelayPos+1)%InitData.nStepsDelay,i)].imag()*
            hMirrorDelayBuf[Index((InitData.DelayPos+1)%InitData.nStepsDelay,i)].imag();
  }
  if(Energy0>0.00000000000000000001)
    DiffEnergy=Energy2/Energy0;
  EDiffOut+=(Energy0-Energy2)*InitData.dt*InitData.dt;
}

void cEquOptic2D::FillBorders() {
  for(int i=0;i<InitData.SizeX;i++) {
    hMirrorDelayBuf[Index(InitData.DelayPos,i)]=hMedia[Index(1,i)].AZp;
  }
  InitData.DelayPos=(InitData.DelayPos+1)%InitData.nStepsDelay;
  if(InitData.UseMirrorDelay) CalcMirrorDelay();
  for(int i=0;i<InitData.SizeX;i++) {
    hMedia[Index(0,i)]=hMedia[Index(1,i)];
    hMedia[Index(0,i)].AZm=BorderAZm(i);
    hMedia[Index(InitData.SizeZ-1,i)]=hMedia[Index(InitData.SizeZ-2,i)];
    hMedia[Index(InitData.SizeZ-1,i)].AZp=BorderAZp(i);
  }
  for(int i=0;i<InitData.SizeZ;i++) {
    hMedia[Index(i,0)]=hMedia[Index(i,1)];
    hMedia[Index(i,0)].AXm=BorderAXm(i);
    hMedia[Index(i,InitData.SizeX-1)]=hMedia[Index(i,InitData.SizeX-2)];
    hMedia[Index(i,InitData.SizeX-1)].AXp=BorderAXp(i);
  }       /*
  for(int i=0;i<InitData.SizeZ;i++) {
    hMedia[i*InitData.SizeX]=hMedia[i*InitData.SizeX+1];
    hMedia[i].AXp=BorderAXp(i);
    hMedia[i*InitData.SizeX+InitData.SizeX-1]=hMedia[i*InitData.SizeX+InitData.SizeX-2];
    hMedia[i*InitData.SizeX+InitData.SizeX-1].AXm=BorderAXm(i);
  }         */
}
/*
void cEquOptic2D::StepRoutine() {
  for(int i=1;i<=InitData.nStepsZ;i++)
    for(int j=1;j<=InitData.nStepsX;j++)
      hMedia[i*InitData.SizeX+j].MakeK(NULL,hMediaK+i*InitData.SizeX+j,0);
  for(int i=1;i<=InitData.nStepsZ;i++)
    for(int j=1;j<=InitData.nStepsX;j++)
      hMedia[i*InitData.SizeX+j].MakeK(hMediaK+i*InitData.SizeX+j,
             hMediaK+InitData.MediaSize+i*InitData.SizeX+j,0.5);
  for(int i=1;i<=InitData.nStepsZ;i++)
    for(int j=1;j<=InitData.nStepsX;j++)
      hMedia[i*InitData.SizeX+j].MakeK(hMediaK+InitData.MediaSize+i*InitData.SizeX+j,
             hMediaK+2*InitData.MediaSize+i*InitData.SizeX+j,0.5);
  for(int i=1;i<=InitData.nStepsZ;i++)
    for(int j=1;j<=InitData.nStepsX;j++)
      hMedia[i*InitData.SizeX+j].MakeK(hMediaK+2*InitData.MediaSize+i*InitData.SizeX+j,
             hMediaK+3*InitData.MediaSize+i*InitData.SizeX+j,1.0);
  for(int i=1;i<=InitData.nStepsZ;i++)
    for(int j=1;j<=InitData.nStepsX;j++) {
      hMediaNew[i*InitData.SizeX+j]=hMedia[i*InitData.SizeX+j]+
      (hMediaK[i*InitData.SizeX+j]+hMediaK[InitData.MediaSize+i*InitData.SizeX+j]*2.0+
       hMediaK[2*InitData.MediaSize+i*InitData.SizeX+j]*2.0+
       hMediaK[3*InitData.MediaSize+i*InitData.SizeX+j])*(1/6.0);

      hMediaNew[i*InitData.SizeX+j].AZp=hMedia[(i-1)*InitData.SizeX+j].AZp+
      (hMediaK[(i-1)*InitData.SizeX+j].AZp+
       hMediaK[InitData.MediaSize+(i-1)*InitData.SizeX+j].AZp*(long double)2.0+
       hMediaK[2*InitData.MediaSize+(i-1)*InitData.SizeX+j].AZp*(long double)2.0+
       hMediaK[3*InitData.MediaSize+(i-1)*InitData.SizeX+j].AZp)*(long double)(1/6.0);

      hMediaNew[i*InitData.SizeX+j].AZm=hMedia[(i+1)*InitData.SizeX+j].AZm+
      (hMediaK[(i+1)*InitData.SizeX+j].AZm+
       hMediaK[InitData.MediaSize+(i+1)*InitData.SizeX+j].AZm*(long double)2.0+
       hMediaK[2*InitData.MediaSize+(i+1)*InitData.SizeX+j].AZm*(long double)2.0+
       hMediaK[3*InitData.MediaSize+(i+1)*InitData.SizeX+j].AZm)*(long double)(1/6.0);

      hMediaNew[i*InitData.SizeX+j].AXp=hMedia[i*InitData.SizeX+j-1].AXp+
      (hMediaK[i*InitData.SizeX+j-1].AXp+
       hMediaK[InitData.MediaSize+i*InitData.SizeX+j-1].AXp*(long double)2.0+
       hMediaK[2*InitData.MediaSize+i*InitData.SizeX+j-1].AXp*(long double)2.0+
       hMediaK[3*InitData.MediaSize+i*InitData.SizeX+j-1].AXp)*(long double)(1/6.0);

      hMediaNew[i*InitData.SizeX+j].AXm=hMedia[i*InitData.SizeX+j+1].AXm+
      (hMediaK[i*InitData.SizeX+j+1].AXm+
       hMediaK[InitData.MediaSize+i*InitData.SizeX+j+1].AXm*(long double)2.0+
       hMediaK[2*InitData.MediaSize+i*InitData.SizeX+j+1].AXm*(long double)2.0+
       hMediaK[3*InitData.MediaSize+i*InitData.SizeX+j+1].AXp)*(long double)(1/6.0);
    }
  memcpy(hMedia,hMediaNew,InitData.MediaSize*sizeof(cOptic2D_Media));
  FillBorders();

  for(int i=0;i<FFTSizeFull-1;i++) hFFTSrc[i]=hFFTSrc[i+1];
  hFFTSrc[FFTSizeFull-1]=hMedia[InitData.SizeZ*InitData.SizeX/3+InitData.SizeX/3].AZp;
}
*/
/*
void cEquOptic2D::StepRoutine() {
  for(int i=1;i<=InitData.nStepsZ;i++)
    for(int j=1;j<=InitData.nStepsX;j++)
      hMedia[i*InitData.SizeX+j]=hMedia[i*InitData.SizeX+j]+hMedia[i*InitData.SizeX+j].RightPart()*InitData.dt;
  for(int i=1;i<=InitData.nStepsZ;i++)
    for(int j=1;j<=InitData.nStepsX;j++) {
      hMedia[i*InitData.SizeX+j].AXm=hMedia[i*InitData.SizeX+j+1].AXm;
      hMedia[i*InitData.SizeX+j].AZm=hMedia[(i+1)*InitData.SizeX+j].AZm;
    }
  for(int i=InitData.nStepsZ;i>=1;i--)
    for(int j=InitData.nStepsX;j>=1;j--) {
      hMedia[i*InitData.SizeX+j].AXp=hMedia[i*InitData.SizeX+j-1].AXp;
      hMedia[i*InitData.SizeX+j].AZp=hMedia[(i-1)*InitData.SizeX+j].AZp;
    }
  FillBorders();

  for(int i=0;i<FFTSizeFull-1;i++) hFFTSrc[i]=hFFTSrc[i+1];
  hFFTSrc[FFTSizeFull-1]=hMedia[InitData.SizeZ*InitData.SizeX/2+InitData.SizeX/2].AZp;
}
*/

void cEquOptic2D::StepRoutine() {
  for(int i=1;i<=InitData.nStepsZ;i++)
    for(int j=1;j<=InitData.nStepsX;j++)
      hMedia[Index(i,j)].MakeNewSeveralSteps(hMediaNew+Index(i,j),i,j);
//      hMediaK[i*InitData.SizeX+j]=hMedia[i*InitData.SizeX+j].RightPart()*InitData.dt;
/*
  for(int i=1;i<=InitData.nStepsZ;i++)
    for(int j=1;j<=InitData.nStepsX;j++)
      hMediaK[InitData.MediaSize+i*InitData.SizeX+j]=
       (hMedia[i*InitData.SizeX+j]+hMediaK[i*InitData.SizeX+j]*0.5).RightPart()*InitData.dt;
  for(int i=1;i<=InitData.nStepsZ;i++)
    for(int j=1;j<=InitData.nStepsX;j++)
      hMediaK[2*InitData.MediaSize+i*InitData.SizeX+j]=
       (hMedia[i*InitData.SizeX+j]+hMediaK[InitData.MediaSize+i*InitData.SizeX+j]*0.5).RightPart()*InitData.dt;
  for(int i=1;i<=InitData.nStepsZ;i++)
    for(int j=1;j<=InitData.nStepsX;j++)
      hMediaK[3*InitData.MediaSize+i*InitData.SizeX+j]=
       (hMedia[i*InitData.SizeX+j]+hMediaK[2*InitData.MediaSize+i*InitData.SizeX+j]).RightPart()*InitData.dt;
*/

  for(int i=1;i<=InitData.nStepsZ;i++)
    for(int j=1;j<=InitData.nStepsX;j++)
      hMedia[Index(i,j)]=hMediaNew[Index(i,j)];
  CalcOutEnergy();
  FillBorders();
  for(int i=0;i<=InitData.nStepsZ;i++)
    for(int j=0;j<=InitData.nStepsX;j++) {
      hMedia[Index(i,j)].AXp=hMedia[Index(i,j+1)].AXp;
      hMedia[Index(i,j)].AZp=hMedia[Index(i+1,j)].AZp;
    }
  for(int i=InitData.nStepsZ+1;i>=1;i--)
    for(int j=InitData.nStepsX+1;j>=1;j--) {
      hMedia[Index(i,j)].AXm=hMedia[Index(i,j-1)].AXm;
      hMedia[Index(i,j)].AZm=hMedia[Index(i-1,j)].AZm;
    }
  for(int i=0;i<FFTSizeFull-1;i++) hFFTSrc[i]=hFFTSrc[i+1];
//  hFFTSrc[FFTSizeFull-1]=hMedia[InitData.SizeZ*InitData.SizeX/2+InitData.SizeX/2].AZp;
  hFFTSrc[FFTSizeFull-1]=hMedia[(InitData.SizeZ-1)*InitData.SizeX+InitData.SizeX/2].AZm;
}

void cEquOptic2D::CalcOutEnergy(){
  long double tmpEnergyXP=0,tmpEnergyXM=0,tmpEnergyZP=0,tmpEnergyZM=0;
  for(int i=InitData.nStepsZ;i>=1;i--) {
    tmpEnergyXM+=abs(hMedia[Index(i,InitData.nStepsX)].AXm)*abs(hMedia[Index(i,InitData.nStepsX)].AXm);
    tmpEnergyXP+=abs(hMedia[Index(i,1)].AXp)*abs(hMedia[Index(i,1)].AXp);
  }
  for(int j=InitData.nStepsX;j>=1;j--) {
    tmpEnergyZM+=abs(hMedia[Index(InitData.nStepsZ,j)].AZm)*abs(hMedia[Index(InitData.nStepsZ,j)].AZm);
    tmpEnergyZP+=abs(hMedia[Index(1,j)].AZp)*abs(hMedia[Index(1,j)].AZp)
                 *(1-InitData.Z0RefCoef)*(1-InitData.Z0RefCoef);
  }

  /*
  for(int i=InitData.nStepsZ;i>=1;i--) {
    tmpEnergyXP+=abs(hMedia[i*InitData.SizeX+InitData.nStepsX+1].AXp)*abs(hMedia[i*InitData.SizeX+InitData.nStepsX+1].AXp);
    tmpEnergyXM+=abs(hMedia[i*InitData.SizeX].AXm)*abs(hMedia[i*InitData.SizeX].AXm);
  }
  for(int j=InitData.nStepsX;j>=1;j--) {
    tmpEnergyZP+=abs(hMedia[(InitData.nStepsZ+1)*InitData.SizeX+j].AZp)*abs(hMedia[(InitData.nStepsZ+1)*InitData.SizeX+j].AZp);
    tmpEnergyZM+=abs(hMedia[j].AZm)*abs(hMedia[j].AZm);
  } */
  EnergyXPOut+=InitData.dt*InitData.dt*tmpEnergyXP;
  EnergyXMOut+=InitData.dt*InitData.dt*tmpEnergyXM;
  EnergyZPOut+=InitData.dt*InitData.dt*tmpEnergyZP;
  EnergyZMOut+=InitData.dt*InitData.dt*tmpEnergyZM;
}

void cEquOptic2D::CalcEnergy(){
  EnergyXP=EnergyXM=EnergyZP=EnergyZM=EnergyR0=EnergyRZ=EnergyRX=EnergyRp=EnergyRm=EDelayBuf=0;
  for(int i=InitData.nStepsZ;i>=1;i--)
    for(int j=InitData.nStepsX;j>=1;j--) {
      EnergyXM+=abs(hMedia[i*InitData.SizeX+j].AXm)*abs(hMedia[i*InitData.SizeX+j].AXm);
      EnergyXP+=abs(hMedia[i*InitData.SizeX+j].AXp)*abs(hMedia[i*InitData.SizeX+j].AXp);
      EnergyZM+=abs(hMedia[i*InitData.SizeX+j].AZm)*abs(hMedia[i*InitData.SizeX+j].AZm);
      EnergyZP+=abs(hMedia[i*InitData.SizeX+j].AZp)*abs(hMedia[i*InitData.SizeX+j].AZp);
      EnergyR0+=hMedia[i*InitData.SizeX+j].R0;
    }
  for(int i=0;i<InitData.nStepsDelay;i++)
    if(i!=(InitData.DelayPos+1)%InitData.nStepsDelay)
      for(int j=InitData.nStepsX;j>=1;j--)
        EDelayBuf+=abs(hMirrorDelayBuf[i*InitData.SizeX+j])*abs(hMirrorDelayBuf[i*InitData.SizeX+j]);
  EDelayBuf*=InitData.dt*InitData.dt;
  EnergyXP*=InitData.dt*InitData.dt;
  EnergyXM*=InitData.dt*InitData.dt;
  EnergyZP*=InitData.dt*InitData.dt;
  EnergyZM*=InitData.dt*InitData.dt;
  EnergyR0*=InitData.dt*InitData.dt*2;
  EnergyRZ*=InitData.dt*InitData.dt*2;
  EnergyRX*=InitData.dt*InitData.dt*2;
  EnergyRp*=InitData.dt*InitData.dt*2;
  EnergyRm*=InitData.dt*InitData.dt*2;
 /* if(EnergyZPlast>0.000000000000000001)
    decZP=logl(sqrtl(EnergyZP/EnergyZPlast))/InitData.dt;
  else decZP=0;*/
  if (t-LastDecTime>0.00000000000000001) {
    if( (ESumLast>0.000000000000000001)&&((EnergyZP+EnergyZM+EnergyXP+EnergyXM+EDelayBuf>0.000000000000000001)))
      decZP=0.5*logl((EnergyZP+EnergyZM+EnergyXP+EnergyXM+EDelayBuf)/ESumLast)/((t-LastDecTime));
    else decZP=0;
    LastDecTime=t;
    EnergyZPlast=EnergyZP;
    ESumLast=EnergyZP+EnergyZM+EnergyXP+EnergyXM+EDelayBuf;
  }
}

void cEquOptic2D::SaveAData(){
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
  ofstream FileOutEt = ofstream("Et.dat",ios::out);
  ofstream FileOutSpectrum = ofstream("Spectrum.dat",ios::out);
  for(int i=0;i<InitData.nStepsZ;i++)
    for(int j=0;j<InitData.nStepsX;j++) {
      FileOutAZm << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hMedia[Index(i+1,j+1)].AZm) << endl;
      FileOutAXm << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hMedia[Index(i+1,j+1)].AXm) << endl;
      FileOutAZp << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hMedia[Index(i+1,j+1)].AZp) << endl;
      FileOutAXp << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hMedia[Index(i+1,j+1)].AXp) << endl;
      FileOutRX << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hMedia[Index(i+1,j+1)].RX) << endl;
      FileOutRZ << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hMedia[Index(i+1,j+1)].RZ) << endl;
      FileOutRp << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hMedia[Index(i+1,j+1)].Rp) << endl;
      FileOutRm << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hMedia[Index(i+1,j+1)].Rm) << endl;
      FileOutR0 << InitData.dt*i << " " << InitData.dt*j
                 << " " << hMedia[Index(i+1,j+1)].R0 << endl;
    }
   for(int i=0;i<nReportNumber;i++) {
     FileOutEt << i*t/nReportNumber << " " << hEnergyData[i] << endl;
   }
   /*
   for(int j=0;j<InitData.nStepsX;j++) {
      FileOutAZp_ph << InitData.dt*j
                 << " " << arg(hMedia[Index(1,j+1)].AZp) << endl;
      FileOutAZm_ph << InitData.dt*j
                 << " " << arg(hMedia[Index(1,j+1)].AZm) << endl;
   }
   */
  FFT.TransForward(hFFTDest,hFFTSrc);
  TCplxLong tmp;
  for(int i=0;i<(FFTSizeFull >> 1);i++) {
    tmp=hFFTDest[i];
    hFFTDest[i]=hFFTDest[i+(FFTSizeFull >> 1)];
    hFFTDest[i+(FFTSizeFull >> 1)]=tmp;
  }
  for(int i=0;i<FFTSizeFull;i++) {
    tmp1=(i-(FFTSizeFull>>1))*2.0*M_PI/FFTSizeFull/InitData.dt;
    FileOutSpectrum << tmp1 << " " << abs(hFFTDest[i]) << endl;
  }
}


