//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquOptic2D_X_01.h"
#include <fstream.h>

//---------------------------------------------------------------------------

#pragma package(smart_init)

void cEquOptic2D_X_base::Prepare2DColor(int Number2D) {
  switch(Number2D) {
    case 0:
      maxAZp=0.000000000000001;
      for(int i=0;i<=InitData.nStepsZ+1;i++)
        for(int j=0;j<=InitData.nStepsX+1;j++) {
          hAbsStore[IndexXZ(j,i)].AZp=abs(hAZp[IndexXZ(j,i)]);
          if(hAbsStore[IndexXZ(j,i)].AZp>maxAZp) maxAZp=hAbsStore[IndexXZ(j,i)].AZp;
        }
    break;
    case 1:
      maxAZm=0.000000000000001;
      for(int i=0;i<=InitData.nStepsZ+1;i++)
        for(int j=0;j<=InitData.nStepsX+1;j++) {
          hAbsStore[IndexXZ(j,i)].AZm=abs(hAZm[IndexXZ(j,i)]);
          if(hAbsStore[IndexXZ(j,i)].AZm>maxAZm) maxAZm=hAbsStore[IndexXZ(j,i)].AZm;
        }
    break;
    case 2:
      maxAXp=0.000000000000001;
      for(int i=0;i<=InitData.nStepsZ+1;i++)
        for(int j=0;j<=InitData.nStepsX+1;j++) {
          hAbsStore[IndexXZ(j,i)].AXp=abs(hAXp[IndexXZ(j,i)]);
          if(hAbsStore[IndexXZ(j,i)].AXp>maxAXp) maxAXp=hAbsStore[IndexXZ(j,i)].AXp;
        }
    break;
    case 3:
      maxAXm=0.000000000000001;
      for(int i=0;i<=InitData.nStepsZ+1;i++)
        for(int j=0;j<=InitData.nStepsX+1;j++) {
          hAbsStore[IndexXZ(j,i)].AXm=abs(hAXm[IndexXZ(j,i)]);
          if(hAbsStore[IndexXZ(j,i)].AXm>maxAXm) maxAXm=hAbsStore[IndexXZ(j,i)].AXm;
        }
    break;
    case 4:
      maxPZp=0.000000000000001;
      for(int i=0;i<=InitData.nStepsZ+1;i++)
        for(int j=0;j<=InitData.nStepsX+1;j++) {
          hAbsStore[IndexXZ(j,i)].PZp=abs(hM[IndexXZ(j,i)].PZp);
          if(hAbsStore[IndexXZ(j,i)].PZp>maxPZp) maxPZp=hAbsStore[IndexXZ(j,i)].PZp;
        }
    case 5:
      maxPZm=0.000000000000001;
      for(int i=0;i<=InitData.nStepsZ+1;i++)
        for(int j=0;j<=InitData.nStepsX+1;j++) {
          hAbsStore[IndexXZ(j,i)].PZm=abs(hM[IndexXZ(j,i)].PZm);
          if(hAbsStore[IndexXZ(j,i)].PZm>maxPZm) maxPZm=hAbsStore[IndexXZ(j,i)].PZm;
        }
    break;
    case 6:
      maxPXp=0.000000000000001;
      for(int i=0;i<=InitData.nStepsZ+1;i++)
        for(int j=0;j<=InitData.nStepsX+1;j++) {
          hAbsStore[IndexXZ(j,i)].PXp=abs(hM[IndexXZ(j,i)].PXp);
          if(hAbsStore[IndexXZ(j,i)].PXp>maxPXp) maxPXp=hAbsStore[IndexXZ(j,i)].PXp;
        }
    break;
    case 7:
      maxPXm=0.000000000000001;
      for(int i=0;i<=InitData.nStepsZ+1;i++)
        for(int j=0;j<=InitData.nStepsX+1;j++) {
          hAbsStore[IndexXZ(j,i)].PXm=abs(hM[IndexXZ(j,i)].PXm);
          if(hAbsStore[IndexXZ(j,i)].PXm>maxPXm) maxPXm=hAbsStore[IndexXZ(j,i)].PXm;
        }
    break;
  }
}

int cEquOptic2D_X_base::Get2DColor(int Number2D, double x, double y) {
  int i0=x*(InitData.nStepsZ-1)+1;
  int i1=i0+1;
  int j0=y*(InitData.nStepsX-1)+1;
  int j1=j0+1;
  double c0=x*(InitData.nStepsZ-1)+1-i0;
  double c1=y*(InitData.nStepsX-1)+1-j0;
  /*
  double c00=c0*c1;
  double c01=c0*(1.0-c1);
  double c10=(1.0-c0)*c1;
  double c11=(1.0-c0)*(1.0-c1);
  */
  double c11=c0*c1;
  double c10=c0*(1.0-c1);
  double c01=(1.0-c0)*c1;
  double c00=(1.0-c0)*(1.0-c1);

  double tmpinv;
  /*
  if((x>1)||(x<0))
    x=0;
  if((y>1)||(y<0))
    y=0;*/
  switch(Number2D) {
    case 0:
      tmpinv=1.0/maxAZp;
      return(tmpinv*255*(hAbsStore[Index_ZX(i0,j0)].AZp*c00+hAbsStore[Index_ZX(i0,j0+1)].AZp*c01+
                  hAbsStore[Index_ZX(i0+1,j0)].AZp*c10+hAbsStore[Index_ZX(i0+1,j0+1)].AZp*c11));
//      return(tmpinv*255*(hAbsStore[Index_ZX(i0,j0)].AZp));
    break;
    case 1:
      tmpinv=1.0/maxAZm;
      return(tmpinv*255*(hAbsStore[Index_ZX(i0,j0)].AZm*c00+hAbsStore[Index_ZX(i0,j0+1)].AZm*c01+
                  hAbsStore[Index_ZX(i0+1,j0)].AZm*c10+hAbsStore[Index_ZX(i0+1,j0+1)].AZm*c11));
    break;
    case 2:
      tmpinv=1.0/maxAXp;
      return(tmpinv*255*(hAbsStore[Index_ZX(i0,j0)].AXp*c00+hAbsStore[Index_ZX(i0,j0+1)].AXp*c01+
                  hAbsStore[Index_ZX(i0+1,j0)].AXp*c10+hAbsStore[Index_ZX(i0+1,j0+1)].AXp*c11));
    break;
    case 3:
      tmpinv=1.0/maxAXm;
      return(tmpinv*255*(hAbsStore[Index_ZX(i0,j0)].AXm*c00+hAbsStore[Index_ZX(i0,j0+1)].AXm*c01+
                  hAbsStore[Index_ZX(i0+1,j0)].AXm*c10+hAbsStore[Index_ZX(i0+1,j0+1)].AXm*c11));
    case 4:
      tmpinv=1.0/maxPZp;
      return(tmpinv*255*(hAbsStore[Index_ZX(i0,j0)].PZp*c00+hAbsStore[Index_ZX(i0,j0+1)].PZp*c01+
                  hAbsStore[Index_ZX(i0+1,j0)].PZp*c10+hAbsStore[Index_ZX(i0+1,j0+1)].PZp*c11));
    break;
    case 5:
      tmpinv=1.0/maxPZm;
      return(tmpinv*255*(hAbsStore[Index_ZX(i0,j0)].PZm*c00+hAbsStore[Index_ZX(i0,j0+1)].PZm*c01+
                  hAbsStore[Index_ZX(i0+1,j0)].PZm*c10+hAbsStore[Index_ZX(i0+1,j0+1)].PZm*c11));
    break;
    case 6:
      tmpinv=1.0/maxPXp;
      return(tmpinv*255*(hAbsStore[Index_ZX(i0,j0)].PXp*c00+hAbsStore[Index_ZX(i0,j0+1)].PXp*c01+
                  hAbsStore[Index_ZX(i0+1,j0)].PXp*c10+hAbsStore[Index_ZX(i0+1,j0+1)].PXp*c11));
    break;
    case 7:
      tmpinv=1.0/maxPXm;
      return(tmpinv*255*(hAbsStore[Index_ZX(i0,j0)].PXm*c00+hAbsStore[Index_ZX(i0,j0+1)].PXm*c01+
                  hAbsStore[Index_ZX(i0+1,j0)].PXm*c10+hAbsStore[Index_ZX(i0+1,j0+1)].PXm*c11));
    break;
    default:
      return(255*x*y);
  }
}

void cEquOptic2D_X_base::PrepareColorData() {
  char *hCD=(char *)ColorData;
  char *hCDGray=(char *)ColorDataGray;
  //hCDGray=(char *)ColorData;
  union {
    int c;
    struct {
      char r,g,b,m;
    } rgb;
  };
  for(int color=0;color<256;color++) {
     hCD[4*color+3]=0;                       //m
     int number_of_colors=255;
     double tt = 6.2831853071795864 * color / (1.484375 * number_of_colors);
     double rv = 1 + sin(tt - 1.9634954084936207);
     double gv = 1 + sin(tt - 0.39269908169872415);
     double bv = 1 + sin(tt + 1.1780972450961724);
     tt = 127.5 * sqrt(sqrt((color + (0.078125 * number_of_colors)) / (0.91796875 * number_of_colors)));
     hCD[4*color+0] = int(tt * rv);          //r
     hCD[4*color+1] = int(tt * gv);          //g
     hCD[4*color+2] = int(tt * bv);          //b
     hCDGray[4*color+0]=hCDGray[4*color+1]=hCDGray[4*color+2]=color;
  }
}

void cEquOptic2D_X_base::InitReport(cDifEquReport *hReport) {
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
void cEquOptic2D_X_base::Report(cDifEquReport *hReport) {
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
  hReport->nVals=1;
  OPTIC_TYPE EnergyA=EnergyXP+EnergyXM+EnergyZP+EnergyZM;
  OPTIC_TYPE dtReport=1;
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
  /*
  hReport->ValText[1]="Energy A "; hReport->Val[1]=EnergyA;
  hReport->ValText[2]="Energy R "; hReport->Val[2]=EnergyR;
  hReport->ValText[hReport->nVals]="Energy wrong ";
  OPTIC_TYPE EOut=EnergyXPOut+EnergyXMOut+EnergyZPOut+EnergyZMOut+EnergyUPOut;
  hReport->Val[hReport->nVals++]=(EnergyR0-EnergyR)+EnergyQ+(EnergyA0-EnergyA)-EOut;
 */
 /*
  hReport->ValText[hReport->nVals]="Energy out X ";
  hReport->Val[hReport->nVals++]=EnergyXPOut+EnergyXMOut;
  hReport->ValText[hReport->nVals]="Energy out Z ";
  hReport->Val[hReport->nVals++]=EnergyZPOut+EnergyZMOut;
  */
  hReport->ValText[hReport->nVals]="Dec A "; hReport->Val[hReport->nVals++]=decA;
  hReport->ValText[hReport->nVals]="Dec ZP "; hReport->Val[hReport->nVals++]=decZP;

  /*
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
  hReport->ValText[hReport->nVals]="Freq ";
  hReport->Val[hReport->nVals++]=Frequency;
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
    if(InitData.f2D1D2DSolidMedia&&(hEnergyData[nStepNumber-1].Out>0.0000000001)) {
      hReport->ValText[hReport->nVals]="Power UP/Out "; hReport->Val[hReport->nVals++]=hEnergyData[nStepNumber-1].OutGenZ/hEnergyData[nStepNumber-1].Out;
    }
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

void cEquOptic2D_X_base::FinalReport(cDifEquReport *hReport) {
  OPTIC_TYPE tmpR=0;
  hReport->nVals=0;
  hReport->ValText[hReport->nVals]="Average R ";
  hReport->Val[hReport->nVals++]=EnergyR/InitData.nStepsX/InitData.nStepsZ/InitData.dt/InitData.dt;
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
    tmpR=hEnergyData[nStepNumber-1].OutXp+hEnergyData[nStepNumber-1].OutXm+
         hEnergyData[nStepNumber-1].OutZp+hEnergyData[nStepNumber-1].OutZm;
    hReport->ValText[hReport->nVals]="OutSum ";
    hReport->Val[hReport->nVals++]=tmpR;
    if(InitData.f2D1D2DSolidMedia&&(hEnergyData[nStepNumber-1].Out>0.0000000001)) {
      hReport->ValText[hReport->nVals]="Power UP/Out "; hReport->Val[hReport->nVals++]=hEnergyData[nStepNumber-1].OutGenZ/hEnergyData[nStepNumber-1].Out;
    }
  }
  hReport->ValText[hReport->nVals]="Dec A "; hReport->Val[hReport->nVals++]=decA;
  hReport->ValText[hReport->nVals]="Dec ZP "; hReport->Val[hReport->nVals++]=decZP;

  hReport->nGraphs=2;
  for(int i=0;i<nStepNumber;i++) hReportBuf[i]=hEnergyData[i].Out;
  hReport->Graph[0].GetData(0,nStepNumber,hReportBuf);
  for(int i=0;i<nStepNumber;i++) hReportBuf[i]=hEnergyData[i].R;
  hReport->Graph[1].GetData(0,nStepNumber,hReportBuf);
/*  hReport->Graph[0].GetDataAbs(0,nStepNumber,hOutA);
  hReport->GraphText[0]="Out A";*/
  hReport->ValText[hReport->nVals]="Max P "; hReport->Val[hReport->nVals++]=hReport->Graph[0].Max;
  if(InitData.Coaxial) {
    hReport->ValText[hReport->nVals]="EoutZ/EGen ";
    hReport->Val[hReport->nVals++]=(hEnergyData[nStepNumber-1].OutZp+hEnergyData[nStepNumber-1].OutZm)/(EnergyZP+EnergyZM)/InitData.AGenCoef;
  }
  FinalSave();
}

void cEquOptic2D_X_base::FinalSave() {
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

     FileOutAZMicro << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hAZm[Index_ZX(i+1,j+1)])+
                           abs(hAZp[Index_ZX(i+1,j+1)])
                 << endl;
/*
      FileOutAZMicro << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hAZm[Index_ZX(i+1,j+1)]-hAZp[Index_ZX(i+1,j+1)])
                 << endl;
      FileOutAXMicro << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hAXm[Index_ZX(i+1,j+1)]-hAXp[Index_ZX(i+1,j+1)])
                 << endl;
*/
      /*
      FileOutAZMicro << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hAZm[Index_ZX(i+1,j+1)])+
                           abs(hAZp[Index_ZX(i+1,j+1)])
                 << endl;
      FileOutAXMicro << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hAXm[Index_ZX(i+1,j+1)])+
                           abs(hAXp[Index_ZX(i+1,j+1)])
                 << endl;
      *//*
      FileOutAZMicro << InitData.dt*i << " " << InitData.dt*j
                 << " " << sqrt(abs(hAZm[Index_ZX(i+1,j+1)])*abs(hAZm[Index_ZX(i+1,j+1)])+
                           abs(hAZp[Index_ZX(i+1,j+1)])*abs(hAZp[Index_ZX(i+1,j+1)]))
                 << endl;
      FileOutAXMicro << InitData.dt*i << " " << InitData.dt*j
                 << " " << sqrt(abs(hAXm[Index_ZX(i+1,j+1)])*abs(hAXm[Index_ZX(i+1,j+1)])+
                           abs(hAXp[Index_ZX(i+1,j+1)])*abs(hAXp[Index_ZX(i+1,j+1)]))
                 << endl;
      *//*
      FileOutAZMicro << InitData.dt*i << " " << InitData.dt*j
                 << " " << sqrt(abs(hAZm[Index_ZX(i+1,j+1)])*abs(hAZm[Index_ZX(i+1,j+1)])+
                           abs(hAZp[Index_ZX(i+1,j+1)])*abs(hAZp[Index_ZX(i+1,j+1)])+
                           (hAZp[Index_ZX(i+1,j+1)]*hAZm[Index_ZX(i+1,j+1)]).real()*2)
                 << endl;
      FileOutAXMicro << InitData.dt*i << " " << InitData.dt*j
                 << " " << sqrt(abs(hAXm[Index_ZX(i+1,j+1)])*abs(hAXm[Index_ZX(i+1,j+1)])+
                           abs(hAXp[Index_ZX(i+1,j+1)])*abs(hAXp[Index_ZX(i+1,j+1)])+
                           (hAXp[Index_ZX(i+1,j+1)]*hAXm[Index_ZX(i+1,j+1)]).real()*2)
                 << endl;
      *//*
      FileOutAZArg << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hAZm[Index_ZX(i+1,j+1)])*abs(hAZm[Index_ZX(i+1,j+1)])+
                           abs(hAZp[Index_ZX(i+1,j+1)])*abs(hAZp[Index_ZX(i+1,j+1)])+
                           (hAZp[Index_ZX(i+1,j+1)]*hAZm[Index_ZX(i+1,j+1)]).real()
                 << endl;
      FileOutAXArg << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hAXm[Index_ZX(i+1,j+1)])*abs(hAXm[Index_ZX(i+1,j+1)])+
                           abs(hAXp[Index_ZX(i+1,j+1)])*abs(hAXp[Index_ZX(i+1,j+1)])+
                           (hAXp[Index_ZX(i+1,j+1)]*hAXm[Index_ZX(i+1,j+1)]).real()
                 << endl;
      */

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
                 InitData.beta*InitData.PCoef*hM[Index_ZX(i+1,j+1)].Rm) << endl;
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

void cEquOptic2D_X_base::LoadInitData(cDifEqu_InitDataRec *hInitData) {
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
  if(newSize!=(InitData.nStepsX+2)*(InitData.nStepsZ+2)) {
    delete [] hReportBuf;
    delete [] hAXp; delete [] hAXm; delete [] hAZp; delete [] hAZm;
    delete [] hM; delete [] hAbsStore;
    hAXp = new complex<OPTIC_TYPE>[newSize*3];
    hAXm = new complex<OPTIC_TYPE>[newSize*3];
    hAZp = new complex<OPTIC_TYPE>[newSize*3];
    hAZm = new complex<OPTIC_TYPE>[newSize*3];
    hM   = new cOptic2D_Media_X[newSize*3];
    hReportBuf = new OPTIC_TYPE[MaxSize];
    hAbsStore = new cOptic2D_AbsStore[ newSize ];
  }

  if ( InitData.f2D1D2DSolidMedia )
    for(int i=0;i<hInit->nThreads+1-nEqus;i++)
      hEqus[i] = new cEquOpticTask_2D1D2DSolidMedia;
  else if ( InitData.f2BraggsSolidMedia )
    for(int i=0;i<hInit->nThreads+1-nEqus;i++)
      hEqus[i] = new cEquOpticTask_2BraggsSolidMedia;
  else
    for(int i=0;i<hInit->nThreads+1-nEqus;i++)
      hEqus[i] = new cEquOpticBaseTask;

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

void cEquOptic2D_X_base::PrepareMedia() {
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
      /*
      hM[IndexXZ(ix,iz)].PXp=polar(InitData.rPolNoise,2.0*M_PI*rand()/RAND_MAX);
      hM[IndexXZ(ix,iz)].PXm=polar(InitData.rPolNoise,2.0*M_PI*rand()/RAND_MAX);
      hM[IndexXZ(ix,iz)].PZp=polar(InitData.rPolNoise,2.0*M_PI*rand()/RAND_MAX);
      hM[IndexXZ(ix,iz)].PZm=polar(InitData.rPolNoise,2.0*M_PI*rand()/RAND_MAX);
      */
      hM[IndexXZ(ix,iz)].PXp=0;//polar(InitData.rPolNoise,2.0*M_PI*rand()/RAND_MAX);
      hM[IndexXZ(ix,iz)].PXm=0;//polar(InitData.rPolNoise,2.0*M_PI*rand()/RAND_MAX);
      hM[IndexXZ(ix,iz)].PZp=0;//polar(InitData.rPolNoise,2.0*M_PI*rand()/RAND_MAX);
      hM[IndexXZ(ix,iz)].PZm=0;//polar(InitData.rPolNoise,2.0*M_PI*rand()/RAND_MAX);

      hM[IndexXZ(ix,iz)].RX=hM[IndexXZ(ix,iz)].RZ=0;
      hM[IndexXZ(ix,iz)].Rp=hM[IndexXZ(ix,iz)].Rm=0;
      hM[IndexXZ(ix,iz)].R0=InitData.StartingInv;
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
        hM[IndexXZ(ix,iz)].PXp=tmpC*0.125;
        tmpC=hAXm[IndexXZ(ix,iz)]*OPTIC_TYPE(4)
            +hAXm[IndexXZ(ix-1,iz)]
            +hAXm[IndexXZ(ix,iz-1)]
            +hAXm[IndexXZ(ix+1,iz)]
            +hAXm[IndexXZ(ix,iz+1)];
        hM[IndexXZ(ix,iz)].PXm=tmpC*0.125;

        tmpC=hAZp[IndexXZ(ix,iz)]*OPTIC_TYPE(4)
            +hAZp[IndexXZ(ix-1,iz)]
            +hAZp[IndexXZ(ix,iz-1)]
            +hAZp[IndexXZ(ix+1,iz)]
            +hAZp[IndexXZ(ix,iz+1)];
        hM[IndexXZ(ix,iz)].PZp=tmpC*0.125;
        tmpC=hAZm[IndexXZ(ix,iz)]*OPTIC_TYPE(4)
            +hAZm[IndexXZ(ix-1,iz)]
            +hAZm[IndexXZ(ix,iz-1)]
            +hAZm[IndexXZ(ix+1,iz)]
            +hAZm[IndexXZ(ix,iz+1)];
        hM[IndexXZ(ix,iz)].PZm=tmpC*0.125;
      }
    for(int ix=1;ix<=InitData.nStepsX;ix++)
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        hAXp[IndexXZ(ix,iz)]=hM[IndexXZ(ix,iz)].PXp;
        hAXm[IndexXZ(ix,iz)]=hM[IndexXZ(ix,iz)].PXm;
        hAZp[IndexXZ(ix,iz)]=hM[IndexXZ(ix,iz)].PZp;
        hAZm[IndexXZ(ix,iz)]=hM[IndexXZ(ix,iz)].PZm;
        hM[IndexXZ(ix,iz)].PXp=hM[IndexXZ(ix,iz)].PXm=0;
        hM[IndexXZ(ix,iz)].PZp=hM[IndexXZ(ix,iz)].PZm=0;
      }
  }
}

OPTIC_TYPE cEquOptic2D_X_base::MicrowaveModePulse() {
  if(!InitData.MicrowaveModeLength) return 0;
  if(t>InitData.MicrowaveModeLength) return 0;
  return InitData.A0*sin(M_PI*t/InitData.MicrowaveModeLength)*
                     sin(M_PI*t/InitData.MicrowaveModeLength);
}

void cEquOptic2D_X_base::AddNoise() {
}

void cEquOptic2D_X_base::PurgeNoise() {
  if(InitData.PurgeNoiseDelay>=0) {
    PurgeNoiseTime++;
    if(PurgeNoiseTime>InitData.PurgeNoiseDelay) {
      PurgeNoiseTime=0;
      PurgeRoutine2();
    }
  }
  if (InitData.DiffractionCoef) {
    int ix=InitData.nStepsX*rand()/RAND_MAX+1;
    int iz=InitData.nStepsZ*rand()/RAND_MAX+1;
    double phase=2*M_PI*rand()/RAND_MAX;
    hAXp[IndexXZ(ix,iz)]+=InitData.DiffractionCoef*polar(InitData.DiffractionCoef,phase);
  }
}

void cEquOptic2D_X_base::PurgeRoutine1() {
  int Size=InitData.SizeX*InitData.SizeZ;
  for(int ix=2;ix<=InitData.nStepsX-1;ix++) {
    hAXp[Size+IndexXZ(ix,1)]=(hAXp[IndexXZ(ix,1)]*OPTIC_TYPE(2)+
               hAXp[IndexXZ(ix-1,1)]+hAXp[IndexXZ(ix+1,1)])*OPTIC_TYPE(0.25);
    hAXm[Size+IndexXZ(ix,InitData.nStepsZ)]=(hAXm[IndexXZ(ix,InitData.nStepsZ)]*OPTIC_TYPE(2)+
               hAXm[IndexXZ(ix-1,InitData.nStepsZ)]+
               hAXm[IndexXZ(ix+1,InitData.nStepsZ)])*OPTIC_TYPE(0.25);
  }
  for(int ix=2;ix<=InitData.nStepsX-1;ix++) {
    hAXp[IndexXZ(ix,1)]=hAXp[Size+IndexXZ(ix,1)];
    hAXm[IndexXZ(ix,InitData.nStepsZ)]=hAXm[Size+IndexXZ(ix,InitData.nStepsZ)];
  }
  for(int iz=2;iz<=InitData.nStepsZ-1;iz++) {
    hAZp[Size+IndexXZ(1,iz)]=(hAZp[IndexXZ(1,iz)]*OPTIC_TYPE(2)+
           hAZp[IndexXZ(1,iz-1)]+hAZp[IndexXZ(1,iz+1)])*OPTIC_TYPE(0.25);
    hAZm[Size+IndexXZ(InitData.nStepsX,iz)]=
           (hAXm[IndexXZ(InitData.nStepsX,iz)]*OPTIC_TYPE(2)+
           hAXm[IndexXZ(InitData.nStepsX,iz-1)]+
           hAXm[IndexXZ(InitData.nStepsX,iz+1)])*OPTIC_TYPE(0.25);
  }
  for(int iz=2;iz<=InitData.nStepsZ-1;iz++) {
    hAZp[IndexXZ(1,iz)]=hAZp[Size+IndexXZ(1,iz)];
    hAZm[IndexXZ(InitData.nStepsX,iz)]=
              hAZm[Size+IndexXZ(InitData.nStepsX,iz)];
  }
}

void cEquOptic2D_X_base::PurgeRoutine2() {
  int Size=InitData.SizeX*InitData.SizeZ;
  complex<OPTIC_TYPE> tmpC;
  for(int ix=2;ix<=InitData.nStepsX-1;ix++)
    for(int iz=2;iz<=InitData.nStepsZ-1;iz++) {
    tmpC=hAXp[IndexXZ(ix,iz)]*OPTIC_TYPE(4)
         +hAXp[IndexXZ(ix-1,iz)]
         +hAXp[IndexXZ(ix,iz-1)]
         +hAXp[IndexXZ(ix+1,iz)]
         +hAXp[IndexXZ(ix,iz+1)];
    hAXp[Size+IndexXZ(ix,iz)]=tmpC*0.125;
    tmpC=hAXm[IndexXZ(ix,iz)]*OPTIC_TYPE(4)
         +hAXm[IndexXZ(ix-1,iz)]
         +hAXm[IndexXZ(ix,iz-1)]
         +hAXm[IndexXZ(ix+1,iz)]
         +hAXm[IndexXZ(ix,iz+1)];
    hAXm[Size+IndexXZ(ix,iz)]=tmpC*0.125;

    tmpC=hAZp[IndexXZ(ix,iz)]*OPTIC_TYPE(4)
         +hAZp[IndexXZ(ix-1,iz)]
         +hAZp[IndexXZ(ix,iz-1)]
         +hAZp[IndexXZ(ix+1,iz)]
         +hAZp[IndexXZ(ix,iz+1)];
    hAZp[Size+IndexXZ(ix,iz)]=tmpC*0.125;
    tmpC=hAZm[IndexXZ(ix,iz)]*OPTIC_TYPE(4)
         +hAZm[IndexXZ(ix-1,iz)]
         +hAZm[IndexXZ(ix,iz-1)]
         +hAZm[IndexXZ(ix+1,iz)]
         +hAZm[IndexXZ(ix,iz+1)];
    hAZm[Size+IndexXZ(ix,iz)]=tmpC*0.125;
  }
  for(int ix=2;ix<=InitData.nStepsX-1;ix++) {
    for(int iz=2;iz<=InitData.nStepsZ-1;iz++) {
      hAXp[IndexXZ(ix,iz)]=hAXp[Size+IndexXZ(ix,iz)];
      hAXm[IndexXZ(ix,iz)]=hAXm[Size+IndexXZ(ix,iz)];
      hAZp[IndexXZ(ix,iz)]=hAZp[Size+IndexXZ(ix,iz)];
      hAZm[IndexXZ(ix,iz)]=hAZm[Size+IndexXZ(ix,iz)];
    }
  }
}

void cEquOptic2D_X_base::ForceLinear() {
  if(InitData.fForceLinear)
    for(int iz=0; iz<=InitData.nStepsZ+1; iz++)
      for(int ix=0; ix<=InitData.nStepsX+1; ix++) {
        hM[IndexXZ(ix,iz)].R0=InitData.Q;
        hM[IndexXZ(ix,iz)].Rm=InitData.QGrid;
        hM[IndexXZ(ix,iz)].Rp=InitData.QGrid;
        hM[IndexXZ(ix,iz)].RX=hM[IndexXZ(ix,iz)].RZ=0;
      }
}

void cEquOptic2D_X_base::FillBorders(int iLayer) {
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
      hAZm[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=InitData.Z0RefCoef*hAZp[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];;
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
        hM[iLayer*Size+IndexXZ(0,iz)].PXp=0;
        hM[iLayer*Size+IndexXZ(0,iz)].PXm=0;
        hM[iLayer*Size+IndexXZ(0,iz)].PZp=hM[iLayer*Size+IndexXZ(1,iz)].PZp;
        hM[iLayer*Size+IndexXZ(0,iz)].PZm=hM[iLayer*Size+IndexXZ(1,iz)].PZm;
        hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PXp=0;
        hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PXm=0;
        hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PZp=hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PZp;
        hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PZm=hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PZm;
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
          hM[iLayer*Size+IndexXZ(ix,0)].PXp=0;//hM[iLayer*Size+IndexXZ(ix,1)].PXp;
          hM[iLayer*Size+IndexXZ(ix,0)].PXm=0;//hM[iLayer*Size+IndexXZ(ix,1)].PXm;
          hM[iLayer*Size+IndexXZ(ix,0)].PZp=0;
          hM[iLayer*Size+IndexXZ(ix,0)].PZm=0;
          hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PXp=0;//hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)].PXp;
          hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PXm=0;//hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)].PXm;
          hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PZp=0;
          hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PZm=0;
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
          hM[iLayer*Size+IndexXZ(0,iz)].PXp=hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PXp;
          hM[iLayer*Size+IndexXZ(0,iz)].PXm=hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PXm;
          hM[iLayer*Size+IndexXZ(0,iz)].PZp=hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PZp;
          hM[iLayer*Size+IndexXZ(0,iz)].PZm=hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PZm;
          hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PXp=hM[iLayer*Size+IndexXZ(1,iz)].PXp;
          hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PXm=hM[iLayer*Size+IndexXZ(1,iz)].PXm;
          hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PZp=hM[iLayer*Size+IndexXZ(1,iz)].PZp;
          hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PZm=hM[iLayer*Size+IndexXZ(1,iz)].PZm;
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
          hM[iLayer*Size+IndexXZ(ix,0)].PXp=hM[iLayer*Size+IndexXZ(ix,1)].PXp;
          hM[iLayer*Size+IndexXZ(ix,0)].PXm=hM[iLayer*Size+IndexXZ(ix,1)].PXm;
          hM[iLayer*Size+IndexXZ(ix,0)].PZp=0;
          hM[iLayer*Size+IndexXZ(ix,0)].PZm=0;
          hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PXp=hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)].PXp;
          hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PXm=hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)].PXm;
          hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PZp=0;
          hM[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)].PZm=0;
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
          hM[iLayer*Size+IndexXZ(0,iz)].PXp=InitData.Z0RefCoef*phase*hM[iLayer*Size+IndexXZ(1,iz)].PXm;
          hM[iLayer*Size+IndexXZ(0,iz)].PXm=InitData.Z0RefCoef*phase*hM[iLayer*Size+IndexXZ(1,iz)].PXp;
          hM[iLayer*Size+IndexXZ(0,iz)].PZp=InitData.Z0RefCoef*phase*hM[iLayer*Size+IndexXZ(1,iz)].PZm;
          hM[iLayer*Size+IndexXZ(0,iz)].PZm=InitData.Z0RefCoef*phase*hM[iLayer*Size+IndexXZ(1,iz)].PZp;
          hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PXp=InitData.Z1RefCoef*phase*hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PXm;
          hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PXm=InitData.Z1RefCoef*phase*hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PXp;
          hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PZp=InitData.Z1RefCoef*phase*hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PZm;
          hM[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)].PZm=InitData.Z1RefCoef*phase*hM[iLayer*Size+IndexXZ(InitData.nStepsX,iz)].PZp;
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

void cEquOptic2D_X_base::StepRoutine() {
/*
  if(InitData.UseFastP)
    StepRoutineEulerNoMedia();
  else
    StepRoutineRunge2NoMedia();
*/
  StepRoutineRunge2NoMedia();

}

void cEquOptic2D_X_base::DoTasks() {
  hTaskCritSection->Acquire();
  ThreadCounter=InitData.nThreads;
  hTaskEvent->ResetEvent();
  for(int i=0;i<InitData.nThreads+1;i++) {
    hEqus[0]->t=t;
  }
  hTaskCritSection->Release();
  for(int i=0;i<InitData.nThreads;i++) {
    hThreads[i]->Resume();
  }
  hEqus[0]->Step();
  if(InitData.nThreads) {
    if (hTaskEvent->WaitFor(THREAD_WAIT_TIME) != wrSignaled) {
      Application->MessageBoxA("Thread wait time is off", "OMG",MB_OK);
    }
  }
}

void cEquOptic2D_X_base::StepRoutineEulerNoMedia() {
  FillBorders(0);
  DoTasks();
  DoTasks();
}

void cEquOptic2D_X_base::StepRoutineRunge2NoMedia() {
  PurgeNoise();
  FillBorders(0);
  ForceLinear();
  DoTasks();
  DoTasks();
  DoTasks();
  CalcOutEnergy();
  //CalcEnergy();
  //hEnergyData[nStepNumber].E=EnergyXP+EnergyXM+EnergyZP+EnergyZM;
  //SaveBitmaps();

  for(int ix=0; ix<InitData.nStepsX; ix++)
    hAZpOutData[nStepNumber*InitData.nStepsX+ix]=hAZp[IndexXZ(ix+1,InitData.nStepsZ)];
  for(int iz=0; iz<InitData.nStepsZ; iz++)
    hAXpOutData[nStepNumber*InitData.nStepsZ+iz]=hAXp[IndexXZ(InitData.nStepsX,iz+1)];

  for(int i=0;i<FFTSizeFull-1;i++) hFFTSrc[i]=hFFTSrc[i+1];
  hFFTSrc[FFTSizeFull-1]=hAZp[InitData.SizeZ*InitData.SizeX/2+InitData.SizeX/2];
//  hFFTSrc[FFTSizeFull-1]=hMedia[(InitData.SizeZ-1)*InitData.SizeX+InitData.SizeX/2].AZm;
}

void cEquOptic2D_X_base::CalcEnergy() {
  int iI;
  EnergyXP=EnergyXM=EnergyZP=EnergyZM=EnergyR=EnergyQ=0;
  for(int ix=1;ix<=InitData.nStepsX;ix++)
    for(int iz=1;iz<=InitData.nStepsZ;iz++) {
      iI=IndexXZ(ix,iz);
      EnergyXP+=hAXp[iI].real()*hAXp[iI].real()+hAXp[iI].imag()*hAXp[iI].imag();
      EnergyXM+=hAXm[iI].real()*hAXm[iI].real()+hAXm[iI].imag()*hAXm[iI].imag();
      EnergyZP+=hAZp[iI].real()*hAZp[iI].real()+hAZp[iI].imag()*hAZp[iI].imag();
      EnergyZM+=hAZm[iI].real()*hAZm[iI].real()+hAZm[iI].imag()*hAZm[iI].imag();
      EnergyR+=hM[iI].R0;
      EnergyQ+=hM[iI].R0Generated;
    }
  EnergyXP*=0.5*InitData.dt*InitData.dt;
  EnergyXM*=0.5*InitData.dt*InitData.dt;
  EnergyZP*=0.5*InitData.dt*InitData.dt;
  EnergyZM*=0.5*InitData.dt*InitData.dt;
  EnergyR*=InitData.dt*InitData.dt;
  EnergyQ*=InitData.dt*InitData.dt;
}

void cEquOptic2D_X_base::CalcOutEnergy() {
  int iI;
  OPTIC_TYPE tmpXP=0,tmpXM=0,tmpZP=0,tmpZM=0,tmpGenZ=0,tmpGenX=0,tmpR=0,tmpGen2D1D2D=0;
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
  if(InitData.f2D1D2DSolidMedia) {
    for(int ix=1;ix<=InitData.nStepsX;ix++) {
      for(int iz=InitData.nStepsZB+1; iz<=InitData.nStepsZ-InitData.nStepsZB2;iz++) {
        tmpGen2D1D2D+=hAZp[IndexXZ(ix,iz)].real()
                *hAZp[IndexXZ(ix,iz)].real()
                +hAZp[IndexXZ(ix,iz)].imag()
                *hAZp[IndexXZ(ix,iz)].imag()
                +hAZm[IndexXZ(ix,iz)].real()
                *hAZm[IndexXZ(ix,iz)].real()
                +hAZm[IndexXZ(ix,iz)].imag()
                *hAZm[IndexXZ(ix,iz)].imag();
      }
    }
  }
  for(int ix=1;ix<=InitData.nStepsX;ix++)
    for(int iz=1;iz<=InitData.nStepsZ;iz++) {
      iI=IndexXZ(ix,iz);
      tmpR+=hM[iI].R0;
    }
  hEnergyData[nStepNumber].R=tmpR*InitData.dt*InitData.dt;
  tmpXP*=0.5*InitData.dt*InitData.dt;  EnergyXPOut+=tmpXP;
  tmpXM*=0.5*InitData.dt*InitData.dt;  EnergyXMOut+=tmpXM;
  tmpZP*=0.5*InitData.dt*InitData.dt;  EnergyZPOut+=tmpZP;
  tmpZM*=0.5*InitData.dt*InitData.dt;  EnergyZMOut+=tmpZM;
  tmpGenZ*=-InitData.dt*InitData.dt*InitData.dt*InitData.AGenZ;
  tmpGenX*=-InitData.dt*InitData.dt*InitData.dt*InitData.AGenX;
  tmpGen2D1D2D*=-InitData.dt*InitData.dt*InitData.dt*InitData.AGenCoef1D;
  EnergyUPOut+=tmpGenZ+tmpGenX+tmpGen2D1D2D;
  hEnergyData[nStepNumber].OutXp=tmpXP/InitData.dt;
  hEnergyData[nStepNumber].OutXm=tmpXM/InitData.dt;
  hEnergyData[nStepNumber].OutZp=tmpZP/InitData.dt;
  hEnergyData[nStepNumber].OutZm=tmpZM/InitData.dt;
  hEnergyData[nStepNumber].OutGenZ=tmpGenZ/InitData.dt+tmpGen2D1D2D/InitData.dt;
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

void cEquOptic2D_X_base::SaveADataMaxE() {
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

  complex<OPTIC_TYPE> tmp;
  OPTIC_TYPE tmpa;
  for(int iz=1;iz<InitData.nStepsZ;iz++) {
    tmp = hAXp[IndexXZ(InitData.nStepsX,iz)];
    tmpa=0;
    if( abs(tmp) > 0.000000001) tmpa = arg(tmp);
    FileOutAXpF << InitData.dt*iz << " " << tmpa << endl;
    
    tmp = hAXm[IndexXZ(1,iz)];
    tmpa=0;
    if( abs(tmp) > 0.000000001) tmpa = arg(tmp);
    FileOutAXmF << InitData.dt*iz << " " << tmpa << endl;
  }
  for(int ix=1;ix<InitData.nStepsX;ix++) {
    tmp = hAZp[IndexXZ(ix,InitData.nStepsZ)];
    tmpa=0;
    if( abs(tmp) > 0.000000001) tmpa = arg(tmp);
    FileOutAZpF << InitData.dt*ix << " " << tmpa << endl;
    
    tmp = hAZp[IndexXZ(ix,InitData.nStepsZ)];
    tmpa=0;
    if( abs(tmp) > 0.000000001) tmpa = arg(tmp);
    FileOutAZmF << InitData.dt*ix << " " << tmpa << endl;
  }
  MaxESaved=1;
}

void cEquOptic2D_X_base::SaveBitmaps() {
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
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      for(int iz=0;iz<InitData.nStepsZ;iz++) {
        hBitmap->Canvas->Pixels[iz][ix]=ColorData[int(255*abs(hAZp[IndexXZ(ix+1,iz+1)])/maxAZp)];
      }
    }

    for(int ix=0;ix<InitData.nStepsX;ix++) {
      for(int iz=0;iz<InitData.nStepsZ;iz++) {
        hBitmap->Canvas->Pixels[iz][InitData.nStepsX+10+ix]=ColorData[int(255*abs(hAXp[IndexXZ(ix+1,iz+1)])/maxAXp)];
      }
    }
    /*
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      for(int iz=0;iz<InitData.nStepsZ;iz++) {
        tmpInt=int(255*abs(hAZp[IndexXZ(ix+1,iz+1)])/normAZp);
        tmpInt = tmpInt > 255 ? 255 : tmpInt;
        hBitmap->Canvas->Pixels[iz][InitData.nStepsX+10+ix]=ColorData[tmpInt];
      }
    }
    */

    hBitmap->Canvas->Pen->Color=clBlack;
    hBitmap->Canvas->MoveTo(0,ComplexBaseLine-InitData.nStepsX);
    hBitmap->Canvas->LineTo(0,ComplexBaseLine);
    hBitmap->Canvas->LineTo(InitData.nStepsZ,ComplexBaseLine);
    int nSt=InitData.tMax/InitData.dt;
    for(int i=0;i<nStepNumber;i++) {
      tmpInt=(hEnergyData[i].OutZp+hEnergyData[i].OutXp)*0.5*InitData.nStepsX/ComplexBitmapNorm;
      hBitmap->Canvas->Pixels[i*InitData.nStepsZ/nSt][ComplexBaseLine-tmpInt]=clGreen;
      tmpInt=(hEnergyData[i].OutZp+hEnergyData[i].OutZp)*0.5*InitData.nStepsX/ComplexBitmapNorm;
      hBitmap->Canvas->Pixels[i*InitData.nStepsZ/nSt][ComplexBaseLine-tmpInt]=clBlue;
    }
    hBitmap->Canvas->Pen->Color=clRed;
    hBitmap->Canvas->MoveTo((nStepNumber-1)*InitData.nStepsZ/nSt,ComplexBaseLine-InitData.nStepsX);
    hBitmap->Canvas->LineTo((nStepNumber-1)*InitData.nStepsZ/nSt,ComplexBaseLine);
    hBitmap->SaveToFile("Bitmaps\\Complex\\"+FileName);
  }

  LastBitmapSaveTime=t;  nBitmap++;
  delete hBitmap;
}

void cEquOptic2D_X_base::PrepareBitmapFolders() {
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

int cEquOptic2D_X_base::TimeToSaveBitmaps() {
  if(!InitData.SaveBitmapInterval) return 0;
  if( (LastBitmapSaveTime!=-1)&&
    (t-LastBitmapSaveTime < InitData.SaveBitmapInterval) ) return 0;
  return 1;
}

