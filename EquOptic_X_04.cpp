//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquOptic_X_04.h"
#include <fstream.h>

//---------------------------------------------------------------------------

#pragma package(smart_init)
void cEquOptic2D_X_Mod2::Prepare2DColor(int Number2D) {
  int s1=(InitData.nStepsZ+2)*(InitData.nStepsX+2);
  int s2=(InitData.nStepsZ+2)*(InitData.nStepsX+2)+(InitData.nStepsZB+2)*(InitData.nStepsX+2);
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
      for(int j=0;j<=InitData.nStepsX+1;j++) {
        for(int i=0;i<=InitData.nStepsZ+1;i++) {
          hAbsStore[IndexXZ(j,i)].PZp=abs(hAZp2[IndexXZ(j,i)]);
          if(hAbsStore[IndexXZ(j,i)].PZp>maxPZp) maxPZp=hAbsStore[IndexXZ(j,i)].PZp;
          hAbsStore[s2+IndexXZ(j,i)].PZp=abs(hAZp1[IndexXZ(j,i)]);
          if(hAbsStore[s2+IndexXZ(j,i)].PZp>maxPZp) maxPZp=hAbsStore[s2+IndexXZ(j,i)].PZp;
        }
        for(int i=0;i<=InitData.nStepsZB+1;i++) {
          hAbsStore[s1+IndexXZ_M(j,i)].PZp=abs(hAZpMedia[IndexXZ_M(j,i)]);
          if(hAbsStore[s1+IndexXZ(j,i)].PZp>maxPZp) maxPZp=hAbsStore[s1+IndexXZ(j,i)].PZp;
        }
      }
    case 5:
      maxPZm=0.000000000000001;
      for(int j=0;j<=InitData.nStepsX+1;j++) {
        for(int i=0;i<=InitData.nStepsZ+1;i++) {
          hAbsStore[IndexXZ(j,i)].PZm=abs(hAZm2[IndexXZ(j,i)]);
          if(hAbsStore[IndexXZ(j,i)].PZm>maxPZm) maxPZm=hAbsStore[IndexXZ(j,i)].PZm;
          hAbsStore[s2+IndexXZ(j,i)].PZm=abs(hAZm1[IndexXZ(j,i)]);
          if(hAbsStore[s2+IndexXZ(j,i)].PZm>maxPZm) maxPZm=hAbsStore[s2+IndexXZ(j,i)].PZm;
        }
        for(int i=0;i<=InitData.nStepsZB+1;i++) {
          hAbsStore[s1+IndexXZ_M(j,i)].PZm=abs(hAZmMedia[IndexXZ_M(j,i)]);
          if(hAbsStore[s1+IndexXZ(j,i)].PZm>maxPZm) maxPZm=hAbsStore[s1+IndexXZ(j,i)].PZm;
        }
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

int cEquOptic2D_X_Mod2::Get2DColor(int Number2D, double x, double y) {
  int nz=InitData.nStepsZ;
  if( (Number2D==4)||(Number2D==5) ) nz=InitData.nStepsZ*2+InitData.nStepsZB+2;
  int i0=x*(nz-1)+1;
  int i1=i0+1;
  int j0=y*(InitData.nStepsX-1)+1;
  int j1=j0+1;
  double c0=x*(nz-1)+1-i0;
  double c1=y*(InitData.nStepsX-1)+1-j0;
  if( (Number2D==4)||(Number2D==5) ) {
    c0=c1; c1=x*(nz-1)+1-i0;
  }
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
      return(tmpinv*255*(hAbsStore[IndexXZ_L(j0,i0)].PZp*c00+hAbsStore[IndexXZ_L(j0,i0+1)].PZp*c01+
                  hAbsStore[IndexXZ_L(j0+1,i0)].PZp*c10+hAbsStore[IndexXZ_L(j0+1,i0+1)].PZp*c11));
    break;
    case 5:
      tmpinv=1.0/maxPZm;
      return(tmpinv*255*(hAbsStore[IndexXZ_L(j0,i0)].PZm*c00+hAbsStore[IndexXZ_L(j0,i0+1)].PZm*c01+
                  hAbsStore[IndexXZ_L(j0+1,i0)].PZm*c10+hAbsStore[IndexXZ_L(j0+1,i0+1)].PZm*c11));
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

int cEquOptic2D_X_Mod2::isBragg(int ix) {
  if(!InitData.nMediaRows) return 1;
  int s2=(InitData.nStepsX-InitData.nMediaRowSize)/InitData.nMediaRows;
  int ss=s2-InitData.nMediaRowSize;
  int sr=(InitData.nStepsX-InitData.nMediaRowSize) % InitData.nMediaRows;
  int dx=ix-1-(s2+1)*sr;
  int res = dx>0 ? (dx % s2 <InitData.nMediaRowSize) : ( (ix-1) % (s2+1) <InitData.nMediaRowSize );
  return res ? 1 : 0;
}

OPTIC_TYPE cEquOptic2D_X_Mod2::ActivePartSize() {
  int c=0; for(int ix=1; ix<=InitData.nStepsX; ix++) if(isBragg(ix)) c++;
  return OPTIC_TYPE(c)/InitData.nStepsX;
}

void cEquOptic2D_X_Mod2::InitReport(cDifEquReport *hReport) {
  int Lines[EQU_MAX_GRAPHS];
  for(int i=0;i<EQU_MAX_GRAPHS;i++) Lines[i]=1;
  Lines[0]=1; Lines[1]=1; Lines[2]=1; Lines[3]=1; Lines[4]=1;
  hReport->LoadText(hReportText);
  hReport->SetGraphSize(EQU_GRAPH_SIZE,Lines);
  hReport->SetMaxTo1();
//  CalcEnergy();
//  nReportNumber++;
//  hEnergyData[nReportNumber-1]=EnergyXP+EnergyXM+EnergyZP+EnergyZM+EDelayBuf;
  LastReportTime=1;
}
void cEquOptic2D_X_Mod2::Report(cDifEquReport *hReport) {
//  nReportNumber++;
  hAXp=hAXp1; hAXm=hAXm1; hAZp=hAZp1; hAZm=hAZm1;
  if(fDebug) {
    hAXp=hAXp2; hAXm=hAXm2; hAZp=hAZp2; hAZm=hAZm2;
  }
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
                                                                                        // --------------------- DEBUG begin
/*
  for(int ix=1; ix<=InitData.nStepsX; ix++)
    hReportBuf[ix-1]=abs(hAZm[IndexXZ(ix,iDebug1)]);
  hReport->Graph[2].GetData(0,InitData.nStepsX,hReportBuf);

  for(int ix=1; ix<=InitData.nStepsX; ix++)
    hReportBuf[ix-1]=abs(hAZm[IndexXZ(ix,iDebug2)]);
  hReport->Graph[3].GetData(0,InitData.nStepsX,hReportBuf);

  if(iDebug5) {
    for(int ix=1; ix<=InitData.nStepsX; ix++)
    hReportBuf[ix-1]=abs(hAZp[IndexXZ(ix,iDebug3)]);
    hReport->Graph[3].GetData(0,InitData.nStepsX,hReportBuf);
  }
*/
                                                                                        // --------------------- DEBUG end

//  hReport->nVals=1;
//  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
  CalcEnergy();
  hReport->nVals=0;
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
//  hReport->ValText[1]="Energy A "; hReport->Val[1]=EnergyA;
  hReport->ValText[hReport->nVals]="Average R ";
  hReport->Val[hReport->nVals++]=EnergyR/InitData.nStepsX/InitData.nStepsZB/InitData.dt/InitData.dt/ActivePartSize();
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
  if( abs(hFFTSrc[FFTSizeFull-2])>0.00000000001 ) {
    hReport->ValText[hReport->nVals]="Freq ";
    Frequency = arg(hFFTSrc[FFTSizeFull-1]/hFFTSrc[FFTSizeFull-2])/InitData.dt;
    hReport->Val[hReport->nVals++]=Frequency;
  }

  int iRx;
  int sR;
  OPTIC_TYPE tmpF;
  if(InitData.nMediaRows) {
    sR =(InitData.nStepsX-InitData.nMediaRowSize)/InitData.nMediaRows;
    for(int iR=0; iR<=InitData.nMediaRows; iR++) {
      iRx = sR*iR+InitData.nMediaRowSize/2;
      if(  (abs(hAZpMedia[IndexXZ_M(iRx,InitData.nStepsZB/2)])>0.00000000001)&&
           (abs(lastAforFreq[iR])>0.00000000001) &&
           ((t-LastReportTime)>0.00000000001) ) {
           tmpF = arg(hAZpMedia[IndexXZ_M(iRx,InitData.nStepsZB/2)]/lastAforFreq[iR])/(t-LastReportTime);
           hReport->ValText[hReport->nVals]="Freq "+AnsiString(iRx)+" ";
           hReport->Val[hReport->nVals++]=tmpF;
           lastAforFreq[iR] = hAZpMedia[IndexXZ_M(iRx,InitData.nStepsZB/2)];
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
  hAXp=hAXp1; hAXm=hAXm1; hAZp=hAZp1; hAZm=hAZm1;
  EnergyZPlast=EnergyZP; ESumLast=EnergyA; LastReportTime=t;
}

void cEquOptic2D_X_Mod2::FinalReport(cDifEquReport *hReport) {
  hReport->nVals=0;
  hReport->ValText[hReport->nVals]="Average R ";
  hReport->Val[hReport->nVals++]=EnergyR/InitData.nStepsX/InitData.nStepsZB/InitData.dt/InitData.dt/ActivePartSize();
  hReport->Val[hReport->nVals++]=Frequency;
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
  hReport->nGraphs=1;

  for(int i=0;i<nStepNumber;i++) hReportBuf[i]=hEnergyData[i].Out;
  hReport->Graph[0].GetData(0,nStepNumber,hReportBuf);
  
/*  hReport->Graph[0].GetDataAbs(0,nStepNumber,hOutA);
  hReport->GraphText[0]="Out A";
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;*/
  FinalSave();
}

void cEquOptic2D_X_Mod2::FinalSave() {
  long double tmp1;
  ofstream FileOutAZp = ofstream("AZp.dat",ios::out);
  //ofstream FileOutAZp_ph = ofstream("AZP_ph.dat",ios::out);
  //ofstream FileOutAZm_ph = ofstream("AZM_ph.dat",ios::out);
  ofstream FileOutAXp = ofstream("AXp.dat",ios::out);
  ofstream FileOutAZm = ofstream("AZm.dat",ios::out);
  ofstream FileOutAXm = ofstream("AXm.dat",ios::out);

  ofstream FileOutAZpL = ofstream("AZpLong.dat",ios::out);
  //ofstream FileOutAXpL = ofstream("AXpLong.dat",ios::out);
  ofstream FileOutAZmL = ofstream("AZmLong.dat",ios::out);
  //ofstream FileOutAXmL = ofstream("AXmLong.dat",ios::out);

  ofstream FileOutAZpt = ofstream("AZpA(t).dat",ios::out);
  ofstream FileOutAZpArgt = ofstream("AZpArg(t).dat",ios::out);

  ofstream FileOutR0 = ofstream("R0.dat",ios::out);
  ofstream FileOutRZ = ofstream("RZ.dat",ios::out);
  ofstream FileOutRX = ofstream("RX.dat",ios::out);
  ofstream FileOutRp = ofstream("Rp.dat",ios::out);
  ofstream FileOutRm = ofstream("Rm.dat",ios::out);
  ofstream FileOutEtOut = ofstream("EtOut.dat",ios::out);
  ofstream FileOutEtOut2= ofstream("EtOut2.dat",ios::out);
  ofstream FileOutAt = ofstream("At.dat",ios::out);
  ofstream FileOutEXp = ofstream("EOutXp.dat",ios::out);
  ofstream FileOutEXm = ofstream("EOutXm.dat",ios::out);
  ofstream FileOutEZp = ofstream("EOutZp.dat",ios::out);
  ofstream FileOutEZm = ofstream("EOutZm.dat",ios::out);
  ofstream FileOutSpectrum = ofstream("Spectrum.dat",ios::out);
  ofstream FileOutAZpFZ = ofstream("AZpFZ.dat",ios::out);
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

      FileOutAZmL << InitData.dt*i << " " << InitData.dt*j
                  << " " << abs(hAZm2[Index_ZX(i+1,j+1)]) << endl;
      FileOutAZpL << InitData.dt*i << " " << InitData.dt*j
                  << " " << abs(hAZp2[Index_ZX(i+1,j+1)]) << endl;

      FileOutAZmL << InitData.dt*(i+InitData.nStepsZ+InitData.nStepsZB) << " " << InitData.dt*j
                  << " " << abs(hAZm1[Index_ZX(i+1,j+1)]) << endl;
      FileOutAZpL << InitData.dt*(i+InitData.nStepsZ+InitData.nStepsZB) << " " << InitData.dt*j
                  << " " << abs(hAZp1[Index_ZX(i+1,j+1)]) << endl;
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

  for(int i=0;i<InitData.nStepsZB;i++)
    for(int j=0;j<InitData.nStepsX;j++) {
      FileOutAZmL << InitData.dt*(i+InitData.nStepsZ) << " " << InitData.dt*j
                  << " " << abs(hAZmMedia[IndexXZ_M(j+1,i+1)]) << endl;
      FileOutAZpL << InitData.dt*(i+InitData.nStepsZ) << " " << InitData.dt*j
                  << " " << abs(hAZpMedia[IndexXZ_M(j+1,i+1)]) << endl;
      FileOutR0 << InitData.dt*(i) << " " << InitData.dt*j
                  << " " << hM2[IndexXZ_M(j+1,i+1)].R0 << endl;
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
    if(InitData.RCoef>0.00000000001)
      FileOutEtOut2 << i*t/nStepNumber << " " << tmpR/InitData.RCoef << endl;
    for(int ix=0; ix<InitData.nStepsX; ix++) {
      FileOutAZpt << InitData.dt*i << " "<< InitData.dt*ix << " " <<
                abs(hAZpOutData[i*InitData.nStepsX+ix]) << endl;
      if( (abs(hAZpOutData[i*InitData.nStepsX+ix])>0.000001)&&(isBragg(ix)) ) {
        FileOutAZpArgt << InitData.dt*i << " "<< InitData.dt*ix << " " <<
                  arg(hAZpOutData[i*InitData.nStepsX+ix])+M_PI << endl;
      } else {
        FileOutAZpArgt << InitData.dt*i << " "<< InitData.dt*ix << " "
                  << 0 <<endl;
      }
    }
  }


  FileOutAZpFZ << InitData.dt << endl;
  FileOutAZpFZ << InitData.nStepsX << endl;
  for(int ix=1; ix<=InitData.nStepsX; ix++) {
    FileOutAZpFZ << hAZp2[IndexXZ(ix,InitData.nStepsZ)].real() << endl;
    FileOutAZpFZ << hAZp2[IndexXZ(ix,InitData.nStepsZ)].imag() << endl;
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

void cEquOptic2D_X_Mod2::LoadInitData(cDifEqu_InitDataRec *hInitData) {
  cDifEqu_InitDataRec_Optic2D_T<OPTIC_TYPE> *hInit=
             (cDifEqu_InitDataRec_Optic2D_T<OPTIC_TYPE> *)hInitData;
  int newSize=(hInit->nStepsX+2)*(hInit->nStepsZ+2);
  int newSizeMedia=(hInit->nStepsX+2)*(hInit->nStepsZB+2);
  int nStepsT=hInit->tMax/hInit->dt+1;
  if(!hTaskCritSection) hTaskCritSection = new TCriticalSection;
  if(!hTaskEvent) hTaskEvent = new TEvent(NULL,true,false,"");
  delete [] hEnergyData; delete [] hAZpOutData; delete [] hReportBuf;
  hEnergyData = new EnergyDataType[hInit->tMax/hInit->dt+1];
  hAZpOutData = new complex<OPTIC_TYPE>[hInit->nStepsX*nStepsT];
  int tmpSize=newSize*2 > nStepsT ? newSize*2 : nStepsT;
  FFTSizeFull = 1 << hInit->FFTSize;
  tmpSize = tmpSize > FFTSizeFull ? tmpSize : FFTSizeFull;
  hReportBuf = new OPTIC_TYPE[tmpSize];
  if(newSize!=(InitData.nStepsX+2)*(InitData.nStepsZ+2)) {
    delete [] hAXp1; delete [] hAXm1; delete [] hAZp1; delete [] hAZm1;
    delete [] hAXp2; delete [] hAXm2; delete [] hAZp2; delete [] hAZm2;
    hAXp1 = new complex<OPTIC_TYPE>[newSize*3];
    hAXm1 = new complex<OPTIC_TYPE>[newSize*3];
    hAZp1 = new complex<OPTIC_TYPE>[newSize*3];
    hAZm1 = new complex<OPTIC_TYPE>[newSize*3];
    hAXp2 = new complex<OPTIC_TYPE>[newSize*3];
    hAXm2 = new complex<OPTIC_TYPE>[newSize*3];
    hAZp2 = new complex<OPTIC_TYPE>[newSize*3];
    hAZm2 = new complex<OPTIC_TYPE>[newSize*3];
  }
  if(newSizeMedia!=(InitData.nStepsX+2)*(InitData.nStepsZB+2)) {
    delete [] hM2;
    delete [] hAZpMedia; delete [] hAZmMedia;
    hM2   = new cOptic2D_Media_X_2Way[newSizeMedia*3];
    hAZpMedia = new complex<OPTIC_TYPE>[newSizeMedia*3];
    hAZmMedia = new complex<OPTIC_TYPE>[newSizeMedia*3];
  }
  delete [] hAbsStore;
  hAbsStore = new cOptic2D_AbsStore[ newSize*2+(hInit->nStepsX+2)*(hInit->nStepsZB+2) ];
  for(int i=0;i<nEqus;i++) delete hEqus[i];
  if(hInit->nMediaRows<0) {
    for(int i=0;i<hInit->nThreads+1-nEqus;i++) {
      hInit->nMediaRows= -hInit->nMediaRows;
      hEqus[i] = new cEquOpticBaseTask_BMB_holes;
      ((cEquOpticBaseTask_BMB_holes *)hEqus[i])->hDeltaArray = DeltaArray;
    }
  } else {
    for(int i=0;i<hInit->nThreads+1-nEqus;i++) {
      hEqus[i] = new cEquOpticBaseTask_BMB;
    }
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
    ((cEquOpticBaseTask_BMB *)hEqus[i])->hAXp1=hAXp1;
    ((cEquOpticBaseTask_BMB *)hEqus[i])->hAXm1=hAXm1;
    ((cEquOpticBaseTask_BMB *)hEqus[i])->hAZp1=hAZp1;
    ((cEquOpticBaseTask_BMB *)hEqus[i])->hAZm1=hAZm1;
    ((cEquOpticBaseTask_BMB *)hEqus[i])->hAXp2=hAXp2;
    ((cEquOpticBaseTask_BMB *)hEqus[i])->hAXm2=hAXm2;
    ((cEquOpticBaseTask_BMB *)hEqus[i])->hAZp2=hAZp2;
    ((cEquOpticBaseTask_BMB *)hEqus[i])->hAZm2=hAZm2;
    ((cEquOpticBaseTask_BMB *)hEqus[i])->hAZpMedia=hAZpMedia;
    ((cEquOpticBaseTask_BMB *)hEqus[i])->hAZmMedia=hAZmMedia;
    ((cEquOpticBaseTask_BMB *)hEqus[i])->hM2=hM2;
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

  for(int i=0; i<128; i++) lastAforFreq[i]=1;

  PrepareMedia();
  FillBorders(0);
  hAXp=hAXp1; hAXm=hAXm1; hAZp=hAZp1; hAZm=hAZm1;
  EnergyXPOut=EnergyXMOut=EnergyZPOut=EnergyZMOut=EnergyQ=0;
  CalcEnergy(); EnergyA0=EnergyXP+EnergyXM+EnergyZP+EnergyZM; EnergyR0=EnergyR;
  ESumLast=EnergyZPlast=0; LastReportTime=0; LastMaxOut=0; EnergyUPOut=0;
  nStepNumber=0; t=0; MaxESaved=0; LastBitmapSaveTime=-1; nBitmap=0;
  iSpectrumX=InitData.nStepsX/3; iSpectrumZ=InitData.nStepsZ/3;
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

void cEquOptic2D_X_Mod2::PrepareMedia() {
  for(int ix=0;ix<=InitData.nStepsX+1;ix++) {
    for(int iz=0;iz<=InitData.nStepsZ+1;iz++) {
      hAXp1[IndexXZ(ix,iz)]=0;
      hAXm1[IndexXZ(ix,iz)]=0;
      hAZp1[IndexXZ(ix,iz)]=0;
      hAZm1[IndexXZ(ix,iz)]=0;
      hAXp2[IndexXZ(ix,iz)]=0;
      hAXm2[IndexXZ(ix,iz)]=0;
      hAZp2[IndexXZ(ix,iz)]=0;
      hAZm2[IndexXZ(ix,iz)]=0;
    }

    for(int iz=0;iz<=InitData.nStepsZB+1;iz++) {
      hAZpMedia[IndexXZ_M(ix,iz)]=sin(M_PI*(iz-1)/(InitData.nStepsZB-1))*
        (InitData.A0+
         InitData.A1*sin(M_PI*(ix-1)/(InitData.nStepsX-1))+
         InitData.A2*sin(2.0*M_PI*(ix-1)/(InitData.nStepsX-1)) )+
        polar(InitData.rPolNoise,2.0*M_PI*rand()/RAND_MAX);
      hAZmMedia[IndexXZ_M(ix,iz)]=0;
      hM2[IndexXZ_M(ix,iz)].PZp=0;
      hM2[IndexXZ_M(ix,iz)].PZm=0;
      hM2[IndexXZ_M(ix,iz)].RZ=0;
      hM2[IndexXZ_M(ix,iz)].R0=InitData.StartingInv;
    }
  }


  for(int i=0; i<InitData.nMediaRows; i++) {
    if(InitData.DeltaCoef<0) {
      DeltaArray[i]=(((OPTIC_TYPE(rand()))/RAND_MAX)-0.5)*InitData.DeltaCoef;
    } else {
      DeltaArray[i] = InitData.DeltaCoef*i+InitData.Delta;
    }
  }
}

void cEquOptic2D_X_Mod2::FillBorders(int iLayer) {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  int SizeMedia=(InitData.nStepsX+2)*(InitData.nStepsZB+2);
  for(int ix=0;ix<=InitData.nStepsX+1;ix++) {
    hAXp1[iLayer*Size+IndexXZ(ix,0)]=0;//hAXp1[iLayer*Size+IndexXZ(ix,1)];
    hAXm1[iLayer*Size+IndexXZ(ix,0)]=0;//hAXm1[iLayer*Size+IndexXZ(ix,1)];
    hAZp1[iLayer*Size+IndexXZ(ix,0)]=hAZpMedia[iLayer*SizeMedia+IndexXZ_M(ix,InitData.nStepsZB)];
    hAZm1[iLayer*Size+IndexXZ(ix,0)]=0;
    hAXp1[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;//hAXp1[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];
    hAXm1[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;//hAXm1[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];
    hAZp1[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;
    hAZm1[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;//hAZp1[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];

    hAXp2[iLayer*Size+IndexXZ(ix,0)]=0;//hAXp2[iLayer*Size+IndexXZ(ix,1)];
    hAXm2[iLayer*Size+IndexXZ(ix,0)]=0;//hAXm2[iLayer*Size+IndexXZ(ix,1)];
    hAZp2[iLayer*Size+IndexXZ(ix,0)]=0;//hAZm2[iLayer*Size+IndexXZ(ix,1)];
    hAZm2[iLayer*Size+IndexXZ(ix,0)]=0;
    hAXp2[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;//hAXp2[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];
    hAXm2[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;//hAXm2[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];
    hAZp2[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;
    hAZm2[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=hAZmMedia[iLayer*SizeMedia+IndexXZ_M(ix,1)];
  }
  for(int iz=0;iz<=InitData.nStepsZ+1;iz++) {
      hAXp1[iLayer*Size+IndexXZ(0,iz)]=0;
      hAXm1[iLayer*Size+IndexXZ(0,iz)]=0;
      hAZp1[iLayer*Size+IndexXZ(0,iz)]=0;//hAZp1[iLayer*Size+IndexXZ(1,iz)];
      hAZm1[iLayer*Size+IndexXZ(0,iz)]=0;//hAZm1[iLayer*Size+IndexXZ(1,iz)];
      hAXp1[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=0;
      hAXm1[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=0;
      hAZp1[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=0;//hAZp1[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
      hAZm1[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=0;//hAZm1[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];

      hAXp2[iLayer*Size+IndexXZ(0,iz)]=0;
      hAXm2[iLayer*Size+IndexXZ(0,iz)]=0;
      hAZp2[iLayer*Size+IndexXZ(0,iz)]=0;//hAZp2[iLayer*Size+IndexXZ(1,iz)];
      hAZm2[iLayer*Size+IndexXZ(0,iz)]=0;//hAZm2[iLayer*Size+IndexXZ(1,iz)];
      hAXp2[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=0;
      hAXm2[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=0;
      hAZp2[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=0;//hAZp2[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
      hAZm2[iLayer*Size+IndexXZ(InitData.nStepsX+1,iz)]=0;//hAZm2[iLayer*Size+IndexXZ(InitData.nStepsX,iz)];
  }

  for(int ix=0;ix<=InitData.nStepsX+1;ix++) {
    hAZpMedia[iLayer*SizeMedia+IndexXZ_M(ix,0)]=hAZp2[iLayer*Size+IndexXZ(ix,InitData.nStepsZ)];
    hAZmMedia[iLayer*SizeMedia+IndexXZ_M(ix,0)]=0;
    hAZpMedia[iLayer*SizeMedia+IndexXZ_M(ix,InitData.nStepsZB+1)]=0;
    hAZmMedia[iLayer*SizeMedia+IndexXZ_M(ix,InitData.nStepsZB+1)]=hAZm1[iLayer*Size+IndexXZ(ix,1)];
    if(InitData.Z0RefCoef) {
      hAZm2[iLayer*Size+IndexXZ(ix,InitData.nStepsZ+1)]=0;
      hAZpMedia[iLayer*SizeMedia+IndexXZ_M(ix,1)]=hAZmMedia[iLayer*SizeMedia+IndexXZ_M(ix,1)]*InitData.Z0RefCoef;
    }
    if(!InitData.fNoMedia) {
      hM2[iLayer*SizeMedia+IndexXZ_M(ix,0)].PZp=hM2[iLayer*SizeMedia+IndexXZ_M(ix,1)].PZm;
      hM2[iLayer*SizeMedia+IndexXZ_M(ix,0)].PZm=0;
      hM2[iLayer*SizeMedia+IndexXZ_M(ix,InitData.nStepsZB+1)].PZp=0;
      hM2[iLayer*SizeMedia+IndexXZ_M(ix,InitData.nStepsZB+1)].PZm=0;
    }
  }
}

void cEquOptic2D_X_Mod2::StepRoutine() {
/*
  if(InitData.UseFastP)
    StepRoutineEulerNoMedia();
  else
    StepRoutineRunge2NoMedia();
*/
  StepRoutineRunge2NoMedia();
}

void cEquOptic2D_X_Mod2::DoTasks() {
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

void cEquOptic2D_X_Mod2::StepRoutineEulerNoMedia() {
  FillBorders(0);
  DoTasks();
  DoTasks();
}

void cEquOptic2D_X_Mod2::StepRoutineRunge2NoMedia() {

  int tPurgeNoiseBack=PurgeNoiseTime;
  hAXp=hAXp1; hAXm=hAXm1; hAZp=hAZp1; hAZm=hAZm1;
  PurgeNoise();
  PurgeNoiseTime=tPurgeNoiseBack;
  hAXp=hAXp2; hAXm=hAXm2; hAZp=hAZp2; hAZm=hAZm2;
  PurgeNoise();
  hAXp=hAXp1; hAXm=hAXm1; hAZp=hAZp1; hAZm=hAZm1;

  FillBorders(0);
  DoTasks();
  DoTasks();
  DoTasks();
  CalcOutEnergy();
  //CalcEnergy();
  //hEnergyData[nStepNumber].E=EnergyXP+EnergyXM+EnergyZP+EnergyZM;
  //SaveBitmaps();
  for(int i=0;i<FFTSizeFull-1;i++) hFFTSrc[i]=hFFTSrc[i+1];
  hFFTSrc[FFTSizeFull-1]=hAXp1[IndexXZ(iSpectrumX,iSpectrumZ)];
//  hFFTSrc[FFTSizeFull-1]=hMedia[(InitData.SizeZ-1)*InitData.SizeX+InitData.SizeX/2].AZm;

  for(int ix=0; ix<InitData.nStepsX; ix++)
    hAZpOutData[nStepNumber*InitData.nStepsX+ix]=hAZp1[IndexXZ(ix+1,InitData.nStepsZ)];
}

void cEquOptic2D_X_Mod2::CalcEnergy() {
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
    if(isBragg(ix))
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

void cEquOptic2D_X_Mod2::CalcOutEnergy() {
  int iI;
  OPTIC_TYPE tmpXP=0,tmpXM=0,tmpZP=0,tmpZM=0,tmpGenZ=0,tmpGenX=0,tmpZ2=0,tmpXP2=0,tmpXM2=0;
  for(int ix=1;ix<=InitData.nStepsX;ix++) {
    tmpZP+= hAZp1[IndexXZ(ix,InitData.nStepsZ)].real()
           *hAZp1[IndexXZ(ix,InitData.nStepsZ)].real()
           +hAZp1[IndexXZ(ix,InitData.nStepsZ)].imag()
           *hAZp1[IndexXZ(ix,InitData.nStepsZ)].imag();
    tmpZM+=  hAZm1[IndexXZ(ix,1)].real()
            *hAZm1[IndexXZ(ix,1)].real()
            +hAZm1[IndexXZ(ix,1)].imag()
            *hAZm1[IndexXZ(ix,1)].imag();
    tmpZ2+=  hAZm2[IndexXZ(ix,1)].real()
            *hAZm2[IndexXZ(ix,1)].real()
            +hAZm2[IndexXZ(ix,1)].imag()
            *hAZm2[IndexXZ(ix,1)].imag();

  }
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    tmpXP+= hAXp1[IndexXZ(InitData.nStepsX,iz)].real()
           *hAXp1[IndexXZ(InitData.nStepsX,iz)].real()
           +hAXp1[IndexXZ(InitData.nStepsX,iz)].imag()
           *hAXp1[IndexXZ(InitData.nStepsX,iz)].imag();
    tmpXM+= hAXm1[IndexXZ(1,iz)].real()
           *hAXm1[IndexXZ(1,iz)].real()
           +hAXm1[IndexXZ(1,iz)].imag()
           *hAXm1[IndexXZ(1,iz)].imag();
    tmpXP2+=hAXp2[IndexXZ(InitData.nStepsX,iz)].real()
           *hAXp2[IndexXZ(InitData.nStepsX,iz)].real()
           +hAXp2[IndexXZ(InitData.nStepsX,iz)].imag()
           *hAXp2[IndexXZ(InitData.nStepsX,iz)].imag();
    tmpXM2+=hAXm2[IndexXZ(1,iz)].real()
           *hAXm2[IndexXZ(1,iz)].real()
           +hAXm2[IndexXZ(1,iz)].imag()
           *hAXm2[IndexXZ(1,iz)].imag();
  }
  if((InitData.AGenZ)||(InitData.AGenX)) {
    for(int ix=1;ix<=InitData.nStepsX;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        tmpGenZ+=hAZp1[IndexXZ(ix,iz)].real()
                *hAZp1[IndexXZ(ix,iz)].real()
                +hAZp1[IndexXZ(ix,iz)].imag()
                *hAZp1[IndexXZ(ix,iz)].imag()
                +hAZm1[IndexXZ(ix,iz)].real()
                *hAZm1[IndexXZ(ix,iz)].real()
                +hAZm1[IndexXZ(ix,iz)].imag()
                *hAZm1[IndexXZ(ix,iz)].imag();
        tmpGenX+=hAXp1[IndexXZ(ix,iz)].real()
                *hAXp1[IndexXZ(ix,iz)].real()
                +hAXp1[IndexXZ(ix,iz)].imag()
                *hAXp1[IndexXZ(ix,iz)].imag()
                +hAXm1[IndexXZ(ix,iz)].real()
                *hAXm1[IndexXZ(ix,iz)].real()
                +hAXm1[IndexXZ(ix,iz)].imag()
                *hAXm1[IndexXZ(ix,iz)].imag();
      }
    }
  }
  tmpXP*=0.5*InitData.dt*InitData.dt;  EnergyXPOut+=tmpXP;
  tmpXM*=0.5*InitData.dt*InitData.dt;  EnergyXMOut+=tmpXM;
  tmpZP*=0.5*InitData.dt*InitData.dt;  EnergyZPOut+=tmpZP;
  tmpZM*=0.5*InitData.dt*InitData.dt;  EnergyZMOut+=tmpZM;
  tmpZ2*=0.5*InitData.dt*InitData.dt;
  tmpXP2*=0.5*InitData.dt*InitData.dt;
  tmpXM2*=0.5*InitData.dt*InitData.dt;
  tmpGenZ*=-InitData.dt*InitData.dt*InitData.dt*InitData.AGenZ;
  tmpGenX*=-InitData.dt*InitData.dt*InitData.dt*InitData.AGenX;
  EnergyUPOut+=tmpGenZ+tmpGenX;
  hEnergyData[nStepNumber].OutXp=tmpXP;
  hEnergyData[nStepNumber].OutXm=tmpXM;
  hEnergyData[nStepNumber].OutZp=tmpZP;
  hEnergyData[nStepNumber].OutZm=tmpZM;
  hEnergyData[nStepNumber].OutZ2=tmpZ2;
  hEnergyData[nStepNumber].OutXp2=tmpXP2;
  hEnergyData[nStepNumber].OutXm2=tmpXM2;
  hEnergyData[nStepNumber].OutGenZ=tmpGenZ;
  hEnergyData[nStepNumber].Out=tmpXP+tmpXM+tmpZP+tmpZM;
  if( (hEnergyData[nStepNumber].Out>InitData.MaxV[0])&&
      (LastMaxOut>hEnergyData[nStepNumber].Out)&&
      (!MaxESaved))
              SaveADataMaxE();
  LastMaxOut=hEnergyData[nStepNumber].Out;
}

