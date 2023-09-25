//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquOptic2D_1DwB_Mod1.h"
#define EXP_EPS 0.0000001

//---------------------------------------------------------------------------

#pragma package(smart_init)

void cEquOptic2D_1DwB_Mod1::InitReport(cDifEquReport *hReport) {
  int Lines[EQU_MAX_GRAPHS];
  for(int i=0;i<EQU_MAX_GRAPHS;i++) Lines[i]=1;
  Lines[0]=1; Lines[1]=1; Lines[2]=1; Lines[3]=1; Lines[4]=1;
  hReport->LoadText(hReportText);
  hReport->SetGraphSize(EQU_GRAPH_SIZE,Lines);
  hReport->SetMaxTo1();
}

void cEquOptic2D_1DwB_Mod1::Report(cDifEquReport *hReport) {
  nReportNumber++;
  for(int iz=1; iz<=InitData.nStepsZ; iz++)
    hReportBuf[iz-1]=abs(hAZp[I(0,iz,0)]);
  hReport->Graph[0].GetData(0,InitData.nStepsZ,hReportBuf);
  for(int iz=1; iz<=InitData.nStepsZ; iz++)
    hReportBuf[iz-1]=abs(hAZm[I(0,iz,0)]);
  hReport->Graph[1].GetData(0,InitData.nStepsZ,hReportBuf);
  for(int iz=1; iz<=InitData.nStepsZ; iz++)
    hReportBuf[iz-1]=abs(hB[iz]);
  hReport->Graph[2].GetData(0,InitData.nStepsZ,hReportBuf);
  for(int iz=1; iz<=InitData.nStepsZ; iz++)
    hReportBuf[iz-1]=abs(hM[iz].R0);
  hReport->Graph[3].GetData(0,InitData.nStepsZ,hReportBuf);
  for(int i=0; i<EQU_GRAPH_SIZE; i++)
    hReportBuf[i]=abs(Ap(InitData.nStepsZ/2,OPTIC_TYPE(i)/EQU_GRAPH_SIZE));
  hReport->Graph[6].GetData(0,EQU_GRAPH_SIZE,hReportBuf);
  for(int i=0; i<EQU_GRAPH_SIZE; i++)
    hReportBuf[i]=0;
  int dx=EQU_GRAPH_SIZE/(InitData.nStepsX+1);
  for(int ix=1; ix<=InitData.nStepsX; ix++)
    hReportBuf[ix*dx]=hEAHarmP[ix-1];
  hReport->Graph[7].GetData(0,EQU_GRAPH_SIZE,hReportBuf);

  for(int ix=0; ix<EQU_GRAPH_SIZE; ix++) {
    hReportBuf[ix]=0;
    for(int ig=0; ig<InitData.nStepsX; ig++) {
      hReportBuf[ix]+=hEAHarmP[ig]*cos(M_PI*ig*ix/(EQU_GRAPH_SIZE-1));
      if(!ix) hReportBuf[ix]/=2.0;
    }
    hReport->Graph[8].GetData(0,EQU_GRAPH_SIZE,hReportBuf);
  }
  for(int ix=0; ix<InitData.nStepsX2; ix++) {
    hReportBuf[ix]=hM[IM(0,InitData.nStepsZ/2,ix)].R0;
  }
  if(InitData.nStepsX2>3) hReport->Graph[9].GetData(0,InitData.nStepsX2,hReportBuf);
  else hReport->Graph[9].GetData(0,10,hReportBuf);


  hReport->nVals=1;
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
  hReport->ValText[hReport->nVals]="Max B ";
  hReport->Val[hReport->nVals++]=hReport->Graph[2].Max;
  CalcEnergy();

  OPTIC_TYPE EnergyA=EnergyZP+EnergyZM+EnergyB;

  hReport->ValText[hReport->nVals]="EA ";
  hReport->Val[hReport->nVals++]=EnergyZP+EnergyZM;
  hReport->ValText[hReport->nVals]="EB ";
  hReport->Val[hReport->nVals++]=EnergyB;
  hReport->ValText[hReport->nVals]="ER ";
  hReport->Val[hReport->nVals++]=EnergyR;
  hReport->ValText[hReport->nVals]="E Wrong ";
  hReport->Val[hReport->nVals++]=EnergyA-EnergyA0;
  hReport->ValText[hReport->nVals]="E Wrong with out ";
  hReport->Val[hReport->nVals++]=EnergyA+EnergyZPOut+EnergyZMOut-EnergyA0;
  hReport->ValText[hReport->nVals]="E Wrong with R ";
  hReport->Val[hReport->nVals++]=EnergyA+EnergyZPOut+EnergyZMOut+EnergyR-EnergyA0-EnergyQ-EnergyR0;
  hReport->ValText[hReport->nVals]="ER Wrong";
  hReport->Val[hReport->nVals++]=EnergyR-EnergyQ-EnergyR0;
  hReport->ValText[hReport->nVals]="P Wrong ";
  hReport->Val[hReport->nVals++]=PGenR+PGenA-POutA;
  hReport->ValText[hReport->nVals]="Tmp Wrong ";
  hReport->Val[hReport->nVals++]=EnergyA+EnergyR-EnergyA0-EnergyQ-EnergyR0;

  hReport->ValText[hReport->nVals]="Inc ZP ";
  OPTIC_TYPE tmpInc=0;
  if(EnergyZPlast > 0.000000001)
    tmpInc=2.0*log(EnergyZP/EnergyZPlast)/InitData.dt;
  hReport->Val[hReport->nVals++]=tmpInc;


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
  /*
  hReport->nVals=3;
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
  hReport->ValText[1]="Max I "; hReport->Val[1]=hReport->Graph[1].Max;
  hReport->ValText[2]="Max Gamma "; hReport->Val[2]=hReport->Graph[2].Max;
  hReport->Graph[2].Max=2; */
  EnergyZPlast=EnergyZP; ESumLast=EnergyA; LastReportTime=t;
}

void cEquOptic2D_1DwB_Mod1::FinalReport(cDifEquReport *hReport) {

  hReport->nGraphs=0;

  hReport->nVals=0; hReport->nGraphs=2;
  for(int i=0;i<nStepNumber;i++) hReportBuf[i]=hEnergyData[i].Out;
  hReport->Graph[0].GetData(0,nStepNumber,hReportBuf);
  for(int i=0;i<nStepNumber;i++) hReportBuf[i]=hEnergyData[i].R;
  hReport->Graph[1].GetData(0,nStepNumber,hReportBuf);

  hReport->nVals=1;
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
  hReport->ValText[hReport->nVals]="Max B ";
  hReport->Val[hReport->nVals++]=hReport->Graph[2].Max;
  CalcEnergy();
  OPTIC_TYPE EnergyA=EnergyZP+EnergyZM+EnergyB;
  hReport->ValText[hReport->nVals]="EA ";
  hReport->Val[hReport->nVals++]=EnergyZP+EnergyZM;
  hReport->ValText[hReport->nVals]="EB ";
  hReport->Val[hReport->nVals++]=EnergyB;
  hReport->ValText[hReport->nVals]="ER ";
  hReport->Val[hReport->nVals++]=EnergyR;
  hReport->ValText[hReport->nVals]="E Wrong ";
  hReport->Val[hReport->nVals++]=EnergyA-EnergyA0;
  hReport->ValText[hReport->nVals]="E Wrong with out ";
  hReport->Val[hReport->nVals++]=EnergyA+EnergyZPOut+EnergyZMOut-EnergyA0;
  hReport->ValText[hReport->nVals]="E Wrong with R ";
  hReport->Val[hReport->nVals++]=EnergyA+EnergyZPOut+EnergyZMOut+EnergyR-EnergyA0-EnergyQ-EnergyR0;
  hReport->ValText[hReport->nVals]="ER Wrong";
  hReport->Val[hReport->nVals++]=EnergyR-EnergyQ-EnergyR0;
  hReport->ValText[hReport->nVals]="P Wrong ";
  hReport->Val[hReport->nVals++]=PGenR+PGenA-POutA;
  hReport->ValText[hReport->nVals]="P Out B ";
  hReport->Val[hReport->nVals++]=PGenA;
  hReport->ValText[hReport->nVals]="P Out A ";
  hReport->Val[hReport->nVals++]=POutA;
  hReport->ValText[hReport->nVals]="P Rel ";
  hReport->Val[hReport->nVals++]=PRel;

  FinalSave();
}

void cEquOptic2D_1DwB_Mod1::FinalSave() {
  long double tmp1;

  //ofstream FileOutAZp = ofstream("AZp.dat",ios::out);

}

void cEquOptic2D_1DwB_Mod1::LoadInitData(cDifEqu_InitDataRec *hInitData) {
  cDifEqu_InitDataRec_Optic2D_T<OPTIC_TYPE> *hInit=
             (cDifEqu_InitDataRec_Optic2D_T<OPTIC_TYPE> *)hInitData;
  int nStepsT=hInit->tMax/hInit->dt+1;
  OPTIC_TYPE MaxSize = hInit->nStepsZ+2 > nStepsT ? hInit->nStepsZ+2 : nStepsT;
  MaxSize = MaxSize >  hInit->nStepsX ? MaxSize : hInit->nStepsX;
  MaxSize = MaxSize > EQU_GRAPH_SIZE ? MaxSize : EQU_GRAPH_SIZE;
  int size = (hInit->nStepsZ+2)*(hInit->nStepsX);
  delete [] hEnergyData;
  hEnergyData = new EnergyDataType[hInit->tMax/hInit->dt+1];
  delete [] hAZpOutData;
  hAZpOutData = new complex<OPTIC_TYPE>[nStepsT];
  delete [] hEAHarmP;
  hEAHarmP = new OPTIC_TYPE[hInit->nStepsX];

  if((InitData.nStepsZ+2)*(InitData.nStepsX)!=size) {
    delete [] hReportBuf;
    delete [] hB; delete [] hAZp; delete [] hAZm;
    delete [] hM; delete [] hAbsStore;
    delete [] hAZpD; delete [] hAZmD; delete [] hBD; delete [] hMD;
    delete [] hDelta; delete [] hDeltaMatrixCoef;
    delete [] hMatrix;
    hB = new  complex<OPTIC_TYPE>[(hInit->nStepsZ+2)*3];
    hAZp = new complex<OPTIC_TYPE>[size*3];
    hAZm = new complex<OPTIC_TYPE>[size*3];
    //hM   = new cOptic2D_Media_2Way_B_Mod1[(hInit->nStepsZ+2)*3];
    hReportBuf = new OPTIC_TYPE[MaxSize];
    hDelta = new complex<OPTIC_TYPE>[hInit->nStepsX];
    hMatrix = new complex<OPTIC_TYPE>[hInit->nStepsX*4];
    hDeltaMatrixCoef = new complex<OPTIC_TYPE>[hInit->nStepsX*4];
  }
  delete [] hM; delete [] htmpM;
  hM = new cOptic2D_Media_2Way_B_Mod1[(hInit->nStepsZ+2)*3*hInit->nStepsX2];
  htmpM = new cOptic2D_Media_2Way_B_Mod1[hInit->nStepsX2];

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
  CalcEnergy();

  EnergyA0=EnergyZP+EnergyZM+EnergyB; EnergyR0=EnergyR;
  ESumLast=EnergyZPlast=0; LastReportTime=0; LastMaxOut=0; EnergyUPOut=0;
  nStepNumber=0; t=0; MaxESaved=0; LastBitmapSaveTime=-1; nBitmap=0;
  maxmaxAZp=maxmaxAZm=0;
  normAZp=normAZm=100;
  MakeGlobalNormBitmap=MakeComplexBitmap=0;
  PurgeNoiseTime=0;
}

void cEquOptic2D_1DwB_Mod1::PrepareMedia() {
  OPTIC_TYPE r1,r2,r3,r4;
  complex<OPTIC_TYPE> one(1,0);
  complex<OPTIC_TYPE> ii(0,1);
  double ddt=InitData.dt;
  int tmpI;
  for(int iz=0;iz<InitData.nStepsZ+2;iz++) {
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      hAZp[I(0,iz,ix)]=0; hAZm[I(0,iz,ix)]=0; hB[iz]=0;
    }
    for(int ix2=0; ix2<InitData.nStepsX2; ix2++) {
      tmpI=IM(0,iz,ix2);
      hM[IM(0,iz,ix2)].R0=InitData.StartingInv;
      hM[IM(0,iz,ix2)].R0Generated=0;
    }
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      hAZp[I(0,iz,ix)]+=InitData.A0*sin(M_PI*iz/2/InitData.nStepsZ)+
                InitData.A1*sin(M_PI*iz/InitData.nStepsZ)+
                InitData.A2*sin(4*M_PI*iz/InitData.nStepsZ);
      hAZm[I(0,iz,ix)]+=InitData.A0*sin(M_PI*iz/InitData.nStepsZ)+
                InitData.A1*sin(2*M_PI*iz/InitData.nStepsZ)+
                InitData.A2*sin(3*M_PI*iz/InitData.nStepsZ);
      if(InitData.rPolNoise) {
        r1=2.0*M_PI*rand()/RAND_MAX;
        r2=2.0*M_PI*rand()/RAND_MAX;
        hAZp[I(0,iz,ix)]+=polar(InitData.rPolNoise,r1);
        hAZm[I(0,iz,ix)]+=polar(InitData.rPolNoise,r2);
      }
    }
  }
  for(int ix=0;ix<InitData.nStepsX;ix++) {
    hDelta[ix]=complex<OPTIC_TYPE>(InitData.AGenZ,-InitData.DeltaCoef*ix*ix);
    /*
    hDeltaMatrixCoef[2*ix]=(one+ii*ddt*hDelta[ix]-exp(ii*ddt*hDelta[ix]))/
                           (ddt*hDelta[ix]*hDelta[ix]);
    hDeltaMatrixCoef[2*ix+1]=( (one-ii*ddt*hDelta[ix])*exp(ii*ddt*hDelta[ix])-one)/
                           (ddt*hDelta[ix]*hDelta[ix]);
    hDeltaMatrixCoef[2*ix+1]*=-ii*InitData.ACoef*exp(ii*ddt*hDelta[ix]);
    */
    if(abs(hDelta[ix])>EXP_EPS) {
      hDeltaMatrixCoef[4*ix]=(exp(hDelta[ix]*ddt)-hDelta[ix]*ddt-one)
         /(hDelta[ix]*hDelta[ix]*ddt);
      hDeltaMatrixCoef[4*ix+1]=(one-exp(hDelta[ix]*ddt)*(one-hDelta[ix]*ddt))
         /(hDelta[ix]*hDelta[ix]*ddt);
      hDeltaMatrixCoef[4*ix+2]=(exp(-hDelta[ix]*ddt)+hDelta[ix]*ddt-one)
         /(hDelta[ix]*hDelta[ix]*ddt);
      hDeltaMatrixCoef[4*ix+3]=(one-exp(-hDelta[ix]*ddt)*(one+hDelta[ix]*ddt))
         /(hDelta[ix]*hDelta[ix]*ddt);
    } else {
      hDeltaMatrixCoef[4*ix]=hDeltaMatrixCoef[4*ix+1]=
        hDeltaMatrixCoef[4*ix+2]=hDeltaMatrixCoef[4*ix+3]=0.5;
    }
  }
}

void cEquOptic2D_1DwB_Mod1::FillBorders() {
  int sz=InitData.SizeZ;
  for(int ix=0;ix<InitData.nStepsX;ix++) {
    hAZp[I(0,0,ix)]=0;
    hAZm[I(0,InitData.nStepsZ+1,ix)]=0;
    hAZp[I(1,0,ix)]=0;
    hAZm[I(1,InitData.nStepsZ+1,ix)]=0;
  }
}

void cEquOptic2D_1DwB_Mod1::StepRoutine() {
   StepRoutineRunge2MatrixExpRx();
}

void cEquOptic2D_1DwB_Mod1::StepRoutineEuler() {
  complex<OPTIC_TYPE> tmpPZp,tmpPZm,tmpPB,tmpA,SumA,SumAP;
  OPTIC_TYPE dt2=InitData.dt;
  complex<OPTIC_TYPE> ii=complex<OPTIC_TYPE>(0,1),tmpAZp,tmpAZm,tmpB;
  cOptic2D_Media_2Way_B tmpM;
  OPTIC_TYPE two=2;
  int sz=InitData.SizeZ;
                        // InitData.AGenCoef - коэффициент связи среды с волной B !
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    SumA=SumAP=0;
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      tmpPZp=hAZp[I(0,iz,ix)]*hM[iz].R0*two*InitData.beta;
      tmpPZm=hAZm[I(0,iz,ix)]*hM[iz].R0*two*InitData.beta;
      SumA+=hAZp[I(0,iz,ix)]+hAZm[I(0,iz,ix)];
      SumAP+=hAZp[I(0,iz,ix)]*conj(tmpPZp)+hAZm[I(0,iz,ix)]*conj(tmpPZm);
      hAZp[I(1,iz,ix)]=dt2*tmpPZp+dt2*(hDelta[ix]+InitData.AGenZ)*hAZp[I(0,iz,ix)]+ii*InitData.ACoef*hB[iz]*dt2;
      hAZm[I(1,iz,ix)]=dt2*tmpPZm+dt2*(hDelta[ix]+InitData.AGenZ)*hAZm[I(0,iz,ix)]+ii*InitData.ACoef*hB[iz]*dt2;
    }
    tmpPB=(hB[iz]*hM[iz].R0*two*InitData.AGenCoef)*InitData.beta;
    hM[iz+sz].R0 =(InitData.Q-InitData.RCoef*hM[iz].R0
                       -(SumAP+hB[iz]*conj(tmpPB)).real() )*dt2;
    hM[iz+sz].R0Generated=(InitData.Q-InitData.RCoef*hM[iz].R0)*dt2;
    hB[iz+sz]=(ii*InitData.ACoef*(SumA)+tmpPB)*dt2;
  }

  
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      hAZm[I(0,iz,ix)]=hAZm[I(0,iz+1,ix)]+hAZm[I(1,iz+1,ix)];
      hAZp[I(0,InitData.nStepsZ+1-iz,ix)]=
           hAZp[I(0,InitData.nStepsZ-iz,ix)]
          +hAZp[I(1,InitData.nStepsZ-iz,ix)];
    }
    hM[iz].M1M2v(hM[iz],hM[iz+sz],1);
    hB[iz]=hB[iz]+hB[iz+sz];
  }

  /*
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      hAZm[I(0,iz,ix)]=hAZm[I(0,iz,ix)]+hAZm[I(1,iz,ix)];
      hAZp[I(0,iz,ix)]=
           hAZp[I(0,iz,ix)]
          +hAZp[I(1,iz,ix)];
    }
    hM[iz].M1M2v(hM[iz],hM[iz+sz],1);
    hB[iz]=hB[iz]+hB[iz+sz];
  }
 */
  FillBorders();
  CalcOutEnergy();
  hAZpOutData[nStepNumber]=hAZp[I(0,InitData.nStepsZ-1,0)];
  for(int i=0;i<FFTSizeFull-1;i++) hFFTSrc[i]=hFFTSrc[i+1];
  hFFTSrc[FFTSizeFull-1]=hAZp[I(0,InitData.nStepsZ/2,0)];
}


void cEquOptic2D_1DwB_Mod1::StepRoutineRunge2MatrixExpRx() {

  complex<OPTIC_TYPE> tmpPZp,tmpPZm,tmpPB,gamma2,e1,bb,tmpPB2,tmpPZp2,tmpPZm2,ee1;
  OPTIC_TYPE dt2=0.5*InitData.dt;
  OPTIC_TYPE ddt=InitData.dt;
  complex<OPTIC_TYPE> ii=complex<OPTIC_TYPE>(0,1),tmpAZp,tmpAZm,tmpB,tAZpM,tAZmM,tBM,SumA,SumAP,SumA2;
  complex<OPTIC_TYPE> BCoef=complex<OPTIC_TYPE>(InitData.AGenX,InitData.Delta);
  complex<OPTIC_TYPE> ia=ii*InitData.ACoef;
  complex<OPTIC_TYPE> q=-ii*InitData.ACoef*dt2;
  OPTIC_TYPE two=2, tmpRx;
  complex<OPTIC_TYPE> one = complex<OPTIC_TYPE>(1,0);
  cOptic2D_Media_2Way_B_Mod1 tmpM;
  int sz=InitData.SizeZ;
  complex<OPTIC_TYPE> q1=complex<OPTIC_TYPE>(1,0);//complex<OPTIC_TYPE>(1,0)/(OPTIC_TYPE(1)-q*q*two);
  /*
  for(int i=0; i<InitData.nStepsX; i++)
    q1+=q*q*two/(one-hDelta[i]*dt2);
  q1=one-BCoef*dt2-q1;
  q1=one/q1;
  */
  for(int i=0; i<InitData.nStepsX; i++)
    q1+=InitData.ACoef*hDelta[i]*K2m(i)*K2p(i);
  q1=one/q1;
                        // InitData.AGenCoef - коэффициент связи среды с волной B !
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    SumA=SumAP=0;
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      tmpPZp=hAZp[I(0,iz,ix)]*hM[iz].R0*two*InitData.beta;
      tmpPZm=hAZm[I(0,iz,ix)]*hM[iz].R0*two*InitData.beta;
      SumAP+=hAZp[I(0,iz,ix)]*conj(tmpPZp)+hAZm[I(0,iz,ix)]*conj(tmpPZm);
      gamma2=hDelta[ix];
      for(int ix2=0; ix2<InitData.nStepsX2; ix2++) gamma2+=hM[IM(0,iz,ix2)].R0*two*InitData.beta*MediaCoef(ix,ix2);
        e1=exp(gamma2*ddt);
      if(abs(gamma2)>EXP_EPS) {
        bb=ia*hB[iz]*(1.0-1.0/e1)/gamma2;
      } else {
        e1=one + gamma2*ddt;
        bb=ia*hB[iz]*ddt;
      }
      hAZp[I(1,iz,ix)]=(hAZp[I(0,iz,ix)]+bb)*e1;
      hAZp[I(1,iz,ix)]=(hAZm[I(1,iz,ix)]+bb)*e1;
    }
    tmpPB=0;
    for(int ix2=0; ix2<InitData.nStepsX2; ix2++) {
      SumA=0;
      for(int ix=0;ix<InitData.nStepsX;ix++) {
        SumA+=(hAZp[I(0,iz,ix)]+hAZm[I(0,iz,ix)])*MediaCoef(ix,ix2);
      }
      tmpPB2=(hB[iz]*hM[IM(0,iz,ix2)].R0*two*InitData.AGenCoef)*InitData.beta;
      hM[IM(1,iz,ix2)].R0 =(InitData.Q-InitData.RCoef*hM[IM(0,iz,ix2)].R0
                         -(SumAP+hB[iz]*conj(tmpPB2)).real() )*dt2;
      hM[IM(1,iz,ix2)].R0Generated=(InitData.Q-InitData.RCoef*hM[IM(0,iz,ix2)].R0)*dt2;
      tmpPB+=tmpPB2;
    }
    hB[iz+sz]=ia*SumA+tmpPB*dt2;
  }

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    SumA=SumAP=SumA2=0;
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      tmpRx=0;
      for(int ix2=0; ix2<InitData.nStepsX2; ix2++) {
        tmpRx+=hM[IM(0,iz,ix2)].R0*MediaCoef(ix,ix2);
      }
      tmpRx=exp(tmpRx*InitData.beta*dt2);
      hAZp[I(2,iz,ix)]=hAZp[I(0,iz-1,ix)]*tmpRx+ii*hDelta[ix]*hB[iz-1]*K1m(ix);
      hAZm[I(2,iz,ix)]=hAZm[I(0,iz+1,ix)]*tmpRx+ii*hDelta[ix]*hB[iz+1]*K1m(ix);
      SumA += -ii*hDelta[ix]*K2m(ix)*(hAZp[I(2,iz,ix)]+hAZm[I(2,iz,ix)]);
      SumA2+= ia*(hAZp[I(0,iz,ix)]+hAZm[I(0,iz,ix)])*K1p(ix);
    }
    hB[iz+2*sz]=q1*(hB[iz]+SumA2-SumA);
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      hAZp[I(2,iz,ix)]=(hAZp[I(2,iz,ix)]-ia*K2p(ix)*hB[iz+2*sz])*exp(hDelta[ix]*ddt);
      hAZm[I(2,iz,ix)]=(hAZm[I(2,iz,ix)]-ia*K2p(ix)*hB[iz+2*sz])*exp(hDelta[ix]*ddt);
    }
  }

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    tmpPB=0;
    tmpB  = hB[iz]+hB[iz+sz]*two;
    for(int ix2=0; ix2<InitData.nStepsX2; ix2++) {
      htmpM[ix2].M1M2v(hM[IM(0,iz,ix2)],hM[IM(1,iz,ix2)],2);
      hM[IM(2,iz,ix2)].R0 = hM[IM(2,iz,ix2)].R0Generated = (InitData.Q-InitData.RCoef*htmpM[ix2].R0)*dt2;
      tmpPB+=tmpB*htmpM[ix2].R0*two*InitData.AGenCoef;
    }
    hB[iz+2*sz]+=(tmpPB+tmpB*BCoef)*dt2;
    for(int ix2=0; ix2<InitData.nStepsX2; ix2++) {
      hM[IM(2,iz,ix2)].R0 -= (tmpB*conj(tmpPB)).real()*dt2;
    }
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      tmpAZp = (hAZp[I(1,iz-1,ix)])*two;
      tmpAZm = (hAZm[I(1,iz+1,ix)])*two;
      e1=0; tmpPZp=tmpPZm=0;
      for(int ix2=0; ix2<InitData.nStepsX2; ix2++) {

        e1+=htmpM[ix2].R0*MediaCoef(ix,ix2);
        tmpPZp2=tmpAZp*htmpM[ix2].R0*MediaCoef(ix,ix2)*two*InitData.beta;
        tmpPZm2=tmpAZm*htmpM[ix2].R0*MediaCoef(ix,ix2)*two*InitData.beta;
        hM[IM(2,iz,ix2)].R0+=- ((tmpAZp*conj(tmpPZp2)).real()+(tmpAZm*conj(tmpPZm2)).real())*dt2;
        tmpPZp+=tmpPZp2;
        tmpPZm+=tmpPZm2;
      }
      e1=exp(e1*InitData.beta*dt2);
      hAZp[iz+2*sz]*=e1;
      hAZm[iz+2*sz]*=e1;
    }
  }

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    for(int ix2=0; ix2<InitData.nStepsX2; ix2++) {
      htmpM[ix2].M1M2v(hM[IM(0,iz,ix2)],hM[IM(1,iz,ix2)],1);
      hM[IM(0,iz,ix2)].M1M2v(htmpM[ix2],hM[IM(2,iz,ix2)],1);
    }
    hB[iz]=hB[iz+2*sz];
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      hAZp[I(0,iz,ix)]=hAZp[I(2,iz,ix)];
      hAZm[I(0,iz,ix)]=hAZm[I(2,iz,ix)];
    }
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


void cEquOptic2D_1DwB_Mod1::StepRoutineRunge2MatrixExp() {
  
  complex<OPTIC_TYPE> tmpPZp,tmpPZm,tmpPB,gamma2,e1,bb;
  OPTIC_TYPE dt2=0.5*InitData.dt;
  OPTIC_TYPE ddt=InitData.dt;
  complex<OPTIC_TYPE> ii=complex<OPTIC_TYPE>(0,1),tmpAZp,tmpAZm,tmpB,tAZpM,tAZmM,tBM,SumA,SumAP,SumA2;
  complex<OPTIC_TYPE> BCoef=complex<OPTIC_TYPE>(InitData.AGenX,InitData.Delta);
  complex<OPTIC_TYPE> ia=ii*InitData.ACoef;
  complex<OPTIC_TYPE> q=-ii*InitData.ACoef*dt2;
  OPTIC_TYPE two=2;
  complex<OPTIC_TYPE> one = complex<OPTIC_TYPE>(1,0);
  cOptic2D_Media_2Way_B_Mod1 tmpM;
  int sz=InitData.SizeZ;
  complex<OPTIC_TYPE> q1=0;//complex<OPTIC_TYPE>(1,0)/(OPTIC_TYPE(1)-q*q*two);
  /*
  for(int i=0; i<InitData.nStepsX; i++)
    q1+=q*q*two/(one-hDelta[i]*dt2);
  q1=one-BCoef*dt2-q1;
  q1=one/q1;
  */
  for(int i=0; i<InitData.nStepsX; i++)
    q1+=InitData.ACoef*hDelta[i]*K2m(i)*K2p(i);
  q1=one/q1;
                        // InitData.AGenCoef - коэффициент связи среды с волной B !
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    SumA=SumAP=0;
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      tmpPZp=hAZp[I(0,iz,ix)]*hM[iz].R0*two*InitData.beta;
      tmpPZm=hAZm[I(0,iz,ix)]*hM[iz].R0*two*InitData.beta;
      SumAP+=hAZp[I(0,iz,ix)]*conj(tmpPZp)+hAZm[I(0,iz,ix)]*conj(tmpPZm);
      gamma2=hDelta[ix]+hM[iz].R0*two*InitData.beta;
      e1=exp(-gamma2*ddt);
      bb=ia*hB[iz]*(1.0-e1)/gamma2;
      SumA+=(hAZp[I(0,iz,ix)]+hAZm[I(0,iz,ix)])*(1.0-1.0/e1)/gamma2;
//      if(!iz) bb/=2.0;
      hAZp[I(1,iz,ix)]=(hAZp[I(0,iz,ix)]+bb)/e1;
      hAZp[I(1,iz,ix)]=(hAZm[I(1,iz,ix)]+bb)/e1;
    }
    tmpPB=(hB[iz]*hM[iz].R0*two*InitData.AGenCoef)*InitData.beta;
    hM[iz+sz].R0 =(InitData.Q-InitData.RCoef*hM[iz].R0
                       -(SumAP+hB[iz]*conj(tmpPB)).real() )*dt2;
    hM[iz+sz].R0Generated=(InitData.Q-InitData.RCoef*hM[iz].R0)*dt2;
    hB[iz+sz]=ia*SumA+tmpPB*dt2;
  }

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    SumA=SumAP=SumA2=0;
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      hAZp[I(2,iz,ix)]=hAZp[I(0,iz-1,ix)]*exp(hM[iz].R0*InitData.beta*dt2)+ii*hDelta[ix]*hB[iz-1]*K1m(ix);
      hAZm[I(2,iz,ix)]=hAZm[I(0,iz+1,ix)]*exp(hM[iz].R0*InitData.beta*dt2)+ii*hDelta[ix]*hB[iz+1]*K1m(ix);
      SumA += -ii*hDelta[ix]*K2m(ix)*(hAZp[I(2,iz,ix)]+hAZm[I(2,iz,ix)]);
      SumA2+= ia*(hAZp[I(0,iz,ix)]+hAZm[I(0,iz,ix)])*K1p(ix);
    }
    hB[iz+2*sz]=q1*(hB[iz]+SumA2-SumA);
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      hAZp[I(2,iz,ix)]=(hAZp[I(2,iz,ix)]-ia*K2p(ix)*hB[iz+2*sz])*exp(hDelta[ix]*ddt);
      hAZm[I(2,iz,ix)]=(hAZm[I(2,iz,ix)]-ia*K2p(ix)*hB[iz+2*sz])*exp(hDelta[ix]*ddt);
    }
  }

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    tmpM.M1M2v(hM[iz],hM[iz+sz],2);
    hM[iz+2*sz].R0 = hM[iz+2*sz].R0Generated=(InitData.Q-InitData.RCoef*tmpM.R0)*dt2;
    e1=exp(tmpM.R0*InitData.beta*dt2);
    tmpB  = hB[iz]+hB[iz+sz]*two;
    tmpPB = tmpB*tmpM.R0*two*InitData.AGenCoef;
    hB[iz+2*sz]+=(tmpPB+tmpB*BCoef)*dt2;
    hM[iz+2*sz].R0 -= (tmpB*conj(tmpPB)).real()*dt2;
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      tmpAZp = (hAZp[I(1,iz-1,ix)])*two;
      tmpAZm = (hAZm[I(1,iz+1,ix)])*two;
      tmpPZp=tmpAZp*tmpM.R0*two*InitData.beta;
      tmpPZm=tmpAZm*tmpM.R0*two*InitData.beta;
      hM[iz+2*sz].R0+=- ((tmpAZp*conj(tmpPZp)).real()+(tmpAZm*conj(tmpPZm)).real())*dt2;
      hAZp[iz+2*sz]*=e1;
      hAZm[iz+2*sz]*=e1;
    }
  }

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    tmpM.M1M2v(hM[iz],hM[iz+sz],1);
    hM[iz].M1M2v(tmpM,hM[iz+2*sz],1);
    hB[iz]=hB[iz+2*sz];
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      hAZp[I(0,iz,ix)]=hAZp[I(2,iz,ix)];
      hAZm[I(0,iz,ix)]=hAZm[I(2,iz,ix)];
    }
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


void cEquOptic2D_1DwB_Mod1::StepRoutineRunge2Matrix() {

  complex<OPTIC_TYPE> tmpPZp,tmpPZm,tmpPB;
  OPTIC_TYPE dt2=0.5*InitData.dt;
  complex<OPTIC_TYPE> ii=complex<OPTIC_TYPE>(0,1),tmpAZp,tmpAZm,tmpB,tAZpM,tAZmM,tBM,SumA,SumAP;
  complex<OPTIC_TYPE> BCoef=complex<OPTIC_TYPE>(InitData.AGenX,InitData.Delta);
  complex<OPTIC_TYPE> ia=ii*InitData.ACoef;
  complex<OPTIC_TYPE> q=-ii*InitData.ACoef*dt2;
  OPTIC_TYPE two=2;
  complex<OPTIC_TYPE> one = complex<OPTIC_TYPE>(1,0);
  cOptic2D_Media_2Way_B_Mod1 tmpM;
  int sz=InitData.SizeZ;
  complex<OPTIC_TYPE> q1=0;//complex<OPTIC_TYPE>(1,0)/(OPTIC_TYPE(1)-q*q*two);
  for(int i=0; i<InitData.nStepsX; i++)
    q1+=q*q*two/(one-hDelta[i]*dt2);
  q1=one-BCoef*dt2-q1;
  q1=one/q1;

                        // InitData.AGenCoef - коэффициент связи среды с волной B !
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    SumA=SumAP=0;
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      tmpPZp=hAZp[I(0,iz,ix)]*hM[iz].R0*two*InitData.beta;
      tmpPZm=hAZm[I(0,iz,ix)]*hM[iz].R0*two*InitData.beta;
      SumA+=hAZp[I(0,iz,ix)]+hAZm[I(0,iz,ix)];
      SumAP+=hAZp[I(0,iz,ix)]*conj(tmpPZp)+hAZm[I(0,iz,ix)]*conj(tmpPZm);
      hAZp[I(1,iz,ix)]=dt2*tmpPZp+dt2*hDelta[ix]*hAZp[I(0,iz,ix)]+ia*hB[iz]*dt2;
      hAZm[I(1,iz,ix)]=dt2*tmpPZm+dt2*hDelta[ix]*hAZm[I(0,iz,ix)]+ia*hB[iz]*dt2;
    }
    tmpPB=(hB[iz]*hM[iz].R0*two*InitData.AGenCoef)*InitData.beta;
    hM[iz+sz].R0 =(InitData.Q-InitData.RCoef*hM[iz].R0
                       -(SumAP+hB[iz]*conj(tmpPB)).real() )*dt2;
    hM[iz+sz].R0Generated=(InitData.Q-InitData.RCoef*hM[iz].R0)*dt2;
    hB[iz+sz]=(ia*SumA+tmpPB)*dt2;
  }

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    SumA=SumAP=0;
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      SumA += (hAZp[I(0,iz-1,ix)]+hAZp[I(1,iz-1,ix)]+hAZm[I(0,iz+1,ix)]+hAZm[I(1,iz+1,ix)])
                        /(complex<OPTIC_TYPE>(1,0)-hDelta[ix]*dt2);
    }
    tmpB   = hB[iz]+hB[iz+sz];
    hB[iz+2*sz]=q1*(tmpB-q*(SumA));
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      hAZp[I(2,iz,ix)]=(hAZp[I(0,iz-1,ix)]+hAZp[I(1,iz-1,ix)])/(one-hDelta[ix]*dt2)-q*hB[iz+2*sz];
      hAZm[I(2,iz,ix)]=(hAZm[I(0,iz+1,ix)]+hAZm[I(1,iz+1,ix)])/(one-hDelta[ix]*dt2)-q*hB[iz+2*sz];
    }
  }

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    tmpM.M1M2v(hM[iz],hM[iz+sz],2);
    hM[iz+2*sz].R0 = hM[iz+2*sz].R0Generated=(InitData.Q-InitData.RCoef*tmpM.R0)*dt2;
    tmpB  = hB[iz]+hB[iz+sz]*two;
    tmpPB = tmpB*tmpM.R0*two*InitData.AGenCoef;
    hB[iz+2*sz]+=(tmpPB+tmpB*BCoef)*dt2;
    hM[iz+2*sz].R0 -= (tmpB*conj(tmpPB)).real()*dt2;
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      tmpAZp = (hAZp[I(0,iz-1,ix)]+hAZp[I(1,iz-1,ix)])*two;
      tmpAZm = (hAZm[I(0,iz+1,ix)]+hAZm[I(1,iz+1,ix)])*two;
      tmpPZp=tmpAZp*tmpM.R0*two*InitData.beta;
      tmpPZm=tmpAZm*tmpM.R0*two*InitData.beta;
      hM[iz+2*sz].R0+=- ((tmpAZp*conj(tmpPZp)).real()+(tmpAZm*conj(tmpPZm)).real())*dt2;
      hAZp[iz+2*sz]+=dt2*tmpPZp;
      hAZm[iz+2*sz]+=dt2*tmpPZm;
    }
  }

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    tmpM.M1M2v(hM[iz],hM[iz+sz],1);
    hM[iz].M1M2v(tmpM,hM[iz+2*sz],1);
    hB[iz]=hB[iz+2*sz];
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      hAZp[I(0,iz,ix)]=hAZp[I(2,iz,ix)];
      hAZm[I(0,iz,ix)]=hAZm[I(2,iz,ix)];
    }
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


void cEquOptic2D_1DwB_Mod1::StepRoutineRunge2() {
  complex<OPTIC_TYPE> tmpPZp,tmpPZm,tmpPB,tmpA,SumA,SumAP;
  OPTIC_TYPE dt2=InitData.dt*0.5;
  complex<OPTIC_TYPE> ii=complex<OPTIC_TYPE>(0,1),tmpAZp,tmpAZm,tmpB;
  cOptic2D_Media_2Way_B_Mod1 tmpM;
  OPTIC_TYPE two=2;
  int sz=InitData.SizeZ;
  int size = (InitData.nStepsZ+2)*(InitData.nStepsX);
  int tt=0;
                        // InitData.AGenCoef - коэффициент связи среды с волной B !
 /*
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    SumA=SumAP=0;
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      tmpPZp=hAZp[I(0,iz,ix)]*hM[iz].R0*two*InitData.beta;
      tmpPZm=hAZm[I(0,iz,ix)]*hM[iz].R0*two*InitData.beta;
      SumA+=hAZp[I(0,iz,ix)]+hAZm[I(0,iz,ix)];
      SumAP+=hAZp[I(0,iz,ix)]*conj(tmpPZp)+hAZm[I(0,iz,ix)]*conj(tmpPZm);
      hAZp[I(1,iz,ix)]=dt2*tmpPZp+dt2*(hDelta[ix]+InitData.AGenZ)*hAZp[I(0,iz,ix)]+ii*InitData.ACoef*hB[iz]*dt2;
      hAZm[I(1,iz,ix)]=dt2*tmpPZm+dt2*(hDelta[ix]+InitData.AGenZ)*hAZm[I(0,iz,ix)]+ii*InitData.ACoef*hB[iz]*dt2;
    }
    tmpPB=(hB[iz]*hM[iz].R0*two*InitData.AGenCoef)*InitData.beta;
    hM[iz+sz].R0 =(InitData.Q-InitData.RCoef*hM[iz].R0
                       -(SumAP+hB[iz]*conj(tmpPB)).real() )*dt2;
    hM[iz+sz].R0Generated=(InitData.Q-InitData.RCoef*hM[iz].R0)*dt2;
    hB[iz+sz]=(ii*InitData.ACoef*(SumA)+tmpPB)*dt2;
  }

  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    tmpM.M1M2v(hM[iz],hM[iz+sz],2);
    SumA=SumAP=0;
    tmpB=hB[iz]+hB[iz+sz]*two;
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      tmpAZp=hAZp[I(0,iz-1,ix)]+hAZp[I(1,iz-1,ix)]*two;
      tmpAZm=hAZm[I(0,iz+1,ix)]+hAZm[I(1,iz+1,ix)]*two;
      tmpPZp=tmpAZp*tmpM.R0*two*InitData.beta;
      tmpPZm=tmpAZm*tmpM.R0*two*InitData.beta;
      SumA+=tmpAZp+tmpAZm;
      SumAP+=tmpAZp*conj(tmpPZp)+tmpAZm*conj(tmpPZm);
      hAZp[I(2,iz,ix)]=dt2*tmpPZp+dt2*(hDelta[ix]+InitData.AGenZ)*tmpAZp+ii*InitData.ACoef*tmpB*dt2;
      hAZm[I(2,iz,ix)]=dt2*tmpPZm+dt2*(hDelta[ix]+InitData.AGenZ)*tmpAZm+ii*InitData.ACoef*tmpB*dt2;
    }
    tmpPB=(tmpB*tmpM.R0*two*InitData.AGenCoef)*InitData.beta;
    hM[iz+2*sz].R0 =(InitData.Q-InitData.RCoef*tmpM.R0
                       -(SumAP+tmpB*conj(tmpPB)).real() )*dt2;
    hM[iz+2*sz].R0Generated=(InitData.Q-InitData.RCoef*tmpM.R0)*dt2;
    hB[iz+2*sz]=(ii*InitData.ACoef*(SumA)+tmpPB)*dt2;
  }


  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      hAZm[I(0,iz,ix)]=hAZm[I(0,iz+1,ix)]+hAZm[I(1,iz+1,ix)]+hAZm[I(2,iz,ix)];
      hAZp[I(0,InitData.nStepsZ+1-iz,ix)]=
           hAZp[I(0,InitData.nStepsZ-iz,ix)]
          +hAZp[I(1,InitData.nStepsZ-iz,ix)];
          +hAZp[I(2,InitData.nStepsZ-iz+1,ix)];
    }
    tmpM.M1M2v(hM[iz],hM[iz+sz],1);
    hM[iz].M1M2v(tmpM,hM[iz+2*sz],1);
    hB[iz]=hB[iz]+hB[iz+sz]+hB[iz+2*sz];
  }
  */
  /*
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    for(int ix=0;ix<InitData.nStepsX;ix++) {
      hAZm[I(0,iz,ix)]=hAZm[I(0,iz,ix)]+hAZm[I(1,iz,ix)];
      hAZp[I(0,iz,ix)]=
           hAZp[I(0,iz,ix)]
          +hAZp[I(1,iz,ix)];
    }
    hM[iz].M1M2v(hM[iz],hM[iz+sz],1);
    hB[iz]=hB[iz]+hB[iz+sz];
  }
  */
  FillBorders();
  CalcOutEnergy();
  hAZpOutData[nStepNumber]=hAZp[I(0,InitData.nStepsZ-1,0)];
  for(int i=0;i<FFTSizeFull-1;i++) hFFTSrc[i]=hFFTSrc[i+1];
  hFFTSrc[FFTSizeFull-1]=hAZp[I(0,InitData.nStepsZ/2,0)];
}


complex<OPTIC_TYPE> cEquOptic2D_1DwB_Mod1::Ap(int iz, double x) {
  complex<OPTIC_TYPE> Res=0;
  for(int ix=0; ix<InitData.nStepsX; ix++)  {
    Res+= hAZp[I(0,iz,ix)]*cos(M_PI*ix*x);
  }
  return Res;
}

void cEquOptic2D_1DwB_Mod1::CalcEnergy() {
  EnergyB=EnergyZP=EnergyZM=EnergyR=EnergyQ=0;
  PGenR=PGenA=PGenZ=PGenX=0;
  for(int ix=0; ix<InitData.nStepsX; ix++) hEAHarmP[ix]=0;
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    EnergyB+=hB[iz].real()*hB[iz].real()+hB[iz].imag()*hB[iz].imag();
    EnergyR+=hM[iz].R0;
    EnergyQ+=hM[iz].R0Generated;
    PGenX+=InitData.AGenX*(hB[iz].real()*hB[iz].real()+hB[iz].imag()*hB[iz].imag());
    for(int ix=0; ix<InitData.nStepsX; ix++)  {
      hEAHarmP[ix]+=hAZp[I(0,iz,ix)].real()*hAZp[I(0,iz,ix)].real()+hAZp[I(0,iz,ix)].imag()*hAZp[I(0,iz,ix)].imag();
      EnergyZP+=hAZp[I(0,iz,ix)].real()*hAZp[I(0,iz,ix)].real()+hAZp[I(0,iz,ix)].imag()*hAZp[I(0,iz,ix)].imag();
      EnergyZM+=hAZm[I(0,iz,ix)].real()*hAZm[I(0,iz,ix)].real()+hAZm[I(0,iz,ix)].imag()*hAZm[I(0,iz,ix)].imag();
      PGenZ+=InitData.AGenZ*(hAZp[I(0,iz,ix)].real()*hAZp[I(0,iz,ix)].real()+hAZp[I(0,iz,ix)].imag()*hAZp[I(0,iz,ix)].imag());
      PGenZ+=InitData.AGenZ*(hAZm[I(0,iz,ix)].real()*hAZm[I(0,iz,ix)].real()+hAZm[I(0,iz,ix)].imag()*hAZm[I(0,iz,ix)].imag());
    }
    PGenR+=InitData.Q-hM[iz].R0*InitData.RCoef;
  }
  PGenX*=InitData.dt;
  PGenZ*=InitData.dt;
  PGenA=PGenX+PGenZ;
  PGenR*=InitData.dt;
  POutA=0;
  for(int ix=0; ix<InitData.nStepsX; ix++)  {
    POutA+=hAZp[I(0,InitData.nStepsZ,ix)].real()*hAZp[I(0,InitData.nStepsZ,ix)].real()+
                hAZp[I(0,InitData.nStepsZ,ix)].imag()*hAZp[I(0,InitData.nStepsZ,ix)].imag();
    POutA+=hAZm[I(0,1,ix)].real()*hAZm[(0,1,ix)].real()+hAZm[(0,1,ix)].imag()*hAZm[(0,1,ix)].imag();
  }
  POutA*=0.5;
  PRel=-PGenR+InitData.Q*InitData.dt*InitData.nStepsZ;

  for(int ix=0; ix<InitData.nStepsX; ix++)  hEAHarmP[ix]*=0.5*InitData.dt;
  EnergyB*=0.5*InitData.dt;
  EnergyZP*=0.5*InitData.dt;
  EnergyZM*=0.5*InitData.dt;
  EnergyR*=InitData.dt;
  EnergyQ*=InitData.dt;
}

void cEquOptic2D_1DwB_Mod1::CalcOutEnergy() {
  OPTIC_TYPE tmpXP=0,tmpXM=0,tmpZP=0,tmpZM=0,tmpGenZ=0,tmpGenX=0,tmpR=0;
  for(int ix=0; ix<InitData.nStepsX; ix++)  {
    tmpZP+= hAZp[I(0,InitData.nStepsZ,ix)].real()*hAZp[I(0,InitData.nStepsZ,ix)].real()
           +hAZp[I(0,InitData.nStepsZ,ix)].imag()*hAZp[I(0,InitData.nStepsZ,ix)].imag();
    tmpZM+=  hAZm[I(0,1,ix)].real()*hAZm[I(0,1,ix)].real()
            +hAZm[I(0,1,ix)].imag()*hAZm[I(0,1,ix)].imag();
  }
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
