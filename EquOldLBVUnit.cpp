//---------------------------------------------------------------------------

#include <vcl.h>
#include <mem.h>
#pragma hdrstop

#include "EquOldLBVUnit.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)


void cDiffurMain::Init(REAL_TYPE pL,REAL_TYPE ptMax,
              int pnSteps,int pnEls,REAL_TYPE prTime,REAL_TYPE pr,
              REAL_TYPE pdelta,REAL_TYPE palpha,
              REAL_TYPE pTInput,REAL_TYPE paInput,
              REAL_TYPE pp,REAL_TYPE pZ0,bool pCalcPower)
{
  if( pnSteps!=nSteps)
  {
    delete [] hIdata; delete [] hIdata1;
    delete [] hAdata; delete [] hAdata1;
    delete [] hInterpolated;
    delete [] hSplineCoefB;
    delete [] hSplineCoefC; delete [] hSplineCoefD;
    hIdata=new complex<REAL_TYPE>[pnSteps+1];
    hIdata1=new complex<REAL_TYPE>[pnSteps+1];
    hAdata=new complex<REAL_TYPE>[pnSteps+1];
    hAdata1=new complex<REAL_TYPE>[pnSteps+1];
    hInterpolated=new complex<REAL_TYPE>[pnSteps+1];
    hSplineCoefB=new REAL_TYPE[pnSteps+1];
    hSplineCoefC=new REAL_TYPE[pnSteps+1];
    hSplineCoefD=new REAL_TYPE[pnSteps+1];
  }
  if( pnEls!=nEls)
  {
    delete [] hEldata; delete [] hElNew;
    delete [] hElk1; delete [] hElk2;
    delete [] hElk3; delete [] hElk4;
    delete [] hElPresent;
    hEldata=new REAL_TYPE[2*pnEls]; hElNew=new REAL_TYPE[2*pnEls];
    hElk1=new REAL_TYPE[2*pnEls]; hElk2=new REAL_TYPE[2*pnEls];
    hElk3=new REAL_TYPE[2*pnEls]; hElk4=new REAL_TYPE[2*pnEls];
    hElPresent = new int[pnEls];
  }
  t=0; L=GammaCoef*pL; tMax=GammaCoef*ptMax; nSteps=pnSteps; nEls=pnEls; tMax_Res=tMax+L;
  rTime=prTime; rValue=pr; delta=pdelta; alpha=palpha;
  dt=L/(REAL_TYPE)nSteps;
  aInput=paInput; if(pTInput>0.000001) TInput_inv=1/pTInput; TInput=pTInput;

  if(pZ0<0) {
    bKappaLinear=true;
    pZ0=-pZ0;
  } else bKappaLinear=false;
  Z0=GammaCoef*pZ0; Z0_inv = (Z0>0.000001) ? 1.0/Z0 : 1.0 ; p=pp;
  if (!bKappaLinear) Ckappa = L*(1+p) / ( (L-Z0)*(1+p)+2.0*Z0 );
    else Ckappa= p>Z0_inv ? L*p/((L-Z0)*(p)+0.5) : L/(L-0.5*p*Z0*Z0);
//  Ckappa = 1.0;

  ElDroppedMaxFrac=0;
  ElDroppedSum=0;
  for(int j=0;j<nSteps+1;j++)
  {
    hIdata[j]=0;
    hAdata[j]=0;  //0.01*complex<REAL_TYPE>(sin(L*j/nSteps),0);
  }

  for(int j=0;j<10;j++) rDebug[j]=-100;
  backOutData=ResIndex=PercentageCounter=nCStep=0; CalcPower=pCalcPower;
}

REAL_TYPE cDiffurMain::CalcMaxA()
{
REAL_TYPE MaxA=0;
REAL_TYPE OutData;
  while(t<tMax)
  {
    Step();
    OutData= (CalcPower) ? hAdata[OutPoint].real()*hAdata[OutPoint].real()+
                           hAdata[OutPoint].imag()*hAdata[OutPoint].imag()
                           : abs(hAdata[OutPoint]);
    OutData/=GammaCoef*GammaCoef*GammaCoef;
    if (OutData>MaxA) MaxA=OutData;
  }
  for(int i=0;i<nSteps+1;i++)
  {
    OutData=(CalcPower) ? hAdata[i].real()*hAdata[i].real()+
                           hAdata[i].imag()*hAdata[i].imag()
                           : abs(hAdata[i]);
    OutData/=GammaCoef*GammaCoef*GammaCoef;
    if (OutData>MaxA) MaxA=OutData;
  }
  return(MaxA);
}

REAL_TYPE cDiffurMain::CalcMaxA_FirstPeak()
{
REAL_TYPE MaxA=0;
REAL_TYPE OutData;
  while(t<tMax)
  {
    Step();
    OutData=(CalcPower) ? hAdata[OutPoint].real()*hAdata[OutPoint].real()+
                           hAdata[OutPoint].imag()*hAdata[OutPoint].imag()
                           : abs(hAdata[OutPoint]);
    OutData/=GammaCoef*GammaCoef*GammaCoef;
    if (OutData>MaxA) MaxA=OutData; else
      if ( (MaxA>1.0)&&(MaxA*0.3>OutData) ) return(MaxA);
  }
  for(int i=0;i<nSteps+1;i++)
  {
    OutData=(CalcPower) ? hAdata[i].real()*hAdata[i].real()+
                           hAdata[i].imag()*hAdata[i].imag()
                           : abs(hAdata[i]);
    OutData/=GammaCoef*GammaCoef*GammaCoef;
    if (OutData>MaxA) MaxA=OutData;
    if (OutData>MaxA) MaxA=OutData; else
      if ( (MaxA>1.0)&&(MaxA*0.3>OutData) ) return(MaxA);
  }
  return(MaxA);
}

void cDiffurMain::PercentStep(ofstream *hOutFile)
{
  REAL_TYPE StopT=tMax*(PercentageCounter+1)*0.01;
  REAL_TYPE ResRm1=1.0/ResRange;
  REAL_TYPE ResT=(tMax_Res)*(ResIndex+1)*ResRm1;
  REAL_TYPE OutData;
  int n0;
  float dx;
  while(t<StopT)
  {
    Step();
    OutData=(CalcPower) ? hAdata[OutPoint].real()*hAdata[OutPoint].real()+
                           hAdata[OutPoint].imag()*hAdata[OutPoint].imag()
                           : abs(hAdata[OutPoint]);
    OutData/=GammaCoef*GammaCoef*GammaCoef;
    if (hOutFile) *hOutFile << t << "  " << OutData << endl;
    while(t>=ResT)
    {
      if(ResIndex<600)
      {
        if(t>=ResT+(tMax_Res)*ResRm1)
        {
          backOutData=backOutData
            +(OutData-backOutData)*(tMax_Res)/((t-ResT)*ResRange+(tMax_Res));
          hResDataA[ResIndex]=backOutData;
        } else hResDataA[ResIndex]=backOutData=OutData;
      }
      ResIndex++;
      ResT=(tMax_Res)*(ResIndex+1)*ResRm1;
    }
  }
  PercentageCounter++;
  if( (PercentageCounter==100) && (tMax_Res>tMax) )
  {
    AfterLastStep();
    for(int i=0;i<nSteps+1;i++)
    {
      OutData=(CalcPower) ? hAdata[i].real()*hAdata[i].real()+
                           hAdata[i].imag()*hAdata[i].imag()
                           : abs(hAdata[i]);
      OutData/=GammaCoef*GammaCoef*GammaCoef;
      hAdata[i]=complex<REAL_TYPE>(OutData,0);
      hIdata[i]=0;
      t+=dt;
      if (hOutFile) *hOutFile << t << "  " << real(hAdata[i]) << endl;
    }
    FillSplineArray(0);
    for(int j=ResIndex;j<ResRange;j++)
    {
      n0=(j-ResIndex)*nSteps/(ResRange-ResIndex);
      dx=(tMax_Res)*(j+1)*ResRm1-tMax-n0*dt;
      hResDataA[j]=real(hAdata[n0])+dx*(hSplineCoefB[n0]
                     +dx*(hSplineCoefC[n0]+dx*hSplineCoefD[n0]));
    }
  }
}

void cDiffurMain_LBV::AfterLastStep()
{
  complex<REAL_TYPE> tmpcplx;
  for(int i=0;i<nSteps/2;i++)
  {
    tmpcplx=hAdata[i];
    hAdata[i]=hAdata[nSteps-i];
    hAdata[nSteps-i]=tmpcplx;
  }
}

REAL_TYPE cDiffurMain::kappa(REAL_TYPE z)
{
  if(Z0<0.000001) return( Ckappa );
  if(!bKappaLinear) return( ( Z0>z) ? Ckappa/sqrt(p*p+z*(1-p*p)*Z0_inv) : Ckappa );
  if(Ckappa*(1.0+p*(z-Z0))<0) return(0);
  return( ( Z0>z) ? Ckappa*(1.0+p*(z-Z0)) : Ckappa );
}

void cDiffurMain::Step()
{
  complex<REAL_TYPE> TempI,TempINew; // for Runge-Kutt's f(new Els)
  complex<REAL_TYPE> TempA; // == (end of space) ? 0 : ...
  REAL_TYPE RungeCoef=1.0/6.0;
  FillSplineArray(0);    // Fill the hInterpolated array
  FillSplineArray(1);
  InitEls();
  j=0;
  int ntmpMovingEls;
  if(hDumpStorage) hDumpStorage->DumpDiffur(hAdata,hIdata,hEldata);
  hIdata1[0]=Calc_I(hEldata);
  for(j=0;j<nSteps;j++)  // The main loop j == z coord
  {
    hAdata1[j]=hAdata[j+1]+(REAL_TYPE)0.5*dt*(hIdata1[j]+hIdata[j+1]);
    DumpRequest.FillReq(this);
   /*       // test - liner interpolation

    if (j<nSteps-1)
     InterpolatedL=(REAL_TYPE)0.5*(hAdata[j+1]+(REAL_TYPE)0.5*dt*hIdata[j+1])
                 +(REAL_TYPE)0.5*(hAdata[j+2]+(REAL_TYPE)0.5*dt*hIdata[j+2]);
     else InterpolatedL=0;
   */

    memcpy(hElNew,hEldata,nEls*2*sizeof(REAL_TYPE));
    TempINew=0;
    ntmpMovingEls=nMovingEls;
    for(int i=0;i<nEls;i++)  // k1 calculating
      if (hElPresent[i]) {
        hElk1[2*i]=dt*RightTH0( (hElNew+2*i) );
        hElk1[2*i+1]=dt*RightTH1( (hElNew+2*i), kappa(dt*((REAL_TYPE)j))*hAdata1[j] );
        hElNew[2*i]+=0.5*hElk1[2*i];
        hElNew[2*i+1]+=0.5*hElk1[2*i+1];
        TempINew+=exp(complex<REAL_TYPE>(0,-hElNew[2*i]));
        CheckIfElStopped(hElNew+2*i,i);
      }

    if (j==(nSteps-1))
    {
      SplineRequest=1;
      FillSplineArray(0);
      FillSplineArray(1);
      SplineRequest=0;
      for(int i=0;i<nEls;i++)
        if (hElPresent[i]) {
          hElk2[2*i]=dt*RightTH0( (hElNew+2*i) );
          hElk2[2*i+1]=dt*RightTH1( (hElNew+2*i), kappa(dt*((REAL_TYPE)j+0.5))*hInterpolated[0] );
          hElNew[2*i]=hEldata[2*i]+0.5*hElk2[2*i];
          hElNew[2*i+1]=hEldata[2*i+1]+0.5*hElk2[2*i+1];
          CheckIfElStopped(hElNew+2*i,i);
        }
      ntmpMovingEls=nMovingEls;
      for(int i=0;i<nEls;i++)
        if (hElPresent[i]) {
          hElk3[2*i]=dt*RightTH0( (hElNew+2*i) );
          hElk3[2*i+1]=dt*RightTH1( (hElNew+2*i),kappa(dt*((REAL_TYPE)j+0.5))*hInterpolated[0] );
          hElNew[2*i]=hEldata[2*i]+hElk3[2*i];
          hElNew[2*i+1]=hEldata[2*i+1]+hElk3[2*i+1];
          CheckIfElStopped(hElNew+2*i,i);
        }
    } else
    {
                             // k2 calculating
      if(ntmpMovingEls)
        TempI=-(REAL_TYPE)2.0*G*kappa(dt*((REAL_TYPE)j+0.5))*TempINew/(REAL_TYPE)ntmpMovingEls;//Calc_I(hElNew);
      else TempI=0;
      TempINew=0;
      ntmpMovingEls=nMovingEls;
      for(int i=0;i<nEls;i++)
        if (hElPresent[i]) {
          hElk2[2*i]=dt*RightTH0( (hElNew+2*i) );
          hElk2[2*i+1]=dt*RightTH1( (hElNew+2*i),
            kappa(dt*((REAL_TYPE)j+0.5))*(hInterpolated[j+1]+(REAL_TYPE)0.5*dt*TempI) );
          hElNew[2*i]=hEldata[2*i]+0.5*hElk2[2*i];
          hElNew[2*i+1]=hEldata[2*i+1]+0.5*hElk2[2*i+1];
          TempINew+=exp(complex<REAL_TYPE>(0,-hElNew[2*i]));
          CheckIfElStopped(hElNew+2*i,i);
        }
                             // k3 calculating
      if(ntmpMovingEls)
        TempI=-(REAL_TYPE)2.0*G*kappa(dt*((REAL_TYPE)j+0.5))*TempINew/(REAL_TYPE)ntmpMovingEls;//Calc_I(hElNew);
      else TempI=0;
      TempINew=0;
      ntmpMovingEls=nMovingEls;
      for(int i=0;i<nEls;i++)
        if (hElPresent[i]) {
          hElk3[2*i]=dt*RightTH0( (hElNew+2*i) );
          hElk3[2*i+1]=dt*RightTH1( (hElNew+2*i),
            kappa(dt*((REAL_TYPE)j+0.5))*(hInterpolated[j+1]+(REAL_TYPE)0.5*dt*TempI) );
          hElNew[2*i]=hEldata[2*i]+hElk3[2*i];
          hElNew[2*i+1]=hEldata[2*i+1]+hElk3[2*i+1];
          TempINew+=exp(complex<REAL_TYPE>(0,-hElNew[2*i]));
          CheckIfElStopped(hElNew+2*i,i);
        }
    }
                             // k4 calculating
    if(ntmpMovingEls)
      TempI=-(REAL_TYPE)2.0*G*kappa(dt*((REAL_TYPE)j+1.0))*TempINew/(REAL_TYPE)ntmpMovingEls;//Calc_I(hElNew);
    else TempI=0;
    TempINew=0;
    TempA= ( j==(nSteps-1) ) ? complex<REAL_TYPE>(0,0)
           : (hAdata[j+2]+(REAL_TYPE)0.5*dt*(TempI+hIdata[j+2]));
    ntmpMovingEls=nMovingEls;
    for(int i=0;i<nEls;i++)
      if (hElPresent[i]) {
        hElk4[2*i]=dt*RightTH0( (hElNew+2*i) );
        hElk4[2*i+1]=dt*RightTH1( (hElNew+2*i), kappa(dt*((REAL_TYPE)j+1.0))*TempA );
        hElNew[2*i]=hEldata[2*i]
          +RungeCoef*(hElk1[2*i]+2*hElk2[2*i]+2*hElk3[2*i]+hElk4[2*i]);
        hElNew[2*i+1]=hEldata[2*i+1]
          +RungeCoef*(hElk1[2*i+1]+2*hElk2[2*i+1]+2*hElk3[2*i+1]+hElk4[2*i+1]);
        TempINew+=exp(complex<REAL_TYPE>(0,-hElNew[2*i]));
        CheckIfElStopped(hElNew+2*i,i);
      }
    if(ntmpMovingEls)
      hIdata1[j+1]=-(REAL_TYPE)2.0*G*kappa(dt*((REAL_TYPE)j+1.0))*TempINew/(REAL_TYPE)ntmpMovingEls;
    else hIdata1[j+1]=0;
    memcpy(hEldata,hElNew,nEls*2*sizeof(REAL_TYPE));
    if(hDumpStorage) hDumpStorage->DumpDiffur((hAdata+j),(hIdata+j),hEldata);
  }
  j=nSteps;
//  hIdata1[nSteps]=-(REAL_TYPE)2.0*G*TempINew/(REAL_TYPE)nEls;//Calc_I(hElNew);
//  hIdata1[nSteps]=Calc_I(hElNew);
  hAdata1[nSteps]=InputField(t+dt);
  memcpy(hIdata,hIdata1,(nSteps+1)*sizeof(complex<REAL_TYPE>));
  memcpy(hAdata,hAdata1,(nSteps+1)*sizeof(complex<REAL_TYPE>));

  REAL_TYPE rAvgGamma=0;
  for(int i=0;i<nEls;i++){ rAvgGamma+=hEldata[2*i+1]; }
  rDebug[0]=rAvgGamma=rAvgGamma/ntmpMovingEls;

  double ElDroppedFrac=(double)(nEls-nMovingEls)/(double)nEls;
  ElDroppedSum+=nEls-nMovingEls;
  if(ElDroppedFrac>ElDroppedMaxFrac) ElDroppedMaxFrac=ElDroppedFrac;
  
  DumpRequest.FillReq(this);
  t+=dt; nCStep++;
}

void cDiffurMain_LBV::Step()
{
  G=1;
  complex<REAL_TYPE> TempI,TempINew; // for Runge-Kutt's f(new Els)
  complex<REAL_TYPE> TempA; // == (end of space) ? 0 : ...
  REAL_TYPE RungeCoef=1.0/6.0;
  FillSplineArray(0);    // Fill the hInterpolated array
  FillSplineArray(1);
  InitEls();
  j=0;
  int ntmpMovingEls;
  if(hDumpStorage) hDumpStorage->DumpDiffur(hAdata,hIdata,hEldata);
  hIdata1[0]=Calc_I(hEldata);
  hAdata1[0]=InputField(t);
  for(j=0;j<nSteps;j++)  // The main loop j == z coord
    {
    DumpRequest.FillReq(this);
    memcpy(hElNew,hEldata,nEls*2*sizeof(REAL_TYPE));
    TempINew=0;
    ntmpMovingEls=nMovingEls;
    for(int i=0;i<nEls;i++)  // k1 calculating
      if (hElPresent[i]) {
        hElk1[2*i]=dt*RightTH0( (hElNew+2*i) );
        hElk1[2*i+1]=dt*RightTH1( (hElNew+2*i), kappa( L-dt*((REAL_TYPE)j))*hAdata1[j] );
        hElNew[2*i]+=0.5*hElk1[2*i];
        hElNew[2*i+1]+=0.5*hElk1[2*i+1];
        TempINew+=exp(complex<REAL_TYPE>(0,-hElNew[2*i]));
        CheckIfElStopped(hElNew+2*i,i);
      }
                             // k2 calculating
    if(ntmpMovingEls)
      TempI=-(REAL_TYPE)2.0*G*kappa(L-dt*((REAL_TYPE)j+0.5))*TempINew/(REAL_TYPE)ntmpMovingEls;//Calc_I(hElNew);
    else TempI=0;
    if (j==0) TempA=InputField(t)+(REAL_TYPE)0.25*dt*(TempI+hIdata[0]);
      else TempA=(hInterpolated[j-1]+(REAL_TYPE)0.5*dt*TempI);
    TempINew=0;
    ntmpMovingEls=nMovingEls;
    for(int i=0;i<nEls;i++)
      if (hElPresent[i]) {
        hElk2[2*i]=dt*RightTH0( (hElNew+2*i) );
        hElk2[2*i+1]=dt*RightTH1( (hElNew+2*i),kappa(L-dt*((REAL_TYPE)j+0.5))*TempA );
        hElNew[2*i]=hEldata[2*i]+0.5*hElk2[2*i];
        hElNew[2*i+1]=hEldata[2*i+1]+0.5*hElk2[2*i+1];
        TempINew+=exp(complex<REAL_TYPE>(0,-hElNew[2*i]));
        CheckIfElStopped(hElNew+2*i,i);
      }
                             // k3 calculating
    if(ntmpMovingEls)
      TempI=-(REAL_TYPE)2.0*G*kappa(L-dt*((REAL_TYPE)j+0.5))*TempINew/(REAL_TYPE)ntmpMovingEls;//Calc_I(hElNew);
    else TempI=0;
    if (j==0) TempA=InputField(t)+(REAL_TYPE)0.25*dt*(TempI+hIdata[0]);
      else TempA=(hInterpolated[j-1]+(REAL_TYPE)0.5*dt*TempI);
    TempINew=0;
    ntmpMovingEls=nMovingEls;
    for(int i=0;i<nEls;i++)
      if (hElPresent[i]) {
        hElk3[2*i]=dt*RightTH0( (hElNew+2*i) );
        hElk3[2*i+1]=dt*RightTH1( (hElNew+2*i),kappa(L-dt*((REAL_TYPE)j+0.5))*TempA );
        hElNew[2*i]=hEldata[2*i]+hElk3[2*i];
        hElNew[2*i+1]=hEldata[2*i+1]+hElk3[2*i+1];
        TempINew+=exp(complex<REAL_TYPE>(0,-hElNew[2*i]));
        CheckIfElStopped(hElNew+2*i,i);
      }
                             // k4 calculating
    if(ntmpMovingEls)
      TempI=-(REAL_TYPE)2.0*G*kappa(L-dt*((REAL_TYPE)j+1.0))*TempINew/(REAL_TYPE)ntmpMovingEls;//Calc_I(hElNew);
    else TempI=0;
    TempINew=0;
    TempA= hAdata[j]+(REAL_TYPE)0.5*dt*(TempI+hIdata[j]);
    ntmpMovingEls=nMovingEls;
    for(int i=0;i<nEls;i++)
      if (hElPresent[i]) {
        hElk4[2*i]=dt*RightTH0( (hElNew+2*i) );
        hElk4[2*i+1]=dt*RightTH1( (hElNew+2*i), kappa(L-dt*((REAL_TYPE)j+1.0))*TempA );
        hElNew[2*i]=hEldata[2*i]
          +RungeCoef*(hElk1[2*i]+2*hElk2[2*i]+2*hElk3[2*i]+hElk4[2*i]);
        hElNew[2*i+1]=hEldata[2*i+1]
          +RungeCoef*(hElk1[2*i+1]+2*hElk2[2*i+1]+2*hElk3[2*i+1]+hElk4[2*i+1]);
        TempINew+=exp(complex<REAL_TYPE>(0,-hElNew[2*i]));
        CheckIfElStopped(hElNew+2*i,i);
      }
    if(ntmpMovingEls)
      hIdata1[j+1]=-(REAL_TYPE)2.0*G*kappa(L-dt*((REAL_TYPE)j+1.0))*TempINew/(REAL_TYPE)ntmpMovingEls;
    else hIdata1[j+1]=0;
    hAdata1[j+1]=hAdata[j]+(REAL_TYPE)0.5*dt*(hIdata1[j+1]+hIdata[j]);
    memcpy(hEldata,hElNew,nEls*2*sizeof(REAL_TYPE));
    if(hDumpStorage) hDumpStorage->DumpDiffur( (hAdata+j),(hIdata+j),hEldata);
  }
  j=nSteps;
//  hIdata1[nSteps]=-(REAL_TYPE)2.0*G*TempINew/(REAL_TYPE)nEls;//Calc_I(hElNew);
//  hIdata1[nSteps]=Calc_I(hElNew);
  memcpy(hIdata,hIdata1,(nSteps+1)*sizeof(complex<REAL_TYPE>));
  memcpy(hAdata,hAdata1,(nSteps+1)*sizeof(complex<REAL_TYPE>));

  double ElDroppedFrac=(double)(nEls-nMovingEls)/(double)nEls;
  ElDroppedSum+=nEls-nMovingEls;
  if(ElDroppedFrac>ElDroppedMaxFrac) ElDroppedMaxFrac=ElDroppedFrac;
  
  DumpRequest.FillReq(this);
  t+=dt; nCStep++;
}

complex<REAL_TYPE> cDiffurMain::Calc_I(REAL_TYPE * hEls)
{
  complex<REAL_TYPE> Res(0,0);
  for(int i=0;i<nEls;i++)
    if (hElPresent[i]) Res+=exp(complex<REAL_TYPE>(0,-hEls[2*i]));
  if(nMovingEls)
    return( -(REAL_TYPE)2.0*G*kappa(j*dt)*Res/(REAL_TYPE)nMovingEls );
  else return( complex<REAL_TYPE>(0,0) );
}

void cDiffurMain::StoreData(REAL_TYPE * hData,int type)
{
  switch(type)
  {
    case 0:
      for(int i=0;i<nSteps+1;i++)
        hData[i]=abs(hAdata[i]);
    break;
    case 1:
      for(int i=0;i<nSteps+1;i++)
        hData[i]=abs(hIdata[i]);
  }
}

void cDiffurMain::FillSplineArray(int ofs)
{
  REAL_TYPE dt1=1.0/dt;
  REAL_TYPE *hData=(REAL_TYPE *)hInterpolated;
  int tmpNSteps=nSteps;
  for(int i=0;i<nSteps+1;i++)
    hData[2*i+ofs]=
       ((REAL_TYPE *)hAdata)[2*i+ofs]+0.5*dt*((REAL_TYPE *)hIdata)[2*i+ofs];
  if(SplineRequest)
  {
    hData[0]=0;
    hData[1]=((REAL_TYPE *)hAdata1)[nSteps-1];
    hData[2]=((REAL_TYPE *)hAdata1)[nSteps-2];
    hData[3]=((REAL_TYPE *)hAdata1)[nSteps-3];
    hData[4]=((REAL_TYPE *)hAdata1)[nSteps-4];
    nSteps=4;
  }
  hSplineCoefD[0]=dt;
  hSplineCoefC[1]=(hData[2+ofs]-hData[ofs])*dt1;
  for(int i=1;i<nSteps;i++)
  {
    hSplineCoefD[i]=dt;
    hSplineCoefB[i]=2.0*(hSplineCoefD[i-1]+hSplineCoefD[i]);
    hSplineCoefC[i+1]=(hData[2*i+2+ofs]-hData[2*i+ofs])*dt1;
    hSplineCoefC[i]=hSplineCoefC[i+1]-hSplineCoefC[i];
  }
  hSplineCoefB[0]=hSplineCoefB[nSteps]=-dt;
  hSplineCoefC[0]=(hSplineCoefC[2]-hSplineCoefC[1])*0.5*0.333333333333;
  hSplineCoefC[nSteps]=(hSplineCoefC[nSteps-1]-hSplineCoefC[nSteps-2])*0.5*0.33333333333;
  for(int i=1;i<nSteps+1;i++)
  {
    REAL_TYPE tmp=dt/hSplineCoefB[i-1];
    hSplineCoefB[i]=hSplineCoefB[i]-tmp*dt;
    hSplineCoefC[i]=hSplineCoefC[i]-tmp*hSplineCoefC[i-1];
  }
  hSplineCoefC[nSteps]=hSplineCoefC[nSteps]/hSplineCoefB[nSteps];
  for(int i=nSteps-1;i>=0;i--)
    hSplineCoefC[i]=(hSplineCoefC[i]-dt*hSplineCoefC[i+1])/hSplineCoefB[i];
  hSplineCoefB[nSteps]=(hData[2*nSteps+ofs]-hData[2*nSteps-2+ofs])*dt1
                 +(hSplineCoefC[nSteps-1]+2.0*hSplineCoefC[nSteps])*dt;
  for(int i=0;i<nSteps;i++)
  {
    hSplineCoefB[i]=(hData[2*i+2+ofs]-hData[2*i+ofs])*dt1
                   -(hSplineCoefC[i+1]+2.0*hSplineCoefC[i])*dt;
    hSplineCoefD[i]=(hSplineCoefC[i+1]-hSplineCoefC[i])*dt1;
    hSplineCoefC[i]=3.0*hSplineCoefC[i];
  }
  hSplineCoefC[nSteps]=3.0*hSplineCoefC[nSteps];
  hSplineCoefD[nSteps]=hSplineCoefD[nSteps-1];

  dt1=0.5*dt;
  for(int i=0;i<nSteps+1;i++)
    hData[2*i+ofs]=hData[2*i+ofs]
          +dt1*(hSplineCoefB[i]+dt1*(hSplineCoefC[i]+dt1*hSplineCoefD[i]));
    nSteps=tmpNSteps;
}

REAL_TYPE cDiffurMain::RightTH0(REAL_TYPE *hEl)
{ return(hEl[1]); }

REAL_TYPE cDiffurMain::RightTH1(REAL_TYPE *hEl, complex<REAL_TYPE> A)
{ return( real(-A*exp(complex<REAL_TYPE>(0,hEl[0]))) ); }


void cDiffurMain::InitEls()
{
  for(int i=0;i<nEls;i++)
    {
      REAL_TYPE Theta0=2.0*M_PI*i/nEls;
      hEldata[2*i]=Theta0+r()*cos(Theta0);
//      hEldata[2*i+1]= (t>=delta) ? -delta-alpha*(t-Z0) : 0;
      hEldata[2*i+1]= (t>=delta) ? -alpha*(t-delta) : 0;
//      hEldata[2*i+1]= -delta*t*t*t*t-alpha*t*t*t;
      hElPresent[i]=1;
    }
  nMovingEls=nEls;
}

REAL_TYPE cDiffurMain_Gamma::RightTH0(REAL_TYPE *hEl)
{
  return( CPierce_inv*(sqrt( hEl[1]*hEl[1]/(hEl[1]*hEl[1]-1))-beta0_inv) );
}

REAL_TYPE cDiffurMain_Gamma::RightTH1(REAL_TYPE *hEl, complex<REAL_TYPE> A)
{ return( -RTH1Coef*real(A*exp(complex<REAL_TYPE>(0,hEl[0]))) ); }


void cDiffurMain_Gamma::InitEls()
{
  for(int i=0;i<nEls;i++)
    {
      REAL_TYPE Theta0=2.0*M_PI*i/nEls;
      hEldata[2*i]=Theta0+r()*cos(Theta0);
      hEldata[2*i+1]=(t>=delta) ? gamma0-RTH1Coef*alpha*(t-delta) : gamma0;
      hElPresent[i]=1;
    }
  nMovingEls=nEls;
}

void cDiffurMain::SecondaryInit_Gamma(REAL_TYPE pCPierce,REAL_TYPE pgamma0)
{
  CPierce=pCPierce; CPierce_inv=1.0/CPierce;
  gamma0=pgamma0;
  REAL_TYPE gamma0_inv=1.0/gamma0;
  beta0=sqrt(1.0-gamma0_inv*gamma0_inv);
  beta0_inv=1.0/beta0;
  RTH1Coef=-CPierce*gamma0*gamma0*gamma0*beta0*beta0*beta0;
}

void cDiffurMain_Gamma::SecondaryInit_Gamma(REAL_TYPE pCPierce,REAL_TYPE pgamma0)
{
  CPierce=pCPierce; CPierce_inv=1.0/CPierce;
  gamma0=pgamma0;
  REAL_TYPE gamma0_inv=1.0/gamma0;
  beta0=sqrt(1.0-gamma0_inv*gamma0_inv);
  beta0_inv=1.0/beta0;
  RTH1Coef=-CPierce*gamma0*gamma0*gamma0*beta0*beta0*beta0;
}

void cDiffurMain_Gamma::CheckIfElStopped(REAL_TYPE *hEl,int i)
{
  if ( hEl[1]< 1.000001 ) { hElPresent[i]=0; nMovingEls--; }
}

REAL_TYPE cDiffurMain_LBV_Gamma::RightTH0(REAL_TYPE *hEl)
{
  return( CPierce_inv*(sqrt( hEl[1]*hEl[1]/(hEl[1]*hEl[1]-1))-beta0_inv) );
}

REAL_TYPE cDiffurMain_LBV_Gamma::RightTH1(REAL_TYPE *hEl, complex<REAL_TYPE> A)
{ return( RTH1Coef*real(A*exp(complex<REAL_TYPE>(0,hEl[0]))) ); }


void cDiffurMain_LBV_Gamma::InitEls()
{
  for(int i=0;i<nEls;i++)
    {
      REAL_TYPE Theta0=2.0*M_PI*i/nEls;
      hEldata[2*i]=Theta0+r()*cos(Theta0);
      hEldata[2*i+1]=hEldata[2*i+1]=(t>=delta) ? gamma0+RTH1Coef*alpha*(t-delta) : gamma0;
      hElPresent[i]=1;
    }
  nMovingEls=nEls;
}

void cDiffurMain_LBV_Gamma::SecondaryInit_Gamma(REAL_TYPE pCPierce,REAL_TYPE pgamma0)
{
  CPierce=pCPierce; CPierce_inv=1.0/CPierce;
  gamma0=pgamma0;
  REAL_TYPE gamma0_inv=1.0/gamma0;
  beta0=sqrt(1.0-gamma0_inv*gamma0_inv);
  beta0_inv=1.0/beta0;
  RTH1Coef=-CPierce*gamma0*gamma0*gamma0*beta0*beta0*beta0;
}

void cDiffurMain_LBV_Gamma::CheckIfElStopped(REAL_TYPE *hEl,int i)
{
  if ( hEl[1]< 1.000001 ) { hElPresent[i]=0; nMovingEls--; }
}


int cDumpRequest::FillReq(cDiffurMain *hDiffur)
{
  int q=0;
  cDumpRequest *tmphNext;
  if (hNext) q=hNext->FillReq(hDiffur);
  if (q==-1)
  {
    tmphNext=hNext;
    hNext=hNext->hNext;
    tmphNext->hNext=NULL;
    delete tmphNext;
  }
  return(Fill(hDiffur));
}

int cDumpRequest_ElU::Fill(cDiffurMain *hDiffur)
{
  if(hDiffur->t>=TReq)
  {
    *hOutFile << hDiffur->L*hDiffur->j/hDiffur->nSteps << "  ";
    for(int i=0;i<hDiffur->nEls;i++)
    {
      *hOutFile << hDiffur->hEldata[2*i+1] << "  ";
    }
    *hOutFile << endl;
    if(hDiffur->j==hDiffur->nSteps) return(-1); else return(0);
  }
  return(0);
}

int cDumpRequest_ThetaPhaseFrame::Fill(cDiffurMain *hDiffur)
{
  if( (hDiffur->t>=TReq)&&
      (hDiffur->j>=zReq/hDiffur->L*(REAL_TYPE)hDiffur->nSteps) )
  {
    for(int i=0;i<hDiffur->nEls;i++)
    {
      *hOutFile << hDiffur->hEldata[2*i] << "  " << hDiffur->hEldata[2*i+1] << endl;
    }
    *hOutFile << hDiffur->hEldata[0]+2*M_PI << "  " << hDiffur->hEldata[1] << endl;
    return(-1);
  }
  return(0);
}

