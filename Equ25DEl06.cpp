//---------------------------------------------------------------------------


#pragma hdrstop

#include "Equ25DEl06.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

//---------------------------------------------------------------------------


#pragma hdrstop

#include "Equ25DEl04.h"
#include <fstream.h>
#include <vcl\math.hpp>

//---------------------------------------------------------------------------

#pragma package(smart_init)

void cEquVerticalASpeedDistr::InitReport(cDifEquReport *hReport) {
  int Lines[EQU_MAX_GRAPHS];
  for(int i=0;i<EQU_MAX_GRAPHS;i++) Lines[i]=1;
  Lines[0]=1; Lines[1]=1; Lines[2]=1; Lines[3]=1;
  hReport->LoadText(hReportText);
  hReport->SetGraphSize(EQU_GRAPH_SIZE,Lines);
  hReport->SetMaxTo1();
}
void cEquVerticalASpeedDistr::Report(cDifEquReport *hReport) {
/*  long double *hGamma = new long double[InitData.nSteps+1];
  for(int i=0;i<=InitData.nSteps;i++) {
    hGamma[i]=0;
    for(int j=0;j<InitData.nEls;j++) {
      if(hElPresent[i*InitData.nEls+j])
        hGamma[i]-=hEl[i*InitData.nEls*2+j*2+1];
    }
    hGamma[i]/=InitData.nEls;
  } */
  hReport->Graph[0].GetDataAbs(0,InitData.nSteps+1,hAdata);
  hReport->Graph[1].GetDataAbs(0,InitData.nSteps+1,hIdata);
  hReport->Graph[2].GetDataAbs(0,InitData.nSteps+1,hIdata);
  hReport->Graph[3].GetDataAbs(0,InitData.nSteps+1,hIdata);
  hReport->Graph[4].GetDataAbs(0,InitData.nSteps+1,hIdata);
/*  hReport->Graph[2].GetData(0,InitData.nSteps+1,hAvgP);
  hReport->Graph[2].Max = hReport->Graph[2].Max > -hReport->Graph[2].Min ?
                          hReport->Graph[2].Max : -hReport->Graph[2].Min;*/
/*  hReport->Graph[3].GetData(0,InitData.ElTime/InitData.dt,hEAData);
  hReport->Graph[4].GetData(0,InitData.ElTime/InitData.dt,hRandomData);*/
  hReport->nVals=2;
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
  hReport->ValText[1]="Max I "; hReport->Val[1]=hReport->Graph[1].Max;
  /*
  hReport->ValText[2]="Max Average P "; hReport->Val[2]=hReport->Graph[2].Max;
  hReport->ValText[3]="Max W "; hReport->Val[3]=maxW;
  hReport->ValText[4]="Max TS "; hReport->Val[4]=maxdeltaTS;
  hReport->ValText[5]="integral A "; hReport->Val[5]=intA;
  hReport->ValText[6]="integral W "; hReport->Val[6]=intW;
  hReport->ValText[7]="integral A-W "; hReport->Val[7]=intA-intW;
  hReport->ValText[8]="max el energy ";
  hReport->ValText[9]="max random";
  hReport->Val[8]=hReport->Graph[3].Max > -hReport->Graph[3].Min ?
                          hReport->Graph[3].Max : hReport->Graph[3].Min;

  hReport->Graph[3].Max = hReport->Graph[3].Max > -hReport->Graph[3].Min ?
                          hReport->Graph[3].Max : -hReport->Graph[3].Min;
  hReport->Val[9]=hReport->Graph[4].Max > -hReport->Graph[4].Min ?
                          hReport->Graph[4].Max : hReport->Graph[4].Min;

  hReport->Graph[4].Max = hReport->Graph[4].Max > -hReport->Graph[4].Min ?
                          hReport->Graph[4].Max : -hReport->Graph[4].Min;
  */
//  hReport->Graph[2].Max=2;
}
void cEquVerticalASpeedDistr::FinalReport(cDifEquReport *hReport) {

  hReport->nVals=1; hReport->nGraphs=1;
  hReport->Graph[0].GetDataAbs(0,nStepNumber,hOutA);
  hReport->GraphText[0]="Out A";
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
/*  hReport->ValText[1]="integral A "; hReport->Val[1]=intA;
  hReport->ValText[2]="integral W "; hReport->Val[2]=intW;
  hReport->ValText[3]="integral A-W "; hReport->Val[3]=intA-intW;

  ofstream FileOut = ofstream("Out.dat",ios::out);
  for(int i=0;i<nStepNumber;i++) {
    FileOut << InitData.dt*i << " " << abs(hOutA[i]) << endl;
  }                                                              */
}

cEquVerticalASpeedDistr::cEquVerticalASpeedDistr() : cEqu25DElModel()
{
  hIdata=hAdata=hIdata1=hAdata1=NULL; hOutA=NULL; hElSpeedFracs=NULL;
  hEldata=hElk1=hElk2=NULL; hElk3=hElk4=NULL; hElNew=NULL;
  hNMovingEls=hElPresent=NULL; InitData.ElDirection=1;
  hSplineCoefBRe=hSplineCoefCRe=hSplineCoefDRe=NULL;
  hSplineCoefBIm=hSplineCoefCIm=hSplineCoefDIm=NULL;
  htmpSplineData=NULL;  hElSpeed=NULL;
  hInterpolated=NULL;
  InitData.nEls=InitData.nSteps=InitData.nSpeedFracs=0;
  InitData.tMax=0;
  nStepNumber=0;
}

cEquVerticalASpeedDistr::~cEquVerticalASpeedDistr()// : ~cDifEqu_DataBankIO()
{
  delete [] hIdata;  delete [] hAdata; delete hOutA;
  delete [] hIdata1;  delete [] hAdata1;  delete [] hElSpeed;
  delete [] hEldata;  delete [] hElk1; delete [] hElk2; delete [] hElNew;
  delete [] hElk3; delete [] hElk4; delete [] hElSpeedFracs;
  delete [] hElPresent; delete [] hNMovingEls; delete [] hElSpeed;
  delete [] hSplineCoefBRe; delete [] hSplineCoefCRe; delete [] hSplineCoefDRe;
  delete [] hSplineCoefBRe; delete [] hSplineCoefCRe; delete [] hSplineCoefDRe;
  delete [] hInterpolated; delete [] htmpSplineData;
  hIdata=hAdata=hIdata1=hAdata1=NULL; hEldata=NULL; hElSpeedFracs=NULL; hElSpeed=NULL;
  hSplineCoefBRe=hSplineCoefCRe=hSplineCoefDRe=NULL;
  hSplineCoefBIm=hSplineCoefCIm=hSplineCoefDIm=NULL;
  hElk1=hElk2=NULL;  hElk3=hElk4=NULL; hElSpeed=NULL; hElNew=NULL;
  hElPresent=hNMovingEls=NULL; hInterpolated=NULL;
  hDump=NULL; htmpSplineData=NULL;         
}
/*
cEquVerticalASpeedDistr::cEquVerticalASpeedDistr() : cEqu25DElModel()
{
  hI=hA=hAk1=hAk2=hABack=hIBack=NULL; hOutA=NULL;
  hEl=hElk1=hElk2=NULL; hAvgP=NULL; hEAData=NULL; hRandomData=NULL;
  hEl2=hElk12=hElk22=NULL;  hdeltaTS=NULL; hNormalP=NULL;
  hNMovingEls=hElPresent=NULL; InitData.ElDirection=1;
  InitData.nEls=InitData.nSteps=0; EATempIData=NULL;
  InitData.tMax=0;
  nStepNumber=0;
}

cEquVerticalASpeedDistr::~cEquVerticalASpeedDistr()// : ~cDifEqu_DataBankIO()
{
  delete [] hI;  delete [] hA; delete [] hOutA;
  delete [] hAk1;  delete [] hAk2; delete [] hIBack;
  delete [] hI;  delete [] hA; delete [] hABack;
  delete [] hEl;  delete [] hElk1; delete [] hElk2;
  delete [] hElPresent; delete [] hNMovingEls;
  delete [] hEl2;  delete [] hElk12; delete [] hElk22;
  delete [] hAvgP; delete [] hdeltaTS; delete [] hNormalP;
  delete [] hEAData;  delete [] EATempIData; delete [] hRandomData;
  hI=hA=hAk1=hAk2=hABack=hIBack=NULL;
  hEl=hElk1=hElk2=NULL; hAvgP=NULL;
  hEl2=hElk12=hElk22=NULL; hdeltaTS=NULL;
  hElPresent=hNMovingEls=NULL;
  hDump=NULL; hEAData=NULL; EATempIData=NULL; hRandomData=NULL;
} */

void cEquVerticalASpeedDistr::MakeNormalP() {
  if(!hNormalP) hNormalP = new long double[InitData.nSteps+1];
  for(int i=0;i<=InitData.nSteps;i++)
    if(InitData.deltaI*InitData.deltaI<0.00000000001) hNormalP[i]=InitData.dt*i;
    else {
      hNormalP[i]=2*sin(0.5*InitData.deltaI*InitData.dt*i)/InitData.deltaI;
    }
}

void cEquVerticalASpeedDistr::LoadInitData(cDifEqu_InitDataRec *hInitData)
{
  cDifEqu_InitDataRec_25DElModel *hInit=(cDifEqu_InitDataRec_25DElModel *)hInitData;
  int nTStepsNew=hInit->tMax/hInit->dt;
  int nTSteps=InitData.tMax/InitData.dt;
  if ( nTSteps!=nTStepsNew ) {
    delete [] hOutA; hOutA=new TCplxLong[nTStepsNew+1];
  }
  if( InitData.nSteps!=hInit->nSteps)  {
    delete [] hIdata;  delete [] hAdata;
    delete [] hIdata1;  delete [] hAdata1;
    delete [] hNMovingEls; delete [] hInterpolated;
    delete [] hSplineCoefBRe; delete [] hSplineCoefCRe; delete [] hSplineCoefDRe;
    delete [] hSplineCoefBRe; delete [] hSplineCoefCRe; delete [] hSplineCoefDRe;
    delete [] htmpSplineData;
    hSplineCoefBRe=new long double[hInit->nSteps+1];
    hSplineCoefCRe=new long double[hInit->nSteps+1];
    hSplineCoefDRe=new long double[hInit->nSteps+1];
    hSplineCoefBIm=new long double[hInit->nSteps+1];
    hSplineCoefCIm=new long double[hInit->nSteps+1];
    hSplineCoefDIm=new long double[hInit->nSteps+1];
    htmpSplineData=new long double[hInit->nSteps+1];
    hIdata=new TCplxLong[hInit->nSteps+1];
    hAdata=new TCplxLong[hInit->nSteps+1];
    hIdata1=new TCplxLong[hInit->nSteps+1];
    hAdata1=new TCplxLong[hInit->nSteps+1];
    hInterpolated = new TCplxLong[hInit->nSteps+1];
    hNMovingEls = new int[hInit->nSteps+1];
  }
  if( (hInit->nEls!=InitData.nEls) || (hInit->nSteps!=InitData.nSteps)
                                   || (hInit->nSpeedFracs!=InitData.nSpeedFracs) )
  {
    delete [] hEldata;    delete [] hElNew;
    delete [] hElk1;  delete [] hElk2;
    delete [] hElk3;  delete [] hElk4;
    delete [] hElPresent; delete [] hElSpeedFracs; delete [] hElSpeed;
    hEldata=new long double[2*hInit->nEls*hInit->nSpeedFracs];
    //hElDataNew=new REAL_TYPE[(hInit->nSteps+1)*2*hInit->nEls];
    hElk1=new long double[2*hInit->nEls*hInit->nSpeedFracs];
    hElk2=new long double[2*hInit->nEls*hInit->nSpeedFracs];
    hElk3=new long double[2*hInit->nEls*hInit->nSpeedFracs];
    hElk4=new long double[2*hInit->nEls*hInit->nSpeedFracs];
    hElNew=new long double[2*hInit->nEls*hInit->nSpeedFracs];
    hElPresent = new int[2*hInit->nEls*hInit->nSpeedFracs];
    hElSpeedFracs = new long double [hInit->nSpeedFracs];
    hElSpeed = new long double [hInit->nSpeedFracs];
  }
  InitData=*hInit; nStepNumber=0;
  t=0; tMax=InitData.tMax; dt=InitData.dt;

  for(int i=0;i<=InitData.nSteps;i++) {
    hAdata[i]=hIdata[i]=0;
  }
  InitKappa();
  MakeSpeedDistibution();
}

long double cEquVerticalASpeedDistr::SimplePhaseDeSync(int nZStep) {
  long double tmp=hNormalP[nZStep]*InitData.AI.real();
  return(tmp*tmp);
}

long double cEquVerticalASpeedDistr::EADesync(int nZStep) {
  long double tmp=(InitData.dt*nZStep-InitData.EADeltaZ0)*InitData.EAkz;
  return(tmp > 0 ? tmp:0);
}

long double cEquVerticalASpeedDistr::SimplePhaseKappa(int nZStep) {
TCplxLong tmp;
long double qq1,qq2;
static nZZ;
  if((InitData.UseSimplePhaseModel)&&
     (InitData.Ai_kappa>0)
     /*&&(InitData.dt*nZStep*abs(InitData.Ai_kappa)<1)*/) {
     /*
     if(InitData.deltaI<0.0000001)
       return(InitData.dt*nZStep*abs(InitData.Ai_kappa));
     else {
       tmp=exp(TCplxLong(0,InitData.deltaI*InitData.dt*nZStep));
       tmp=(tmp-TCplxLong(1,0))*InitData.Ai_kappa/InitData.deltaI;
       qq1=abs(tmp); qq2=InitData.dt*nZStep*abs(InitData.Ai_kappa);
       if(nZZ!=nZStep) {
         nZZ=nZStep;
       }
       return( abs(tmp)<1 ? abs(tmp) : 1 );
     }
     */
     tmp=hNormalP[nZStep]*abs(InitData.Ai_kappa);
     return( abs(tmp)<1 ? abs(tmp) : 1 );
  } else return(1);

//  return( (InitData.UseSimplePhaseModel)&&(InitData.Ai_kappa>0)&&(InitData.dt*nZStep*abs(InitData.Ai_kappa)<1) ?
//           InitData.dt*nZStep*abs(InitData.Ai_kappa) : 1);
}

void cEquVerticalASpeedDistr::InitEls3(int rMul) {
  long double r=rMul*InitData.r;
  long double rnd1=((long double)rand()*M_PI*InitData.rndCoef)/RAND_MAX;
  long double rnd2=RandG(0,InitData.rndCoef);
  for(int i=0;i<InitData.nEls;i++) {
      long double Theta0=2.0*M_PI*i/InitData.nEls;
      for(int j=0;j<InitData.nSpeedFracs;j++) {
        rnd2=RandG(0,InitData.rndCoef);
        hEldata[2*(j*InitData.nEls+i)]=Theta0+r*cos(Theta0+rnd2);
//      hEldata[2*i+1]= (t>=delta) ? -delta-alpha*(t-Z0) : 0;
        hEldata[2*(j*InitData.nEls+i)+1]=hElSpeed[j];
        // (t>=delta) ? -alpha*(t-delta) : 0;
//      hEldata[2*i+1]= -delta*t*t*t*t-alpha*t*t*t;
        hElPresent[j*InitData.nEls+i]=1;
      }
  }
  nMovingEls=InitData.nEls;
}

void cEquVerticalASpeedDistr::SavePData(long double *hPdata) {
  ofstream FileOut = ofstream("P.dat",ios::out);
  for(int i=0;i<InitData.nSteps;i++) {
    FileOut << InitData.dt*i << " " << hPdata[i] << endl;
  }
}

void cEquVerticalASpeedDistr::SaveDataArray(long double x1, long double x2, long double *hData, int n, AnsiString filename){
  ofstream FileOut = ofstream(filename.c_str(),ios::out);
  for(int i=0;i<n;i++) {
    FileOut << x1+(long double)i*(x2-x1)/(n-1) << " " << hData[i] << endl;
  }
}

void cEquVerticalASpeedDistr::LoadEAData(){
  ifstream File("EAData.dat",ios::in);
  long double tmpData[5000];
  long double tmp;
  int fin=0,dataSize,Size=InitData.ElTime/InitData.dt;
  /*
  for(int i=0;i<InitData.ElTime/InitData.dt;i++) {
    File >> tmp;
    File >> tmpData[i];
    if(tmpData[i]*tmpData[i]>0.000000001) fin=1;
    if((tmpData[i]*tmpData[i]<0.000000001)&&(fin)) break;
  } */
  for(dataSize=0;!File.eof();dataSize++) {
    File >> tmp;
    tmpData[dataSize]=0;
    File >> tmpData[dataSize];
    /*
    if(abs(tmpData[dataSize])>0.0000001)
      fin=1;
    if((abs(tmpData[dataSize])<0.0000001)&&(fin)) break;
    */
  }
  double c1,c2;
  int    i1,i2;
  for(int i=0;i<Size;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    hEAData[i]=c1*tmpData[i1]+c2*tmpData[i2];
  }
}



long double cEquVerticalASpeedDistr::CalcEnergy(){
  long double xbase=hEAData[nStepNumber-1];
  long double xm1=xbase-InitData.EADelta;
  long double xm2=xbase-2*InitData.EADelta;
  long double xp1=xbase+InitData.EADelta;
  long double xp2=xbase+2*InitData.EADelta;
  long double ym1,ybase,ym2,yp1,yp2,k,a,b,c;
  long double delta;   /*
  for(int i=0;i<InitData.EAnSteps;i++) {
    xm1=xbase-InitData.EADelta;
    xm2=xbase-2*InitData.EADelta;
    xp1=xbase+InitData.EADelta;
    xp2=xbase+2*InitData.EADelta;
    ym1=EACalcMaxIPos(xm1);
    ym2=EACalcMaxIPos(xm2);
    yp1=EACalcMaxIPos(xp1);
    yp2=EACalcMaxIPos(xp2);
    ybase=(yp1+yp2+ym1+ym2)*0.25;
    k=(yp1+yp2-ym1-ym2)/(6.0*InitData.EADelta);
    if (k*k>0.0000000001) {
      delta=MaxAPos-ybase;
      delta=delta/k;
      xbase+=delta;
    }
  }                      */
  return(xbase);
}

void cEquVerticalASpeedDistr::MakeSpeedDistibution() {
  long double tmp=0;
  for(int i=0;i<InitData.nSpeedFracs;i++) {
    hElSpeedFracs[i]=sin(M_PI*(i+1)/(InitData.nSpeedFracs+1));
    hElSpeedFracs[i]*=hElSpeedFracs[i];
    tmp+=hElSpeedFracs[i];
  }
  for(int i=0;i<InitData.nSpeedFracs;i++) {
    hElSpeedFracs[i]/=tmp;
    hElSpeed[i]=((long double)(i+1))/(InitData.nSpeedFracs+1)*InitData.SpeedDistribFactor*2.0
                -InitData.SpeedDistribFactor;
  }
}
long double cEquVerticalASpeedDistr::Right0()
{ return(hElNew[2*workElNum+1]); }

long double cEquVerticalASpeedDistr::Right1()
{ return( real(-workA*exp(complex<long double>(0,hElNew[2*workElNum]))) ); }

TCplxLong cEquVerticalASpeedDistr::Calc_I(long double * hEls)
{
  TCplxLong Res(0,0),Res1(0,0);
  for(int i=0;i<InitData.nSpeedFracs;i++) {
    Res1=0;
    for(int j=0;j<InitData.nEls;j++)
      if (hElPresent[i*InitData.nEls+j])
        Res1+=exp(TCplxLong(0,-hEls[2*(i*InitData.nEls+j)]));
    Res+=Res1*hElSpeedFracs[i];
  }
  if(nMovingEls)
    return( (long double)2.0*Res/(long double)nMovingEls );
  else return( TCplxLong(0,0) );
}

void cEquVerticalASpeedDistr::StepLOV1()
{
  complex<long double> TempI,TempINew; // for Runge-Kutt's f(new Els)
  complex<long double> TempA; // == (end of space) ? 0 : ...
  long double RungeCoef=1.0/6.0;
  long double aMax=0,tmpAbsA;
  if (InitData.OutMaxA) {
    for(int i=0;i<InitData.nSteps;i++) {
      tmpAbsA=abs(hAdata[i]);
      if(tmpAbsA>aMax) aMax=tmpAbsA;
    }
    hOutA[nStepNumber]=aMax;
  } else hOutA[nStepNumber]=hAdata[InitData.nSteps];
  FillSplineArray();    // Fill the hInterpolated array
//  for(int i=0;i<InitData.nSteps;i++)
//    hInterpolated[i]=(hAdata[i]+hIdata[i]*dt*((long double)0.5)+
//                      hAdata[i+1]+hIdata[i+1]*dt*((long double)0.5))*((long double)0.5);
  InitEls3(1);
  int j=0;
  int ntmpMovingEls;
//  if(hDumpStorage) hDumpStorage->DumpDiffur(hAdata,hIdata,hEldata);
  hIdata1[0]=Kappa(0)*Calc_I(hEldata);
  for(j=0;j<InitData.nSteps;j++)  // The main loop j == z coord
  {
    hAdata1[j]=hAdata[j]+(long double)0.5*dt*(hIdata1[j]+hIdata[j]);
   /*       // test - linear interpolation

    if (j<nSteps-1)
     InterpolatedL=(long double)0.5*(hAdata[j+1]+(long double)0.5*dt*hIdata[j+1])
                 +(long double)0.5*(hAdata[j+2]+(long double)0.5*dt*hIdata[j+2]);
     else InterpolatedL=0;
   */

    memcpy(hElNew,hEldata,InitData.nEls*InitData.nSpeedFracs*2*sizeof(long double));
    TempINew=0;
    ntmpMovingEls=nMovingEls;
    workA=Kappa(dt*((long double)j))*hAdata1[j];
    for(workElNum=0;workElNum<InitData.nEls*InitData.nSpeedFracs;workElNum++)  // k1 calculating
      if (hElPresent[workElNum]) {
        hElk1[2*workElNum]=dt*Right0();
        hElk1[2*workElNum+1]=dt*Right1();
        hElNew[2*workElNum]+=0.5*hElk1[2*workElNum];
        hElNew[2*workElNum+1]+=0.5*hElk1[2*workElNum+1];
//        TempINew+=exp(complex<long double>(0,-hElNew[2*workElNum]));
//        CheckIfElStopped(hElNew+2*workElNum,workElNum);
      }
                                 // k2 calculating
//    if(ntmpMovingEls)
//      TempI=(long double)2.0*kappa(dt*((long double)j+0.5))*TempINew/(long double)ntmpMovingEls;
//    else TempI=0;
    TempI=Kappa(dt*((long double)j))*Calc_I(hElNew);
    TempINew=0;
    ntmpMovingEls=nMovingEls;
    workA=Kappa(dt*((long double)j+0.5))*
          (hInterpolated[j]+(long double)0.5*dt*(TempI));
    for(workElNum=0;workElNum<InitData.nEls*InitData.nSpeedFracs;workElNum++)
      if (hElPresent[workElNum]) {
        hElk2[2*workElNum]=dt*Right0();
        hElk2[2*workElNum+1]=dt*Right1();
        hElNew[2*workElNum]=hEldata[2*workElNum]+0.5*hElk2[2*workElNum];
        hElNew[2*workElNum+1]=hEldata[2*workElNum+1]+0.5*hElk2[2*workElNum+1];
//        TempINew+=exp(complex<long double>(0,-hElNew[2*workElNum]));
//          CheckIfElStopped(hElNew+2*workElNum,workElNum);
      }
                             // k3 calculating
//    if(ntmpMovingEls)
//      TempI=(long double)2.0*kappa(dt*((long double)j+0.5))*TempINew/(long double)ntmpMovingEls;
//    else TempI=0;
    TempI=Kappa(dt*((long double)j+0.5))*Calc_I(hElNew);
    TempINew=0;
    ntmpMovingEls=nMovingEls;
    workA=Kappa(dt*((long double)j+0.5))*
          (hInterpolated[j]+(long double)0.5*dt*(TempI));
    for(workElNum=0;workElNum<InitData.nEls*InitData.nSpeedFracs;workElNum++)
      if (hElPresent[workElNum]) {
        hElk3[2*workElNum]=dt*Right0();
        hElk3[2*workElNum+1]=dt*Right1();
        hElNew[2*workElNum]=hEldata[2*workElNum]+hElk3[2*workElNum];
        hElNew[2*workElNum+1]=hEldata[2*workElNum+1]+hElk3[2*workElNum+1];
        TempINew+=exp(complex<long double>(0,-hElNew[2*workElNum]));
//          CheckIfElStopped(hElNew+2*workElNum,workElNum);
      }

                                 // k4 calculating
//    if(ntmpMovingEls)
//      TempI=(long double)2.0*kappa(dt*((long double)j+1.0))*TempINew/(long double)ntmpMovingEls;
//    else TempI=0;
    TempI=Kappa(dt*((long double)j+0.5))*Calc_I(hElNew);
    TempINew=0;
    workA=Kappa(dt*((long double)j+1.0))*(hAdata[j+1]+(long double)0.5*dt*(TempI+hIdata[j+1]));
    ntmpMovingEls=nMovingEls;
    for(workElNum=0;workElNum<InitData.nEls*InitData.nSpeedFracs;workElNum++)
      if (hElPresent[workElNum]) {
        hElk4[2*workElNum]=dt*Right0();
        hElk4[2*workElNum+1]=dt*Right1();
        hElNew[2*workElNum]=hEldata[2*workElNum]
          +RungeCoef*(hElk1[2*workElNum]+2*hElk2[2*workElNum]+2*hElk3[2*workElNum]+hElk4[2*workElNum]);
        hElNew[2*workElNum+1]=hEldata[2*workElNum+1]
          +RungeCoef*(hElk1[2*workElNum+1]+2*hElk2[2*workElNum+1]+2*hElk3[2*workElNum+1]+hElk4[2*workElNum+1]);
//        TempINew+=exp(complex<long double>(0,-hElNew[2*workElNum]));
//        CheckIfElStopped(hElNew+2*workElNum,workElNum);
      }
    /*
    if(ntmpMovingEls)
      hIdata1[j+1]=(long double)2.0*kappa(dt*((long double)j+1.0))*TempINew/(long double)ntmpMovingEls;
    else hIdata1[j+1]=0;*/
    hIdata1[j+1]=Kappa(dt*((long double)j+1.0))*Calc_I(hElNew);
    memcpy(hEldata,hElNew,InitData.nEls*InitData.nSpeedFracs*2*sizeof(long double));
//    if(hDumpStorage) hDumpStorage->DumpDiffur((hAdata+j),(hIdata+j),hEldata);
  }
  j=InitData.nSteps;
  hAdata1[j]=hAdata[j]+(long double)0.5*dt*(hIdata1[j]+hIdata[j]);
//  hIdata1[nSteps]=-(long double)2.0*G*TempINew/(long double)nEls;//Calc_I(hElNew);
//  hIdata1[nSteps]=Calc_I(hElNew);
  memcpy(hIdata,hIdata1,(InitData.nSteps+1)*sizeof(complex<long double>));
  if(InitData.MoveAData) {
    memcpy(hAdata+1,hAdata1,(InitData.nSteps)*sizeof(complex<long double>));
    hAdata[0]=0;
  } else {
    memcpy(hAdata,hAdata1,(InitData.nSteps+1)*sizeof(complex<long double>));
  }

//  long double rAvgGamma=0;
//  for(int i=0;i<nEls;i++){ rAvgGamma+=hEldata[2*i+1]; }
//  rDebug[0]=rAvgGamma=rAvgGamma/ntmpMovingEls;

//  double ElDroppedFrac=(double)(nEls-nMovingEls)/(double)nEls;
//  ElDroppedSum+=nEls-nMovingEls;
//  if(ElDroppedFrac>ElDroppedMaxFrac) ElDroppedMaxFrac=ElDroppedFrac;

//  DumpRequest.FillReq(this);
//  t+=dt; nCStep++;
}

void cEquVerticalASpeedDistr::FillSplineCoefArray(long double *hData,
          long double *hSplineCoefB,long double *hSplineCoefC,long double *hSplineCoefD) {
  long double dt1=1.0/InitData.dt;
  int nSteps=InitData.nSteps;
  hSplineCoefD[0]=dt;
  hSplineCoefC[1]=(hData[1]-hData[0])*dt1;
  for(int i=1;i<nSteps;i++)
  {
    hSplineCoefD[i]=dt;
    hSplineCoefB[i]=2.0*(hSplineCoefD[i-1]+hSplineCoefD[i]);
    hSplineCoefC[i+1]=(hData[i+1]-hData[i])*dt1;
    hSplineCoefC[i]=hSplineCoefC[i+1]-hSplineCoefC[i];
  }
  hSplineCoefB[0]=hSplineCoefB[nSteps]=-dt;
  hSplineCoefC[0]=(hSplineCoefC[2]-hSplineCoefC[1])*0.5*0.333333333333;
  hSplineCoefC[nSteps]=(hSplineCoefC[nSteps-1]-hSplineCoefC[nSteps-2])*0.5*0.33333333333;
  for(int i=1;i<nSteps+1;i++)
  {
    long double tmp=dt/hSplineCoefB[i-1];
    hSplineCoefB[i]=hSplineCoefB[i]-tmp*dt;
    hSplineCoefC[i]=hSplineCoefC[i]-tmp*hSplineCoefC[i-1];
  }
  hSplineCoefC[nSteps]=hSplineCoefC[nSteps]/hSplineCoefB[nSteps];
  for(int i=nSteps-1;i>=0;i--)
    hSplineCoefC[i]=(hSplineCoefC[i]-dt*hSplineCoefC[i+1])/hSplineCoefB[i];
  hSplineCoefB[nSteps]=(hData[nSteps]-hData[nSteps-1])*dt1
                 +(hSplineCoefC[nSteps-1]+2.0*hSplineCoefC[nSteps])*dt;
  for(int i=0;i<nSteps;i++)
  {
    hSplineCoefB[i]=(hData[i+1]-hData[i])*dt1
                   -(hSplineCoefC[i+1]+2.0*hSplineCoefC[i])*dt;
    hSplineCoefD[i]=(hSplineCoefC[i+1]-hSplineCoefC[i])*dt1;
    hSplineCoefC[i]=3.0*hSplineCoefC[i];
  }
  hSplineCoefC[nSteps]=3.0*hSplineCoefC[nSteps];
  hSplineCoefD[nSteps]=hSplineCoefD[nSteps-1];
}

void cEquVerticalASpeedDistr::FillSplineArray() {
  long double dt1=1.0/InitData.dt;
  long double *hData=(long double *)hInterpolated;
  int tmpNSteps=InitData.nSteps;
  int nSteps=InitData.nSteps;
  for(int i=0;i<nSteps+1;i++)
    htmpSplineData[i]=
       hAdata[i].real()+hIdata[i].real()*dt*((long double)0.5);
  FillSplineCoefArray(htmpSplineData,hSplineCoefBRe,hSplineCoefCRe,hSplineCoefDRe);
  for(int i=0;i<nSteps+1;i++)
    htmpSplineData[i]=
       hAdata[i].imag()+hIdata[i].imag()*dt*((long double)0.5);
  FillSplineCoefArray(htmpSplineData,hSplineCoefBIm,hSplineCoefCIm,hSplineCoefDIm);

  dt1=0.5*dt;
  for(int i=0;i<nSteps+1;i++)
    hInterpolated[i]=TCplxLong(
       hAdata[i].real()+hIdata[i].real()*dt*((long double)0.5)
          +dt1*(hSplineCoefBRe[i]+dt1*(hSplineCoefCRe[i]+dt1*hSplineCoefDRe[i])),
       hAdata[i].imag()+hIdata[i].imag()*dt*((long double)0.5)
          +dt1*(hSplineCoefBIm[i]+dt1*(hSplineCoefCIm[i]+dt1*hSplineCoefDIm[i])));
    nSteps=tmpNSteps;
}

void cEquVerticalASpeedDistr::StepRoutine() {
  StepLOV1();
  //StepRoutineRunge2();
}


