//---------------------------------------------------------------------------


#pragma hdrstop

#include <fstream.h>
#include <vcl\math.hpp>
#include "Equ25DEl03.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)
cEqu2PhaseElectron cEqu2PhaseElectron::operator+(cEqu2PhaseElectron El) {
  cEqu2PhaseElectron tmp;
  tmp.tI=tI+El.tI; tmp.tS=tS+El.tS;
  tmp.P=P+El.P; tmp.W=W+El.W; tmp.P1=P1+El.P1; tmp.P2=P2+El.P2;
  tmp.Present=Present|El.Present;
  return(tmp);
}
cEqu2PhaseElectron cEqu2PhaseElectron::operator*(long double mul) {
  cEqu2PhaseElectron tmp;
  tmp.tI=tI*mul; tmp.tS=tS*mul; tmp.P=P*mul; tmp.W=W*mul; tmp.P1=P1*mul; tmp.P2=P2*mul;
  tmp.Present=Present;
  return(tmp);
}

void cEqu2PhaseElModel::InitReport(cDifEquReport *hReport) {
  int Lines[EQU_MAX_GRAPHS];
  for(int i=0;i<EQU_MAX_GRAPHS;i++) Lines[i]=1;
  Lines[0]=1; Lines[1]=1; Lines[2]=1; Lines[3]=1;
  hReport->LoadText(hReportText);
  hReport->SetGraphSize(EQU_GRAPH_SIZE,Lines);
  hReport->SetMaxTo1();
}
void cEqu2PhaseElModel::Report(cDifEquReport *hReport) {
/*  long double *hGamma = new long double[InitData.nSteps+1];
  for(int i=0;i<=InitData.nSteps;i++) {
    hGamma[i]=0;
    for(int j=0;j<InitData.nEls;j++) {
      if(hElPresent[i*InitData.nEls+j])
        hGamma[i]-=hEl[i*InitData.nEls*2+j*2+1];
    }
    hGamma[i]/=InitData.nEls;
  } */
  hReport->Graph[0].GetDataAbs(0,InitData.nSteps+1,hA);
  hReport->Graph[1].GetDataAbs(0,InitData.nSteps+1,hI);
  hReport->Graph[2].GetData(0,InitData.nSteps+1,hAvgP);
  hReport->Graph[2].Max = hReport->Graph[2].Max > -hReport->Graph[2].Min ?
                          hReport->Graph[2].Max : -hReport->Graph[2].Min;
  hReport->Graph[3].GetData(0,InitData.ElTime/InitData.dt,hEAData);
  hReport->Graph[4].GetData(0,InitData.ElTime/InitData.dt,hRandomData);
  hReport->nVals=10;
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
  hReport->ValText[1]="Max I "; hReport->Val[1]=hReport->Graph[1].Max;
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
//  hReport->Graph[2].Max=2;
}
void cEqu2PhaseElModel::FinalReport(cDifEquReport *hReport) {
  hReport->nVals=4; hReport->nGraphs=1;
  hReport->Graph[0].GetDataAbs(0,nStepNumber,hOutA);
  hReport->GraphText[0]="Out A";
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
  hReport->ValText[1]="integral A "; hReport->Val[1]=intA;
  hReport->ValText[2]="integral W "; hReport->Val[2]=intW;
  hReport->ValText[3]="integral A-W "; hReport->Val[3]=intA-intW;
  
  ofstream FileOut = ofstream("Out.dat",ios::out);
  for(int i=0;i<nStepNumber;i++) {
    FileOut << InitData.dt*i << " " << abs(hOutA[i]) << endl;
  }
  if(InitData.UseEnergyAdapt)
    SaveDataArray(0,InitData.ElTime,hEAData,InitData.ElTime/InitData.dt,"EAData.dat");
}

cEqu2PhaseElModel::cEqu2PhaseElModel() : cEqu25DElModel()
{
  hI=hA=hAk1=hAk2=hABack=hIBack=NULL; hOutA=NULL;
  hEl=hElk1=hElk2=NULL; hAvgP=NULL; hEAData=NULL; hRandomData=NULL;
  hEl2=hElk12=hElk22=NULL;  hdeltaTS=NULL; hNormalP=NULL;
  hNMovingEls=hElPresent=NULL; InitData.ElDirection=1;
  InitData.nEls=InitData.nSteps=0; EATempIData=NULL;
  InitData.tMax=0;
  nStepNumber=0;
}

cEqu2PhaseElModel::~cEqu2PhaseElModel()// : ~cDifEqu_DataBankIO()
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
}

void cEqu2PhaseElModel::MakeNormalP() {
  if(!hNormalP) hNormalP = new long double[InitData.nSteps+1];
  for(int i=0;i<=InitData.nSteps;i++)
    if(InitData.deltaI*InitData.deltaI<0.00000000001) hNormalP[i]=InitData.dt*i;
    else {
      hNormalP[i]=2*sin(0.5*InitData.deltaI*InitData.dt*i)/InitData.deltaI;
    }
}

void cEqu2PhaseElModel::LoadInitData(cDifEqu_InitDataRec *hInitData)
{
  cDifEqu_InitDataRec_25DElModel *hInit=(cDifEqu_InitDataRec_25DElModel *)hInitData;
  int nTStepsNew=hInit->tMax/hInit->dt;
  int nTSteps=InitData.tMax/InitData.dt;
  if ( nTSteps!=nTStepsNew ) {
    delete [] hOutA; hOutA=new TCplxLong[nTStepsNew+1];
    delete [] hEAData; hEAData=new long double[nTStepsNew+1];
    delete [] hRandomData; hRandomData=new long double[nTStepsNew+1];
  }
  if( InitData.nSteps!=hInit->nSteps)  {
    delete [] hI;  delete [] hA; delete [] hAvgP;
    delete [] hAk1;  delete [] hAk2; delete [] EATempIData;
    delete [] hABack; delete [] hIBack;
    delete [] hNormalP; hNormalP=NULL;
    EATempIData=new long double[hInit->nSteps+1];
    hI=new TCplxLong[hInit->nSteps+1];
    hA=new TCplxLong[hInit->nSteps+1];
    hAvgP=new long double[hInit->nSteps+1];
    hABack=new TCplxLong[hInit->nSteps+1];
    hIBack=new TCplxLong[hInit->nSteps+1];
    hAk1=new TCplxLong[hInit->nSteps+1];
    hAk2=new TCplxLong[hInit->nSteps+1];
  }
  if( (hInit->nEls!=InitData.nEls) || (hInit->nSteps!=InitData.nSteps) )
  {
    delete [] hEl2;// delete [] hElNew;
    delete [] hElk12;  delete [] hElk22; delete [] hdeltaTS;
    if (hInit->UseSimplePhaseModel) {
      hEl2=new cEqu2PhaseElectron[hInit->nEls];
      hElk12=new cEqu2PhaseElectron[hInit->nEls];
      hElk22=new cEqu2PhaseElectron[hInit->nEls];
      hdeltaTS = new long double[hInit->nEls];
    } else {
      hEl2=new cEqu2PhaseElectron[hInit->nEls*hInit->nEls];
      hElk12=new cEqu2PhaseElectron[hInit->nEls*hInit->nEls];
      hElk22=new cEqu2PhaseElectron[hInit->nEls*hInit->nEls];
      hdeltaTS = new long double[hInit->nEls*hInit->nEls];
    }
  }
  InitData=*hInit; nStepNumber=0; maxW=maxdeltaTS=intA=intW=0;
  t=0; tMax=InitData.tMax; dt=InitData.dt; MaxAPos=0;

  MakeNormalP();
  for(int i=0;i<InitData.nSteps+1;i++) { hA[i]=0; hI[i]=0; }
  for(int i=0;i<nTStepsNew;i++) { hEAData[i]=0; }
  for(int i=0;i<nTStepsNew;i++) { hRandomData[i]=0; }
  if((InitData.EAnSteps<0)&&(InitData.UseEnergyAdapt)) {
    LoadEAData();
    InitData.UseEnergyAdapt=false;
  }
  lastRandom=0;
  lastRandom1=RandG(0,0.5)*InitData.rndCoef;
  lastRandom2=0;
  randomize();
}

long double cEqu2PhaseElModel::SimplePhaseDeSync(int nZStep) {
  long double tmp=hNormalP[nZStep]*InitData.AI.real();
  return(tmp*tmp);
}

long double cEqu2PhaseElModel::EADesync(int nZStep) {
  long double tmp=(InitData.dt*nZStep-InitData.EADeltaZ0)*InitData.EAkz;
  return(tmp > 0 ? tmp:0);
}

long double cEqu2PhaseElModel::SimplePhaseKappa(int nZStep) {
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

void cEqu2PhaseElModel::InitEls2(int rMul)
{
  long double r=rMul*InitData.r;
  TCplxLong TempI=0;
  int nElColumns = InitData.UseSimplePhaseModel ? 1 : InitData.nEls;
  long double rnd1=((long double)rand()*M_PI*InitData.rndCoef)/RAND_MAX;
  long double rnd2=((long double)rand()*M_PI*InitData.rndCoef)/RAND_MAX;
  rnd1=RandomRange(0,1000000000);
  rnd2=RandomRange(0,1000000000);
  rnd1*=M_PI*InitData.rndCoef/1000000000;
  rnd2*=M_PI*InitData.rndCoef/1000000000;
//  rnd2=rnd1=RandG(lastRandom,InitData.rndCoef);
//  rnd2=rnd1=sin(InitData.rndCoef*t);
  //rnd2=rnd1=lastRandom;
  if(InitData.NoRandom) {                     
    rnd1=rnd2=(sin(1379*t)*sin(1379*t)+sin(573*t)*sin(573*t));
  }
  if(t<InitData.ElTime) {
    for(int i=0;i<nElColumns;i++)
    {
      long double ThetaI0=2.0*M_PI*i/InitData.nEls;
      for(int j=0;j<InitData.nEls;j++) {
        long double ThetaS0=2.0*M_PI*j/InitData.nEls;

        hEl2[InitData.nEls*i+j].tI=ThetaI0+r*cos(ThetaI0/*+rnd1*/);//*M_PI*InitData.rndCoef);
        hEl2[InitData.nEls*i+j].tS=ThetaS0+r*cos(ThetaS0+rnd2);//*M_PI*InitData.rndCoef);
        hEl2[InitData.nEls*i+j].P=InitData.UseSimplePhaseModel ? 1 : 0;
        hEl2[InitData.nEls*i+j].P1=0;
        hEl2[InitData.nEls*i+j].P2=0;
        hEl2[InitData.nEls*i+j].W=0;
        hEl2[InitData.nEls*i+j].Present=1;
        TempI+=(hEl2[InitData.nEls*i+j].P)*exp(TCplxLong(0,-hEl2[InitData.nEls*i+j].tS));
      }
    }
    NMovingEls=InitData.nEls;
  } else {
    for(int i=0;i<nElColumns*InitData.nEls;i++)
      hEl2[i].Present=0;
    NMovingEls=0;
  }
  TempI*=InitData.ICoef*4*M_PI/(nElColumns*InitData.nEls);
  hI[0]=TempI;
//  lastRandom=rnd1;
  hRandomData[nStepNumber]=lastRandom;
  lastRandom+=lastRandom1*InitData.dt;
//  lastRandom1+=lastRandom2*InitData.dt;
  lastRandom1+=RandG(0,InitData.rndCoef)*InitData.dt-InitData.rndCoef*lastRandom*0.1*InitData.dt;
}

cEqu2PhaseElectron cEqu2PhaseElModel::RightPartEl2(long double El_TI,long double El_TS,
                             TCplxLong El_P,long double El_W,TCplxLong A,int nZStep)
{
  cEqu2PhaseElectron tmp;
  long double qq1,qq2,qq3;
  TCplxLong qqC;
  if(InitData.UseSimplePhaseModel) {
    tmp.tS=El_W+SimplePhaseDeSync(nZStep)+hEAData[nStepNumber]-EADesync(nZStep);
    tmp.tI=0;
    tmp.P=0;
    tmp.P1=InitData.AI*exp(TCplxLong(0,El_TI));
    tmp.W=SimplePhaseKappa(nZStep)*( -A*exp(TCplxLong(0,El_TS)) ).real();
    tmp.Present=0;
  } else {
    tmp.tS=InitData.deltaS+El_W+El_P.real()*El_P.real()+El_P.imag()*El_P.imag();
    tmp.tI=InitData.deltaI-El_W;
    tmp.P=-A*exp(TCplxLong(0,El_TS))-InitData.AI*exp(TCplxLong(0,El_TI));
    tmp.P1=-InitData.AI*exp(TCplxLong(0,El_TI));
    tmp.P2=-InitData.AI*exp(TCplxLong(0,InitData.deltaI*InitData.dt*nZStep));
    tmp.W=( conj(El_P)*A*exp(TCplxLong(0,El_TS))*((long double)2.0) ).real();
    /*
    qqC=El_P*exp(TCplxLong(0,-El_TS))*InitData.dt;
    qq2=abs(A+qqC)*abs(A+qqC)-abs(A)*abs(A);
    qq3=InitData.dt*tmp.W; */
    tmp.Present=0;
  }
  return(tmp);
}

void cEqu2PhaseElModel::SavePData(long double *hPdata) {
  ofstream FileOut = ofstream("P.dat",ios::out);
  for(int i=0;i<InitData.nSteps;i++) {
    FileOut << InitData.dt*i << " " << hPdata[i] << endl;
  }
}

void cEqu2PhaseElModel::SaveDataArray(long double x1, long double x2, long double *hData, int n, AnsiString filename){
  ofstream FileOut = ofstream(filename.c_str(),ios::out);
  for(int i=0;i<n;i++) {
    FileOut << x1+(long double)i*(x2-x1)/(n-1) << " " << hData[i] << endl;
  }
}

void cEqu2PhaseElModel::LoadEAData(){
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

void cEqu2PhaseElModel::StepRoutine() {
  StepRoutineRunge2();
}

void cEqu2PhaseElModel::StepRoutineRunge2() {
//  if(hDump) hDump->Dump(hA,hI,hEl,hElPresent);
  TCplxLong TempI;
  long double qq1,qq2,qq3,avgW;
  if(hDump) hDump->Dump(hA,hI,NULL,NULL);
  int nElColumns = InitData.UseSimplePhaseModel ? 1 : InitData.nEls;

  if(InitData.UseSimplePhaseModel&&InitData.UseEnergyAdapt&&(t>InitData.EAStartTime)){
                        // ---------- Energy Adaptation System

    if((t<InitData.EAEndTime)) {
      hEAData[nStepNumber]=CalcEnergy();
      EACoef=(hEAData[nStepNumber-1]-hEAData[nStepNumber-6])/5;
    }
    else {
      hEAData[nStepNumber]=hEAData[nStepNumber-1]+EACoef;
    }
    /*
    if((t>=InitData.EAEndTime)&&(t<InitData.EAEndTime+InitData.dt))
      SaveEnergyAdaptCalibrationData();*/
  }

  InitEls2(1);
  int SaveP = (nStepNumber==0) ? 1 : 0;
//  if( (t>=2.0)&&(t<2.0+InitData.dt) ) SaveP=1;
//  maxW=0;
  for(int i=0;i<nElColumns*InitData.nEls;i++) hdeltaTS[i]=0;
  memcpy(hABack,hA,(InitData.nSteps+1)*sizeof(TCplxLong));
  memcpy(hIBack,hI,(InitData.nSteps+1)*sizeof(TCplxLong));
  hOutA[nStepNumber]=hA[InitData.nSteps];//*hA[InitData.nSteps];
  long double *hPdata;
  if(SaveP) hPdata = new long double[InitData.nSteps+1];
  hA[0]=GetAInput();
  hAvgP[0]=0;
  if(SaveP) hPdata[0]=0;

  MaxAPos=0;
  for(int i=0;i<InitData.nSteps;i++) {
//    if(SaveP) hPdata[i]=abs(hEl2[0].P)*abs(hEl2[0].P);
    hAvgP[i+1]=0;
    TempI=0;
    if(NMovingEls) {
      for(int j=0;j<nElColumns*InitData.nEls;j++) {
        if(hEl2[j].Present) {
          hElk12[j]=RightPartEl2(hEl2[j].tI,hEl2[j].tS,hEl2[j].P,hEl2[j].W,hA[i],i)*InitData.dt;
          TempI+=(hEl2[j].P+hElk12[j].P)*exp(TCplxLong(0,-hEl2[j].tS-hElk12[j].tS));
          qq1=abs(hEl2[j].P);
          qq2=abs(hEl2[j].P1);//abs(InitData.AI)*sin(InitData.deltaI*i*InitData.dt)/InitData.deltaI;
          hdeltaTS[j]+=InitData.dt*(qq1*qq1-qq2*qq2);
          //hdeltaTS[j]+=InitData.dt*hEl2[j].W;
        }
      }
    }
    TempI*=InitData.ICoef*4*M_PI/(nElColumns*InitData.nEls);
    hAk1[i]=(TempI+hIBack[i])*SimplePhaseKappa(i)*InitData.dt*((long double)0.5);
    TempI=0;
    if(NMovingEls) {
      for(int j=0;j<nElColumns*InitData.nEls;j++) {
        if(hEl2[j].Present) {
          hElk22[j]=RightPartEl2(hEl2[j].tI+hElk12[j].tI,hEl2[j].tS+hElk12[j].tS,
                                 hEl2[j].P+hElk12[j].P,hEl2[j].W+hElk12[j].W,
                                 hABack[i]+hAk1[i],i)*InitData.dt;
          hEl2[j]=hEl2[j]+(hElk12[j]+hElk22[j])*((long double)0.5);
          TempI+=(hEl2[j].P)*exp(TCplxLong(0,-hEl2[j].tS));
          hAvgP[i+1]+=abs(hEl2[j].P2)*abs(hEl2[j].P2)-hNormalP[i+1]*hNormalP[i+1];
//          hAvgP[i+1]+=hNormalP[i+1]*hNormalP[i+1];
//          hAvgP[i+1]+=abs(hEl2[j].P2)*abs(hEl2[j].P2);
          //abs(hEl2[j].P)*abs(hEl2[j].P);//-abs(hEl2[j].P2)*abs(hEl2[j].P2);
          maxW = abs(hEl2[j].W)>abs(maxW) ? hEl2[j].W : maxW;
        }
      }
    }
    if(NMovingEls) {
      for(int j=0;j<nElColumns*InitData.nEls;j++) {
        if (abs(hdeltaTS[j])>abs(maxdeltaTS)) maxdeltaTS=hdeltaTS[j];
      }
    }
    TempI*=InitData.ICoef*4*M_PI/(nElColumns*InitData.nEls);
    hAvgP[i+1]/=(nElColumns*InitData.nEls);
    if(SaveP) hPdata[i+1]=hAvgP[i+1];
    hAk2[i]=(TempI+hIBack[i])*InitData.dt*SimplePhaseKappa(i)*((long double)0.5);
    hA[i+1]=hABack[i]+(hAk1[i]+hAk2[i])*((long double)0.5);
    if(abs(hA[i+1])>abs(hA[MaxAPos])) MaxAPos=i+1;
    qq1=0;
    if(NMovingEls) {
      for(int j=0;j<nElColumns*InitData.nEls;j++) {
        qq1+=(hElk12[j].W+hElk22[j].W)*((long double)0.5);
      }
    }
    qq1*=4*M_PI*InitData.dt/(nElColumns*InitData.nEls);
    qq2=abs(hA[i+1])*abs(hA[i+1])-abs(hABack[i])*abs(hABack[i]);
    hI[i+1]=TempI;
    int qq=hEl2[0].tS/2.0/M_PI;
    //hAvgP[i+1]=hEl2[0].tS-2.0*M_PI*qq;
  }
  avgW=0;
  if(NMovingEls) {
    for(int j=0;j<nElColumns*InitData.nEls;j++) {
      avgW+=hEl2[j].W;
    }
  }
  avgW*=4*M_PI*InitData.dt/(nElColumns*InitData.nEls);
  intW+=avgW;
  intA+=abs(hA[InitData.nSteps-1])*abs(hA[InitData.nSteps-1])*InitData.dt;

  //hOutA[nStepNumber]=hA[InitData.nSteps];//*hA[InitData.nSteps];
  if(SaveP) { SavePData(hPdata); delete [] hPdata; }
}

long double cEqu2PhaseElModel::CalcEnergy(){
  long double xbase=hEAData[nStepNumber-1];
  long double xm1=xbase-InitData.EADelta;
  long double xm2=xbase-2*InitData.EADelta;
  long double xp1=xbase+InitData.EADelta;
  long double xp2=xbase+2*InitData.EADelta;
  long double ym1,ybase,ym2,yp1,yp2,k,a,b,c;
  long double delta;
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
  }
  return(xbase);
}

int cEqu2PhaseElModel::EACalcMaxIPos(long double delta){
  long double tmp=hEAData[nStepNumber];
  long double tmp1=InitData.rndCoef;
  hEAData[nStepNumber]=delta;
  InitData.rndCoef=0;
  int nElColumns = InitData.UseSimplePhaseModel ? 1 : InitData.nEls;
  TCplxLong TempI;
  int EAMaxIPos;
  long double EAMaxI;
  EAMaxIPos=0; EAMaxI=0;
  InitEls2(1);
  for(int i=0;i<InitData.nSteps;i++) {
    TempI=0;
    if(NMovingEls) {
      for(int j=0;j<nElColumns*InitData.nEls;j++) {
        if(hEl2[j].Present) {
          hEl2[j]=hEl2[j]+RightPartEl2(hEl2[j].tI,hEl2[j].tS,hEl2[j].P,hEl2[j].W,hA[i],i)*InitData.dt;
          TempI+=(hEl2[j].P)*exp(TCplxLong(0,-hEl2[j].tS));
        }
      }
    }
    EATempIData[i]=abs(TempI);
    if(EATempIData[i]>EAMaxI){
      EAMaxI=EATempIData[i];
      EAMaxIPos=i;
    }
  }
  for(int i=0;i<InitData.nSteps-1;i++)
    if((EATempIData[i]>EATempIData[i+1])&&(EATempIData[i]>EAMaxI*0.8)) {
      EAMaxIPos=i;
      break;
    }
  hEAData[nStepNumber]=tmp;
  InitData.rndCoef=tmp1;
  return(EAMaxIPos);
}

void cEqu2PhaseElModel::SaveEnergyAdaptCalibrationData(){
  int nElColumns = InitData.UseSimplePhaseModel ? 1 : InitData.nEls;
  long double *EACalibration = new long double [801];
  for(int k=0;k<801;k++) {
    EACalibration[k]=EACalcMaxIPos(4.0*double(k-200)/100)*InitData.dt;
  }
  SaveDataArray(-8,8,EACalibration,801,"EA.dat");
//  hEAData[nStepNumber]=0;
  delete [] EACalibration;
}










void cEqu2PhaseElModel::StepRoutineRunge4() {
/*
//  if(hDump) hDump->Dump(hA,hI,hEl,hElPresent);
  long double RungeCoef[5];
  RungeCoef[0]=0; RungeCoef[1]=0; RungeCoef[2]=0; RungeCoef[3]=0; RungeCoef[4]=0;
  TCplxLong TempI;
  if(hDump) hDump->Dump(hA,hI,NULL,NULL);
  memcpy(hABack,hA,(InitData.nSteps+1)*sizeof(TCplxLong));
  memcpy(hIBack,hI,(InitData.nSteps+1)*sizeof(TCplxLong));
  hOutA[nStepNumber]=hA[InitData.nSteps];//*hA[InitData.nSteps];
  InitEls2(1);
  hA[0]=GetAInput();
  for(int i=0;i<InitData.nSteps;i++) {
    TempI=0;
    if(NMovingEls) {
      for(int j=0;j<InitData.nEls*InitData.nEls;j++) {
        if(hEl2[j].Present) {
          hElk2[j]+=RightPartEl2(hEl2[j].tI,hEl2[j].tS,hEl2[j].P,hEl2[j].W,hA[i],i)*InitData.dt;
          TempI+=(hEl2[j].P)*exp(TCplxLong(0,-hEl2[j].tS));
        }
      }
    }
    TempI*=InitData.ICoef*4*M_PI/(InitData.nEls*InitData.nEls);
    hAk1[i]=(TempI+hIBack[i])*InitData.dt*((long double)0.5);
    TempI=0;
    if(NMovingEls) {
      for(int j=0;j<InitData.nEls*InitData.nEls;j++) {
        if(hEl2[j].Present) {
          hElk22[j]=RightPartEl2(hEl2[j].tI+hElk12[j].tI,hEl2[j].tS+hElk12[j].tS,
                                 hEl2[j].P+hElk12[j].P,hEl2[j].W+hElk12[j].W,
                                 hABack[i]+hAk1[i],i)*InitData.dt;
          hEl2[j]=hEl2[j]+(hElk12[j]+hElk22[j])*((long double)0.5);
          TempI+=(hEl2[j].P)*exp(TCplxLong(0,-hEl2[j].tS));
        }
      }
    }
    TempI*=InitData.ICoef*4*M_PI/(InitData.nEls*InitData.nEls);
    hAk2[i]=(TempI+hIBack[i])*InitData.dt*((long double)0.5);
    hA[i+1]=hABack[i]+(hAk1[i]+hAk2[i])*((long double)0.5);
    hI[i+1]=TempI;
  }
  //hOutA[nStepNumber]=hA[InitData.nSteps];//*hA[InitData.nSteps];*/
}

