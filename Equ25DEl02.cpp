//---------------------------------------------------------------------------


#pragma hdrstop

#include "Equ25DEl02.h"
#include <fstream.h>



//---------------------------------------------------------------------------

#pragma package(smart_init)

void cEqu25DElModel::InitReport(cDifEquReport *hReport) {
  int Lines[EQU_MAX_GRAPHS];
  for(int i=0;i<EQU_MAX_GRAPHS;i++) Lines[i]=1;
  Lines[0]=1; Lines[1]=1; Lines[2]=1; Lines[3]=1;
  hReport->LoadText(hReportText);
  hReport->SetGraphSize(EQU_GRAPH_SIZE,Lines);
  hReport->SetMaxTo1();
}
void cEqu25DElModel::Report(cDifEquReport *hReport) {
  long double *hGamma = new long double[InitData.nSteps+1];
  for(int i=0;i<=InitData.nSteps;i++) {
    hGamma[i]=0;
    for(int j=0;j<InitData.nEls;j++) {
      if(hElPresent[i*InitData.nEls+j])
        hGamma[i]-=hEl[i*InitData.nEls*2+j*2+1];
    }
    hGamma[i]/=InitData.nEls;
  }
  hReport->Graph[0].GetDataAbs(0,InitData.nSteps+1,hA);
  hReport->Graph[1].GetDataAbs(0,InitData.nSteps+1,hI);
  hReport->Graph[2].GetData(0,InitData.nSteps+1,hGamma);
  hReport->nVals=3;
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
  hReport->ValText[1]="Max I "; hReport->Val[1]=hReport->Graph[1].Max;
  hReport->ValText[2]="Max Gamma "; hReport->Val[2]=hReport->Graph[2].Max;
  hReport->Graph[2].Max=2;
  delete [] hGamma;
}
void cEqu25DElModel::FinalReport(cDifEquReport *hReport) {
  hReport->nVals=1; hReport->nGraphs=1;
  hReport->Graph[0].GetDataAbs(0,nStepNumber,hOutA);
  hReport->GraphText[0]="Out A";
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
}

cEqu25DElModel::cEqu25DElModel() : cDifEqu_DataBankIO()
{
  hI=hA=NULL; hOutA=NULL;
  hEl=hElk1=hElk2=NULL; EATempIData=NULL;
  hNMovingEls=hElPresent=NULL; InitData.ElDirection=1;
  InitData.nEls=InitData.nSteps=0;
  InitData.tMax=0;
  nStepNumber=0;  Got1stMax=0;  InitData.fDeleteDeltaZData=InitData.fDeleteKappaData=0;
}

cEqu25DElModel::~cEqu25DElModel()// : ~cDifEqu_DataBankIO()
{
  if(InitData.fDeleteKappaData) {
    delete [] InitData.hKappaData;
    InitData.hKappaData=NULL;
    InitData.fDeleteKappaData=0;
  }
  if(InitData.fDeleteDeltaData) {
    delete [] InitData.hDeltaData;
    InitData.hDeltaData=NULL;
    InitData.fDeleteDeltaData=0;
  }
  if(InitData.fDeleteDeltaZData) {
    delete [] InitData.hDeltaZData;
    InitData.hDeltaZData=NULL;
    InitData.fDeleteDeltaZData=0;
  }
  delete [] hI;  delete [] hA; delete hOutA;
  delete [] hEl;  delete [] hElk1; delete [] hElk2;
  delete [] hElPresent; delete [] hNMovingEls; delete [] EATempIData;
  hI=hA=NULL; hEl=NULL;
  hElk1=hElk2=NULL;  EATempIData=NULL;
  hElPresent=hNMovingEls=NULL;
  hDump=NULL;         
}

void cEqu25DElModel::LoadInitData(cDifEqu_InitDataRec *hInitData)
{
  cDifEqu_InitDataRec_25DElModel *hInit=(cDifEqu_InitDataRec_25DElModel *)hInitData;
  int nTStepsNew=hInit->tMax/hInit->dt;
  int nTSteps=InitData.tMax/InitData.dt;
    if ( nTSteps!=nTStepsNew ) {
    delete [] hOutA; hOutA=new TCplxLong[nTStepsNew+1];
  }
  if(hInit->UseEnergyAdapt) hInit->DeltaMode=1;
  if(hInit->fDeleteKappaData) {
    delete [] hInit->hKappaData;
    hInit->hKappaData=NULL;
  }
  if(hInit->fDeleteDeltaData) {
    delete [] hInit->hDeltaData;
    hInit->hDeltaData=NULL;
  }
  if(hInit->DeltaMode==1) {
    hInit->hDeltaData=new long double[nTStepsNew+1];
    for(int i=0;i<nTStepsNew+1;i++) hInit->hDeltaData[i]=0;
    hInit->fDeleteDeltaData=1;
  }
  if( InitData.nSteps!=hInit->nSteps)  {
    delete [] EATempIData;
    delete [] hI;  delete [] hA;
    delete [] hNMovingEls;
    hI=new TCplxLong[hInit->nSteps+1];
    hA=new TCplxLong[hInit->nSteps+1];
    hNMovingEls = new int[hInit->nSteps+1];
    EATempIData = new long double[hInit->nSteps+1];
  }
  if( (hInit->nEls!=InitData.nEls) || (hInit->nSteps!=InitData.nSteps) )
  {
    delete [] hEl;// delete [] hElNew;
    delete [] hElk1;  delete [] hElk2;
    delete [] hElPresent;
    hEl=new long double[(hInit->nSteps+1)*2*hInit->nEls];
    //hElDataNew=new REAL_TYPE[(hInit->nSteps+1)*2*hInit->nEls];
    hElk1=new long double[(hInit->nSteps+1)*2*hInit->nEls];
    hElk2=new long double[(hInit->nSteps+1)*2*hInit->nEls];
    hElPresent = new int[(hInit->nSteps+1)*hInit->nEls];
  }
  InitData=*hInit; nStepNumber=0;
  t=0; tMax=InitData.tMax; dt=InitData.dt;

  double tmp=InitData.ElTime;
  if(!InitData.ElsOnStart) InitData.ElTime=0;
  for(int i=0;i<InitData.nSteps+1;i++) InitEls(i,0);
  InitData.ElTime=tmp; nZStepNumber=0;
  InitKappa();

  CalcI();
  FillAData();
}

long double cEqu25DElModel::Delta() {
  if(InitData.DeltaMode==1) return(InitData.hDeltaData[nStepNumber]);
  return( InitData.deltaV1+t*InitData.deltaV2 );
//  return( (InitData.deltaV1+t*InitData.deltaV2)>0 ? InitData.deltaV1+t*InitData.deltaV2 : 0  );
}

long double cEqu25DElModel::DeltaZ() {
  long double midZ,c1,c2;
  int i1,i2,imid;
  if(InitData.DeltaZMode==2) {
    if(z==zBack) return zDeltaZBack;
    i1=0; i2=InitData.nDeltaZSteps-1;
    while(i2-i1>1) {
      imid=(i1+i2)/2;
      if(InitData.hDeltaZData[2*imid]<z) i1=imid;
        else i2=imid;
    }
    c1=(InitData.hDeltaZData[2*i2]-z)/(InitData.hDeltaZData[2*i2]-InitData.hDeltaZData[2*i1]);
    c2=(z-InitData.hDeltaZData[2*i1])/(InitData.hDeltaZData[2*i2]-InitData.hDeltaZData[2*i1]);
    zBack=z; zDeltaZBack=InitData.hDeltaZData[2*i1+1]*c1+InitData.hDeltaZData[2*i2+1]*c2;
    return zDeltaZBack;
  }
  if(InitData.fUseLOVModel) {
    if(z>InitData.L-InitData.EADeltaZ0) return 0;
    return((InitData.L-InitData.EADeltaZ0-z)*InitData.EAkz);
  }
  if(InitData.EADeltaZ0<0) return InitData.EAkz;
  if(z<InitData.EADeltaZ0) return 0;
  return((z-InitData.EADeltaZ0)*InitData.EAkz);
}

long double cEqu25DElModel::DeltaZEXT(long double z,long double t,cDifEqu_InitDataRec_25DElModel *hInitData) {
  long double midZ,c1,c2;
  int i1,i2,imid;
  if(hInitData->DeltaZMode==2) {
    i1=0; i2=hInitData->nDeltaZSteps-1;
    while(i2-i1>1) {
      imid=(i1+i2)/2;
      if(hInitData->hDeltaZData[2*imid]<z) i1=imid;
        else i2=imid;
    }
    c1=(hInitData->hDeltaZData[2*i2]-z)/(hInitData->hDeltaZData[2*i2]-hInitData->hDeltaZData[2*i1]);
    c2=(z-hInitData->hDeltaZData[2*i1])/(hInitData->hDeltaZData[2*i2]-hInitData->hDeltaZData[2*i1]);
    return hInitData->hDeltaZData[2*i1+1]*c1+hInitData->hDeltaZData[2*i2+1]*c2;
  }
  if(hInitData->fUseLOVModel) {
    if(z>hInitData->L-hInitData->EADeltaZ0) return 0;
    return((hInitData->L-hInitData->EADeltaZ0-z)*hInitData->EAkz);
  }
  if(InitData.EADeltaZ0<0) return InitData.EAkz;
  if(z<hInitData->EADeltaZ0) return 0;
  return((z-hInitData->EADeltaZ0)*hInitData->EAkz);
}

long double cEqu25DElModel::AInput() {
long double tmp;
  if(InitData.AInputMode==0) {
    return(InitData.AInputV1);
  }
  if(InitData.AInputMode==1) {
    tmp=(t-InitData.AInputV2)*(t-InitData.AInputV2)/
                    (InitData.AInputV3*InitData.AInputV3);
    return(InitData.AInputV1*exp(-tmp));
  }
  if(InitData.AInputMode==2) {
    if(t<InitData.AInputV3*3)
      tmp=(t-InitData.AInputV3*3)*(t-InitData.AInputV3*3)/
                    (InitData.AInputV3*InitData.AInputV3);
    if( (t>=InitData.AInputV3*3)&&(t<=InitData.AInputV3*3+InitData.AInputV2) )
      tmp=0;
    if( (t>InitData.AInputV3*3+InitData.AInputV2) )
      tmp=(t-InitData.AInputV3*3-InitData.AInputV2)*
          (t-InitData.AInputV3*3-InitData.AInputV2)/
                    (InitData.AInputV3*InitData.AInputV3);
    return(InitData.AInputV1*exp(-tmp));
  }
  return 0;
}

void cEqu25DElModel::InitKappa() {
  long double Z0=InitData.KappaV1, p=InitData.KappaV2, L=InitData.dt*InitData.nSteps;
  KappaZInv = (Z0>0.000001) ? 1.0/Z0 : 1.0 ;
  if (InitData.KappaMode==KAPPAMODE_TOMSK) {
/*
    KappaC =
     InitData.dt*InitData.nSteps*(1.0+InitData.KappaV2) /
       ( (InitData.dt*InitData.nSteps-InitData.KappaV1)*(1.0+InitData.KappaV2)+2.0*InitData.KappaV1 );*/
    if(InitData.KappaV4<0) KappaC = L*(1.0+p) / ( (L-Z0)*(1.0+p)+2.0*Z0 );
      else KappaC=InitData.KappaV4;
  }
  if (InitData.KappaMode==KAPPAMODE_LINEAR)
    if(InitData.KappaV4<0) KappaC= p>KappaZInv ? L*p/((L-Z0)*(p)+0.5) : L/(L-0.5*p*Z0*Z0);
       else KappaC=InitData.KappaV4;
  if(InitData.fLoadKappaData) {
    LoadKappaData();
  }
}


long double cEqu25DElModel::Kappa(long double z,long double t) {
  if (InitData.KappaMode==KAPPAMODE_SIMPLE) return(InitData.KappaV4>0? InitData.KappaV4 :1);
  long double Z0=InitData.KappaV1, p=InitData.KappaV2,tmp1,tmp2;
  long double tmp=0;
  if(InitData.KappaMode==KAPPAMODE_SIN2) {
    tmp=sin(M_PI*z/InitData.dt/InitData.nSteps);
    return tmp*tmp;
  }
  if(InitData.KappaMode==KAPPAMODE_MOVEA) {
    if( ((z-t)<InitData.KappaV1)&&((z-t)>0) ) return 1; else return 0;
  }
  if(InitData.KappaMode==KAPPAMODE_MOVEASIN2) {
    tmp=sin(M_PI*(z-t)/InitData.KappaV1);
    if( ((z-t)<InitData.KappaV1)&&((z-t)>0) ) return tmp*tmp; else return 0;
  }
  if (InitData.KappaMode==KAPPAMODE_TOMSK) {
    if(Z0<0.000001) return( KappaC );
    if(z>Z0) return KappaC;
    tmp1=p*p+z*(1.0-p*p)*KappaZInv;
    if(tmp1<0.1) {
      tmp2=0;
    }
    return( KappaC/sqrt(tmp1));
    return( ( Z0>z) ? KappaC/sqrt(p*p+z*(1.0-p*p)*KappaZInv) : KappaC );
  }
  if (InitData.KappaMode==KAPPAMODE_LINEAR) {
    if(Z0<0.000001) return( KappaC );
    if(KappaC*(1.0+p*(z-Z0))<0) return(0);
    return( ( Z0>z) ? KappaC*(1.0+p*(z-Z0)) : KappaC );
  }
  long double c1,c2;
  int i1,i2;
  if (InitData.KappaMode==KAPPAMODE_LINEARFRACS) {
    dzKappa=InitData.dt*(InitData.nSteps-1)/(InitData.nKappaSteps-1);
    c1=z/dzKappa;
    i1=c1;
    i2=i1+1;
    i2= i2 >= InitData.nKappaSteps ? InitData.nKappaSteps-1 : i2;
    i1= i1 >= InitData.nKappaSteps ? InitData.nKappaSteps-1 : i1;
    c1=c1-i1;
    c2=1-c1;
    /*
    while (i1>=InitData.nKappaSteps) i1--;
    if( i1>=InitData.nKappaSteps )
      i1=InitData.nKappaSteps-1;
    return InitData.hKappaData[i1]+InitData.hKappaData[i2];*/
    return(InitData.hKappaData[i1]*c2+InitData.hKappaData[i2]*c1);
  }
  if (InitData.KappaMode==KAPPAMODE_CONSTFRACS) {
    dzKappa=InitData.dt*(InitData.nSteps-1)/(InitData.nKappaSteps-1);
    i1=z/dzKappa;
    i1= i1 >= InitData.nKappaSteps ? InitData.nKappaSteps-1 : i1;
    return(InitData.hKappaData[i1]);
  }
  long double zFrac=z/(InitData.nKappaSteps-1);
  if (InitData.KappaMode==KAPPAMODE_FOURIER) {
    for(int i=0; i<InitData.nKappaSteps>>1;i++) {
      tmp+=1.0+InitData.hKappaData[2*i]*sin(zFrac*i*M_PI)+
               InitData.hKappaData[2*i+1]*cos(zFrac*i*M_PI);
    }
    return tmp;
  }
  return(1);
}

long double cEqu25DElModel::KappaEXT(long double z,long double t,cDifEqu_InitDataRec_25DElModel *hInitData) {
  if(!hInitData) hInitData=&InitData;
  if (hInitData->KappaMode==KAPPAMODE_SIMPLE)
    return(hInitData->KappaV4>0? hInitData->KappaV4 :1);
  long double Z0=hInitData->KappaV1, p=hInitData->KappaV2,tmp1,tmp2,tmpKappaZInv,tmpKappaC,L=hInitData->dt*hInitData->nSteps;
  long double tmp=0;
  if(hInitData->KappaMode==KAPPAMODE_SIN2) {
    tmp=sin(M_PI*z/hInitData->dt/hInitData->nSteps);
    return tmp*tmp;
  }
  if(hInitData->KappaMode==KAPPAMODE_MOVEA) {
    if( ((z-t)<hInitData->KappaV1)&&((z-t)>0) ) return 1; else return 0;
  }
  if(InitData.KappaMode==KAPPAMODE_MOVEASIN2) {
    tmp=sin(M_PI*(z-t)/hInitData->KappaV1);
    if( ((z-t)<hInitData->KappaV1)&&((z-t)>0) ) return tmp*tmp; else return 0;
  }
  if (hInitData->KappaMode==KAPPAMODE_TOMSK) {
    tmpKappaZInv = (Z0>0.000001) ? 1.0/Z0 : 1.0 ;
    if(hInitData->KappaV4<0) tmpKappaC = L*(1.0+p) / ( (L-Z0)*(1.0+p)+2.0*Z0 );
      else tmpKappaC=hInitData->KappaV4;
    if(Z0<0.000001) return( tmpKappaC );
    if(z>Z0) return tmpKappaC;
    tmp1=p*p+z*(1.0-p*p)*tmpKappaZInv;
    if(tmp1<0.1) {
      tmp2=0;
    }
    return( tmpKappaC/sqrt(tmp1));
    return( ( Z0>z) ? KappaC/sqrt(p*p+z*(1-p*p)*tmpKappaZInv) : tmpKappaC );
  }
  if (hInitData->KappaMode==KAPPAMODE_LINEAR) {
    if(Z0<0.000001) return( KappaC );
    if(KappaC*(1.0+p*(z-Z0))<0) return(0);
    return( ( Z0>z) ? KappaC*(1.0+p*(z-Z0)) : KappaC );
  }
  long double c1,c2;
  int i1,i2;
  long double tmpdzKappa=hInitData->dt*(hInitData->nSteps-1)/(hInitData->nKappaSteps-1);
  if (hInitData->KappaMode==KAPPAMODE_LINEARFRACS) {
    c1=z/tmpdzKappa;
    i1=c1;
    i2=i1+1;
    i2= i2 >= hInitData->nKappaSteps ? hInitData->nKappaSteps-1 : i2;
    i1= i1 >= hInitData->nKappaSteps ? hInitData->nKappaSteps-1 : i1;
    c1=c1-i1;
    c2=1-c1;   /*
    if( i1>=hInitData->nKappaSteps )
      i1=hInitData->nKappaSteps-1;
    return 1;
    while (i1>=hInitData->nKappaSteps) i1--;*/
    return(hInitData->hKappaData[i1]*c2+hInitData->hKappaData[i2]*c1);
  }
  if (hInitData->KappaMode==KAPPAMODE_CONSTFRACS) {
    dzKappa=hInitData->dt*(hInitData->nSteps-1)/(((long double)hInitData->nKappaSteps)+0.001);
    i1=z/dzKappa;
    i1= i1 >= hInitData->nKappaSteps ? hInitData->nKappaSteps-1 : i1;
    return(hInitData->hKappaData[i1]);
  }
  long double zFrac=z/(hInitData->nKappaSteps-1);
  if (hInitData->KappaMode==KAPPAMODE_FOURIER) {
    for(int i=0; i<(hInitData->nKappaSteps>>1);i++) {
      tmp+=1.0+hInitData->hKappaData[2*i]*sin(zFrac*i*M_PI)+
               hInitData->hKappaData[2*i+1]*cos(zFrac*i*M_PI);
    }
    return tmp;
  }
  return(1);
}

long double cEqu25DElModel::AvgKappa(long double t,cDifEqu_InitDataRec_25DElModel *hInitData) {
  long double tmp=-KappaEXT(0,t,hInitData)-KappaEXT(hInitData->dt*hInitData->nSteps,t,hInitData);
  for(int i=0;i<hInitData->nSteps+1;i++)
    tmp+=KappaEXT(hInitData->dt*i,t,hInitData);
  tmp/=hInitData->nSteps;
  return tmp;
}

int cEqu25DElModel::LoadKappaData(cDifEqu_InitDataRec_25DElModel *hData) {
  long double tmpData[10000];
  int nKappaSteps=0;
  if(!hData) hData=&InitData;
  if (hData->KappaMode==KAPPAMODE_LINEARFRACS) {
    ifstream KappaStream("Kappa.dat",ios::in);
    if(hData->fDeleteKappaData) {
      delete [] hData->hKappaData;
      hData->hKappaData=NULL;
    }
    while(1){
      tmpData[nKappaSteps]=-10987654321;
      KappaStream >> tmpData[nKappaSteps];
      KappaStream >> tmpData[nKappaSteps];
      if(tmpData[nKappaSteps]==-10987654321) break;
      nKappaSteps++;
    }
    hData->hKappaData=new long double[nKappaSteps];
    for(int i=0;i<nKappaSteps;i++) hData->hKappaData[i]=tmpData[i];
    hData->nKappaSteps=nKappaSteps;
    InitData.fDeleteKappaData=1;
  }
  return 0;
}

int cEqu25DElModel::LoadDeltaZData(cDifEqu_InitDataRec_25DElModel *hData) {
  long double tmpData[10000];
  int nDeltaZSteps=0;
  if(!hData) hData=&InitData;
  if (hData->DeltaZMode==2) {
    ifstream DeltaZStream("DeltaZ.dat",ios::in);
    if(hData->fDeleteDeltaZData) {
      delete [] hData->hDeltaZData;
      hData->hDeltaZData=NULL;
    }
    while(1){
      tmpData[2*nDeltaZSteps]=-10987654321;
      DeltaZStream >> tmpData[2*nDeltaZSteps];
      DeltaZStream >> tmpData[2*nDeltaZSteps+1];
      if(tmpData[2*nDeltaZSteps]==-10987654321) break;
      nDeltaZSteps++;
    }
    hData->hDeltaZData=new long double[2*nDeltaZSteps];
    for(int i=0;i<2*nDeltaZSteps;i++) hData->hDeltaZData[i]=tmpData[i];
    hData->nDeltaZSteps=nDeltaZSteps;
    InitData.fDeleteDeltaZData=1;
  }
  return 0;
}

void cEqu25DElModel::InitEls(int nElFrac,int rMul)
{
  long double r=rMul*InitData.r;
  if(t<InitData.ElTime) {
    for(int i=0;i<InitData.nEls;i++)
    {
      long double Theta0=2.0*M_PI*i/InitData.nEls;
      hEl[nElFrac*2*InitData.nEls+2*i]=Theta0+r*cos(Theta0);
      hEl[nElFrac*2*InitData.nEls+2*i+1]=0;
      hElPresent[nElFrac*InitData.nEls+i]=1;
    }
    hNMovingEls[nElFrac]=InitData.nEls;
  } else {
    for(int i=0;i<InitData.nEls;i++)
      hElPresent[nElFrac*InitData.nEls+i]=0;
    hNMovingEls[nElFrac]=0;
  }
}

void cEqu25DElModel::FillAData()
{
  TCplxLong ACoef=TCplxLong(0,1);
  ACoef=-sqrt(ACoef/(long double)M_PI);
/*  ACoef=-sqrt(ACoef*((long double)M_PI));
  ACoef=TCplxLong(1,0)/ACoef; */
//  ACoef=1;

  hA[0]=GetAInput();
//  A[1]=hP[0]/sqrt(dt);
  hA[1]=sqrt(InitData.dt)*(hI[0]+hI[1]*(long double)2.0)*(long double)(2.0/3.0);
  hA[1]*=ACoef;
  hA[1]+=GetAInput();
  for(int i=2;i<InitData.nSteps+1;i++)
  {
    hA[i]=0;
    for(int i2=0;i2<i;i2++)
      hA[i]+= ( (i2==0)|(i2==i-1) ) ?
                 hI[i2]*((long double)0.5)/sqrt(InitData.dt*(i-i2)) :
                 hI[i2]/sqrt(InitData.dt*(i-i2));
    hA[i]*=InitData.dt;
    hA[i]+=sqrt(InitData.dt)*(hI[i-1]+hI[i]*((long double)2.0))*(long double)(2.0/3.0);
//    A[i]+=hP[0]/sqrt(dz*i));
    hA[i]*=ACoef;
//    A[i]+=GetAInput();//complex<REAL_TYPE>(0,aInput);
    hA[1]+=GetAInput();
  }
}

void cEqu25DElModel::StepRoutine() {
  if(hDump) hDump->Dump(hA,hI,hEl,hElPresent);
  hOutA[nStepNumber]=hA[InitData.nSteps];//*hA[InitData.nSteps];
  for(int i=0;i<InitData.nSteps+1;i++) {
    nZStepNumber=i;
    for(int j=0;j<InitData.nEls;j++) {
      if(hElPresent[InitData.nEls*i+j]) {
        hElk1[2*InitData.nEls*i+2*j]=
          RightPartEl0(hEl[2*InitData.nEls*i+2*j],hEl[2*InitData.nEls*i+2*j+1],hA[i])*InitData.dt;
        hElk1[2*InitData.nEls*i+2*j+1]=
          RightPartEl1(hEl[2*InitData.nEls*i+2*j],hEl[2*InitData.nEls*i+2*j+1],hA[i])*InitData.dt;
      }
    }

//    hI[i]+=-hA[i]*InitData.dt;
  }
  for(int i=0;i<InitData.nSteps+1;i++) {
    nZStepNumber=i;
    for(int j=0;j<InitData.nEls;j++) {
      if(hElPresent[InitData.nEls*i+j]) {
        hElk2[2*InitData.nEls*i+2*j]=
          RightPartEl0(hEl[2*InitData.nEls*i+2*j]+hElk1[2*InitData.nEls*i+2*j],
                       hEl[2*InitData.nEls*i+2*j+1]+hElk1[2*InitData.nEls*i+2*j+1],
                       hA[i])*InitData.dt;
        hElk2[2*InitData.nEls*i+2*j+1]=
          RightPartEl1(hEl[2*InitData.nEls*i+2*j]+hElk1[2*InitData.nEls*i+2*j],
                       hEl[2*InitData.nEls*i+2*j+1]+hElk1[2*InitData.nEls*i+2*j+1],
                       hA[i])*InitData.dt;
        hEl[2*InitData.nEls*i+2*j]+=
          (hElk1[2*InitData.nEls*i+2*j]+hElk2[2*InitData.nEls*i+2*j])*((long double)0.5);
        hEl[2*InitData.nEls*i+2*j+1]+=(hElk1[2*InitData.nEls*i+2*j+1]+hElk2[2*InitData.nEls*i+2*j+1])*((long double)0.5);
      }
    }
  }
  CalcI();
  StepDataMove();
  CalcI();
  FillAData();
  //hOutA[nStepNumber]=hA[InitData.nSteps];//*hA[InitData.nSteps];
}

/*
long double cEqu25DElModel::Kappa(int zstep, int tstep) {
  if (!InitData.KappaMode) return(1);
  //long double z=zstep*InitData.dt;
  if (InitData.KappaMode==1)
    return( ((InitData.KappaV1-InitData.dt*(zstep-tstep))>0)&&
            (zstep>tstep)? 1 : 0);
  return(1);
} */

void cEqu25DElModel::CalcI() {
TCplxLong TempI=0;
  for(int i=0;i<InitData.nSteps+1;i++) {
    TempI=0;
    for(int j=0;j<InitData.nEls;j++) {
      if(hElPresent[InitData.nEls*i+j]) {
        TempI+=exp(TCplxLong(0,-hEl[2*i*InitData.nEls+2*j]));
      }
    }
    hI[i]=(long double)2.0*Kappa(i,nStepNumber)*TempI/(long double)InitData.nEls;
  }
}

long double cEqu25DElModel::RightPartEl0(long double El0,long double El1,TCplxLong A)
{ return(El1); }

long double cEqu25DElModel::RightPartEl1(long double El0,long double El1,TCplxLong A)
//{ return( imag(A*exp(TCplxLong(0,El0))) ); }
{ return( real(A*exp(TCplxLong(0,El0))*Kappa(nZStepNumber,nStepNumber)) ); }



void cEqu25DElModel::StepDataMove()
{
/*
  AvgGamma=0;
  for(int i=0;i<nEls;i++) if(h25DElPresent[nSteps*nEls+i])
    AvgGamma+=h25DElData[2*nSteps*nEls+2*i+1];
  if(h25DNMovingEls[nSteps]) AvgGamma/=h25DNMovingEls[nSteps];*/
  if(InitData.ElDirection>0) {
    for(int i=InitData.nSteps;i>0;i--) {
      memcpy(hEl+InitData.nEls*i*2,hEl+InitData.nEls*(i-1)*2,InitData.nEls*2*sizeof(long double));
      memcpy(hElPresent+InitData.nEls*i,hElPresent+InitData.nEls*(i-1),InitData.nEls*sizeof(int));
      hNMovingEls[i]=hNMovingEls[i-1];
    }
    InitEls(0,1);
  } else {
    for(int i=0;i<InitData.nSteps;i++) {
      memcpy(hEl+InitData.nEls*i*2,hEl+InitData.nEls*(i+1)*2,InitData.nEls*2*sizeof(long double));
      memcpy(hElPresent+InitData.nEls*i,hElPresent+InitData.nEls*(i+1),InitData.nEls*sizeof(int));
      hNMovingEls[i]=hNMovingEls[i+1];
    }
    InitEls(InitData.nSteps,1);
  }
}

long double cEqu25DElModel::GetMaxA()
{
  long double Res=abs(hA[0]),tmp;
  for(int i=1;i<=InitData.nSteps;i++) {
    tmp=abs(hA[i]);
    Res = Res > tmp ? Res : tmp;
  }
  return( Res );
}

long double cEqu25DElModel::CalcEnergy(){
  long double xbase=InitData.hDeltaData[nStepNumber-1];
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
//    ym2=EACalcMaxIPos(xm2);
    yp1=EACalcMaxIPos(xp1);
//    yp2=EACalcMaxIPos(xp2);
//    ybase=(yp1+yp2+ym1+ym2)*0.25;
    ybase=(yp1+ym1)*0.5;
    k=(yp1-ym1)/(2.0*InitData.EADelta);
    if (k*k>0.0000000001) {
      delta=MaxAPos-ybase;
      delta=delta/k;
      xbase+=delta;
    }
  }
  return(xbase);
}

int cEqu25DElModel::EACalcMaxIPos(long double delta){
  return 0;
}

void cEqu25DElModel::EAMakeDelta() {
  int n;
  if(InitData.UseEnergyAdapt&&(t>InitData.EAStartTime)){
                        // ---------- Energy Adaptation System

    if((t<InitData.EAEndTime)) {
      InitData.hDeltaData[nStepNumber]=CalcEnergy();
      n=((long double)0.5)/dt;
      EACoef=(InitData.hDeltaData[nStepNumber-1]-InitData.hDeltaData[nStepNumber-1-n])/n;
    }
    else {
      InitData.hDeltaData[nStepNumber]=InitData.hDeltaData[nStepNumber-1]+EACoef;
    }
    /*
    if((t>=InitData.EAEndTime)&&(t<InitData.EAEndTime+InitData.dt))
      SaveEnergyAdaptCalibrationData();*/
  }
}
