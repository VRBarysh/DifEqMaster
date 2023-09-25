//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "OpticParaBasicEqu.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TFormOpticBasic *FormOpticBasic;
//---------------------------------------------------------------------------
__fastcall TFormOpticBasic::TFormOpticBasic(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------
cEquOpticBasic::cEquOpticBasic() : cDifEqu_DataBankIO() {
  A=NULL; P=Pk1=Pk2=NULL; R=Rk1=Rk2=NULL; InitData.nSteps=0;
}

cEquOpticBasic::~cEquOpticBasic() {
  if(A) delete [] A;
  if(P) delete [] P; if(Pk1) delete [] Pk1; if(Pk2) delete [] Pk2;
  if(R) delete [] R; if(Rk1) delete [] Rk1; if(Rk2) delete [] Rk2;
}


void cEquOpticBasic::InitReport(cDifEquReport *hReport) {
  int Lines[EQU_MAX_GRAPHS];
  for(int i=0;i<EQU_MAX_GRAPHS;i++) Lines[i]=1;
  Lines[0]=1; Lines[1]=1; Lines[2]=1; Lines[3]=1;
  hReport->LoadText(hReportText);
  hReport->SetGraphSize(EQU_GRAPH_SIZE,Lines);
  hReport->SetMaxTo1();
}
void cEquOpticBasic::Report(cDifEquReport *hReport) {
  hReport->Graph[0].GetDataAbs(0,InitData.nSteps,A);
  hReport->Graph[1].GetDataAbs(0,InitData.nSteps,P);
  hReport->Graph[2].GetData(0,InitData.nSteps,R);
  hReport->nVals=3;
  hReport->ValText[0]="Max A "; hReport->Val[0]=hReport->Graph[0].Max;
  hReport->ValText[1]="Max P "; hReport->Val[1]=hReport->Graph[1].Max;
  hReport->ValText[2]="Min R "; hReport->Val[2]=hReport->Graph[2].Min;
}
void cEquOpticBasic::FinalReport(cDifEquReport *hReport) {
  hReport->nVals=0;
}

void cEquOpticBasic::LoadInitData(cDifEqu_InitDataRec *hInitData) {
  cDifEqu_InitDataRec_OptiBasic *hInit=(cDifEqu_InitDataRec_OptiBasic *)hInitData;
  tMax=hInitData->tMax; dt=hInitData->dt; t=0;
  if(InitData.nSteps!=hInit->nSteps) {
    if(A) delete [] A; A = new TCplxLong[hInit->nSteps];
    if(P) delete [] P; P = new TCplxLong[hInit->nSteps];
    if(Pk1) delete [] Pk1; Pk1 = new TCplxLong[hInit->nSteps];
    if(Pk2) delete [] Pk2; Pk2 = new TCplxLong[hInit->nSteps];
//    if(PNew) delete [] PNew; PNew = new TCplxLong[hInit->nSteps];
    if(R) delete [] R; R = new long double[hInit->nSteps];
    if(Rk1) delete [] Rk1; Rk1 = new long double[hInit->nSteps];
    if(Rk2) delete [] Rk2; Rk2 = new long double[hInit->nSteps];
//    if(RNew) delete [] RNew; RNew = new long double[hInit->nSteps];
  }
  InitData=*hInit;
  for(int i=0;i<InitData.nSteps;i++) {
    A[i]=0; P[i]=InitData.rPolNoise; R[i]=1;
  }
}

void cEquOpticBasic::FillAData(TCplxLong *hP)
{
  TCplxLong ACoef=TCplxLong(0,1);
  ACoef=sqrt(-ACoef/(long double)M_PI);
  A[0]=0;
//  A[1]=hP[0]/sqrt(dt);
  A[1]=sqrt(InitData.dz)*(hP[0]+hP[1]*(long double)2.0)*(long double)(2.0/3.0);
  A[1]*=ACoef;
  for(int i=2;i<InitData.nSteps;i++)
  {
    A[i]=0;
    for(int i2=0;i2<i;i2++)
      A[i]+= ( (i2==0)|(i2==i-1) ) ?
                 hP[i2]*((long double)0.5)/sqrt(InitData.dz*(i-i2)) :
                 hP[i2]/sqrt(InitData.dz*(i-i2));
    A[i]*=InitData.dz;
    A[i]+=sqrt(InitData.dz)*(hP[i-1]+hP[i]*((long double)2.0))*(long double)(2.0/3.0);
//    A[i]+=hP[0]/sqrt(dz*i));
    A[i]*=ACoef;
//    A[i]+=GetAInput();//complex<REAL_TYPE>(0,aInput);
  }
}

TCplxLong cEquOpticBasic::RightPartP(TCplxLong P,long double R, TCplxLong A) {
  return(-P*InitData.alpha+A*TCplxLong(0,1)*R*InitData.beta);
}

long double cEquOpticBasic::RightPartR(TCplxLong P,long double R, TCplxLong A){
  return(-(R-1)*InitData.alpha1+(A*conj(P)).imag());
}

void cEquOpticBasic::StepRoutine() {
  for(int i=0;i<InitData.nSteps;i++) {
    Pk1[i]=RightPartP(P[i],R[i],A[i])*InitData.dz;
    Rk1[i]=RightPartR(P[i],R[i],A[i])*InitData.dz;
  }
  for(int i=0;i<InitData.nSteps;i++) {
    Pk2[i]=RightPartP(P[i]+Pk1[i],R[i]+Rk1[i],A[i])*InitData.dz;
    Rk2[i]=RightPartR(P[i]+Pk1[i],R[i]+Rk1[i],A[i])*InitData.dz;
    P[i]+=(Pk1[i]+Pk2[i])*((long double)0.5);
    R[i]+=(Rk1[i]+Rk2[i])*((long double)0.5);
  }
  FillAData(P);
}

//---------------------------------------------------------------------------

cEquOpticRealTimeAxis::cEquOpticRealTimeAxis() : cEquOpticBasic() {
  PStorage=NULL; iPStorBegin=iPStorEnd=0;
}

cEquOpticRealTimeAxis::~cEquOpticRealTimeAxis() {
  if(A) delete [] A;
  if(PStorage) delete [] PStorage; if(Pk1) delete [] Pk1; if(Pk2) delete [] Pk2;
  if(R) delete [] R; if(Rk1) delete [] Rk1; if(Rk2) delete [] Rk2;
}

void cEquOpticRealTimeAxis::LoadInitData(cDifEqu_InitDataRec *hInitData) {
  cDifEqu_InitDataRec_OptiBasic *hInit=(cDifEqu_InitDataRec_OptiBasic *)hInitData;
  tMax=hInitData->tMax; dt=hInitData->dt; t=0;
  if(InitData.nSteps!=hInit->nSteps) {
    if(A) delete [] A; A = new TCplxLong[hInit->nSteps];
    if(PStorage) delete [] PStorage; PStorage = new TCplxLong[hInit->nSteps*(hInit->nSteps+1)];
    if(Pk1) delete [] Pk1; Pk1 = new TCplxLong[hInit->nSteps];
    if(Pk2) delete [] Pk2; Pk2 = new TCplxLong[hInit->nSteps];
//    if(PNew) delete [] PNew; PNew = new TCplxLong[hInit->nSteps];
    if(R) delete [] R; R = new long double[hInit->nSteps];
    if(Rk1) delete [] Rk1; Rk1 = new long double[hInit->nSteps];
    if(Rk2) delete [] Rk2; Rk2 = new long double[hInit->nSteps];
//    if(RNew) delete [] RNew; RNew = new long double[hInit->nSteps];
  }
  InitData=*hInit;
  iPStorBegin=iPStorEnd=0;
  P=PStorage;
  for(int i=0;i<InitData.nSteps;i++) {
    A[i]=0; P[i]=InitData.rPolNoise; R[i]=1;
  }
}

void cEquOpticRealTimeAxis::FillAData(TCplxLong *hP)
{
  TCplxLong ACoef=TCplxLong(0,1);
  int iEnd=0;
  ACoef=sqrt(-ACoef/(long double)M_PI);
  for(int i=0;i<InitData.nSteps;i++) {
    iEnd=i-( i>StoredLines() ? StoredLines() : i );
    A[i]=0;
    if(iEnd==i) continue;
    if(iEnd-i>1)
      for(int i2=iEnd;i2<i;i2++)
        A[i]+= ( (i2==iEnd)|(i2==i-1) ) ?
                 PStor(i-i2,i2)*((long double)0.5)/sqrt(InitData.dz*(i-i2)) :
                 PStor(i-i2,i2)/sqrt(InitData.dz*(i-i2));
    A[i]+=sqrt(InitData.dz)*(PStor(1,i-1)+PStor(0,i)*(long double)2.0)*(long double)(2.0/3.0);
    A[i]*=ACoef;
  }
}

void cEquOpticRealTimeAxis::StepRoutine() {
TCplxLong *P1=PStorage+InitData.nSteps*
            ( (iPStorBegin+InitData.nSteps-1)%InitData.nSteps );
  for(int i=0;i<InitData.nSteps;i++) {
    Pk1[i]=RightPartP(P[i],R[i],A[i])*InitData.dz;
    Rk1[i]=RightPartR(P[i],R[i],A[i])*InitData.dz;
  }
  for(int i=0;i<InitData.nSteps;i++) {
    Pk2[i]=RightPartP(P[i]+Pk1[i],R[i]+Rk1[i],A[i])*InitData.dz;
    Rk2[i]=RightPartR(P[i]+Pk1[i],R[i]+Rk1[i],A[i])*InitData.dz;
    P1[i]=P[i]+(Pk1[i]+Pk2[i])*((long double)0.5);
    R[i]+=(Rk1[i]+Rk2[i])*((long double)0.5);
  }
  P=P1;
  iPStorBegin=(iPStorBegin+InitData.nSteps-1) % InitData.nSteps;
  if(iPStorBegin==iPStorEnd)
    iPStorEnd=(iPStorEnd+InitData.nSteps-1) % InitData.nSteps;
  FillAData(P);
}

//---------------------------------------------------------------------------

cMasterOpticBasic::cMasterOpticBasic() : cEquMaster_Controls() {
  hForm=FormOpticBasic;
  ReportText.nGraphs=3;
  ReportText.GraphText[0]="abs A";
  ReportText.GraphText[1]="abs P";
  ReportText.GraphText[2]="R";
}

void cMasterOpticBasic::KillEquation() {
  if(hEquation) delete hEquation;
}
void cMasterOpticBasic::PrepareEquation() {
  cDifEqu_InitDataRec_OptiBasic InitData;
  InitData.tMax=FormOpticBasic->Edit_tMax->Text.ToDouble();
  InitData.dt=FormOpticBasic->Edit_dt->Text.ToDouble();
  InitData.nSteps=FormOpticBasic->Edit_nSteps->Text.ToInt();
  InitData.rPolNoise=FormOpticBasic->Edit_PNoise->Text.ToDouble();
  InitData.dz=FormOpticBasic->Edit_dz->Text.ToDouble();
  InitData.alpha=FormOpticBasic->Edit_alpha->Text.ToDouble();
  InitData.alpha1=FormOpticBasic->Edit_alpha1->Text.ToDouble();
  InitData.beta=FormOpticBasic->Edit_beta->Text.ToDouble();
  if(FormOpticBasic->CheckBox1->Checked) hEquation = new cEquOpticRealTimeAxis;
    else hEquation = new cEquOpticBasic;
  hEquation->LoadInitData(&InitData);
  hEquation->SetHReportText(&ReportText);
}
//---------------------------------------------------------------------------
