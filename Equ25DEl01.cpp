//---------------------------------------------------------------------------

#include <vcl.h>
#include <fstream.h>
#pragma hdrstop

#include "Equ25DEl01.h"
#include "Equ25DEl04.h"
#include "Equ25DEl05.h"
#include "Equ2WaveElModel.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TFormEqu25DEl01 *FormEqu25DEl01;
//---------------------------------------------------------------------------
__fastcall TFormEqu25DEl01::TFormEqu25DEl01(TComponent* Owner)
        : TForm(Owner)
{

}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------

cEqu25DEl_optimizer::cEqu25DEl_optimizer() : cDifEqu_DataBankIO()
{
  hEquation=NULL; hVariable=NULL;
}

cEqu25DEl_optimizer::~cEqu25DEl_optimizer()// : ~cDifEqu_DataBankIO()
{
  delete hEquation;  hEquation=NULL;
}

void cEqu25DEl_optimizer::InitReport(cDifEquReport *hReport) {
  int Lines[EQU_MAX_GRAPHS];
  for(int i=0;i<EQU_MAX_GRAPHS;i++) Lines[i]=1;
  Lines[0]=1; Lines[1]=1; Lines[2]=1; Lines[3]=1;
  hReport->nVals=0; hReport->nGraphs=0;
//  hReport->LoadText(hReportText);
  hReport->SetGraphSize(EQU_GRAPH_SIZE,Lines);
  hReport->SetMaxTo1();
}

void cEqu25DEl_optimizer::Report(cDifEquReport *hReport) {
}

void cEqu25DEl_optimizer::FinalReport(cDifEquReport *hReport) {
  hReport->nGraphs=2;
  hReport->GraphText[0]="OutGraph";
  hReport->Graph[0].GetData(0,Var.nSteps,MaxOut);
  hReport->GraphText[1]="MaxTotalGraph";
  hReport->Graph[1].GetData(0,Var.nSteps,MaxTotal);
  hReport->nVals=4;
  hReport->ValText[0]="Max Out A = ";  hReport->Val[0]=MaxOutMax;
  hReport->ValText[1]="Optimum Value = ";  hReport->Val[1]=BestOutValue;
  hReport->ValText[2]="Max A = ";  hReport->Val[2]=MaxTotalMax;
  hReport->ValText[3]="Optimum Value = ";  hReport->Val[3]=BestTotalValue;

  ofstream FileOut = ofstream("MaxOut.dat",ios::out);
  ofstream FileTotal = ofstream("MaxTotal.dat",ios::out);
  long double value;
  for(int i=0;i<Var.nSteps;i++) {
    value=Var.V1+(Var.V2-Var.V1)*i/(Var.nSteps-1);
    FileOut << value << " " << MaxOut[i] << endl;
    FileTotal << value << " " << MaxTotal[i] << endl;
  }
}

void cEqu25DEl_optimizer::LoadInitData(cDifEqu_InitDataRec *hInitData)
{
  cDifEqu_InitDataRec_25DElModel *hInit=(cDifEqu_InitDataRec_25DElModel *)hInitData;
  InitData=*hInit;
//  hVariable=&InitData.ElTime;
//  hVariable=&InitData.deltaI;
//  hVariable=&InitData.EAkz;
//  hVariable=&InitData.deltaV1;
//  hVariable=&InitData.AInputV3;
  hVariable=&InitData.deltaV2;
  nTSteps=(InitData.tMax-0.00000001)/InitData.dt;
  t=0; tMax=nTSteps*Var.nSteps-0.5; dt=1;
  tStep=0; vStep=0;
  *hVariable=Var.V1; BestOutValue=BestTotalValue=Var.V1;
  BestOutTime=BestTotalTime=0;
  MaxOut[0]=MaxTotal[0]=0; MaxOutMax=MaxTotalMax=0;
  hEquation->LoadInitData(&InitData);
}

void cEqu25DEl_optimizer::StepRoutine() {
  hEquation->Step();
  tStep++;
  if(hEquation->GetOutA()>MaxOut[vStep]) {
    MaxOut[vStep]=hEquation->GetOutA();
    if(MaxOut[vStep]>MaxOutMax) {
      MaxOutMax=MaxOut[vStep];
      BestOutValue=*hVariable;
      BestOutTime=tStep*InitData.dt;
    }
  }
  long double tmpMax=hEquation->GetMaxA();
  if(tmpMax>MaxTotal[vStep]) {
    MaxTotal[vStep]=tmpMax;
    if(MaxTotal[vStep]>MaxTotalMax) {
      MaxTotalMax=MaxTotal[vStep];
      BestTotalValue=*hVariable;
      BestTotalTime=tStep*InitData.dt;
    }
  }
  if(tStep==nTSteps) {
    tStep=0; vStep++; MaxTotal[vStep]=MaxOut[vStep]=0;
    *hVariable=Var.V1+(Var.V2-Var.V1)*vStep/(Var.nSteps-1);
    hEquation->LoadInitData(&InitData);
  }
}

//---------------------------------------------------------------------------

cMaster25DBasic::cMaster25DBasic() : cEquMaster_Controls() {
  hForm=FormEqu25DEl01;
  hDumpForm=Dump25DEl;
  UsesPSO=1;
  ReportText.nGraphs=5;
  ReportText.GraphText[0]="abs A";
  ReportText.GraphText[1]="abs I";
  ReportText.GraphText[2]="Average Gamma";
  ReportText.GraphText[3]="Electron Energy";
  ReportText.GraphText[4]="Starting I phase";
  hEquation=NULL;
}

void cMaster25DBasic::KillEquation() {
  if(hEquation) delete hEquation;
}
void cMaster25DBasic::PrepareEquation() {
  cEqu25DElModel *tmpEqu;
  cDifEqu_InitDataRec_25DElModel InitData;
  InitData.tMax=FormEqu25DEl01->Edit_tMax->Text.ToDouble();
  InitData.L=FormEqu25DEl01->Edit_L->Text.ToDouble();
  InitData.nSteps=FormEqu25DEl01->Edit_nSteps->Text.ToInt();
  InitData.r=FormEqu25DEl01->Edit_r->Text.ToDouble();
  InitData.dt=InitData.L/InitData.nSteps;
  InitData.rTime=FormEqu25DEl01->Edit_rTime->Text.ToDouble();
  InitData.nEls=FormEqu25DEl01->Edit_nEls->Text.ToDouble();
  InitData.ICoef=FormEqu25DEl01->Edit_ICoef->Text.ToDouble();
  InitData.AI=FormEqu25DEl01->Edit_AI->Text.ToDouble();
  InitData.deltaI=FormEqu25DEl01->Edit_deltaI->Text.ToDouble();                 
  InitData.deltaS=FormEqu25DEl01->Edit_deltaS->Text.ToDouble();
  InitData.rndCoef=FormEqu25DEl01->Edit_rndCoef->Text.ToDouble();
  InitData.Ai_kappa=FormEqu25DEl01->Edit_Ai_kappa->Text.ToDouble();
  InitData.nSpeedFracs=FormEqu25DEl01->Edit_nSpeedFracs->Text.ToInt();
  InitData.SpeedDistribFactor=FormEqu25DEl01->Edit_SpeedDistrFactor->Text.ToDouble();
  InitData.MoveAData = FormEqu25DEl01->Check_MoveAData->Checked ? 1 : 0;
  InitData.OutMaxA = FormEqu25DEl01->Check_maxA->Checked ? 1 : 0;
  InitData.NoRandom = InitData.rndCoef>=0 ? 0 : 1;
  InitData.rndCoef =  InitData.rndCoef>=0 ? InitData.rndCoef : -InitData.rndCoef;
  InitData.ElMode = FormEqu25DEl01->Edit_ElMode->Text.ToInt();
  InitData.DeltaMode = 0;
  InitData.DeltaZMode = FormEqu25DEl01->Edit_DeltaZMode->Text.ToInt();
  InitData.nDeltaZSteps = FormEqu25DEl01->Edit_nDeltaZSteps->Text.ToInt();
  InitData.DeltaZMax = FormEqu25DEl01->Edit_DeltaZMax->Text.ToDouble();
  InitData.KappaMode = FormEqu25DEl01->Edit_KappaMode->Text.ToInt();
  InitData.KappaV1=FormEqu25DEl01->Edit_KappaV1->Text.ToDouble();
  InitData.KappaV2=FormEqu25DEl01->Edit_KappaV2->Text.ToDouble();
  InitData.KappaV3=FormEqu25DEl01->Edit_KappaV3->Text.ToDouble();
  InitData.KappaV4=FormEqu25DEl01->Edit_KappaV4->Text.ToDouble();
  InitData.KappaV5=FormEqu25DEl01->Edit_KappaV5->Text.ToDouble();
  InitData.KappaV6=FormEqu25DEl01->Edit_KappaV6->Text.ToDouble();
  InitData.KappaV7=FormEqu25DEl01->Edit_KappaV7->Text.ToDouble();
  InitData.KappaV8=FormEqu25DEl01->Edit_KappaV8->Text.ToDouble();
  InitData.KappaV9=FormEqu25DEl01->Edit_KappaV9->Text.ToDouble();
  InitData.KappaV10=FormEqu25DEl01->Edit_KappaV10->Text.ToDouble();
  InitData.KappaV11=FormEqu25DEl01->Edit_KappaV11->Text.ToDouble();
  InitData.deltaV1=FormEqu25DEl01->Edit_DeltaV1->Text.ToDouble();
  InitData.deltaV2=FormEqu25DEl01->Edit_DeltaV2->Text.ToDouble();
  InitData.AInputMode=FormEqu25DEl01->Edit_AInputMode->Text.ToInt();
  InitData.AInputV1=FormEqu25DEl01->Edit_AInputV1->Text.ToDouble();
  InitData.AInputV2=FormEqu25DEl01->Edit_AInputV2->Text.ToDouble();
  InitData.AInputV3=FormEqu25DEl01->Edit_AInputV3->Text.ToDouble();
  InitData.Nu=FormEqu25DEl01->Edit_Nu->Text.ToDouble();
  InitData.ElDirection=1;
  InitData.ElTime=FormEqu25DEl01->Edit_ElTime->Text.ToDouble();
  InitData.ElsOnStart= FormEqu25DEl01->Check_ElsOnStart->Checked ? 1 : 0;
  InitData.UseSimplePhaseModel=FormEqu25DEl01->Check_Simple2Phase->Checked;
  InitData.UseEnergyAdapt=FormEqu25DEl01->Check_UseEA->Checked;
  InitData.fLoadKappaData=FormEqu25DEl01->Check_LoadKappaData->Checked;
  InitData.fLoadDeltaZData=FormEqu25DEl01->Check_LoadDeltaZData->Checked;
  InitData.fUseLOVModel=FormEqu25DEl01->Check_LOVModel->Checked;
  InitData.fNeed1stPeak=FormEqu25DEl01->Check_1stPeak->Checked;
  InitData.EAStartTime=FormEqu25DEl01->Edit_EAStartTime->Text.ToDouble();
  InitData.EAEndTime=FormEqu25DEl01->Edit_EAEndTime->Text.ToDouble();
  InitData.EAnSteps=FormEqu25DEl01->Edit_EAnSteps->Text.ToInt();
  InitData.EADelta=FormEqu25DEl01->Edit_EADelta->Text.ToDouble();
  InitData.EADeltaZ0=FormEqu25DEl01->Edit_EA_z0->Text.ToDouble();
  InitData.EAkz=FormEqu25DEl01->Edit_EAkz->Text.ToDouble();
  InitData.TSliceSaveTime=FormEqu25DEl01->Edit_TSliceSave->Text.ToDouble();
  InitData.TSliceSaveZ=FormEqu25DEl01->Edit_TSliceSaveZ->Text.ToDouble();
  InitData.PSOMode=PSOForm->Edit_PSOMode->Text.ToInt();
  InitData.PSOData.nSteps=PSOForm->Edit_nSteps->Text.ToInt();
  InitData.PSOData.nPoints=PSOForm->Edit_nPoints->Text.ToInt();
  InitData.PSOData.Omega=PSOForm->Edit_Omega->Text.ToDouble();
  InitData.PSOData.VCoef=PSOForm->Edit_VCoef->Text.ToDouble();
  InitData.PSOData.c1=PSOForm->Edit_c1->Text.ToDouble();
  InitData.PSOData.c2=PSOForm->Edit_c2->Text.ToDouble();
  InitData.PSOData.c3m=PSOForm->Edit_c3m->Text.ToDouble();
  InitData.PSOData.nParked=PSOForm->Edit_nParked->Text.ToInt();
  InitData.PSOData.MaxLifeTime=PSOForm->Edit_LifeTime->Text.ToInt();
  cOptimizingVariation OptiVar;
  OptiVar.V1=FormEqu25DEl01->Edit_oV1->Text.ToDouble();
  OptiVar.V2=FormEqu25DEl01->Edit_oV2->Text.ToDouble();
  OptiVar.nSteps=FormEqu25DEl01->Edit_onSteps->Text.ToInt();
  if(FormEqu25DEl01->Check_PSO->Checked) {
    if (hDump) delete hDump;
    hDump = NULL;
    hEquation = new cEqu25DEl_PSOOptimizer;
    //tmpEqu = new cEqu25DElModel;
    tmpEqu = new cEquVerticalASpeedDistr;
    ((cEqu25DEl_PSOOptimizer *)hEquation)->SetHEquation((cEqu25DElModel *)tmpEqu);
    hEquation->LoadInitData(&InitData);
    return;
  }
  if(!FormEqu25DEl01->Check_optim->Checked) {
    if (hDump) delete hDump; hDump=NULL;
    if(FormEqu25DEl01->Check_Dump->Checked) {
      hDump= new cDump25DEl;
      ((cDump25DEl *)hDump)->LoadInitData(&InitData);
    }
    ((TDump25DEl *)hDumpForm)->SetHDump((cDump25DEl *)hDump);
    if(FormEqu25DEl01->Check_VertA->Checked) {
      hEquation = new cEquVerticalASpeedDistr;
      hEquation->LoadInitData(&InitData);
      hEquation->SetHReportText(&ReportText);
      ((cEqu25DElModel *)hEquation)->SetHDump((cDump25DEl *)hDump);
      return;
    }
    if(!FormEqu25DEl01->Check_2Phase->Checked) {
      hEquation = new cEqu25DElModel;
      hEquation->LoadInitData(&InitData);
      hEquation->SetHReportText(&ReportText);
      ((cEqu25DElModel *)hEquation)->SetHDump((cDump25DEl *)hDump);
    } else {
      hEquation = new cEqu2PhaseElModel;
      hEquation->LoadInitData(&InitData);
      hEquation->SetHReportText(&ReportText);
      ((cEqu2PhaseElModel *)hEquation)->SetHDump((cDump25DEl *)hDump);
    }
  } else {
    if (hDump) delete hDump;
    hDump = NULL;
    hEquation = new cEqu25DEl_optimizer;
    //tmpEqu = new cEqu25DElModel;
    //tmpEqu = new cEqu2PhaseElModel;
    tmpEqu = new cEquVerticalASpeedDistr;
    ((cEqu25DEl_optimizer *)hEquation)->SetHEquation((cEqu25DElModel *)tmpEqu);
    ((cEqu25DEl_optimizer *)hEquation)->LoadOptimizationData(OptiVar);
    hEquation->LoadInitData(&InitData);
  }
}

void cMaster25DBasic::AutoSave() {
//  cDump25DEl *hD=(cDump25DEl *)hDump;
  if(!hDump) return;
  if(hDump) ((cDump25DEl *)hDump)->SaveCURR("");
/*  TCplxLong A0,A1;
  int zs=10,xs=0,ts=10;
  hD->MakeXProfileCoef(0.1,100);
  for(int i=0;i<20;i++) {
    zs=i+10;
    zs=1; ts=1; xs=0;
    A0=hD->hA[(hD->Init.nSteps+1)*ts+zs];
    A1=hD->GetAVal(zs,xs,ts);
  } */
}
//---------------------------------------------------------------------------

cMaster2WaveModel::cMaster2WaveModel() : cEquMaster_Controls() {
  hForm=FormEqu25DEl01;
  hDumpForm=Dump25DEl;
  UsesPSO=1;
  ReportText.nGraphs=7;
  ReportText.GraphText[0]="A+(z)";
  ReportText.GraphText[1]="A-(z)";
  ReportText.GraphText[2]="A-(t,z=L)";
  ReportText.GraphText[3]="I(z)";
  ReportText.GraphText[4]="Arg(A+(z))";
  ReportText.GraphText[5]="Arg(A-(t,z=L))";
  ReportText.GraphText[6]="A+(t,z=L)";
  /*
  ReportText.GraphText[4]="Spline Re";
  ReportText.GraphText[5]="Spline Im";
  ReportText.GraphText[6]="Line Re";
  ReportText.GraphText[7]="Line Im";
  */
  hEquation=NULL;
}

void cMaster2WaveModel::KillEquation() {
  if(hEquation) delete hEquation;
}
void cMaster2WaveModel::PrepareEquation() {
  cEqu25DElModel *tmpEqu;
  cDifEqu_InitDataRec_2WaveElModel InitData;
  InitData.tMax=FormEqu25DEl01->Edit_tMax->Text.ToDouble();
  InitData.L=FormEqu25DEl01->Edit_L->Text.ToDouble();
  InitData.nSteps=FormEqu25DEl01->Edit_nSteps->Text.ToInt();
  InitData.RandSeed=FormEqu25DEl01->Edit_RandSeed->Text.ToInt();
  InitData.r=FormEqu25DEl01->Edit_r->Text.ToDouble();
  InitData.dt=InitData.L/InitData.nSteps;
  InitData.rTime=FormEqu25DEl01->Edit_rTime->Text.ToDouble();
  InitData.nEls=FormEqu25DEl01->Edit_nEls->Text.ToDouble();
  InitData.ICoef=FormEqu25DEl01->Edit_ICoef->Text.ToDouble();
  InitData.AI=FormEqu25DEl01->Edit_AI->Text.ToDouble();
  InitData.deltaI=FormEqu25DEl01->Edit_deltaI->Text.ToDouble();                 
  InitData.deltaS=FormEqu25DEl01->Edit_deltaS->Text.ToDouble();
  InitData.rndCoef=FormEqu25DEl01->Edit_rndCoef->Text.ToDouble();
  InitData.Ai_kappa=FormEqu25DEl01->Edit_Ai_kappa->Text.ToDouble();
  InitData.nSpeedFracs=FormEqu25DEl01->Edit_nSpeedFracs->Text.ToInt();
  InitData.SpeedDistribFactor=FormEqu25DEl01->Edit_SpeedDistrFactor->Text.ToDouble();
  InitData.MoveAData = FormEqu25DEl01->Check_MoveAData->Checked ? 1 : 0;
  InitData.OutMaxA = FormEqu25DEl01->Check_maxA->Checked ? 1 : 0;
  InitData.NoRandom = InitData.rndCoef>=0 ? 0 : 1;
  InitData.rndCoef =  InitData.rndCoef>=0 ? InitData.rndCoef : -InitData.rndCoef;
  InitData.ElMode = FormEqu25DEl01->Edit_ElMode->Text.ToInt();
  InitData.DeltaMode = 0;
  InitData.DeltaZMode = FormEqu25DEl01->Edit_DeltaZMode->Text.ToInt();
  InitData.nDeltaZSteps = FormEqu25DEl01->Edit_nDeltaZSteps->Text.ToInt();
  InitData.DeltaZMax = FormEqu25DEl01->Edit_DeltaZMax->Text.ToDouble();
  InitData.KappaMode = FormEqu25DEl01->Edit_KappaMode->Text.ToInt();
  InitData.KappaV1=FormEqu25DEl01->Edit_KappaV1->Text.ToDouble();
  InitData.KappaV2=FormEqu25DEl01->Edit_KappaV2->Text.ToDouble();
  InitData.KappaV3=FormEqu25DEl01->Edit_KappaV3->Text.ToDouble();
  InitData.KappaV4=FormEqu25DEl01->Edit_KappaV4->Text.ToDouble();
  InitData.KappaV5=FormEqu25DEl01->Edit_KappaV5->Text.ToDouble();
  InitData.KappaV6=FormEqu25DEl01->Edit_KappaV6->Text.ToDouble();
  InitData.KappaV7=FormEqu25DEl01->Edit_KappaV7->Text.ToDouble();
  InitData.KappaV8=FormEqu25DEl01->Edit_KappaV8->Text.ToDouble();
  InitData.KappaV9=FormEqu25DEl01->Edit_KappaV9->Text.ToDouble();
  InitData.KappaV10=FormEqu25DEl01->Edit_KappaV10->Text.ToDouble();
  InitData.KappaV11=FormEqu25DEl01->Edit_KappaV11->Text.ToDouble();
  InitData.deltaV1=FormEqu25DEl01->Edit_DeltaV1->Text.ToDouble();
  InitData.deltaV2=FormEqu25DEl01->Edit_DeltaV2->Text.ToDouble();
  InitData.AInputMode=FormEqu25DEl01->Edit_AInputMode->Text.ToInt();
  InitData.AInputV1=FormEqu25DEl01->Edit_AInputV1->Text.ToDouble();
  InitData.AInputV2=FormEqu25DEl01->Edit_AInputV2->Text.ToDouble();
  InitData.AInputV3=FormEqu25DEl01->Edit_AInputV3->Text.ToDouble();
  InitData.Nu=FormEqu25DEl01->Edit_Nu->Text.ToDouble();
  InitData.ElDirection=1;
  InitData.ElTime=FormEqu25DEl01->Edit_ElTime->Text.ToDouble();
  InitData.ElsOnStart= FormEqu25DEl01->Check_ElsOnStart->Checked ? 1 : 0;
  InitData.UseSimplePhaseModel=FormEqu25DEl01->Check_Simple2Phase->Checked;
  InitData.UseEnergyAdapt=FormEqu25DEl01->Check_UseEA->Checked;
  InitData.fLoadKappaData=FormEqu25DEl01->Check_LoadKappaData->Checked;
  InitData.fLoadDeltaZData=FormEqu25DEl01->Check_LoadDeltaZData->Checked;
  InitData.fUseLOVModel=FormEqu25DEl01->Check_LOVModel->Checked;
  InitData.fNeed1stPeak=FormEqu25DEl01->Check_1stPeak->Checked;
  InitData.EAStartTime=FormEqu25DEl01->Edit_EAStartTime->Text.ToDouble();
  InitData.EAEndTime=FormEqu25DEl01->Edit_EAEndTime->Text.ToDouble();
  InitData.EAnSteps=FormEqu25DEl01->Edit_EAnSteps->Text.ToInt();
  InitData.EADelta=FormEqu25DEl01->Edit_EADelta->Text.ToDouble();
  InitData.EADeltaZ0=FormEqu25DEl01->Edit_EA_z0->Text.ToDouble();
  InitData.EAkz=FormEqu25DEl01->Edit_EAkz->Text.ToDouble();
  InitData.TSliceSaveTime=FormEqu25DEl01->Edit_TSliceSave->Text.ToDouble();
  InitData.TSliceSaveZ=FormEqu25DEl01->Edit_TSliceSaveZ->Text.ToDouble();
  InitData.nThreads=FormEqu25DEl01->Edit_nThreads->Text.ToInt();
  InitData.ACoef=FormEqu25DEl01->Edit_ACoef->Text.ToDouble();
  
  InitData.PSOMode=PSOForm->Edit_PSOMode->Text.ToInt();
  InitData.PSOData.nSteps=PSOForm->Edit_nSteps->Text.ToInt();
  InitData.PSOData.nPoints=PSOForm->Edit_nPoints->Text.ToInt();
  InitData.PSOData.Omega=PSOForm->Edit_Omega->Text.ToDouble();
  InitData.PSOData.VCoef=PSOForm->Edit_VCoef->Text.ToDouble();
  InitData.PSOData.c1=PSOForm->Edit_c1->Text.ToDouble();
  InitData.PSOData.c2=PSOForm->Edit_c2->Text.ToDouble();
  InitData.PSOData.c3m=PSOForm->Edit_c3m->Text.ToDouble();
  InitData.PSOData.nParked=PSOForm->Edit_nParked->Text.ToInt();
  InitData.PSOData.MaxLifeTime=PSOForm->Edit_LifeTime->Text.ToInt();
  cOptimizingVariation OptiVar;
  OptiVar.V1=FormEqu25DEl01->Edit_oV1->Text.ToDouble();
  OptiVar.V2=FormEqu25DEl01->Edit_oV2->Text.ToDouble();
  OptiVar.nSteps=FormEqu25DEl01->Edit_onSteps->Text.ToInt();

  hEquation = new cEqu2WaveElModel;
  hEquation->LoadInitData(&InitData);
  hEquation->SetHReportText(&ReportText);
}

void cMaster2WaveModel::AutoSave() {
}


