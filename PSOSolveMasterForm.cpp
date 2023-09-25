//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "PSOSolveMasterForm.h"
#include "EquPSO2DBraggSolver.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TPSOSolverForm *PSOSolverForm;
//---------------------------------------------------------------------------
__fastcall TPSOSolverForm::TPSOSolverForm(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

cMasterPSOSOlver::cMasterPSOSOlver() : cEquMaster_Controls() {
  hForm=PSOSolverForm;
  hDumpForm=NULL;
  UsesPSO=1;
  ReportText.nGraphs=2;
  ReportText.GraphText[0]="Аналитич. зависимость";
  ReportText.GraphText[1]="Результат поиска";
  hEquation=NULL;
}

void cMasterPSOSOlver::KillEquation() {
  if(hEquation) delete hEquation;
}
void cMasterPSOSOlver::PrepareEquation() {
  cEquPSO2DBraggSolver *tmpEqu;
  cDifEqu_InitDataRec_PSOSolver InitData;
  InitData.deltaReMin=PSOSolverForm->Edit_deltaReMin->Text.ToDouble();
  InitData.deltaReMax=PSOSolverForm->Edit_deltaReMax->Text.ToDouble();
  InitData.deltaImMin=PSOSolverForm->Edit_deltaImMin->Text.ToDouble();
  InitData.deltaImMax=PSOSolverForm->Edit_deltaImMax->Text.ToDouble();
  InitData.gammaReMin=PSOSolverForm->Edit_gammaReMin->Text.ToDouble();
  InitData.gammaReMax=PSOSolverForm->Edit_gammaReMax->Text.ToDouble();
  InitData.gammaImMin=PSOSolverForm->Edit_gammaImMin->Text.ToDouble();
  InitData.gammaImMax=PSOSolverForm->Edit_gammaImMax->Text.ToDouble();
  InitData.epsilon=PSOSolverForm->Edit_epsilon->Text.ToDouble();
  InitData.Lx=PSOSolverForm->Edit_lx->Text.ToDouble();
  InitData.Lz=PSOSolverForm->Edit_lz->Text.ToDouble();
  InitData.alpha=PSOSolverForm->Edit_alpha->Text.ToDouble();
  InitData.LzMin=PSOSolverForm->Edit_LzMin->Text.ToDouble();
  InitData.nScaleSteps=PSOSolverForm->Edit_ScaleSteps->Text.ToInt();
  InitData.mode=PSOSolverForm->Edit_mode->Text.ToInt();
  InitData.SaveBitmapsInt=PSOSolverForm->Edit_SaveBitmapsInt->Text.ToInt();
  InitData.tnx=PSOSolverForm->Edit_nx->Text.ToInt();
  InitData.tnz=PSOSolverForm->Edit_nz->Text.ToInt();
  InitData.PSOData.nSteps=PSOForm->Edit_nSteps->Text.ToInt();
  InitData.PSOData.nPoints=PSOForm->Edit_nPoints->Text.ToInt();
  InitData.PSOData.Omega=PSOForm->Edit_Omega->Text.ToDouble();
  InitData.PSOData.VCoef=PSOForm->Edit_VCoef->Text.ToDouble();
  InitData.PSOData.c1=PSOForm->Edit_c1->Text.ToDouble();
  InitData.PSOData.c2=PSOForm->Edit_c2->Text.ToDouble();
  InitData.PSOData.c3m=PSOForm->Edit_c3m->Text.ToDouble();
  InitData.PSOData.SpotSize=PSOForm->Edit_SpotSize->Text.ToDouble();
  InitData.PSOData.BlackSpotSize=PSOForm->Edit_BlackSpot->Text.ToDouble();
  InitData.PSOData.BlackListValue=PSOForm->Edit_BlackVal->Text.ToDouble();
  InitData.PSOData.nGroups=PSOForm->Edit_nGroups->Text.ToInt();
  InitData.PSOData.nParked=PSOForm->Edit_nParked->Text.ToInt();
  InitData.PSOData.MaxLifeTime=PSOForm->Edit_LifeTime->Text.ToInt();
  InitData.PSOData.BlackTimeInterval=PSOForm->Edit_BlackTime->Text.ToInt();
  if (hDump) delete hDump;
  hDump = NULL;
  hEquation = new cEquPSO2DBraggSolver;
  hEquation->LoadInitData(&InitData);
}

void cMasterPSOSOlver::AutoSave() {
//  cDump25DEl *hD=(cDump25DEl *)hDump;
  if(!hDump) return;
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
