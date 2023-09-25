//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "Equation02.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TDummyForm2 *DummyForm2;
//---------------------------------------------------------------------------
__fastcall TDummyForm2::TDummyForm2(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------
/*
void cDummyEquation::InitReport(cDifEquReport *hReport) {
  hReport->nGraphs=0;
}
void cDummyEquation::Report(cDifEquReport *hReport) {
  hReport->nGraphs=0;
}
void cDummyEquation::FinalReport(cDifEquReport *hReport) {
  hReport->nGraphs=0;
}
  */
//---------------------------------------------------------------------------

cDummyMaster2::cDummyMaster2() : cDummyMaster() {
  hForm=DummyForm2;
  ReportText.nGraphs=4;
  ReportText.GraphText[0]="Sin";
  ReportText.GraphText[1]="Cos";
  ReportText.GraphText[2]="random1";
  ReportText.GraphText[3]="random2";
}

void cDummyMaster2::KillEquation() {
  if(hEquation) delete hEquation;
}
void cDummyMaster2::PrepareEquation() {
  cDifEqu_InitDataRec InitData;
  InitData.tMax=100;
  InitData.dt=0.0000001; 
  hEquation = new cDummyEquation;
  hEquation->LoadInitData(&InitData);
  hEquation->SetHReportText(&ReportText);
}


