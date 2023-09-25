//---------------------------------------------------------------------------

#include <vcl.h>
#include <math.h>
#pragma hdrstop

#include "Equation01.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TDummyForm *DummyForm;
//---------------------------------------------------------------------------
__fastcall TDummyForm::TDummyForm(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------

void cDummyEquation::InitReport(cDifEquReport *hReport) {
  int Lines[EQU_MAX_GRAPHS];
  for(int i=0;i<EQU_MAX_GRAPHS;i++) Lines[i]=1;
  Lines[0]=1; Lines[1]=2; Lines[2]=1; Lines[3]=1;
  hReport->LoadText(hReportText);
  hReport->SetGraphSize(EQU_GRAPH_SIZE,Lines);
  hReport->SetMaxTo1();
}
void cDummyEquation::Report(cDifEquReport *hReport) {
  double Data[300];
  for(int i=0;i<300;i++) Data[i]=sin(t+2*M_PI*i/299);
  hReport->Graph[0].GetData(0,300,Data);
  hReport->Graph[1].GetData(1,300,Data);
  for(int i=0;i<300;i++) Data[i]=cos(t+2*M_PI*i/299);
  hReport->Graph[1].GetData(0,300,Data);
  for(int i=0;i<300;i++) Data[i]=random(65535)*cos(t+2*M_PI*i/299);
  hReport->Graph[2].GetData(0,300,Data);
  for(int i=0;i<300;i++) Data[i]=random(65535)*cos(t+2*M_PI*i/299);
  hReport->Graph[3].GetData(0,300,Data);
}
void cDummyEquation::FinalReport(cDifEquReport *hReport) {
}

//---------------------------------------------------------------------------

cDummyMaster::cDummyMaster() : cEquMaster_Controls() {
  hForm=DummyForm;
}

void cDummyMaster::KillEquation() {
}
void cDummyMaster::PrepareEquation() {
}
//---------------------------------------------------------------------------

