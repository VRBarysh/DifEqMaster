//---------------------------------------------------------------------------

#include <vcl.h>
#include <CheckLst.hpp>
#pragma hdrstop

#include "EquBaseThread.h"
#pragma package(smart_init)
//---------------------------------------------------------------------------
//   Important: Methods and properties of objects in VCL can only be
//   used in a method called using Synchronize, for example:
//
//      Synchronize(UpdateCaption);
//
//   where UpdateCaption could look like:
//
//      void __fastcall cEquBaseThread::UpdateCaption()
//      {
//        Form1->Caption = "Updated in a thread";
//      }
//---------------------------------------------------------------------------

__fastcall cEquBaseThread::cEquBaseThread(bool CreateSuspended)
        : TThread(CreateSuspended)
{
}
//---------------------------------------------------------------------------
void __fastcall cEquBaseThread::Execute()
{
        //---- Place thread code here ----
}
//---------------------------------------------------------------------------


__fastcall cEquGraphThread::cEquGraphThread(bool CreateSuspended)
        : cEquBaseThread(CreateSuspended)
{
}
//---------------------------------------------------------------------------

void __fastcall cEquGraphThread::BeforeStepsCall() {
  LastStep=0;
  Draw(hReport);
}

void __fastcall cEquGraphThread::StepCall() {
  Draw(hReport);
  Handlers.hProgressBar->Position=hEquation->GetDoneFrac()*1000;
  SloMo=Handlers.hSloMo->Checked;
}

void __fastcall cEquGraphThread::AfterLastStepCall() {
  LastStep=1;
  Draw(hReport);
}

void __fastcall cEquGraphThread::SaveBitmaps() {
  hEquation->SaveBitmaps();
}

void __fastcall cEquGraphThread::Execute()
{
  cDifEquReport Report;
  hReport=&Report;
  hEquation->InitReport(&Report);
  Randomize(); randomize();
  Synchronize(BeforeStepsCall);
  while((hEquation->GetDoneFrac()<1.0)&&(!Terminated)) {
    if(!SloMo) hEquation->Step();
    if(hEquation->TimeToSaveBitmaps()) Synchronize(SaveBitmaps);
    if(*(Handlers.hTimerValue)!=OldTimer) {
      OldTimer=*(Handlers.hTimerValue);
      if(SloMo) hEquation->Step();          
      hEquation->Report(&Report);
      Synchronize(StepCall);
    }
  }
  hEquation->FinalReport(&Report);
  Synchronize(AfterLastStepCall);
  PostMessage(Handlers.hMainFormWindow,WM_USER+1,0,0);
}

void cEquGraphThread::DrawGraph(cEquIOGraph *hGraph,int n,int nMax) {
  int BaseLine=Handlers.hBitmap->Height/nMax*(0.5+n);
  double Scale= (hGraph->Max<0.0000001) ? 1 : (0.5*Handlers.hBitmap->Height/nMax)/(hGraph->Max);
  Handlers.hBitmap->Canvas->Pen->Color=clBlack;
  Handlers.hBitmap->Canvas->MoveTo(Handlers.hBitmap->Width,BaseLine);
  Handlers.hBitmap->Canvas->LineTo(0,BaseLine);
  for(int j=0;j<hGraph->nLines;j++) {
    Handlers.hBitmap->Canvas->MoveTo(0,BaseLine);
    for(int i=0;i<hGraph->Size;i++)
      Handlers.hBitmap->Canvas->LineTo(i,BaseLine-Scale*hGraph->Value(j,i));
  }
}

void cEquGraphThread::DrawGraph2D(cEquIOGraph *hGraph,int n,int nMax) {
  int StartLine=Handlers.hBitmap->Height/nMax*n+5;
  int EndLine=Handlers.hBitmap->Height/nMax*(n+1)-5;
  union {
    int c;
    struct {
      char r,g,b,m;
    } rgb;
  };
  int color;
  hEquation->Prepare2DColor(hGraph->Number2D);
  double invi=1.0/double(hGraph->Size-1);
  double invj=1.0/double(EndLine-StartLine);
  for(int j=0;j<=EndLine-StartLine;j++) {
    for(int i=0;i<hGraph->Size;i++) {
      rgb.m=0;
      color=hEquation->Get2DColor(hGraph->Number2D,double(i)*invi,double(j)*invj);
      int number_of_colors=255;
      double t = 6.2831853071795864 * color / (1.484375 * number_of_colors);
      double rv = 1 + sin(t - 1.9634954084936207);
      double gv = 1 + sin(t - 0.39269908169872415);
      double bv = 1 + sin(t + 1.1780972450961724);
      t = 127.5 * sqrt(sqrt((color + (0.078125 * number_of_colors)) / (0.91796875 * number_of_colors)));
      rgb.r = int(t * rv);
      rgb.g = int(t * gv);
      rgb.b = int(t * bv);
      Handlers.hBitmap->Canvas->Pixels[i][StartLine+j]=c;
    }
  }
}

void cEquGraphThread::Draw(cDifEquReport *hReport) {
  int nChecked=0;
  int Checked[EQU_MAX_GRAPH_CHECKS];
  Handlers.hBitmap->Canvas->FillRect(TRect(0,0,Handlers.hImage->ClientWidth,Handlers.hImage->ClientHeight));
  if(!LastStep) {
    for(int i=0;i<EQU_MAX_GRAPH_CHECKS;i++)  Checked[i]=0;
    for(int i=0;i<hReport->nGraphs;i++)
      if(Handlers.hGraphSelector->Checked[i]) Checked[nChecked++]=i;
  } else {
    nChecked=hReport->nGraphs;
    for(int i=0;i<hReport->nGraphs;i++) Checked[i]=i;
  }
  if(nChecked) for(int i=0;i<nChecked;i++)
    if((hReport->Graph+Checked[i])->Number2D==-1)
      DrawGraph(hReport->Graph+Checked[i],i,nChecked);
    else
      DrawGraph2D(hReport->Graph+Checked[i],i,nChecked);
  for(int i=0;i<hReport->nVals;i++)
    Handlers.hBitmap->Canvas->TextOutA(10,20*i+20,hReport->ValText[i]+AnsiString(hReport->Val[i]));
//  Handlers.hImage->Canvas->CopyRect(TRect(0,0,Handlers.hImage->ClientWidth,Handlers.hImage->ClientHeight),
//           Handlers.hBitmap->Canvas,TRect(0,0,Handlers.hImage->ClientWidth,Handlers.hImage->ClientHeight));
  Handlers.hMainCanvas->Draw(0,0,Handlers.hBitmap);
}

