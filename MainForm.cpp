//---------------------------------------------------------------------------

#include <vcl.h>
#include <vcl\syncobjs.hpp>
#include <math.h>
#include <dir.h>
#pragma hdrstop

#include "MainForm.h"
#include "EquGraph.h"
#include "EquBaseThread.h"
#include "Equation01.h"
#include "Equation02.h"
#include "OpticParaBasicEqu.h"
#include "Equ25DEl01.h"
#include "EquOptic2DMasterForm.h"
#include "PSOSolveMasterForm.h"
#include "CephesLib.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"

using namespace std;

TForm1 *Form1;
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{
  TimerValue=0;
  /*
  long double dx=3,dt=0.1;
  long double sqPI=sqrt(M_PI);
  long double sqdt=sqrt(dt),sq2=sqrt((long double)2);
  long double tmpRe=0,tmpIm=0;
  long double x,x3,a,cs,sn,d1,d2,sqd1,sqd2,cs1,cs2,sn1,sn2;
  for(int j=0;j<100;j++) {
    x=dx*j; x3=x*x*x; a=-0.25*x*x/dt; sn=sin(a); cs=cos(a);
    double tmpFC,tmpFS,tmpFC2,tmpFS2;
    fresnl(sq2*x/(2*sqPI*sqdt),&tmpFS,&tmpFC);
    tmpRe=-0.5*x*sq2*sqPI+2*sqdt*cs+x*sq2*sqPI*tmpFS;
    tmpIm=-0.5*x*sq2*sqPI+2*sqdt*sn+x*sq2*sqPI*tmpFC;
//    hXProfileA[j]=TCplxLong(tmpRe,tmpIm);
    tmpRe=-0.5*x3*sq2*sqPI+4*sqdt*sqdt*sqdt*cs+2*x*x*sqdt*sn+x3*sq2*sqPI*tmpFC;
    tmpIm=0.5*x3*sq2*sqPI+4*sqdt*sqdt*sqdt*sn-2*x*x*sqdt*cs-x3*sq2*sqPI*tmpFS;
    tmpIm/=6;  tmpRe/=6;
    for(int i=1;i<100;i++) {
      d1=(i+1)*dt; sqd1=sqrt(d1); d2=i*dt; sqd2=sqrt(d2);
      a=-0.25*x*x/d1; cs1=cos(a); sn1=sin(a);
      a=0.25*x*x/d2; cs2=cos(a); sn2=sin(a);
      fresnl(sq2*x/(2*sqPI*sqd1),&tmpFS,&tmpFC);
      fresnl(sq2*x/(2*sqPI*sqd2),&tmpFS2,&tmpFC2);
      tmpRe=2*sqd1*cs1+x*sq2*sqPI*tmpFS-2*sqd2*cs2-x*sq2*sqPI*tmpFS2;
      tmpIm=2*sqd1*sn1+x*sq2*sqPI*tmpFC+2*sqd2*sn2-x*sq2*sqPI*tmpFC2;
      tmpRe=4*sqd1*sqd1*sqd1*cs1+2*x*x*sqd1*sn1+x3*sq2*sqPI*tmpFC
           -4*sqd2*sqd2*sqd2*cs2+2*x*x*sqd2*sn2-x3*sq2*sqPI*tmpFC2;
      tmpIm=4*sqd2*sqd2*sqd2*sn2+2*x*x*sqd2*cs2+x3*sq2*sqPI*tmpFS2
           +4*sqd1*sqd1*sqd1*sn1-2*x*x*sqd1*cs1-x3*sq2*sqPI*tmpFS;
      tmpIm/=6;  tmpRe/=6;
      x=x;
    }
  } */
}
//---------------------------------------------------------------------------

void __fastcall TForm1::GraphSelectorClickCheck(TObject *Sender)
{
  int nGraphChecks=0;
  for(int i=0; i<GraphSelector->Count;i++)
    if(GraphSelector->Checked[i]) nGraphChecks++;
  if(nGraphChecks>=EQU_MAX_GRAPH_CHECKS)
    for(int i=0; i<GraphSelector->Count;i++)
      if(!GraphSelector->Checked[i]) GraphSelector->ItemEnabled[i]=false;
  if(nGraphChecks<EQU_MAX_GRAPH_CHECKS)
    for(int i=0; i<GraphSelector->Count;i++) GraphSelector->ItemEnabled[i]=true;
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Timer1Timer(TObject *Sender)
{
  TimerValue++;
}
//---------------------------------------------------------------------------

void TForm1::ThreadExitMessageHandler(TMessage &Message)
{
  //Application->MessageBox("Done!","Job is Done!",MB_OK);
  if(hSelectedMaster) hSelectedMaster->AutoSave();     
}
void __fastcall TForm1::FormClose(TObject *Sender, TCloseAction &Action)
{
  for(int i=0;i<nOfMasters;i++) delete hEquationMaster[i];
  delete hBitmap;
}               
//---------------------------------------------------------------------------

void __fastcall TForm1::EquationSelectorClick(TObject *Sender)
{
  if(hSelectedMaster!=hEquationMaster[EquationSelector->ItemIndex]) {
    if(hSelectedMaster) hSelectedMaster->Deactivate();
    hSelectedMaster=hEquationMaster[EquationSelector->ItemIndex];
    hSelectedMaster->Activate();
  } else { hSelectedMaster->ShowForm(); }
}
//---------------------------------------------------------------------------

void __fastcall TForm1::FormShow(TObject *Sender)
{
  char Dir[1000];
  cEquVCLHandlers Handlers;
  Handlers.hImage=Image1;
  Handlers.hGraphSelector=GraphSelector;
  Handlers.hProgressBar=ProgressBar1;
  Handlers.hMainFormWindow=Handle;
  Handlers.hTimerValue=&TimerValue;
  Handlers.hSloMo=Check_SloMo;
  Handlers.hBitmap=new Graphics::TBitmap;
  Handlers.hBitmap->Width=Handlers.hImage->ClientWidth;
  Handlers.hBitmap->Height=Handlers.hImage->ClientHeight;
  hBitmap=Handlers.hBitmap;
  Handlers.hMainCanvas=Canvas;
  int disk = getdisk();
  Dir[0]=disk+'A';  Dir[1]=':'; Dir[2]='\\';
  getcurdir(0,Dir+3);
  Form1->Caption=Dir;
  nOfMasters=6;
  hEquationMaster[0] = new cDummyMaster;
  hEquationMaster[0]->SetVCLHandlers(Handlers);
  hEquationMaster[1] = new cDummyMaster2;
  hEquationMaster[1]->SetVCLHandlers(Handlers);
  hEquationMaster[2] = new cMasterOpticBasic;
  hEquationMaster[2]->SetVCLHandlers(Handlers);
  hEquationMaster[3] = new cMaster25DBasic;
  hEquationMaster[3]->SetVCLHandlers(Handlers);
  hEquationMaster[4] = new cMasterOptic2D;
  hEquationMaster[4]->SetVCLHandlers(Handlers);
  hEquationMaster[5] = new cMasterPSOSOlver;
  hEquationMaster[5]->SetVCLHandlers(Handlers);
  hEquationMaster[6] = new cMaster2WaveModel;
  hEquationMaster[6]->SetVCLHandlers(Handlers);
  hWorkingMaster=hSelectedMaster=NULL;
  EquationSelector->Items->Add("Dummy Equation");
  EquationSelector->Items->Add("Dummy Equation2");
  EquationSelector->Items->Add("Optic basic");
  EquationSelector->Items->Add("2.5D Electron model");
  EquationSelector->Items->Add("2D Optics");
  EquationSelector->Items->Add("PSO Solver");
  EquationSelector->Items->Add("Задача Екатерины Рудольфовны");
  Randomize();
  randomize();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::StartButtonClick(TObject *Sender)
{
  if(hSelectedMaster) {
    hSelectedMaster->Start();
  }
}
//---------------------------------------------------------------------------


void __fastcall TForm1::FormPaint(TObject *Sender)
{
  Canvas->Draw(0,0,hBitmap);
}
//---------------------------------------------------------------------------

void __fastcall TForm1::DumpButtonClick(TObject *Sender)
{
  if(hSelectedMaster) hSelectedMaster->ShowDumpForm();        
}
//---------------------------------------------------------------------------


