//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquBaseMaster.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)
cEquMaster_Controls::~cEquMaster_Controls() {
  KillThread();
  delete hEquation;
  delete hDump;
}

void cEquMaster_Controls::KillThread() {
  if(hThread) {
    hThread->Terminate();
    hThread->Resume();
    hThread->WaitFor();
    delete hThread;
    hThread=NULL;
  }
}

void cEquMaster_Controls::InitGraphSelector() {
  Handlers.hGraphSelector->Items->Clear();
  for(int i=0;i<ReportText.nGraphs;i++)
    Handlers.hGraphSelector->Items->Add(ReportText.GraphText[i]);
}

void cEquMaster_Controls::Activate() {
  InitGraphSelector();
  if(hForm) hForm->Show();
  if (UsesPSO) PSOForm->Show();
}

void cEquMaster_Controls::Deactivate() {
  KillThread();
  KillEquation();
  Handlers.hGraphSelector->Items->Clear();
  if(hForm) hForm->Hide();
}

void cEquMaster_Controls::Start() {
  KillThread();
  KillEquation();
  PrepareEquation();
  hThread = new cEquGraphThread(true);
  hThread->Priority=tpIdle;
  hThread->SetVCLHandlers(Handlers);
  hThread->SetHEquation(hEquation);
  hThread->Resume();
}
