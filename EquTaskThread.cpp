//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "EquTaskThread.h"
#pragma package(smart_init)

//---------------------------------------------------------------------------

//   Important: Methods and properties of objects in VCL can only be
//   used in a method called using Synchronize, for example:
//
//      Synchronize(UpdateCaption);
//
//   where UpdateCaption could look like:
//
//      void __fastcall TTaskThread::UpdateCaption()
//      {
//        Form1->Caption = "Updated in a thread";
//      }
//---------------------------------------------------------------------------

__fastcall TTaskThread::TTaskThread(bool CreateSuspended)
        : TThread(CreateSuspended)
{
}
//---------------------------------------------------------------------------
void __fastcall TTaskThread::Execute()
{
  while (!hEquation->Finished()) {
    hEquation->Step();
    hCritSection->Acquire();
//    if( !(--(*hCounter)) ) hEvent->SetEvent();
    (*hCounter)--;
    if(!(*hCounter)) {
      hEvent->SetEvent();
    }
    hCritSection->Release();
    Suspend();
  }
}
//---------------------------------------------------------------------------
void TTaskThread::Init(cEquBaseTask *hEq, TEvent *hEv,
                       TCriticalSection *hCr, int *hCount) {
  hEquation=hEq; hEvent=hEv; hCritSection=hCr; hCounter=hCount; 
}
