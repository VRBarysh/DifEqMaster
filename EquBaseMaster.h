//---------------------------------------------------------------------------

#ifndef EquBaseMasterH
#define EquBaseMasterH

#include <vcl.h>
#include "EquBaseThread.h"
#include "EquBaseDump.h"
#include "PSOFormUnit.h"

//---------------------------------------------------------------------------
class cEquBaseMaster
{
public:
  virtual void Activate()=0;
  virtual void Deactivate()=0;
  virtual void Start()=0;
  virtual void Pause()=0;
  virtual void Resume()=0;
protected:
};

class cEquMaster_Controls : public cEquBaseMaster
{
public:
  cEquMaster_Controls() {hForm=NULL; hThread=NULL; hEquation=NULL; hDump=NULL; hDumpForm=NULL; UsesPSO=0;}
  ~cEquMaster_Controls();
  void SetVCLHandlers(cEquVCLHandlers pHandlers) { Handlers=pHandlers; }
  virtual void Activate();
  void ShowForm() { if(hForm) hForm->Show(); if (UsesPSO) PSOForm->Show();}
  void ShowDumpForm() { if(hDumpForm) hDumpForm->Show(); }
  void InitGraphSelector();
  void Deactivate();
  void Pause() {if(hThread) hThread->Suspend();}
  void Resume() {if(hThread) hThread->Resume();}
  void Start();
  void KillThread();

  virtual void AutoSave() {}
  virtual void KillEquation()=0;
  virtual void PrepareEquation()=0;
protected:
cEquVCLHandlers Handlers;



TForm           *hForm,*hDumpForm;
cEquBaseDump    *hDump;
cEquGraphThread *hThread;
cDifEqu_DataBankIO *hEquation;
cDifEquReportText ReportText;
int UsesPSO;
};
#endif
