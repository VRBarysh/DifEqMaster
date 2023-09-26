//---------------------------------------------------------------------------

#ifndef Equation01H
#define Equation01H
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Math.h>

#include "EquBaseMaster.h"

//---------------------------------------------------------------------------
class TDummyForm : public TForm
{
__published:	// IDE-managed Components
private:	// User declarations
public:		// User declarations
        __fastcall TDummyForm(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TDummyForm *DummyForm;



//---------------------------------------------------------------------------

class cDummyEquation : public cDifEqu_DataBankIO
{
public:
  ~cDummyEquation() {}
  void StepRoutine() {}
//    { double q; for(int i=0;i<2000000000;i++) q=sqrt(i); }

  void InitReport(cDifEquReport *hReport);
  void Report(cDifEquReport *hReport);
  void FinalReport(cDifEquReport *hReport);

  void LoadInitData(cDifEqu_InitDataRec *hInitData)
    { tMax=hInitData->tMax; dt=hInitData->dt; t=0; }
protected:
};

class cDummyMaster : public cEquMaster_Controls
{
public:
  cDummyMaster();
  virtual ~cDummyMaster() { KillThread(); KillEquation(); }

  void KillEquation();
  void PrepareEquation();
protected:
};

//---------------------------------------------------------------------------
#endif
