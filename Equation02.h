//---------------------------------------------------------------------------

#ifndef Equation02H
#define Equation02H
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include "Equation01.h"
//---------------------------------------------------------------------------
class TDummyForm2 : public TForm
{
__published:	// IDE-managed Components
private:	// User declarations
public:		// User declarations
        __fastcall TDummyForm2(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TDummyForm2 *DummyForm2;
/*
class cDummyEquation : public cDifEqu_DataBankIO
{
public:
  ~cDummyEquation() {}
  void StepRoutine() {}

  void InitReport(cDifEquReport *hReport);
  void Report(cDifEquReport *hReport);
  void FinalReport(cDifEquReport *hReport);

  void LoadInitData(cDifEqu_InitDataRec *hInitData) {}
protected:
};
*/
class cDummyMaster2 : public cDummyMaster
{
public:
  cDummyMaster2();
  virtual ~cDummyMaster2() { KillThread(); KillEquation(); }

  void KillEquation();
  void PrepareEquation();
protected:
};

//---------------------------------------------------------------------------
#endif
