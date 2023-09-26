//---------------------------------------------------------------------------

#ifndef EquBaseDumpH
#define EquBaseDumpH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ComCtrls.hpp>
#include <ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TDumpForm : public TForm
{
__published:	// IDE-managed Components
        TTrackBar *TrackBar1;
        TLabeledEdit *LabeledEdit1;
private:	// User declarations
public:		// User declarations
        __fastcall TDumpForm(TComponent* Owner);
};

class cEquBaseDump
{
public:
  cEquBaseDump();
  virtual void Show()=0;
  int Enabled() { return(DumpEnabled); }
  int Ready() { return(DumpReady); }
protected:
int DumpEnabled, DumpReady;
TDumpForm *hForm;
};

//---------------------------------------------------------------------------
extern PACKAGE TDumpForm *DumpForm;
//---------------------------------------------------------------------------
#endif
