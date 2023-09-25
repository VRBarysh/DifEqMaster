//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "EquBaseDump.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TDumpForm *DumpForm;
//---------------------------------------------------------------------------
__fastcall TDumpForm::TDumpForm(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------

cEquBaseDump::cEquBaseDump() {
  DumpEnabled=DumpReady=0; hForm=DumpForm;
}
