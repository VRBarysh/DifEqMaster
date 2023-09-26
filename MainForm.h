//---------------------------------------------------------------------------

#ifndef MainFormH
#define MainFormH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <CheckLst.hpp>
#include <ComCtrls.hpp>

#include "EquBaseMaster.h"
#include "EquBaseDump.h"

#define MAX_EQUATIONS          50

//---------------------------------------------------------------------------
class TForm1 : public TForm
{
__published:	// IDE-managed Components
        TImage *Image1;
        TButton *StartButton;
        TCheckListBox *GraphSelector;
        TProgressBar *ProgressBar1;
        TTimer *Timer1;
        TListBox *EquationSelector;
        TCheckBox *Check_SloMo;
        TButton *DumpButton;
        void __fastcall GraphSelectorClickCheck(TObject *Sender);
        void __fastcall Timer1Timer(TObject *Sender);                                          
        void __fastcall FormClose(TObject *Sender, TCloseAction &Action);
        void __fastcall EquationSelectorClick(TObject *Sender);
        void __fastcall FormShow(TObject *Sender);
        void __fastcall StartButtonClick(TObject *Sender);
        void __fastcall FormPaint(TObject *Sender);
        void __fastcall DumpButtonClick(TObject *Sender);
private:	// User declarations
        cEquMaster_Controls *hEquationMaster[MAX_EQUATIONS];
        cEquMaster_Controls *hWorkingMaster;
        cEquMaster_Controls *hSelectedMaster;
        int                 nOfMasters;
public:		// User declarations
        int TryTimer() { return(TimerValue); }
        void ClearTimer() { TimerValue=0; }
        int  TimerValue;
        Graphics::TBitmap *hBitmap;
        __fastcall TForm1(TComponent* Owner);
        void ThreadExitMessageHandler(TMessage &Message);
protected:
        BEGIN_MESSAGE_MAP
                VCL_MESSAGE_HANDLER(WM_USER+1, TMessage, ThreadExitMessageHandler)
        END_MESSAGE_MAP(TForm)
};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------
#endif
