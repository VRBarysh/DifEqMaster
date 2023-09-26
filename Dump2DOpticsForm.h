//---------------------------------------------------------------------------

#ifndef Dump25DElFormH
#define Dump25DElFormH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ComCtrls.hpp>
#include <ExtCtrls.hpp>
#include "Gr3DPlane.h"
#include "Dump25DElModel.h"

//---------------------------------------------------------------------------
class TDump2DOpticsForm : public TForm
{
__published:	// IDE-managed Components
        TTrackBar *TrackBar1;
        TLabeledEdit *Edit_t;
        TLabeledEdit *Edit_nXSteps;
        TLabeledEdit *Edit_dx;
        TButton *Button_XZSlice;
        TButton *Button_XTSlice;
        TLabeledEdit *Edit_pos1;
        TLabeledEdit *Edit_scale1;
        TLabeledEdit *Edit_pos2;
        TLabeledEdit *Edit_scale2;
        TLabeledEdit *Edit_scale3;
        TLabeledEdit *Edit_angle1;
        TLabeledEdit *Edit_angle2;
        TButton *Button_ZTSlice;
        TCheckBox *Check_DoubleX;
        TCheckBox *Check_DrawLines;
        TLabeledEdit *Edit_3DSkip;
        TButton *Button_Save;
        void __fastcall FormPaint(TObject *Sender);
        void __fastcall FormCreate(TObject *Sender);
        void __fastcall Button_ZTSliceClick(TObject *Sender);
        void __fastcall FormDestroy(TObject *Sender);
        void __fastcall Button_XTSliceClick(TObject *Sender);
        void __fastcall Button_XZSliceClick(TObject *Sender);
        void __fastcall TrackBar1Change(TObject *Sender);
        void __fastcall Button_SaveClick(TObject *Sender);
private:	// User declarations
        Graphics::TBitmap *hBitmap;
        cDump25DEl *hDump;
        cGr3DInitData Gr3DInit;
        DrawDumpInitData DrawDumpInit;
        int SliceType;

        void FillInitData();
        void DrawZTSlice();
        void DrawXTSlice();
        void DrawXZSlice();
public:		// User declarations
        void SetHDump(cDump25DEl *phDump) {hDump=phDump;}
        __fastcall TDump2DOpticsForm(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TDump2DOpticsForm *Dump2DOpticsForm;
//---------------------------------------------------------------------------
#endif
