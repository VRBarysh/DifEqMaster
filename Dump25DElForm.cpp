//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "Dump25DElForm.h"
#include "Gr3DPlane.h"

#define SLICE_ZT       1
#define SLICE_XZ       2
#define SLICE_XT       3
#define SLICE_Z        4
#define SLICE_NONE     0

//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TDump25DEl *Dump25DEl;
//---------------------------------------------------------------------------
__fastcall TDump25DEl::TDump25DEl(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall TDump25DEl::FormPaint(TObject *Sender)
{
  Canvas->Draw(0,0,hBitmap);        
}
//---------------------------------------------------------------------------
void __fastcall TDump25DEl::FormCreate(TObject *Sender)
{
  hBitmap=new Graphics::TBitmap;
  hBitmap->Width=640;
  hBitmap->Height=480;
  SliceType=SLICE_NONE;
}
//---------------------------------------------------------------------------
void TDump25DEl::FillInitData() {
  Gr3DInit.PosX=Edit_pos1->Text.ToDouble();
  Gr3DInit.PosY=Edit_pos2->Text.ToDouble();
  Gr3DInit.ScaleX=Edit_scale1->Text.ToDouble();
  Gr3DInit.ScaleY=Edit_scale2->Text.ToDouble();
  Gr3DInit.ScaleZ=Edit_scale3->Text.ToDouble();
  Gr3DInit.Angle1=Edit_angle1->Text.ToDouble();
  Gr3DInit.Angle2=Edit_angle2->Text.ToDouble();
  DrawDumpInit.t=Edit_t->Text.ToDouble();
  DrawDumpInit.dx=Edit_dx->Text.ToDouble();
  DrawDumpInit.nXSteps=Edit_nXSteps->Text.ToInt();
  DrawDumpInit.DoubleSided=Check_DoubleX->Checked;
  DrawDumpInit.DrawLines=Check_DrawLines->Checked;
  DrawDumpInit.Skip=Edit_3DSkip->Text.ToInt();
}
//---------------------------------------------------------------------------
void TDump25DEl::DrawZTSlice() {
  if(!hDump) return;
  FillInitData();
  hDump->DrawZTSlice(hBitmap,Gr3DInit,DrawDumpInit);
  Canvas->Draw(0,0,hBitmap);
}
//---------------------------------------------------------------------------
void TDump25DEl::DrawXTSlice() {
  if(!hDump) return;
  FillInitData();
  hDump->DrawXTSlice(hBitmap,Gr3DInit,DrawDumpInit);
  Canvas->Draw(0,0,hBitmap);
}
//---------------------------------------------------------------------------
void TDump25DEl::DrawXZSlice() {
  if(!hDump) return;
  FillInitData();
  hDump->DrawXZSlice(hBitmap,Gr3DInit,DrawDumpInit);
  Canvas->Draw(0,0,hBitmap);
}
//---------------------------------------------------------------------------
void TDump25DEl::DrawZSlice() {
  if(!hDump) return;
  FillInitData();
  hDump->DrawZSlice(hBitmap,Gr3DInit,DrawDumpInit);
  Canvas->Draw(0,0,hBitmap);
}
//---------------------------------------------------------------------------
void __fastcall TDump25DEl::Button_ZTSliceClick(TObject *Sender)
{
  if(!hDump) return;
  SliceType=SLICE_ZT;
  DrawZTSlice();
}
//---------------------------------------------------------------------------

void __fastcall TDump25DEl::FormDestroy(TObject *Sender)
{
  delete hBitmap;
}
//---------------------------------------------------------------------------

void __fastcall TDump25DEl::Button_XTSliceClick(TObject *Sender)
{
  if(!hDump) return;
  SliceType=SLICE_XT;
  DrawXTSlice();
}
//---------------------------------------------------------------------------

void __fastcall TDump25DEl::Button_XZSliceClick(TObject *Sender)
{
  if(!hDump) return;
  SliceType=SLICE_XZ;
  DrawXZSlice();
}
//---------------------------------------------------------------------------

void __fastcall TDump25DEl::TrackBar1Change(TObject *Sender)
{
  if(!hDump) return;
  double frac=(double)TrackBar1->Position/(double)TrackBar1->Max;;
  switch(SliceType) {
    case SLICE_XZ:
      Edit_t->Text=AnsiString((double)frac*hDump->GetMaxTime());
      DrawXZSlice();
    break;
    case SLICE_XT:
      Edit_t->Text=AnsiString((double)frac*hDump->GetMaxZ()*frac);
      DrawXTSlice();
    break;
    case SLICE_Z:
      Edit_t->Text=AnsiString((double)frac*hDump->GetMaxTime()*frac);
      DrawZSlice();
    break;
  }
}
//---------------------------------------------------------------------------

void __fastcall TDump25DEl::Button_SaveClick(TObject *Sender)
{
  if(!hDump) return;
  FillInitData();
  switch(SliceType) {
    case SLICE_XZ:
      hDump->SaveXZSlice(DrawDumpInit);
    break;
    case SLICE_XT:
      hDump->SaveXTSlice(DrawDumpInit);
    break;
    case SLICE_ZT:
      hDump->SaveZTSlice(DrawDumpInit);
    break;
  }
}
//---------------------------------------------------------------------------

void __fastcall TDump25DEl::Button_ZSliceClick(TObject *Sender)
{
  if(!hDump) return;
  SliceType=SLICE_Z;
  DrawZSlice();
}
//---------------------------------------------------------------------------

