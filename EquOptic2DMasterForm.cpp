//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "EquOptic2DMasterForm.h"
#include "EquOptic2D_X_03.h"
#include "EquOptic2D04.h"
#include "EquOptic2D_X_01.h"
#include "EquOptic2D_1DwB.h"
#include "EquOptic2D_1DwB_Mod1.h"
#include "EquOptic_X_04.h"
#include "EquOptic2D_X_MultiMedia.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TEquOptic2DForm *EquOptic2DForm;
//---------------------------------------------------------------------------
__fastcall TEquOptic2DForm::TEquOptic2DForm(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------

cMasterOptic2D::cMasterOptic2D() : cEquMaster_Controls() {
  hForm=EquOptic2DForm;
//  hDumpForm=Dump25DEl;
  hDump=NULL;
  hEquation=NULL;
  ReportText.nGraphs=11;
  ReportText.GraphText[0]="AZp - middle Z";
  ReportText.GraphText[1]="AZp - middle X";
  ReportText.GraphText[2]="AXp - middle Z";
  ReportText.GraphText[3]="AXp - middle X";
  ReportText.GraphText[4]="Spectrum AZp";
  ReportText.GraphText[5]="AZp(t)";
  ReportText.GraphText[6]="AZp 2D";
  ReportText.GraphText[7]="AZm 2D";
  ReportText.GraphText[8]="AXp 2D";
  ReportText.GraphText[9]="AXm 2D";
  ReportText.GraphText[10]="Energy(t)";
  ReportText00=ReportText;
  ReportText01.nGraphs=10;
  ReportText01.GraphText[0]="AZm(Z)";
  ReportText01.GraphText[1]="AZp(Z)";
  ReportText01.GraphText[2]="AZmBragg(Z)";
  ReportText01.GraphText[3]="AZpBragg(Z)";
  ReportText01.GraphText[6]="AZm(X)";
  ReportText01.GraphText[7]="AZp(X)";
  ReportText01.GraphText[8]="AZmBragg(X)";
  ReportText01.GraphText[9]="AZpBragg(X)";
  ReportText01.GraphText[4]="Spectrum AZp";
  ReportText01.GraphText[5]="AZp(t)";
  ReportText02.nGraphs=14;
  ReportText02.GraphText[0]="AZp - middle Z";
  ReportText02.GraphText[1]="AZp - middle X";
  ReportText02.GraphText[2]="AXp - middle Z";
  ReportText02.GraphText[3]="AXp - middle X";
  ReportText02.GraphText[4]="Спектр";
  ReportText02.GraphText[5]="F(t)";
  ReportText02.GraphText[6]="AZp 2D";
  ReportText02.GraphText[7]="AZm 2D";
  ReportText02.GraphText[8]="AXp 2D";
  ReportText02.GraphText[9]="AXm 2D";
  ReportText02.GraphText[10]="PZp 2D";
  ReportText02.GraphText[11]="PZm 2D";
  ReportText02.GraphText[12]="PXp 2D";
  ReportText02.GraphText[13]="PXm 2D";

  ReportText03.nGraphs=6;
  ReportText03.GraphText[0]="AZp";
  ReportText03.GraphText[1]="AZm";
  ReportText03.GraphText[2]="B";
  ReportText03.GraphText[3]="R0";
  ReportText03.GraphText[4]="Спектр";
  ReportText03.GraphText[5]="F(t)";

  ReportText04.nGraphs=10;
  ReportText04.GraphText[0]="AZp";
  ReportText04.GraphText[1]="AZm";
  ReportText04.GraphText[2]="B";
  ReportText04.GraphText[3]="R0";
  ReportText04.GraphText[4]="Спектр";
  ReportText04.GraphText[5]="F(t)";
  ReportText04.GraphText[6]="A(z=L/2,x)";
  ReportText04.GraphText[7]="Амплитуды гармоник A";
  ReportText04.GraphText[8]="Поперечная структура A";
  ReportText04.GraphText[9]="Поперечная структура R(z=Lz/2)";
  nOfEquation=0;
  nOfReport=0;
}

void cMasterOptic2D::KillEquation() {
  if(hEquation) delete hEquation;
}

void cMasterOptic2D::AutoSave() {
}

void cMasterOptic2D::PrepareEquation() {
  cEquOptic2D *tmpEqu;
  cDump2DOpticModel *hEquDump;
  cDifEqu_InitDataRec_Optic2D InitData;
  cDifEqu_InitDataRec_Optic2D_T<OPTIC_TYPE> InitDataT;
  InitData.tMax=EquOptic2DForm->Edit_tMax->Text.ToDouble();
  InitData.dt=EquOptic2DForm->Edit_dt->Text.ToDouble();
  InitData.RCoef=EquOptic2DForm->Edit_TimeR->Text.ToDouble();
  InitData.RCoef2=EquOptic2DForm->Edit_TimeR2->Text.ToDouble();
  InitData.PCoef=EquOptic2DForm->Edit_TimeP->Text.ToDouble();
  InitData.beta=EquOptic2DForm->Edit_beta->Text.ToDouble();
  InitData.nStepsZ=EquOptic2DForm->Edit_NStepsZ->Text.ToInt();
  InitData.nStepsZB=EquOptic2DForm->Edit_nStepsZB->Text.ToInt();
  InitData.nStepsZB2=EquOptic2DForm->Edit_nStepsZB2->Text.ToInt();
  InitData.nStepsDelay=EquOptic2DForm->Edit_nStepsZT->Text.ToInt();
  InitData.nStepsX=EquOptic2DForm->Edit_NStepsX->Text.ToInt();
  InitData.nStepsX2=EquOptic2DForm->Edit_nStepsX2->Text.ToInt();
  InitData.nTSubSteps=EquOptic2DForm->Edit_nTSubSteps->Text.ToInt();
  InitData.FFTSize=EquOptic2DForm->Edit_FFTSize->Text.ToInt();
  InitData.DumpInt=EquOptic2DForm->Edit_DumpInt->Text.ToInt();
  InitData.FFTScale=EquOptic2DForm->Edit_FFTScale->Text.ToInt();
  InitData.DumpReq=EquOptic2DForm->Check_DumpReq->Checked ? 1 : 0;
  InitData.A0=EquOptic2DForm->Edit_A0->Text.ToDouble();
  InitData.A1=EquOptic2DForm->Edit_A1->Text.ToDouble();
  InitData.A2=EquOptic2DForm->Edit_A2->Text.ToDouble();
  InitData.ACoef=EquOptic2DForm->Edit_ACoef->Text.ToDouble();
  InitData.ACoef2=EquOptic2DForm->Edit_ACoef2->Text.ToDouble();
  InitData.ACoef1D=EquOptic2DForm->Edit_ACoef1D->Text.ToDouble();
  InitData.APhase1D=EquOptic2DForm->Edit_APhase1D->Text.ToDouble();
  InitData.AGenCoef1D=EquOptic2DForm->Edit_AGenCoef1D->Text.ToDouble();
  InitData.APhase=EquOptic2DForm->Edit_APhase->Text.ToDouble();
  InitData.AGenCoef=EquOptic2DForm->Edit_AGenCoef->Text.ToDouble();
  InitData.AGenZ=EquOptic2DForm->Edit_AGenZ->Text.ToDouble();
  InitData.AGenX=EquOptic2DForm->Edit_AGenX->Text.ToDouble();
  InitData.ABackCoef=EquOptic2DForm->Edit_ABackCoef->Text.ToDouble();
  InitData.rPolNoise=EquOptic2DForm->Edit_PolNoise->Text.ToDouble();
  InitData.Q=EquOptic2DForm->Edit_Q->Text.ToDouble();
  InitData.QGrid=EquOptic2DForm->Edit_QGrid->Text.ToDouble();
  InitData.StartingInv=EquOptic2DForm->Edit_StartingInv->Text.ToDouble();
  InitData.Z0RefCoef=EquOptic2DForm->Edit_Z0RefCoef->Text.ToDouble();
  InitData.Z1RefCoef=EquOptic2DForm->Edit_Z1RefCoef->Text.ToDouble();
  InitData.Z0RefPhase=EquOptic2DForm->Edit_ZRefPhase->Text.ToDouble();
  InitData.DiffractionCoef=EquOptic2DForm->Edit_DiffractionCoef->Text.ToDouble();
  InitData.MicrowaveModeLength=EquOptic2DForm->Edit_MicrowaveMode->Text.ToDouble();
  InitData.MaxV[0]=EquOptic2DForm->Edit_MaxE->Text.ToDouble();
  InitData.UseSimpleEC=EquOptic2DForm->Check_UseSimpleEC->Checked ? 1 : 0;
  InitData.UseMirrorDelay=EquOptic2DForm->Check_UseMirrorDelay->Checked ? 1 : 0;
  InitData.UseFastP=EquOptic2DForm->Check_UseFastP->Checked ? 1 : 0;
  InitData.LockX=EquOptic2DForm->Check_LockX->Checked ? 1 : 0;
  InitData.Coaxial=EquOptic2DForm->Check_Coaxial->Checked ? 1 : 0;
  InitData.fForceLinear=EquOptic2DForm->Check_ForceLinear->Checked ? 1 : 0;
  InitDataT.SaveBitmapInterval=EquOptic2DForm->Edit_SaveBitmapInterval->Text.ToDouble();
  InitData.BitmapNorm=EquOptic2DForm->Edit_BitmapNorm->Text.ToDouble();
  InitData.Delta=EquOptic2DForm->Edit_Delta->Text.ToDouble();
  InitData.DeltaCoef=EquOptic2DForm->Edit_DeltaCoef->Text.ToDouble();
  InitData.fNoMedia=EquOptic2DForm->Check_NoMedia->Checked ? 1 : 0;
  InitData.fTE=EquOptic2DForm->Check_TE->Checked ? 1 : 0;
  InitData.f2BraggsSolidMedia = EquOptic2DForm->Check_2braggs->Checked ? 1 : 0;
  InitData.f2D1D2DSolidMedia = EquOptic2DForm->Check_3braggs->Checked ? 1 : 0;
  InitData.nThreads=EquOptic2DForm->Edit_nThreads->Text.ToInt();
  InitData.nMediaRows=EquOptic2DForm->Edit_nMediaRows->Text.ToInt();
  InitData.nMediaRowSize=EquOptic2DForm->Edit_nMediaRowSize->Text.ToInt();
  InitData.PurgeNoiseDelay=EquOptic2DForm->Edit_PurgeNoiseDelay->Text.ToInt();
  InitData.nMedia=EquOptic2DForm->Edit_nMedia->Text.ToInt();
  InitData.MediaGaussCoef=EquOptic2DForm->Edit_nMediaGaussCoef->Text.ToDouble();

  InitDataT.tMax=EquOptic2DForm->Edit_tMax->Text.ToDouble();
  InitDataT.dt=EquOptic2DForm->Edit_dt->Text.ToDouble();
  InitDataT.RCoef=EquOptic2DForm->Edit_TimeR->Text.ToDouble();
  InitDataT.RCoef2=EquOptic2DForm->Edit_TimeR2->Text.ToDouble();
  InitDataT.PCoef=EquOptic2DForm->Edit_TimeP->Text.ToDouble();
  InitDataT.beta=EquOptic2DForm->Edit_beta->Text.ToDouble();
  InitDataT.nStepsZ=EquOptic2DForm->Edit_NStepsZ->Text.ToInt();
  InitDataT.nStepsZB=EquOptic2DForm->Edit_nStepsZB->Text.ToInt();
  InitDataT.nStepsZB2=EquOptic2DForm->Edit_nStepsZB2->Text.ToInt();
  InitDataT.nStepsDelay=EquOptic2DForm->Edit_nStepsZT->Text.ToInt();
  InitDataT.nStepsX=EquOptic2DForm->Edit_NStepsX->Text.ToInt();
  InitDataT.nStepsX2=EquOptic2DForm->Edit_nStepsX2->Text.ToInt();
  InitDataT.nTSubSteps=EquOptic2DForm->Edit_nTSubSteps->Text.ToInt();
  InitDataT.FFTSize=EquOptic2DForm->Edit_FFTSize->Text.ToInt();
  InitDataT.DumpInt=EquOptic2DForm->Edit_DumpInt->Text.ToInt();
  InitDataT.FFTScale=EquOptic2DForm->Edit_FFTScale->Text.ToInt();
  InitDataT.DumpReq=EquOptic2DForm->Check_DumpReq->Checked ? 1 : 0;
  InitDataT.A0=EquOptic2DForm->Edit_A0->Text.ToDouble();
  InitDataT.A1=EquOptic2DForm->Edit_A1->Text.ToDouble();
  InitDataT.A2=EquOptic2DForm->Edit_A2->Text.ToDouble();
  InitDataT.ACoef=EquOptic2DForm->Edit_ACoef->Text.ToDouble();
  InitDataT.ACoef2=EquOptic2DForm->Edit_ACoef2->Text.ToDouble();
  InitDataT.ACoef1D=EquOptic2DForm->Edit_ACoef1D->Text.ToDouble();
  InitDataT.APhase1D=EquOptic2DForm->Edit_APhase1D->Text.ToDouble();
  InitDataT.AGenCoef1D=EquOptic2DForm->Edit_AGenCoef1D->Text.ToDouble();
  InitDataT.APhase=EquOptic2DForm->Edit_APhase->Text.ToDouble();
  InitDataT.AGenZ=EquOptic2DForm->Edit_AGenZ->Text.ToDouble();
  InitDataT.AGenX=EquOptic2DForm->Edit_AGenX->Text.ToDouble();
  InitDataT.ABackCoef=EquOptic2DForm->Edit_ABackCoef->Text.ToDouble();
  InitDataT.AGenCoef=EquOptic2DForm->Edit_AGenCoef->Text.ToDouble();
  InitDataT.rPolNoise=EquOptic2DForm->Edit_PolNoise->Text.ToDouble();
  InitDataT.Q=EquOptic2DForm->Edit_Q->Text.ToDouble();
  InitDataT.QGrid=EquOptic2DForm->Edit_QGrid->Text.ToDouble();
  InitDataT.StartingInv=EquOptic2DForm->Edit_StartingInv->Text.ToDouble();
  InitDataT.Z0RefCoef=EquOptic2DForm->Edit_Z0RefCoef->Text.ToDouble();
  InitDataT.Z1RefCoef=EquOptic2DForm->Edit_Z1RefCoef->Text.ToDouble();
  InitDataT.Z0RefPhase=EquOptic2DForm->Edit_ZRefPhase->Text.ToDouble();
  InitDataT.DiffractionCoef=EquOptic2DForm->Edit_DiffractionCoef->Text.ToDouble();
  InitDataT.MicrowaveModeLength=EquOptic2DForm->Edit_MicrowaveMode->Text.ToDouble();
  InitDataT.MaxV[0]=EquOptic2DForm->Edit_MaxE->Text.ToDouble();
  InitDataT.UseSimpleEC=EquOptic2DForm->Check_UseSimpleEC->Checked ? 1 : 0;
  InitDataT.UseMirrorDelay=EquOptic2DForm->Check_UseMirrorDelay->Checked ? 1 : 0;
  InitDataT.UseFastP=EquOptic2DForm->Check_UseFastP->Checked ? 1 : 0;
  InitDataT.LockX=EquOptic2DForm->Check_LockX->Checked ? 1 : 0;
  InitDataT.Coaxial=EquOptic2DForm->Check_Coaxial->Checked ? 1 : 0;
  InitDataT.fForceLinear=EquOptic2DForm->Check_ForceLinear->Checked ? 1 : 0;
  InitDataT.SaveBitmapInterval=EquOptic2DForm->Edit_SaveBitmapInterval->Text.ToDouble();
  InitDataT.BitmapNorm=EquOptic2DForm->Edit_BitmapNorm->Text.ToDouble();
  InitDataT.Delta=EquOptic2DForm->Edit_Delta->Text.ToDouble();
  InitDataT.DeltaCoef=EquOptic2DForm->Edit_DeltaCoef->Text.ToDouble();
  InitDataT.fNoMedia=EquOptic2DForm->Check_NoMedia->Checked ? 1 : 0;
  InitDataT.fTE=EquOptic2DForm->Check_TE->Checked ? 1 : 0;
  InitDataT.f2BraggsSolidMedia = EquOptic2DForm->Check_2braggs->Checked ? 1 : 0;
  InitDataT.f2D1D2DSolidMedia = EquOptic2DForm->Check_3braggs->Checked ? 1 : 0;
  InitDataT.nThreads=EquOptic2DForm->Edit_nThreads->Text.ToInt();
  InitDataT.nMediaRows=EquOptic2DForm->Edit_nMediaRows->Text.ToInt();
  InitDataT.nMediaRowSize=EquOptic2DForm->Edit_nMediaRowSize->Text.ToInt();
  InitDataT.PurgeNoiseDelay=EquOptic2DForm->Edit_PurgeNoiseDelay->Text.ToInt();
  InitDataT.nMedia=EquOptic2DForm->Edit_nMedia->Text.ToInt();
  InitDataT.MediaGaussCoef=EquOptic2DForm->Edit_nMediaGaussCoef->Text.ToDouble();

  if(EquOptic2DForm->Check_2wB->Checked) {
    hEquation = new cEquOptic2D_1DwB;
    hEquation->LoadInitData(&InitDataT);
    hEquation->SetHReportText(&ReportText);
    if(nOfReport!=3) {
      ReportText=ReportText03;
      InitGraphSelector();
      nOfReport=3;
    }
  } else if(EquOptic2DForm->Check_1DwBMod1->Checked) {
    hEquation = new cEquOptic2D_1DwB_Mod1;
    hEquation->LoadInitData(&InitDataT);
    hEquation->SetHReportText(&ReportText);
    if(nOfReport!=4) {
      ReportText=ReportText04;
      InitGraphSelector();
      nOfReport=4;
    }
  } else {
    if(EquOptic2DForm->Check_Equ01->Checked)
      hEquation = new cEquOptic2D_X_Mod2;
    else
    if(InitDataT.nMedia)
      hEquation = new cEquOptic2D_X_MultiMedia;
    else
      hEquation = new cEquOptic2D_X_base;
    hEquation->LoadInitData(&InitDataT);
    hEquation->SetHReportText(&ReportText);
    if(nOfReport!=2) {
      ReportText=ReportText02;
      InitGraphSelector();
      nOfReport=2;
    }
  }
  /*
  if(EquOptic2DForm->Check_Equ01->Checked) {
    hEquation = new cEquOptic2D_01;
    nOfEquation=1;
    if(nOfReport!=1) {
      ReportText=ReportText01; nOfReport=1;
      InitGraphSelector();
    }
  } else {
    hEquation = new cEquOptic2D;
    nOfEquation=0;
    if(nOfReport!=0) {
      ReportText=ReportText00; nOfReport=0;
      InitGraphSelector();
    }
  }
  hEquation->LoadInitData(&InitData);
  hEquation->SetHReportText(&ReportText);
  if (InitData.DumpReq) {
    hEquDump = new cDump2DOpticModel;
    ((cEquOptic2D *)hEquation)->SetHDump(hEquDump);
    hDump=(cEquBaseDump *)hEquDump;
  }
  */
}

//---------------------------------------------------------------------------



