//---------------------------------------------------------------------------

#ifndef EquOptic2DMasterFormH
#define EquOptic2DMasterFormH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include "EquOptic2D.h"
#include <ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TEquOptic2DForm : public TForm
{
__published:	// IDE-managed Components
        TLabeledEdit *Edit_NStepsZ;
        TLabeledEdit *Edit_NStepsX;
        TLabeledEdit *Edit_TimeP;
        TLabeledEdit *Edit_TimeR;
        TLabeledEdit *Edit_PolNoise;
        TLabeledEdit *Edit_tMax;
        TLabeledEdit *Edit_A0;
        TLabeledEdit *Edit_dt;
        TLabeledEdit *Edit_ACoef;
        TLabeledEdit *Edit_FFTSize;
        TLabeledEdit *Edit_DumpInt;
        TCheckBox *Check_DumpReq;
        TLabeledEdit *Edit_nTSubSteps;
        TCheckBox *Check_UseSimpleEC;
        TLabeledEdit *Edit_AGenCoef;
        TCheckBox *Check_UseFastP;
        TLabeledEdit *Edit_beta;
        TCheckBox *Check_LockX;
        TLabeledEdit *Edit_Q;
        TLabeledEdit *Edit_Z0RefCoef;
        TLabeledEdit *Edit_A1;
        TLabeledEdit *Edit_nStepsZB;
        TLabeledEdit *Edit_nStepsZT;
        TCheckBox *Check_Equ01;
        TLabeledEdit *Edit_FFTScale;
        TLabeledEdit *Edit_ZRefPhase;
        TLabeledEdit *Edit_Z1RefCoef;
        TLabeledEdit *Edit_DiffractionCoef;
        TCheckBox *Check_UseMirrorDelay;
        TLabeledEdit *Edit_A2;
        TLabeledEdit *Edit_nThreads;
        TCheckBox *Check_NoMedia;
        TCheckBox *Check_Coaxial;
        TLabeledEdit *Edit_MaxE;
        TLabeledEdit *Edit_AGenZ;
        TLabeledEdit *Edit_AGenX;
        TLabeledEdit *Edit_SaveBitmapInterval;
        TLabeledEdit *Edit_BitmapNorm;
        TLabeledEdit *Edit_StartingInv;
        TLabeledEdit *Edit_ABackCoef;
        TLabeledEdit *Edit_MicrowaveMode;
        TCheckBox *Check_TE;
        TLabeledEdit *Edit_PurgeNoiseDelay;
        TCheckBox *Check_2wB;
        TLabeledEdit *Edit_Delta;
        TCheckBox *Check_1DwBMod1;
        TLabeledEdit *Edit_DeltaCoef;
        TLabeledEdit *Edit_nStepsX2;
        TLabeledEdit *Edit_QGrid;
        TCheckBox *Check_ForceLinear;
        TLabeledEdit *Edit_ACoef2;
        TLabeledEdit *Edit_nMediaRows;
        TLabeledEdit *Edit_nMediaRowSize;
        TLabeledEdit *Edit_APhase;
        TLabeledEdit *Edit_nMedia;
        TLabeledEdit *Edit_nMediaGaussCoef;
        TLabeledEdit *Edit_TimeR2;
        TCheckBox *Check_2braggs;
        TLabeledEdit *Edit_nStepsZB2;
        TCheckBox *Check_3braggs;
        TLabeledEdit *Edit_ACoef1D;
        TLabeledEdit *Edit_APhase1D;
        TLabeledEdit *Edit_AGenCoef1D;
private:	// User declarations
public:		// User declarations
        __fastcall TEquOptic2DForm(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TEquOptic2DForm *EquOptic2DForm;
//---------------------------------------------------------------------------

class cMasterOptic2D : public cEquMaster_Controls
{
public:
  cMasterOptic2D();
  virtual ~cMasterOptic2D() { KillThread(); KillEquation(); delete hDump;}

  void KillEquation();
  void PrepareEquation();
  void AutoSave();
protected:
  cDifEquReportText ReportText00, ReportText01, ReportText02, ReportText03, ReportText04;
  int nOfEquation,nOfReport;
};

#endif
