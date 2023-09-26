//---------------------------------------------------------------------------

#ifndef Equ25DEl01H
#define Equ25DEl01H
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>

#include "EquBaseMaster.h"
#include "Dump25DElModel.h"
#include "Dump25DElForm.h"
#include "Equ25DEl02.h"
#include "Equ25DEl03.h"
#include <ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TFormEqu25DEl01 : public TForm
{
__published:	// IDE-managed Components
        TLabeledEdit *Edit_nSteps;
        TLabeledEdit *Edit_nEls;
        TLabeledEdit *Edit_tMax;
        TLabeledEdit *Edit_L;
        TLabeledEdit *Edit_r;
        TLabeledEdit *Edit_rTime;
        TLabeledEdit *Edit_ElTime;
        TLabeledEdit *Edit_onSteps;
        TLabeledEdit *Edit_oV1;
        TLabeledEdit *Edit_oV2;
        TCheckBox *Check_optim;
        TCheckBox *Check_ElsOnStart;
        TLabeledEdit *Edit_ICoef;
        TLabeledEdit *Edit_deltaI;
        TLabeledEdit *Edit_deltaS;
        TCheckBox *Check_2Phase;
        TLabeledEdit *Edit_AI;
        TCheckBox *Check_Dump;
        TCheckBox *Check_Simple2Phase;
        TLabeledEdit *Edit_rndCoef;
        TLabeledEdit *Edit_Ai_kappa;
        TCheckBox *Check_UseEA;
        TLabeledEdit *Edit_EAStartTime;
        TLabeledEdit *Edit_EAnSteps;
        TLabeledEdit *Edit_EADelta;
        TLabeledEdit *Edit_EAEndTime;
        TLabeledEdit *Edit_EA_z0;
        TLabeledEdit *Edit_EAkz;
        TLabeledEdit *Edit_nSpeedFracs;
        TLabeledEdit *Edit_SpeedDistrFactor;
        TCheckBox *Check_VertA;
        TCheckBox *Check_MoveAData;
        TLabeledEdit *Edit_KappaMode;
        TLabeledEdit *Edit_KappaV1;
        TLabeledEdit *Edit_KappaV2;
        TLabeledEdit *Edit_KappaV3;
        TCheckBox *Check_maxA;
        TCheckBox *Check_PSO;
        TCheckBox *Check_LOVModel;
        TLabeledEdit *Edit_KappaV4;
        TLabeledEdit *Edit_KappaV5;
        TCheckBox *Check_LoadKappaData;
        TLabeledEdit *Edit_KappaV6;
        TLabeledEdit *Edit_KappaV7;
        TLabeledEdit *Edit_KappaV8;
        TLabeledEdit *Edit_KappaV9;
        TLabeledEdit *Edit_KappaV10;
        TLabeledEdit *Edit_KappaV11;
        TCheckBox *Check_1stPeak;
        TLabeledEdit *Edit_Nu;
        TLabeledEdit *Edit_ElMode;
        TLabeledEdit *Edit_DeltaV1;
        TLabeledEdit *Edit_DeltaV2;
        TLabeledEdit *Edit_AInputMode;
        TLabeledEdit *Edit_AInputV1;
        TLabeledEdit *Edit_AInputV3;
        TLabeledEdit *Edit_AInputV2;
        TLabeledEdit *Edit_DeltaZMode;
        TLabeledEdit *Edit_nDeltaZSteps;
        TLabeledEdit *Edit_DeltaZMax;
        TLabeledEdit *LabeledEdit4;
        TLabeledEdit *LabeledEdit5;
        TLabeledEdit *LabeledEdit6;
        TCheckBox *Check_LoadDeltaZData;
        TLabeledEdit *Edit_TSliceSave;
        TLabeledEdit *Edit_TSliceSaveZ;
        TLabeledEdit *Edit_nThreads;
        TLabeledEdit *Edit_ACoef;
        TLabeledEdit *Edit_RandSeed;
private:	// User declarations
public:		// User declarations
        __fastcall TFormEqu25DEl01(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TFormEqu25DEl01 *FormEqu25DEl01;
//---------------------------------------------------------------------------

class cOptimizingVariation {
public:
long double V1,V2;
int nSteps;                                      
};

class cEqu25DEl_optimizer : public cDifEqu_DataBankIO
{
public:
  cEqu25DEl_optimizer();
  ~cEqu25DEl_optimizer();

  void InitReport(cDifEquReport *hReport);
  void Report(cDifEquReport *hReport);
  void FinalReport(cDifEquReport *hReport);
  void StepRoutine();

  void LoadInitData(cDifEqu_InitDataRec *hInitData);
  void LoadOptimizationData( cOptimizingVariation pVar) { Var=pVar; }
  
  void SetHEquation(cEqu25DElModel *phEquation)
    { hEquation=phEquation; }
protected:
cDifEqu_InitDataRec_25DElModel    InitData;
cOptimizingVariation              Var;
long double                       *hVariable,BestOutValue,BestTotalValue;
long double                       BestOutTime,BestTotalTime;
cEqu25DElModel                    *hEquation;
long double                       MaxOut[100000],MaxTotal[100000],MaxOutMax,MaxTotalMax;
int                               tStep,vStep,nTSteps;
};

class cMaster25DBasic : public cEquMaster_Controls
{
public:
  cMaster25DBasic();
  virtual ~cMaster25DBasic() { KillThread(); KillEquation(); delete hDump; }

  void KillEquation();
  void PrepareEquation();
  void AutoSave();
protected:
int EquType;
cDifEquReportText ReportText00,ReportText01,ReportText02,ReportText03;
};

class cMaster2WaveModel : public cEquMaster_Controls
{
public:
  cMaster2WaveModel();
  virtual ~cMaster2WaveModel() { KillThread(); KillEquation(); }

  void KillEquation();
  void PrepareEquation();
  void AutoSave();
protected:
int EquType;
cDifEquReportText ReportText00,ReportText01,ReportText02,ReportText03;
};

#endif
