//---------------------------------------------------------------------------

#ifndef PSOSolveMasterFormH
#define PSOSolveMasterFormH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include "EquBaseMaster.h"
//---------------------------------------------------------------------------
class TPSOSolverForm : public TForm
{
__published:	// IDE-managed Components
        TLabeledEdit *Edit_deltaReMin;
        TLabeledEdit *Edit_deltaImMin;
        TLabeledEdit *Edit_deltaImMax;
        TLabeledEdit *Edit_gammaReMin;
        TLabeledEdit *Edit_gammaReMax;
        TLabeledEdit *Edit_gammaImMin;
        TLabeledEdit *Edit_gammaImMax;
        TLabeledEdit *Edit_deltaReMax;
        TLabeledEdit *Edit_lx;
        TLabeledEdit *Edit_lz;
        TLabeledEdit *Edit_alpha;
        TLabeledEdit *Edit_epsilon;
        TLabeledEdit *Edit_nx;
        TLabeledEdit *Edit_nz;
        TLabeledEdit *Edit_mode;
        TLabeledEdit *Edit_LzMin;
        TLabeledEdit *Edit_ScaleSteps;
        TLabeledEdit *Edit_SaveBitmapsInt;
private:	// User declarations
public:		// User declarations
        __fastcall TPSOSolverForm(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TPSOSolverForm *PSOSolverForm;
//---------------------------------------------------------------------------


class cMasterPSOSOlver : public cEquMaster_Controls
{
public:
  cMasterPSOSOlver();
  virtual ~cMasterPSOSOlver() { KillThread(); KillEquation(); }

  void KillEquation();
  void PrepareEquation();
  void AutoSave();
protected:
int EquType;
};

#endif
