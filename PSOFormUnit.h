//---------------------------------------------------------------------------

#ifndef PSOFormUnitH
#define PSOFormUnitH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TPSOForm : public TForm
{
__published:	// IDE-managed Components
        TLabeledEdit *Edit_nSteps;
        TLabeledEdit *Edit_nPoints;
        TLabeledEdit *Edit_Omega;
        TLabeledEdit *Edit_c1;
        TLabeledEdit *Edit_c2;
        TLabeledEdit *Edit_VCoef;
        TLabeledEdit *Edit_nGroups;
        TLabeledEdit *Edit_SpotSize;
        TLabeledEdit *Edit_LifeTime;
        TLabeledEdit *Edit_nParked;
        TLabeledEdit *Edit_c3m;
        TCheckBox *Check_NoLocalsInSpots;
        TLabeledEdit *Edit_BlackTime;
        TLabeledEdit *Edit_BlackSpot;
        TLabeledEdit *Edit_BlackVal;
        TLabeledEdit *Edit_PSOMode;
private:	// User declarations
public:		// User declarations
        __fastcall TPSOForm(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TPSOForm *PSOForm;
//---------------------------------------------------------------------------

#endif
