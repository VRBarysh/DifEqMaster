//---------------------------------------------------------------------------

#ifndef OpticParaBasicEquH
#define OpticParaBasicEquH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>

#include "EquBaseMaster.h"
#include <ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TFormOpticBasic : public TForm
{
__published:	// IDE-managed Components
        TLabeledEdit *Edit_nSteps;
        TLabeledEdit *Edit_tMax;
        TLabeledEdit *Edit_dt;
        TLabeledEdit *Edit_dz;
        TLabeledEdit *Edit_alpha;
        TLabeledEdit *Edit_alpha1;
        TLabeledEdit *Edit_beta;
        TLabeledEdit *Edit_PNoise;
        TCheckBox *CheckBox1;
private:	// User declarations
public:		// User declarations
        __fastcall TFormOpticBasic(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TFormOpticBasic *FormOpticBasic;
//---------------------------------------------------------------------------

class cDifEqu_InitDataRec_OptiBasic : public cDifEqu_InitDataRec
{
public:
int nSteps;
long double rPolNoise,dz,alpha,alpha1,beta;
};

class cEquOpticBasic : public cDifEqu_DataBankIO
{
public:
  cEquOpticBasic();
  ~cEquOpticBasic();
  void StepRoutine();

  void InitReport(cDifEquReport *hReport);
  void Report(cDifEquReport *hReport);
  void FinalReport(cDifEquReport *hReport);

  virtual void LoadInitData(cDifEqu_InitDataRec *hInitData);

  virtual void FillAData( TCplxLong *hP );

  TCplxLong RightPartP(TCplxLong P,long double R, TCplxLong A);
  long double RightPartR(TCplxLong P,long double R, TCplxLong A);
protected:
TCplxLong *A,*P,*Pk1,*Pk2,*Pnew;
long double *R,*Rk1,*Rk2,*Rnew;
cDifEqu_InitDataRec_OptiBasic InitData;
};

class cEquOpticRealTimeAxis : public cEquOpticBasic
{
public:
  cEquOpticRealTimeAxis();
  ~cEquOpticRealTimeAxis();
  void StepRoutine();

  TCplxLong PStor(int line, int i)
    { return(PStorage[InitData.nSteps*((iPStorBegin+line)%InitData.nSteps)+i]); }
  int StoredLines() {return( (iPStorEnd+InitData.nSteps-iPStorBegin) % InitData.nSteps);}

  void FillAData( TCplxLong *hP );
  
  virtual void LoadInitData(cDifEqu_InitDataRec *hInitData);
protected:
TCplxLong *PStorage;
int iPStorBegin,iPStorEnd;
};

class cMasterOpticBasic : public cEquMaster_Controls
{
public:
  cMasterOpticBasic();
  virtual ~cMasterOpticBasic() { KillThread(); KillEquation(); }

  void KillEquation();
  void PrepareEquation();
protected:
};


#endif
