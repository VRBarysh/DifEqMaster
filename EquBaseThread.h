//---------------------------------------------------------------------------

#ifndef EquBaseThreadH
#define EquBaseThreadH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <CheckLst.hpp>
#include <ComCtrls.hpp>
#include "EquBaseAbstract.h"

#define EQU_MAX_GRAPH_CHECKS  3
//---------------------------------------------------------------------------
class cEquBaseThread : public TThread
{
private:
protected:
        void __fastcall Execute();
        virtual void __fastcall BeforeStepsCall()=0;
        virtual void __fastcall StepCall()=0;
        virtual void __fastcall AfterLastStepCall()=0;

        virtual void Draw(cDifEquReport *hReport)=0;

cDifEqu_ThreadIO        *hEquation;
public:
        __fastcall cEquBaseThread(bool CreateSuspended);
        void SetHEquation(cDifEqu_ThreadIO *phEquation)
          { hEquation=phEquation; }
};

class cEquVCLHandlers
{
public:
  cEquVCLHandlers() {hImage=NULL; hGraphSelector=NULL;
                     hProgressBar=NULL; hMainFormWindow=NULL;
                     hTimerValue=NULL; }

TImage                  *hImage;
Graphics::TBitmap       *hBitmap; 
TCheckListBox           *hGraphSelector;
TProgressBar            *hProgressBar;
TCheckBox               *hSloMo;
HWND                    hMainFormWindow;
TCanvas                 *hMainCanvas;
int *hTimerValue;
};

class cEquGraphThread : public cEquBaseThread
{
private:
protected:
        void __fastcall Execute();
        void __fastcall BeforeStepsCall();
        void __fastcall StepCall();
        void __fastcall AfterLastStepCall();
        void __fastcall SaveBitmaps();

        void Draw(cDifEquReport *hReport);
        virtual void DrawGraph(cEquIOGraph *hGraph,int n,int nMax);
        virtual void DrawGraph2D(cEquIOGraph *hGraph,int n,int nMax);
        int OldTimer;
        bool SloMo;
cEquVCLHandlers         Handlers;
cDifEquReport           *hReport;
int                     LastStep;
public:
        __fastcall cEquGraphThread(bool CreateSuspended);

        void SetVCLHandlers(cEquVCLHandlers NewHandlers)
          { Handlers=NewHandlers; }
};
//---------------------------------------------------------------------------
#endif
