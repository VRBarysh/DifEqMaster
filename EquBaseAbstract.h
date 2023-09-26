//---------------------------------------------------------------------------

#ifndef EquBaseAbstractH
#define EquBaseAbstractH

#include <vcl.h>
#include <complex.h>
#include "EquGraph.h"

#define EQU_MAX_GRAPHS   20
#define EQU_MAX_VALS     50
//---------------------------------------------------------------------------

typedef complex<long double> TCplxLong;
typedef complex<double> TCplx;
typedef complex<float> TCplxShort;

class cDifEquation
{
public:
  int Step() { StepRoutine(); t+=dt; return(++nStepNumber); };
  virtual void StepRoutine()=0;
  double GetDoneFrac() { return(t/tMax); }

protected:
long double t,dt,tMax;
int         nStepNumber,nReportNumber;
int         EquType;
};


class cDifEquReportText
{
public:
  cDifEquReportText() {nGraphs=nVals=0;} 
int nGraphs,nVals;
AnsiString GraphText[EQU_MAX_GRAPHS],ValText[EQU_MAX_VALS];
};

class cDifEquReport
{
public:
int nVals,nGraphs;
double      Val[EQU_MAX_VALS];
AnsiString  ValText[EQU_MAX_VALS],GraphText[EQU_MAX_GRAPHS];
cEquIOGraph Graph[EQU_MAX_GRAPHS];

  cDifEquReport() { nVals=nGraphs=0; }
  void LoadText(cDifEquReportText *hReportText);
  void SetMaxTo1() { for(int i=0;i<EQU_MAX_GRAPHS;i++) Graph[i].Max=1; }
  void SetGraphSize(int GraphSize)
    { for(int i=0;i<EQU_MAX_GRAPHS;i++) Graph[i].Init(1,GraphSize); }
  void SetGraphSize(int GraphSize,int *hLines);

protected:
};

class cDifEqu_ThreadIO : public cDifEquation
{
public:
  virtual void InitReport(cDifEquReport *hReport)=0;
  virtual void Report(cDifEquReport *hReport)=0;
  virtual void FinalReport(cDifEquReport *hReport)=0;
  virtual int TimeToSaveBitmaps() { return 0; }
  virtual void SaveBitmaps() {}
  void SetHReportText(cDifEquReportText *phReportText)
    { hReportText=phReportText; }
  virtual void FillReportText() {}
  virtual int Get2DColor(int Number2D, double x, double y) {return 255*x*y;}
  virtual void Prepare2DColor(int Number2D) {}
  
protected:
cDifEquReportText *hReportText;
};

class cDifEqu_InitDataRec {
public:
long double t,dt,tMax;
};

class cDifEqu_DataBankIO : public cDifEqu_ThreadIO
{
public:
  virtual void LoadInitData(cDifEqu_InitDataRec *hInitData)=0;
protected:
};
#endif
