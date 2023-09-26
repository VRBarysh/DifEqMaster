//---------------------------------------------------------------------------

#ifndef EquPSO2DBraggSolverH
#define EquPSO2DBraggSolverH
//---------------------------------------------------------------------------

#include "EquBaseMaster.h"
#include "EquBaseDump.h"
#include "Gr3DPlane.h"
#include "MultiMaxPSO.h"
#include "MTempPSO.h"
//---------------------------------------------------------------------------

class cDifEqu_InitDataRec_PSOSolver : public cDifEqu_InitDataRec {
public:
  long double Lx,Lz,alpha,epsilon,LzMin,Lz0,Lx0;
  long double deltaReMin,deltaReMax,gammaReMin,gammaReMax,
              deltaImMin,deltaImMax,gammaImMin,gammaImMax;
  int mode,tnx,tnz,nScaleSteps,SaveBitmapsInt;
  CPSOInitData PSOData;
};

class cEquPSO2DBraggSolver : public cDifEqu_DataBankIO
{
public:
  cEquPSO2DBraggSolver();
  ~cEquPSO2DBraggSolver();
/*
                        mode =
                                0 - Normal 4D
                                1 - 4D in epsilon around the analytic result
                                2 - 2D in epsilon around the analytic result
                                3 - Save 2D Plane
                                4 - Repeat 2D with scaling down
*/

  void InitReport(cDifEquReport *hReport);
  void Report(cDifEquReport *hReport);
  void FinalReport(cDifEquReport *hReport);

  void SaveDeltaImGammaImPlane();
  void SaveDeltaReIm();
  void SaveGammaReIm();
  void SaveResults();

  void SaveBitmaps();
  void PrepareBitmapFolders();

  void StepRoutine();

  void StepRoutineRepeat();

  virtual long double MainFunction();
  inline TCplxLong LambdaX();
  inline TCplxLong LambdaZ();

  TCplxLong GammaAnal(int mode,int nx,int nz);
  TCplxLong DeltaAnal(int mode,int nx,int nz);

  void LoadInitData(cDifEqu_InitDataRec *hInitData);

protected:
cDifEqu_InitDataRec_PSOSolver     InitData;
CPSOMultiMax                      PSOOptimizer;
CPSOInitComp                      *hPointData;
long double                       *hVariable,BestOutValue,BestTotalValue;
long double                       BestOutTime,BestTotalTime;
long double                       MaxOut,*hResult,*hAnalResult;
long double                       curDeltaRe,curDeltaIm,curGammaRe,curGammaIm;
TCplxLong                         delta,gamma;
int                               tStep,vStep,nTSteps,fError,iScaleStep,iPSOStep;
int                               nBitmap,SaveBitmapInterval,LastBitmapSaveTime;
};


#endif
