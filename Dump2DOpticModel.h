//---------------------------------------------------------------------------

#ifndef Dump2DOpticModelH
#define Dump2DOpticModelH

#include <vcl.h>
#include "EquBaseMaster.h"
#include "EquBaseDump.h"
#include "Gr3DPlane.h"
#include "EquOptic2D01.h"
//---------------------------------------------------------------------------
 /*
class DrawDumpInitData {
public:
int nXSteps;
long double t,dx;
bool DoubleSided,DrawLines;
int Skip;
};
*/
class cDump2DOpticModel : public cEquBaseDump
{
public:
  cDump2DOpticModel();
  ~cDump2DOpticModel() { delete [] hMediaDump; }

  void LoadInitData(cDifEqu_InitDataRec_Optic2D *hInit);
  void Dump(cOptic2D_Media *hMedia);

  void Show() {}
  /*
  void SaveCURR(char *hFileName);
  void SaveXTSlice(DrawDumpInitData DrawInit);
  void SaveXZSlice(DrawDumpInitData DrawInit);
  void SaveZTSlice(DrawDumpInitData DrawInit);

  void MakeXProfileCoef(long double pdx, int pnXSteps);
  void DrawZTSlice(Graphics::TBitmap *hBitmap,cGr3DInitData Gr3DInit,DrawDumpInitData DrawInit);
  void DrawXTSlice(Graphics::TBitmap *hBitmap,cGr3DInitData Gr3DInit,DrawDumpInitData DrawInit);
  void DrawXZSlice(Graphics::TBitmap *hBitmap,cGr3DInitData Gr3DInit,DrawDumpInitData DrawInit);
  TCplxLong GetAVal(int ZStep,int XStep,int TStep);
  double GetMaxTime() {return(Step*Init.dt);}
  double GetMaxZ() {return(Init.nSteps*Init.dt);}

   */
cDifEqu_InitDataRec_Optic2D Init;
int DumpEls,nTSteps,Step,DumpSizeX,DumpSizeZ,DumpSizeT;
cOptic2D_Media_Dump *hMediaDump;
int *hElPresent;

int DumpChanged;
protected:
};


/*
class cDumpStorage
{
public:
  cDumpStorage(){ hADump=hIDump=NULL; hElDump=hMaxADump=NULL;
                  nEls=nSteps=nStepsT=nOptiSteps1=nOptiSteps2=nOptiSteps3=0;
                  Dumping=DUMPING_NONE; DumpReady=DumpRequired=DumpElectrons=false; }
  ~cDumpStorage(){ delete [] hADump; delete [] hIDump;
                   delete [] hElDump; delete [] hMaxADump; }

  void Init_Diffur(int pnSteps, int pnEls, REAL_TYPE pL, REAL_TYPE ptMax);
  void Init_Opti(int Type,
            int pnOptiSteps1,int pnOptiSteps2=1,int pnOptiSteps3=1);

  void DumpOpti(REAL_TYPE MaxA);
  void DumpDiffur( complex<REAL_TYPE> *hA,
                   complex<REAL_TYPE> *hI, REAL_TYPE *hEls);
  bool Ready() { return(DumpReady); }

  void Save_AUDump(char *hFileName,float time);
  void Save_AIDump(char *hFileName,float time);
  void Save_AIUDump(char *hFileName,float time);

  void SaveSurfaceZT(char *hFileName);
  void SaveCURR(char *hFileName);

  void Draw_AUDump(TImage *hImage,float time);
  void Draw_AIUDump(TImage *hImage,float time);

  complex<float> GetAValue(REAL_TYPE z,REAL_TYPE t);
  complex<float> GetIValue(REAL_TYPE z,REAL_TYPE t);

  float GetEndTime() {return(EndTime);}
  float GetL() { return(L); }
  double GetMaxA() { return(MaxA); }
  double GetMaxI() { return(MaxI); }

  void SetRequiredState(bool pDumpRequired,bool pDumpElectrons)
    { DumpRequired=pDumpRequired; DumpElectrons=pDumpElectrons; }

protected:
complex<float>        *hADump,*hIDump;
float                 *hElDump;
float                 *hMaxADump;
REAL_TYPE             EndTime,L;
int                   nSteps,nStepsT,nEls;
int                   nOptiSteps1,nOptiSteps2,nOptiSteps3;
int                   Step,StepT,OptiStep1,OptiStep2,OptiStep3;
int                   Dumping;
bool                  DumpReady,DumpRequired,DumpElectrons;
double                MaxI,MaxA,MinU,MaxU;
};
*/

#endif
