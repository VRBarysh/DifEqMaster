//---------------------------------------------------------------------------

#ifndef EquOpticTask_Wide1DBraggH
#define EquOpticTask_Wide1DBraggH
//---------------------------------------------------------------------------

#include <complex.h>

#include "EquOptic2D01.h"
#include "EquOptic2D_X_02.h"
#include "EquOpticTask.h"
#include "FFTUnit.h"

class cOptic2D_Media_X_1D {
public:
  complex<OPTIC_TYPE> PZp,PZm,RZ; // r2z, r2x, rz+x, rz-x
  OPTIC_TYPE R0,R0Generated;   // r0
};

class cEquOpticTask_Wide1DBragg : public cEquBaseTask {
public:

cEquOpticTask_Wide1DBragg() { fFinished=0; RungeStep=MediaRungeStep=0; fDebug=0;}

inline int IndexXZ(int ix,int iz) { return((InitData.nStepsZ+2)*ix+iz); }
inline int IsBragg(int ix) { return( (((InitData.nStepsX>>1)-ix)*((InitData.nStepsX>>1)-ix)<=(InitData.nStepsXBragg)*(InitData.nStepsXBragg)) ? 1 : 0); }

virtual void Step();
//virtual void StepNoMedia();
//virtual void StepNoMediaMatrix();
//virtual void StepNoMediaMatrixMMode();
//virtual void StepNoMediaEuler();
//virtual void StepNoMediaRunge2();

virtual void StepRunge4FastP();
//virtual void StepMediaRunge2NoR2();
//virtual void StepMediaRunge2TE();
//virtual void StepMediaRunge2FastP();
//virtual void StepMediaRunge2FastPTE();
//virtual void StepMediaRunge2FixedNonlinear();
int Finished() {return(fFinished);}

cDifEqu_InitDataRec_Optic2D_T<OPTIC_TYPE> InitData;
complex<OPTIC_TYPE> *hAp,*hAm;                // +0 - base
cOptic2D_Media_X_2Way    *hM;                 // +Size - k1
                                              // +2*Size - k2
                                              // +3*Size - k3
                                              // +4*Size - k4
complex<OPTIC_TYPE> *htmpApReal,*hApReal,*htmpAmReal,*hAmReal;

CFFTrans<OPTIC_TYPE> FFT;

OPTIC_TYPE          dtT,t;

int RungeStep,MediaRungeStep; // 0 - make k1, 1 - make k2, 2 - make new
int fFinished,fDebug,ix1st,ixLast;
};

#endif
