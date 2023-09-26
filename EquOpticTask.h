//---------------------------------------------------------------------------

#ifndef EquOpticTaskH
#define EquOpticTaskH
//---------------------------------------------------------------------------

#include <complex.h>

#include "EquOptic2D01.h"
#include "EquOptic2D_X_02.h"
#include "EquBaseTask.h"

class cEquOpticBaseTask : public cEquBaseTask {
public:

cEquOpticBaseTask() { fFinished=0; RungeStep=MediaRungeStep=0; fDebug=0;}

virtual void FillMatrixes();

inline int IndexXZ(int ix,int iz) { return((InitData.nStepsZ+2)*ix+iz); }

virtual void Step();
virtual void StepNoMedia();
virtual void StepNoMediaMatrix();
virtual void StepNoMediaMatrixMMode();
virtual void StepNoMediaEuler();
virtual void StepNoMediaRunge2();

virtual void StepMediaRunge2();
virtual void StepMediaRunge2NoR2();
virtual void StepMediaRunge2TE();
virtual void StepMediaRunge2FastP();
virtual void StepMediaRunge2FastPTE();
virtual void StepMediaRunge2FixedNonlinear();
int Finished() {return(fFinished);}

cDifEqu_InitDataRec_Optic2D_T<OPTIC_TYPE> InitData;
complex<OPTIC_TYPE> *hAXp,*hAXm,*hAZp,*hAZm;  // +0 - base
cOptic2D_Media_X    *hM;                      // +Size - k1
                                              // +2*Size - k2
OPTIC_TYPE          dtT,t;
complex<OPTIC_TYPE> Matrix4[16], Matrix3[9], Matrix2[4];

int RungeStep,MediaRungeStep; // 0 - make k1, 1 - make k2, 2 - make new
int fFinished,fDebug,ix1st,ixLast;
};

#endif
