//---------------------------------------------------------------------------

#ifndef EquOpticTask_2WayH
#define EquOpticTask_2WayH
//---------------------------------------------------------------------------

#include "EquOptic2D_X_01.h"

//---------------------------------------------------  mirror - media - bragg
class cEquOpticBaseTask_2Way : public cEquOpticBaseTask {
public:

cEquOpticBaseTask_2Way() { fFinished=0; RungeStep=MediaRungeStep=0; fDebug=0; }

inline int IndexXZ(int ix,int iz) { return((InitData.nStepsZ+2)*ix+iz); }
inline int IndexXZ_M(int ix,int iz) { return((InitData.nStepsZB+2)*ix+iz); }

virtual void Step();
// virtual void StepNoMedia();
virtual void StepNoMediaMatrix();
virtual void StepNoMediaEuler();
virtual void StepNoMediaRunge2();

virtual void StepMediaRunge2();
virtual void StepMediaRunge2FastP();
int Finished() {return(fFinished);}

OPTIC_TYPE BetaF(int ix);

complex<OPTIC_TYPE> *hAZpMedia,*hAZmMedia;    // +0 - base
cOptic2D_Media_X_2Way    *hM2;                      // +Size - k1
                                                    // +2*Size - k2
};

#endif
