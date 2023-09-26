//---------------------------------------------------------------------------

#ifndef EquOpticTask_MultiMediaH
#define EquOpticTask_MultiMediaH
//---------------------------------------------------------------------------
#include "EquOpticTask.h"

class cEquOpticTask_MultiMedia : public cEquOpticBaseTask {
public:

cEquOpticTask_MultiMedia() { fFinished=0; RungeStep=MediaRungeStep=0; fDebug=0;}

inline int IndexXZM(int ix,int iz,int im) { return( ((InitData.nStepsZ+2)*ix+iz)*InitData.nMedia+im ); }

OPTIC_TYPE MediaDistr(int im);
OPTIC_TYPE MediaDelta(int im);
void MakeMediaDistrNormCoef();

virtual void StepMediaRunge2();
//virtual void StepMediaRunge2NoR2();
//virtual void StepMediaRunge2TE();
//virtual void StepMediaRunge2FastP();
//virtual void StepMediaRunge2FastPTE();
//virtual void StepMediaRunge2FixedNonlinear();




OPTIC_TYPE MediaDistNormCoef;

};

#endif
