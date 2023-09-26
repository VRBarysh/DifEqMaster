//---------------------------------------------------------------------------

#ifndef EquOpticTask_BMB_holesH
#define EquOpticTask_BMB_holesH

#include "EquOpticTask_BMB.h"
//---------------------------------------------------------------------------

//---------------------------------------------------  bragg - media - bragg with holes
class cEquOpticBaseTask_BMB_holes : public cEquOpticBaseTask_BMB {
public:

cEquOpticBaseTask_BMB_holes() { fFinished=0; RungeStep=MediaRungeStep=0; fDebug=0; hDeltaArray=NULL; }

int isBragg(int ix);
OPTIC_TYPE Delta(int ix);

virtual void StepNoMediaMatrix();
virtual void StepMediaRunge2FastP();
virtual void StepMediaRunge2();

complex<OPTIC_TYPE> *hAXp1,*hAXm1,*hAZp1,*hAZm1;
complex<OPTIC_TYPE> *hAXp2,*hAXm2,*hAZp2,*hAZm2;

OPTIC_TYPE *hDeltaArray;

};
#endif
