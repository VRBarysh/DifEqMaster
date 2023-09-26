//---------------------------------------------------------------------------

#ifndef EquOpticTask_BMBH
#define EquOpticTask_BMBH
//--------------------------------------------------------------------------

#include "EquOpticTask_2Way.h"

//---------------------------------------------------  bragg - media - bragg
class cEquOpticBaseTask_BMB : public cEquOpticBaseTask_2Way {
public:

cEquOpticBaseTask_BMB() { fFinished=0; RungeStep=MediaRungeStep=0; fDebug=0;}

virtual void Step();

complex<OPTIC_TYPE> *hAXp1,*hAXm1,*hAZp1,*hAZm1;
complex<OPTIC_TYPE> *hAXp2,*hAXm2,*hAZp2,*hAZm2;

};

#endif
