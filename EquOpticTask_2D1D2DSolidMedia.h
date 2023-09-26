//---------------------------------------------------------------------------

#ifndef EquOpticTask_2D1D2DSolidMediaH
#define EquOpticTask_2D1D2DSolidMediaH
//---------------------------------------------------------------------------

#include "EquOpticTask_2BraggsSolidMedia.h"

class cEquOpticTask_2D1D2DSolidMedia : public cEquOpticTask_2BraggsSolidMedia {
public:

cEquOpticTask_2D1D2DSolidMedia() { fFinished=0; RungeStep=MediaRungeStep=0; fDebug=0;}

virtual void StepNoMediaMatrix();

};

#endif
