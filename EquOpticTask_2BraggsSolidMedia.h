//---------------------------------------------------------------------------

#ifndef EquOpticTask_2BraggsSolidMediaH
#define EquOpticTask_2BraggsSolidMediaH
//---------------------------------------------------------------------------

#include "EquOpticTask.h"

class cEquOpticTask_2BraggsSolidMedia : public cEquOpticBaseTask {
public:

cEquOpticTask_2BraggsSolidMedia() { fFinished=0; RungeStep=MediaRungeStep=0; fDebug=0;}

virtual void FillMatrixes();
//void FillMatrix(complex<OPTIC_TYPE> Alpha, complex<OPTIC_TYPE> *hM4, complex<OPTIC_TYPE> *hM3, complex<OPTIC_TYPE> *hM2);

//inline int IndexXZ(int ix,int iz) { return((InitData.nStepsZ+2)*ix+iz); }

//virtual void Step();
//virtual void StepNoMedia();
virtual void StepNoMediaMatrix();
//virtual void StepNoMediaMatrixMMode();
//virtual void StepNoMediaEuler();
//virtual void StepNoMediaRunge2();

//virtual void StepMediaRunge2();
//virtual void StepMediaRunge2NoR2();
//virtual void StepMediaRunge2TE();
//virtual void StepMediaRunge2FastP();
//virtual void StepMediaRunge2FastPTE();
//virtual void StepMediaRunge2FixedNonlinear();
//int Finished() {return(fFinished);}


complex<OPTIC_TYPE> Matrix4V2[16], Matrix3V2[9], Matrix2V2[4];

};

#endif
