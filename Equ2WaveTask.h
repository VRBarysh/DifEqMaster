//---------------------------------------------------------------------------

#ifndef Equ2WaveTaskH
#define Equ2WaveTaskH

#include "Dump2WaveElModel.h"
#include "EquBaseTask.h"

//---------------------------------------------------------------------------

class cEqu2WaveBaseTask : public cEquBaseTask {
public:

cEqu2WaveBaseTask() { fFinished=0; RungeStep=0;}

virtual void Step();
inline void Step01();

inline REAL_TYPE RightPart0(REAL_TYPE El1);
inline REAL_TYPE RightPart1(complex<REAL_TYPE> Ap, REAL_TYPE El0);

inline int iE0(int iz,int ie) { return iz*InitData.nEls*2+2*ie; }
inline int iE1(int iz,int ie) { return iz*InitData.nEls*2+2*ie+1; }

int Finished() {return(fFinished);}

cDifEqu_InitDataRec_2WaveElModel InitData;
complex<REAL_TYPE> *hKappaData;               // 0.0, 0.5, 1.0, 1.5, 2.0...
complex<REAL_TYPE> *hAp,*hAm,*hI;             // +0 - base
REAL_TYPE          *hEl;                      // +Size - k1
                                              // +2*Size - k2
complex<REAL_TYPE> *htmpAp;                                              
REAL_TYPE          dtT;

int RungeStep;                  // 0 - make k1, 1 - make k2, 2 - make new
int fFinished,iz1st,izLast,ElSize;
};

#endif
 