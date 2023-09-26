//---------------------------------------------------------------------------

#ifndef EquOptic2D03H
#define EquOptic2D03H

#include "EquOptic2D01.h"
//---------------------------------------------------------------------------

class cOptic2D_Media_2Waves
{
public:
  cOptic2D_Media_2Waves() {}
  ~cOptic2D_Media_2Waves() {}

  cOptic2D_Media_2Waves operator+(const cOptic2D_Media_2Waves &m1);
  void operator+=(const cOptic2D_Media_2Waves &m1);
  void operator+=(const cOptic2D_Media_2Waves m1);
  cOptic2D_Media_2Waves operator*(long double m);

  void MakeNewSeveralSteps(cOptic2D_Media_2Waves *hNew,const int i,const int j);
  void MakeNewSimple(cOptic2D_Media_2Waves *hNew,const int i,const int j);
  void MakeNewSimpleEC(cOptic2D_Media_2Waves *hNew,const int i,const int j);
  void MakeNewSimpleEC1(cOptic2D_Media_2Waves *hNew,const int i,const int j);

cOptic2D_Media_2Waves RightPart();
cOptic2D_Media_2Waves RightPartMedia();
inline TCplxLong RightAZp();
inline TCplxLong RightAZm();
inline TCplxLong RightPZp();
inline TCplxLong RightPZm();
inline TCplxLong RightPZpShort();
inline TCplxLong RightPZmShort();
inline long double RightR0();
inline TCplxLong RightRZ();

inline TCplxLong RightPZpSimple();
inline TCplxLong RightPZmSimple();
inline TCplxLong RightPZpShortSimple();
inline TCplxLong RightPZmShortSimple();
inline long double RightR0Simple();
inline TCplxLong RightRZSimple();

  TCplxLong AZp,AZm,PZp,PZm,RZ; // r2z, r2x, rz+x, rz-x
  long double R0;   // r0
  cDifEqu_InitDataRec_Optic2D *hInitData;
};

#endif
