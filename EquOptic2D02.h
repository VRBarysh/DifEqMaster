//---------------------------------------------------------------------------

#ifndef EquOptic2D02H
#define EquOptic2D02H

#include "EquOptic2D01.h"
//---------------------------------------------------------------------------

class cOptic2D_Media_BraggOnly
{
public:
  cOptic2D_Media_BraggOnly() {}
  ~cOptic2D_Media_BraggOnly() {}

//  cOptic2D_Media operator+(cOptic2D_Media &m1);
  cOptic2D_Media_BraggOnly operator+(const cOptic2D_Media_BraggOnly &m1);
  void operator+=(const cOptic2D_Media_BraggOnly &m1);
  void operator+=(const cOptic2D_Media_BraggOnly m1);
  cOptic2D_Media_BraggOnly operator*(long double m);

  void MakeNewSeveralSteps(cOptic2D_Media_BraggOnly *hNew,const int i,const int j);
  void MakeNewSimple(cOptic2D_Media_BraggOnly *hNew,const int i,const int j);
  void MakeNewSimpleEC1(cOptic2D_Media_BraggOnly *hNew,const int i,const int j);

cOptic2D_Media_BraggOnly RightPart();
cOptic2D_Media_BraggOnly RightPartMedia();
inline TCplxLong RightAZp();
inline TCplxLong RightAZm();
inline TCplxLong RightAXp();
inline TCplxLong RightAXm();

  TCplxLong AZp,AZm,AXp,AXm;
  cDifEqu_InitDataRec_Optic2D *hInitData;
};

#endif
