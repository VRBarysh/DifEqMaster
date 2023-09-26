//---------------------------------------------------------------------------

#ifndef EquOptic2D_X_02H
#define EquOptic2D_X_02H
//---------------------------------------------------------------------------
#include <complex.h>
#include "EquOptic2D01.h"

class cOptic2D_Media_X {
public:
  /*
  inline cOptic2D_Media_X operator+(cOptic2D_Media_X m);
  inline cOptic2D_Media_X operator+=(cOptic2D_Media_X m);
  inline cOptic2D_Media_X operator*(OPTIC_TYPE v);
  inline cOptic2D_Media_X operator*=(OPTIC_TYPE v);
  inline void Add(cOptic2D_Media_X &m);
  inline void Mul(OPTIC_TYPE v);
  */
  inline void M1M2v(cOptic2D_Media_X &m1,cOptic2D_Media_X &m2,OPTIC_TYPE v) {
    PXp=m1.PXp+m2.PXp*v; PXm=m1.PXm+m2.PXm*v;
    PZp=m1.PZp+m2.PZp*v; PZm=m1.PZm+m2.PZm*v;
    R0=m1.R0+m2.R0*v; RX=m1.RX+m2.RX*v; RZ=m1.RZ+m2.RZ*v;
    Rp=m1.Rp+m2.Rp*v; Rm=m1.Rm+m2.Rm*v; R0Generated=m1.R0Generated+m2.R0Generated*v;
  }

  complex<OPTIC_TYPE> PZp,PZm,PXp,PXm,RZ,RX,Rp,Rm; // r2z, r2x, rz+x, rz-x
  OPTIC_TYPE R0,R0Generated;   // r0
};

class cOptic2D_Media_X_2Way {
public:
  /*
  inline cOptic2D_Media_X operator+(cOptic2D_Media_X m);
  inline cOptic2D_Media_X operator+=(cOptic2D_Media_X m);
  inline cOptic2D_Media_X operator*(OPTIC_TYPE v);
  inline cOptic2D_Media_X operator*=(OPTIC_TYPE v);
  inline void Add(cOptic2D_Media_X &m);
  inline void Mul(OPTIC_TYPE v);
  */
  inline void M1M2v(cOptic2D_Media_X_2Way &m1,cOptic2D_Media_X_2Way &m2,OPTIC_TYPE v) {
    PZp=m1.PZp+m2.PZp*v; PZm=m1.PZm+m2.PZm*v;
    R0=m1.R0+m2.R0*v; RZ=m1.RZ+m2.RZ*v;
    R0Generated=m1.R0Generated+m2.R0Generated*v;
  }

  complex<OPTIC_TYPE> PZp,PZm,RZ; // r2z, r2x, rz+x, rz-x
  OPTIC_TYPE R0,R0Generated;   // r0
};

#endif
