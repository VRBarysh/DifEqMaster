//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquOptic2D_X_02.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

/*
cOptic2D_Media_X cOptic2D_Media_X::operator+(cOptic2D_Media_X m) {
  cOptic2D_Media_X tmp=(*this);
  tmp.PXp+=m.PXp; tmp.PXm+=m.PXm; tmp.PZp+=m.PZp; tmp.PZm+=m.PZm;
  tmp.R0+=m.R0; tmp.RX+=m.RX; tmp.RZ+=m.RZ; tmp.Rp+=m.Rp; tmp.Rm+=m.Rm;
  return(tmp);
}

cOptic2D_Media_X cOptic2D_Media_X::operator+=(cOptic2D_Media_X m) {
  PXp+=m.PXp; PXm+=m.PXm; PZp+=m.PZp; PZm+=m.PZm;
  R0+=m.R0; RX+=m.RX; RZ+=m.RZ; Rp+=m.Rp; Rm+=m.Rm;
  return(*this);
}

cOptic2D_Media_X cOptic2D_Media_X::operator*(OPTIC_TYPE v) {
  cOptic2D_Media_X tmp=(*this);
  tmp.PXp*=v; tmp.PXm*=v; tmp.PZp*=v; tmp.PZm*=v;
  tmp.R0*=v; tmp.RX*=v; tmp.RZ*=v; tmp.Rp*=v; tmp.Rm*=v;
  return(tmp);
}

cOptic2D_Media_X cOptic2D_Media_X::operator*=(OPTIC_TYPE v) {
  PXp*=v; PXm*=v; PZp*=v; PZm*=v;
  R0*=v; RX*=v; RZ*=v; Rp*=v; Rm*=v;
  return(*this);
}

void cOptic2D_Media_X::Add(cOptic2D_Media_X &m) {
  PXp+=m.PXp; PXm+=m.PXm; PZp+=m.PZp; PZm+=m.PZm;
  R0+=m.R0; RX+=m.RX; RZ+=m.RZ; Rp+=m.Rp; Rm+=m.Rm;
}

void cOptic2D_Media_X::Mul(OPTIC_TYPE v) {
  PXp*=v; PXm*=v; PZp*=v; PZm*=v;
  R0*=v; RX*=v; RZ*=v; Rp*=v; Rm*=v;
}
*/
/*
void cOptic2D_Media_X::M1M2v(cOptic2D_Media_X &m1,cOptic2D_Media_X &m2,OPTIC_TYPE v) {
  PXp=m1.PXp+m2.PXp*v; PXm=m1.PXm+m2.PXm*v;
  PZp=m1.PZp+m2.PZp*v; PXm=m1.PZm+m2.PZm*v;
  R0=m1.R0+m1.R0*v; RX=m1.RX+m2.RX*v; RZ=m1.RZ+m2.RZ*v;
  Rp=m1.Rp+m2.Rp*v; Rm=m1.Rm+m2.Rm*v;
} */
