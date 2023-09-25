//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquOptic2D03.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

cOptic2D_Media_2Waves cOptic2D_Media_2Waves::operator+(const cOptic2D_Media_2Waves &m1) {
  cOptic2D_Media_2Waves Res;
  Res.AZp=AZp+m1.AZp; Res.AZm=AZm+m1.AZm;
  Res.PZp=PZp+m1.PZp; Res.PZm=PZm+m1.PZm;
  Res.R0=R0+m1.R0; Res.RZ=RZ+m1.RZ;
  Res.hInitData=hInitData;
  return(Res);
}

void cOptic2D_Media_2Waves::operator+=(const cOptic2D_Media_2Waves &m1) {
  AZp+=m1.AZp; AZm+=m1.AZm;
  PZp+=m1.PZp; PZm+=m1.PZm;
  R0+=m1.R0; RZ+=m1.RZ;
}

void cOptic2D_Media_2Waves::operator+=(const cOptic2D_Media_2Waves m1) {
  AZp+=m1.AZp; AZm+=m1.AZm;
  PZp+=m1.PZp; PZm+=m1.PZm;
  R0+=m1.R0; RZ+=m1.RZ;
}

cOptic2D_Media_2Waves cOptic2D_Media_2Waves::operator*(long double m) {
  cOptic2D_Media_2Waves Res;
  Res.AZp=AZp*m; Res.AZm=AZm*m;
  Res.PZp=PZp*m; Res.PZm=PZm*m;
  Res.R0=R0*m; Res.RZ=RZ*m;
  Res.hInitData=hInitData;
  return(Res);
}

TCplxLong cOptic2D_Media_2Waves::RightAZp()
    {return( PZp+hInitData->AGenCoef*AZp );}
TCplxLong cOptic2D_Media_2Waves::RightAZm()
    {return( PZm+hInitData->AGenCoef*AZm );}
TCplxLong cOptic2D_Media_2Waves::RightPZp()
    {return((AZp*R0*((long double)2.0)+AZm*RZ)*hInitData->beta - PZp*hInitData->PCoef);}
TCplxLong cOptic2D_Media_2Waves::RightPZm()
    {return((AZm*R0*((long double)2.0)+AZp*conj(RZ))*hInitData->beta  - PZm*hInitData->PCoef);}
TCplxLong cOptic2D_Media_2Waves::RightPZpShort()
  {return((AZp*R0*((long double)2.0)+AZm*RZ)*hInitData->beta);}
TCplxLong cOptic2D_Media_2Waves::RightPZmShort()
  {return((AZm*R0*((long double)2.0)+AZp*conj(RZ))*hInitData->beta);}
long double cOptic2D_Media_2Waves::RightR0() {
  return(hInitData->Q-hInitData->RCoef*R0-(AZp*conj(PZp)+AZm*conj(PZm)).real() );
}
TCplxLong cOptic2D_Media_2Waves::RightRZ() {
  return(-RZ*hInitData->RCoef -(AZp*conj(PZm)+PZp*conj(AZm)) );
}

TCplxLong cOptic2D_Media_2Waves::RightPZpSimple()
    {return((AZp*R0*((long double)2.0))*hInitData->beta - PZp*hInitData->PCoef);}
TCplxLong cOptic2D_Media_2Waves::RightPZmSimple()
    {return((AZm*R0*((long double)2.0))*hInitData->beta  - PZm*hInitData->PCoef);}
TCplxLong cOptic2D_Media_2Waves::RightPZpShortSimple()
  {return((AZp*R0*((long double)2.0))*hInitData->beta);}
TCplxLong cOptic2D_Media_2Waves::RightPZmShortSimple()
  {return((AZm*R0*((long double)2.0))*hInitData->beta);}
long double cOptic2D_Media_2Waves::RightR0Simple() {
  return(hInitData->Q-hInitData->RCoef*R0-(AZp*conj(PZp)+AZm*conj(PZm)).real() );
}
TCplxLong cOptic2D_Media_2Waves::RightRZSimple() { return(0); }

cOptic2D_Media_2Waves cOptic2D_Media_2Waves::RightPartMedia() {
 cOptic2D_Media_2Waves Res;
 Res.AZp=Res.AZm=0;
 Res.PZp=RightPZp(); Res.PZm=RightPZm();
 Res.R0=RightR0(); Res.RZ=RightRZ();
 return(Res);
}

cOptic2D_Media_2Waves cOptic2D_Media_2Waves::RightPart() {
  cOptic2D_Media_2Waves Res;

  if(hInitData->UseFastP) {
    PZp=hInitData->PCoefFast*RightPZpShort();
    PZm=hInitData->PCoefFast*RightPZmShort();
  }
  Res.AZp=RightAZp(); Res.AZm=RightAZm();
  Res.PZp=RightPZp(); Res.PZm=RightPZm();
  Res.R0=RightR0(); Res.RZ=RightRZ();
/*
  if(hInitData->UseFastP) {
    PZp=hInitData->PCoefFast*RightPZpShortSimple();
    PZm=hInitData->PCoefFast*RightPZmShortSimple();
  }
  Res.AZp=RightAZp(); Res.AZm=RightAZm();
  Res.PZp=RightPZpSimple(); Res.PZm=RightPZmSimple();
  Res.R0=RightR0Simple(); Res.RZ=RightRZSimple();*/
  return(Res);
}

void cOptic2D_Media_2Waves::MakeNewSimple(cOptic2D_Media_2Waves *hNew,const int i,const int j){
  cOptic2D_Media_2Waves k1,k2,k3,k4;
  k1=RightPart()*hInitData->dt;
  k2=(*this+k1*0.5).RightPart()*hInitData->dt;
  k3=(*this+k2).RightPart()*hInitData->dt;
  k4=(*this+k3).RightPart()*hInitData->dt;
  *hNew=*this+(k1+k2*2.0+k3*2.0+k4)*(((long double)1.0)/6.0);
}
void cOptic2D_Media_2Waves::MakeNewSeveralSteps(cOptic2D_Media_2Waves *hNew,const int i,const int j){
  cOptic2D_Media_2Waves tmpNew=*this;
  long double tmpdt=hInitData->dt;
  hInitData->dt=tmpdt/hInitData->nTSubSteps;
  for(int i=0;i<hInitData->nTSubSteps;i++){
    tmpNew.MakeNewSimple(hNew,i,j);
    tmpNew=*hNew;
  }
  hInitData->dt=tmpdt;
}
