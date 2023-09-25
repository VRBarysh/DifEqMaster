//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquOptic2D02.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

cOptic2D_Media_BraggOnly
       cOptic2D_Media_BraggOnly::operator+(const cOptic2D_Media_BraggOnly &m1) {
  cOptic2D_Media_BraggOnly Res;
  Res.AZp=AZp+m1.AZp; Res.AZm=AZm+m1.AZm; Res.AXp=AXp+m1.AXp; Res.AXm=AXm+m1.AXm;
  Res.hInitData=hInitData;
  return(Res);
}

void cOptic2D_Media_BraggOnly::operator+=(const cOptic2D_Media_BraggOnly &m1) {
  AZp+=m1.AZp; AZm+=m1.AZm; AXp+=m1.AXp; AXm+=m1.AXm;
}

void cOptic2D_Media_BraggOnly::operator+=(const cOptic2D_Media_BraggOnly m1) {
  AZp+=m1.AZp; AZm+=m1.AZm; AXp+=m1.AXp; AXm+=m1.AXm;
}

cOptic2D_Media_BraggOnly cOptic2D_Media_BraggOnly::operator*(long double m) {
  cOptic2D_Media_BraggOnly Res;
  Res.AZp=AZp*m; Res.AZm=AZm*m; Res.AXp=AXp*m; Res.AXm=AXm*m;
  Res.hInitData=hInitData;
  return(Res);
}

TCplxLong cOptic2D_Media_BraggOnly::RightAZp()
    {return( TCplxLong(0,-1)*hInitData->ACoef*(AXp+AXm)+hInitData->AGenCoef*AZp );}
TCplxLong cOptic2D_Media_BraggOnly::RightAZm()
    {return( TCplxLong(0,-1)*hInitData->ACoef*(AXp+AXm)+hInitData->AGenCoef*AZm );}
TCplxLong cOptic2D_Media_BraggOnly::RightAXp()
    {return( TCplxLong(0,-1)*hInitData->ACoef*(AZp+AZm)+hInitData->AGenCoef*AXp );}
TCplxLong cOptic2D_Media_BraggOnly::RightAXm()
    {return( TCplxLong(0,-1)*hInitData->ACoef*(AZp+AZm)+hInitData->AGenCoef*AXm );}

cOptic2D_Media_BraggOnly cOptic2D_Media_BraggOnly::RightPart() {
  cOptic2D_Media_BraggOnly Res;
  Res.AZp=RightAZp(); Res.AZm=RightAZm();
  Res.AXp=RightAXp(); Res.AXm=RightAXm();
  Res.hInitData=hInitData;
  return(Res);
}

void cOptic2D_Media_BraggOnly::
        MakeNewSimple(cOptic2D_Media_BraggOnly *hNew,const int i,const int j){
  cOptic2D_Media_BraggOnly k1,k2,k3,k4,tmp;
  /*
  tmp=RightPart();
  k1=tmp*hInitData->dt;
  tmp=*this+k1*0.5;
  tmp=tmp.RightPart;
  k2=tmp*hInitData->dt;
  tmp=*this+k2;
  tmp=tmp.RightPart;
  k3=tmp*hInitData->dt;
  tmp=*this+k3;
  tmp=tmp.RightPart;
  k4=tmp*hInitData->dt;
  *hNew=*this+(k1+k2*2.0+k3*2.0+k4)*(((long double)1.0)/6.0);
  */
//  long double eAZp=abs(AZp)*abs(AZp);
//  long double eAZm=abs(AZm)*abs(AZm);
//  long double eAXp=abs(AXp)*abs(AXp);
//  long double eAXm=abs(AXm)*abs(AXm);
//  long double eSum=eAXp+eAXm+eAZp+eAZm;
  k1=RightPart()*hInitData->dt;
//  long double eAZpk1=abs(AZp+k1.AZp)*abs(AZp+k1.AZp);
//  long double eAZmk1=abs(AZm+k1.AZm)*abs(AZm+k1.AZm);
//  long double eAXpk1=abs(AXp+k1.AXp)*abs(AXp+k1.AXp);
//  long double eAXmk1=abs(AXm+k1.AXm)*abs(AXm+k1.AXm);
//  long double eSumk1=eAXpk1+eAXmk1+eAZpk1+eAZmk1;
  k2=(*this+k1*0.5).RightPart()*hInitData->dt;
//  long double eAZpk2=abs(AZp+k2.AZp)*abs(AZp+k2.AZp);
//  long double eSumk2=eAZpk2*4.0;
  k3=(*this+k2).RightPart()*hInitData->dt;
//  long double eAZpk3=abs(AZp+k3.AZp)*abs(AZp+k3.AZp);
//  long double eSumk3=eAZpk3*4.0;
  k4=(*this+k3).RightPart()*hInitData->dt;
//  long double eAZpk4=abs(AZp+k4.AZp)*abs(AZp+k4.AZp);
//  long double eSumk4=eAZpk4*4.0;
  *hNew=*this+(k1+k2*2.0+k3*2.0+k4)*(((long double)1.0)/6.0);
//  long double neAZp=abs(hNew->AZp)*abs(hNew->AZp);
//  long double neAZm=abs(hNew->AZm)*abs(hNew->AZm);
//  long double neAXp=abs(hNew->AXp)*abs(hNew->AXp);
//  long double neAXm=abs(hNew->AXm)*abs(hNew->AXm);
//  long double neSum=neAXp+neAXm+neAZp+neAZm;
}

void cOptic2D_Media_BraggOnly::
        MakeNewSeveralSteps(cOptic2D_Media_BraggOnly *hNew,const int i,const int j){
  cOptic2D_Media_BraggOnly tmpNew=*this;
  long double tmpdt=hInitData->dt;
  hInitData->dt=tmpdt/hInitData->nTSubSteps;
  for(int i=0;i<hInitData->nTSubSteps;i++){
    if(hInitData->UseSimpleEC) tmpNew.MakeNewSimpleEC1(hNew,i,j);
      else tmpNew.MakeNewSimple(hNew,i,j);
    tmpNew=*hNew;
  }
  hInitData->dt=tmpdt;
}

void cOptic2D_Media_BraggOnly::MakeNewSimpleEC1(cOptic2D_Media_BraggOnly *hNew,const int i,const int j){
  cOptic2D_Media_BraggOnly k1,k2,k3,k4,middle;
  long double energyXP=AXp.real()*AXp.real()+AXp.imag()*AXp.imag();
  long double energyXM=AXm.real()*AXm.real()+AXm.imag()*AXm.imag();
  long double energyZP=AZp.real()*AZp.real()+AZp.imag()*AZp.imag();
  long double energyZM=AZm.real()*AZm.real()+AZm.imag()*AZm.imag();
  k1=RightPart()*hInitData->dt;
  k2=(*this+k1*0.5).RightPart()*hInitData->dt;
  k3=(*this+k2).RightPart()*hInitData->dt;
  k4=(*this+k3).RightPart()*hInitData->dt;
  *hNew=*this+(k1+k2*2.0+k3*2.0+k4);

  long double energyXP1=hNew->AXp.real()*hNew->AXp.real()+hNew->AXp.imag()*hNew->AXp.imag();
  long double energyXM1=hNew->AXm.real()*hNew->AXm.real()+hNew->AXm.imag()*hNew->AXm.imag();
  long double energyZP1=hNew->AZp.real()*hNew->AZp.real()+hNew->AZp.imag()*hNew->AZp.imag();
  long double energyZM1=hNew->AZm.real()*hNew->AZm.real()+hNew->AZm.imag()*hNew->AZm.imag();
  long double tmpCoef;
  if(energyXP1+energyXM1+energyZP1+energyZM1>0.00000000001) {
    tmpCoef=sqrt((energyXP+energyXM+energyZP+energyZM)/(energyXP1+energyXM1+energyZP1+energyZM1));
    hNew->AXp*=tmpCoef;
    hNew->AXm*=tmpCoef;
    hNew->AZp*=tmpCoef;
    hNew->AZm*=tmpCoef;
  }
}

