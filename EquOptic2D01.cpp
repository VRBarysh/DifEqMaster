//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquOptic2D01.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)
/*
cOptic2D_Media cOptic2D_Media::operator+(cOptic2D_Media &m1) {
  cOptic2D_Media Res;
  Res.AZp=AZp+m1.AZp; Res.AZm=AZm+m1.AZm; Res.AXp=AXp+m1.AXp; Res.AXm=AXm+m1.AXm;
  Res.PZp=PZp+m1.PZp; Res.PZm=PZm+m1.PZm; Res.PXp=PXp+m1.PXp; Res.PXm=PXm+m1.PXm;
  Res.R0=R0+m1.R0; Res.RZ=RZ+m1.RZ; Res.RX=RX+m1.RX; Res.Rp=Rp+m1.Rp; Res.Rm=Rm+m1.Rm;
  Res.hInitData=hInitData;
  return(Res);
} */

cOptic2D_Media cOptic2D_Media::operator+(const cOptic2D_Media &m1) {
  cOptic2D_Media Res;
  Res.AZp=AZp+m1.AZp; Res.AZm=AZm+m1.AZm; Res.AXp=AXp+m1.AXp; Res.AXm=AXm+m1.AXm;
  Res.PZp=PZp+m1.PZp; Res.PZm=PZm+m1.PZm; Res.PXp=PXp+m1.PXp; Res.PXm=PXm+m1.PXm;
  Res.R0=R0+m1.R0; Res.RZ=RZ+m1.RZ; Res.RX=RX+m1.RX; Res.Rp=Rp+m1.Rp; Res.Rm=Rm+m1.Rm;
  Res.hInitData=hInitData;
  return(Res);
}

void cOptic2D_Media::operator+=(const cOptic2D_Media &m1) {
  AZp+=m1.AZp; AZm+=m1.AZm; AXp+=m1.AXp; AXm+=m1.AXm;
  PZp+=m1.PZp; PZm+=m1.PZm; PXp+=m1.PXp; PXm+=m1.PXm;
  R0+=m1.R0; RZ+=m1.RZ; RX+=m1.RX; Rp+=m1.Rp; Rm+=m1.Rm;
}

void cOptic2D_Media::operator+=(const cOptic2D_Media m1) {
  AZp+=m1.AZp; AZm+=m1.AZm; AXp+=m1.AXp; AXm+=m1.AXm;
  PZp+=m1.PZp; PZm+=m1.PZm; PXp+=m1.PXp; PXm+=m1.PXm;
  R0+=m1.R0; RZ+=m1.RZ; RX+=m1.RX; Rp+=m1.Rp; Rm+=m1.Rm;
}

cOptic2D_Media cOptic2D_Media::operator*(long double m) {
  cOptic2D_Media Res;
  Res.AZp=AZp*m; Res.AZm=AZm*m; Res.AXp=AXp*m; Res.AXm=AXm*m;
  Res.PZp=PZp*m; Res.PZm=PZm*m; Res.PXp=PXp*m; Res.PXm=PXm*m;
  Res.R0=R0*m; Res.RZ=RZ*m; Res.RX=RX*m; Res.Rp=Rp*m; Res.Rm=Rm*m;
  Res.hInitData=hInitData;
  return(Res);
}

TCplxLong cOptic2D_Media::RightAZp()
    {return( TCplxLong(0,-1)*hInitData->ACoef*(AXp+AXm)+PZp+hInitData->AGenCoef*AZp );}
TCplxLong cOptic2D_Media::RightAZm()
    {return( TCplxLong(0,-1)*hInitData->ACoef*(AXp+AXm)+PZm+hInitData->AGenCoef*AZm );}
TCplxLong cOptic2D_Media::RightAXp()
    {return( TCplxLong(0,-1)*hInitData->ACoef*(AZp+AZm)+PXp+hInitData->AGenCoef*AXp );}
TCplxLong cOptic2D_Media::RightAXm()
    {return( TCplxLong(0,-1)*hInitData->ACoef*(AZp+AZm)+PXm+hInitData->AGenCoef*AXm );}
TCplxLong cOptic2D_Media::RightPZp()
    {return((AZp*R0*((long double)2.0)+AZm*RZ+AXp*Rm+AXm*Rp)*hInitData->beta -
               PZp*hInitData->PCoef);}
TCplxLong cOptic2D_Media::RightPZm()
    {return((AZm*R0*((long double)2.0)+AZp*conj(RZ)+AXp*conj(Rp)+AXm*conj(Rm))*hInitData->beta -
               PZm*hInitData->PCoef);}
TCplxLong cOptic2D_Media::RightPXp()
    {return((AXp*R0*((long double)2.0)+AXm*RX+AZp*conj(Rm)+AZm*Rp)*hInitData->beta -
               PXp*hInitData->PCoef);}
TCplxLong cOptic2D_Media::RightPXm()
    {return((AXm*R0*((long double)2.0)+AXp*RZ+AZp*conj(Rp)+AZm*Rm)*hInitData->beta -
               PXm*hInitData->PCoef);}
TCplxLong cOptic2D_Media::RightPZpShort()
  {return((AZp*R0*((long double)2.0)+AZm*RZ+AXp*Rm+AXm*Rp)*hInitData->beta);}
TCplxLong cOptic2D_Media::RightPZmShort()
  {return((AZm*R0*((long double)2.0)+AZp*conj(RZ)+AXp*conj(Rp)+AXm*conj(Rm))*hInitData->beta);}
TCplxLong cOptic2D_Media::RightPXpShort()
  {return((AXp*R0*((long double)2.0)+AXm*RX+AZp*conj(Rm)+AZm*Rp)*hInitData->beta);}
TCplxLong cOptic2D_Media::RightPXmShort()
  {return((AXm*R0*((long double)2.0)+AXp*RZ+AZp*conj(Rp)+AZm*Rm)*hInitData->beta);}
long double cOptic2D_Media::RightR0() {
  return(hInitData->Q-hInitData->RCoef*R0-(AZp*conj(PZp)+AXp*conj(PXp)+AZm*conj(PZm)+AXm*conj(PXm)).real() );
}


TCplxLong cOptic2D_Media::RightRZ() {
  return(-RZ*hInitData->RCoef -(AZp*conj(PZm)+PZp*conj(AZm)) );
}
TCplxLong cOptic2D_Media::RightRX() {
  return(-RX*hInitData->RCoef -(AXp*conj(PXm)+PXp*conj(AXm)) );
}
TCplxLong cOptic2D_Media::RightRp() {
  return(-Rp*hInitData->RCoef -(PZp*conj(AXm)+AZp*conj(PXm)+PXp*conj(AZm)+AXp*conj(PZm)) );
}
TCplxLong cOptic2D_Media::RightRm() {
  return(-Rm*hInitData->RCoef -(PXm*conj(AZm)+AZp*conj(PXp)+AXm*conj(PZm)+PZp*conj(AXp)) );
}

/*
TCplxLong cOptic2D_Media::RightRZ() {
  return(0);
}
TCplxLong cOptic2D_Media::RightRX() {
  return(0);
}
TCplxLong cOptic2D_Media::RightRp() {
  return(0);
}
TCplxLong cOptic2D_Media::RightRm() {
  return(0);
} */
/*
TCplxLong cOptic2D_Media::RightRZ() {
  return(-RZ*hInitData->RCoef -(AZp*conj(PZm)+PZp*conj(AZm)) );
}
TCplxLong cOptic2D_Media::RightRX() {
  return(-RX*hInitData->RCoef -(AXp*conj(PXm)+PXp*conj(AXm)) );
}
TCplxLong cOptic2D_Media::RightRp() {
  return(-Rp*hInitData->RCoef -(PZp*conj(AXm)+AZp*conj(PXm)+PXp*conj(AZm)+AXp*conj(PZm)) );
}
TCplxLong cOptic2D_Media::RightRm() {
  return(-Rm*hInitData->RCoef -(PXm*conj(AZm)+AZp*conj(PXp)+AXm*conj(PZm)+PZp*conj(AXp)) );
}
*/

TCplxLong cOptic2D_Media::RightPZpSimple()
    {return((AZp*R0*((long double)2.0))*hInitData->beta - PZp*hInitData->PCoef);}
TCplxLong cOptic2D_Media::RightPZmSimple()
    {return((AZm*R0*((long double)2.0))*hInitData->beta - PZm*hInitData->PCoef);}
TCplxLong cOptic2D_Media::RightPXpSimple()
    {return((AXp*R0*((long double)2.0))*hInitData->beta - PXp*hInitData->PCoef);}
TCplxLong cOptic2D_Media::RightPXmSimple()
    {return((AXm*R0*((long double)2.0))*hInitData->beta - PXm*hInitData->PCoef);}
TCplxLong cOptic2D_Media::RightPZpShortSimple()
  {return((AZp*R0*((long double)2.0))*hInitData->beta);}
TCplxLong cOptic2D_Media::RightPZmShortSimple()
  {return((AZm*R0*((long double)2.0))*hInitData->beta);}
TCplxLong cOptic2D_Media::RightPXpShortSimple()
  {return((AXp*R0*((long double)2.0))*hInitData->beta);}
TCplxLong cOptic2D_Media::RightPXmShortSimple()
  {return((AXm*R0*((long double)2.0))*hInitData->beta);}
long double cOptic2D_Media::RightR0Simple() {
  return(hInitData->Q-hInitData->RCoef*R0-(AZp*conj(PZp)+AXp*conj(PXp)+AZm*conj(PZm)+AXm*conj(PXm)).real() );
}
TCplxLong cOptic2D_Media::RightRZSimple() { return(0); }
TCplxLong cOptic2D_Media::RightRXSimple() { return(0); }
TCplxLong cOptic2D_Media::RightRpSimple() { return(0); }
TCplxLong cOptic2D_Media::RightRmSimple() { return(0); }

cOptic2D_Media cOptic2D_Media::RightPartMedia() {
 cOptic2D_Media Res;
 Res.AZp=Res.AZm=Res.AXp=Res.AXm=0;
 Res.PZp=RightPZp(); Res.PZm=RightPZm();
 Res.PXp=RightPXp(); Res.PXm=RightPXm();
 Res.R0=RightR0();
 Res.RZ=RightRZ(); Res.RX=RightRX();
 Res.Rp=RightRp(); Res.Rm=RightRm();
 return(Res);
}

cOptic2D_Media cOptic2D_Media::RightPart() {
  cOptic2D_Media Res;
  if(hInitData->UseFastP) {
    PXp=hInitData->PCoefFast*RightPXpShort();
    PXm=hInitData->PCoefFast*RightPXmShort();
    PZp=hInitData->PCoefFast*RightPZpShort();
    PZm=hInitData->PCoefFast*RightPZmShort();
  }
  Res.AZp=RightAZp(); Res.AZm=RightAZm();
  Res.AXp=RightAXp(); Res.AXm=RightAXm();
  Res.PZp=RightPZp(); Res.PZm=RightPZm();
  Res.PXp=RightPXp(); Res.PXm=RightPXm();
  Res.R0=RightR0();
  Res.RZ=RightRZ(); Res.RX=RightRX();
  Res.Rp=RightRp(); Res.Rm=RightRm();

  /*if(hInitData->UseFastP) {
    PXp=hInitData->PCoefFast*RightPXpShortSimple();
    PXm=hInitData->PCoefFast*RightPXmShortSimple();    // Only R0
    PZp=hInitData->PCoefFast*RightPZpShortSimple();
    PZm=hInitData->PCoefFast*RightPZmShortSimple();
  }
  Res.AZp=RightAZp(); Res.AZm=RightAZm();
  Res.AXp=RightAXp(); Res.AXm=RightAXm();
  Res.PZp=RightPZpSimple(); Res.PZm=RightPZmSimple();
  Res.PXp=RightPXpSimple(); Res.PXm=RightPXmSimple();
  Res.R0=RightR0Simple();
  Res.RZ=RightRZSimple(); Res.RX=RightRXSimple();
  Res.Rp=RightRpSimple(); Res.Rm=RightRmSimple();*/
  return(Res);
}

void cOptic2D_Media::MakeNewSimple(cOptic2D_Media *hNew,const int i,const int j){
  cOptic2D_Media k1,k2,k3,k4;
  k1=RightPart()*hInitData->dt;
  k2=(*this+k1*0.5).RightPart()*hInitData->dt;
  k3=(*this+k2).RightPart()*hInitData->dt;
  k4=(*this+k3).RightPart()*hInitData->dt;
  *hNew=*this+(k1+k2*2.0+k3*2.0+k4)*(((long double)1.0)/6.0);
}
void cOptic2D_Media::MakeNewSeveralSteps(cOptic2D_Media *hNew,const int i,const int j){
  cOptic2D_Media tmpNew=*this;
  long double tmpEnergy=AZp.real()*AZp.real()+AZp.imag()*AZp.imag()+
                        AZm.real()*AZm.real()+AZm.imag()*AZm.imag()+
                        AXp.real()*AXp.real()+AXp.imag()*AXp.imag()+
                        AXm.real()*AXm.real()+AXm.imag()*AXm.imag();
  long double tmpdt=hInitData->dt;
  hInitData->dt=tmpdt/hInitData->nTSubSteps;
  for(int i=0;i<hInitData->nTSubSteps;i++){
    if(hInitData->UseSimpleEC) tmpNew.MakeNewSimpleEC1(hNew,i,j);
      else tmpNew.MakeNewSimple(hNew,i,j);
    tmpNew=*hNew;
  }
  hInitData->dt=tmpdt;
  long double tmpEnergyNew=hNew->AZp.real()*hNew->AZp.real()+hNew->AZp.imag()*hNew->AZp.imag()+
           hNew->AZm.real()*hNew->AZm.real()+hNew->AZm.imag()*hNew->AZm.imag()+
           hNew->AXp.real()*hNew->AXp.real()+hNew->AXp.imag()*hNew->AXp.imag()+
           hNew->AXm.real()*hNew->AXm.real()+hNew->AXm.imag()*hNew->AXm.imag();
  hInitData->EWrong+=(tmpEnergy-tmpEnergyNew)*tmpdt*tmpdt;
}

void cOptic2D_Media::MakeNewSimpleEC1(cOptic2D_Media *hNew,const int i,const int j){
  cOptic2D_Media k1,k2,k3,k4,middle;
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

void cOptic2D_Media::MakeNewSimpleEC(cOptic2D_Media *hNew,const int i,const int j){
  cOptic2D_Media k1,k2,k3,k4,middle;
  long double energyXP=AXp.real()*AXp.real()+AXp.imag()*AXp.imag();
  long double energyXM=AXm.real()*AXm.real()+AXm.imag()*AXm.imag();
  long double energyZP=AZp.real()*AZp.real()+AZp.imag()*AZp.imag();
  long double energyZM=AZm.real()*AZm.real()+AZm.imag()*AZm.imag();
  k1=RightPart()*hInitData->dt;
  k2=(*this+k1*0.5).RightPart()*hInitData->dt;
  k3=(*this+k2).RightPart()*hInitData->dt;
  k4=(*this+k3).RightPart()*hInitData->dt;
  *hNew=*this+(k1+k2*2.0+k3*2.0+k4);

  long double energyXP1=hNew->AXp.real()*AXp.real()+hNew->AXp.imag()*AXp.imag();
  long double energyXM1=hNew->AXm.real()*AXm.real()+hNew->AXm.imag()*AXm.imag();
  long double energyZP1=hNew->AZp.real()*AZp.real()+hNew->AZp.imag()*AZp.imag();
  long double energyZM1=hNew->AZm.real()*AZm.real()+hNew->AZm.imag()*AZm.imag();
  middle=(*this+*hNew)*(long double)0.5;
  long double DeltaXpZm=(TCplxLong(0,1)*conj(middle.AXp)*middle.AZm).real()*(2*hInitData->ACoef);
  long double DeltaXpZp=(TCplxLong(0,1)*conj(middle.AXp)*middle.AZp).real()*(2*hInitData->ACoef);
  long double DeltaXmZm=(TCplxLong(0,1)*conj(middle.AXm)*middle.AZm).real()*(2*hInitData->ACoef);
  long double DeltaXmZp=(TCplxLong(0,1)*conj(middle.AXm)*middle.AZp).real()*(2*hInitData->ACoef);
  long double tmpCoef;
  long double E[4],E1[4],dE[4];
  int Ez[4]; Ez[0]=Ez[1]=Ez[2]=Ez[3]=0;
  E[0]=AXp.real()*AXp.real()+AXp.imag()*AXp.imag();
  E[2]=AXm.real()*AXm.real()+AXm.imag()*AXm.imag();
  E[1]=AZp.real()*AZp.real()+AZp.imag()*AZp.imag();
  E[3]=AZm.real()*AZm.real()+AZm.imag()*AZm.imag();
  dE[0]=DeltaXpZm;
  dE[1]=-DeltaXmZm;          //E1[i]=E[i]+dE[i]-dE[(i+3)%4]
  dE[2]=DeltaXmZm;
  dE[3]=-DeltaXpZm;
  for(int i=0;i<4;i++) {
    if(E[i]+dE[i]-dE[(i+3)%4]<0){
      if((!Ez[i+1])||(!Ez[(i+3)%4])) {
        tmpCoef=-E[i]/(dE[i]-dE[(i+3)%4]);
        dE[i]*=tmpCoef;  dE[(i+3)%4]*=tmpCoef;
        Ez[i]=0;
      } else {

      }
    }
  }
  if( (energyXP+DeltaXpZm+DeltaXpZp<0) ) {
    AXp=0; energyXP1=0;
    tmpCoef=energyXP/(-DeltaXpZm-DeltaXpZp);
    DeltaXpZm*=tmpCoef;  DeltaXpZp*=tmpCoef;
  }
  if( (energyXM+DeltaXmZm+DeltaXmZp<0) ) {
    AXm=0; energyXM1=0;
    tmpCoef=energyXM/(-DeltaXmZm-DeltaXmZp);
    DeltaXmZm*=tmpCoef;  DeltaXmZp*=tmpCoef;
  }
  if( (energyZP-DeltaXpZp-DeltaXmZp<0) ) {
    AZp=0; energyZP1=0;
    tmpCoef=energyZP/(DeltaXpZp+DeltaXmZp);
    DeltaXpZp*=tmpCoef;  DeltaXmZp*=tmpCoef;
  }
  if( (energyZM-DeltaXmZm-DeltaXpZm<0) ) {
    AZm=0; energyZM1=0;
    tmpCoef=energyZM/(DeltaXmZm+DeltaXpZm);
    DeltaXmZm*=tmpCoef;  DeltaXpZm*=tmpCoef;
  }
  if(energyXP1>0.000000001) {
    hNew->AXp*=sqrt((energyXP+DeltaXpZm+DeltaXpZp)/energyXP1);
  }
  if(energyXM1>0.000000001) {
    hNew->AXm*=sqrt((energyXM+DeltaXmZm+DeltaXmZp)/energyXM1);
  }
  if(energyZP1>0.000000001) {
    hNew->AZp*=sqrt((energyZP-DeltaXmZp-DeltaXpZp)/energyZP1);
  }
  if(energyZM1>0.000000001) {
    hNew->AZm*=sqrt((energyZM-DeltaXpZm-DeltaXmZm)/energyZM1);
  }
}

void cOptic2D_Media::MakeK(cOptic2D_Media *hk,cOptic2D_Media *hkNew,long double coef) {
  cOptic2D_Media tmp=*this;                               // --- Media
  //if(hk) tmp+=(*hk)*coef;
  if(hk) tmp+=(*hk)*coef;
  tmp.AZp=AZp*(1-coef)+(this-hInitData->SizeX)->AZp*coef;
  tmp.AZm=AZm*(1-coef)+(this+hInitData->SizeX)->AZm*coef;
  tmp.AXp=AXp*(1-coef)+(this-1)->AXp*coef;
  tmp.AXm=AXm*(1-coef)+(this+1)->AXm*coef;
  if(hk) {
   tmp.AZp=(hk->AZp*(1-coef)+(hk-hInitData->SizeX)->AZp*coef)*coef;
   tmp.AZm=(hk->AZm*(1-coef)+(hk+hInitData->SizeX)->AZm*coef)*coef;
   tmp.AXp=(hk->AXp*(1-coef)+(hk-1)->AXp*coef)*coef;
   tmp.AXm=(hk->AXm*(1-coef)+(hk+1)->AXm*coef)*coef;
  }
  cOptic2D_Media Res=tmp.RightPartMedia();

  tmp=(*this)*(1-coef);                                   // --- AZp, wrong AZm!
  tmp+=(*(this+hInitData->SizeX))*coef;
  tmp.AZp=AZp;
  tmp.AXp=AXp*(1-coef)+(this+hInitData->SizeX-1)->AXp*coef;
  tmp.AXm=AXm*(1-coef)+(this+hInitData->SizeX+1)->AXm*coef;
  if (hk) {
    tmp+=(*hk)*(1-coef)*coef;
    tmp+=(*(hk+hInitData->SizeX))*coef;
    tmp.AZp=hk->AZp*coef;
    tmp.AXp=(hk->AXp*(1-coef)+(hk+hInitData->SizeX-1)->AXp*coef)*coef;
    tmp.AXm=(hk->AXm*(1-coef)+(hk+hInitData->SizeX+1)->AXm*coef)*coef;
  }
  Res.AZp=tmp.RightAZp();

  tmp=(*this)*(1-coef);                                   // --- AZm, wrong AZp!
  tmp+=(*(this-hInitData->SizeX))*coef;
  tmp.AZm=AZm;
  tmp.AXp=AXp*(1-coef)+(this-hInitData->SizeX-1)->AXp*coef;
  tmp.AXm=AXm*(1-coef)+(this-hInitData->SizeX+1)->AXm*coef;
  if (hk) {
    tmp+=(*hk)*(1-coef)*coef;
    tmp+=(*(hk-hInitData->SizeX))*coef;
    tmp.AZm=hk->AZm*coef;
    tmp.AXp=(hk->AXp*(1-coef)+(hk-hInitData->SizeX-1)->AXp*coef)*coef;
    tmp.AXm=(hk->AXm*(1-coef)+(hk-hInitData->SizeX+1)->AXm*coef)*coef;
  }
  Res.AZm=tmp.RightAZm();

  tmp=(*this)*(1-coef);                                   // --- AXp, wrong AXm!
  tmp+=(*(this+1))*coef;
  tmp.AXp=AXp;
  tmp.AZp=AZp*(1-coef)+(this-hInitData->SizeX+1)->AZp*coef;
  tmp.AZm=AZm*(1-coef)+(this+hInitData->SizeX+1)->AZm*coef;
  if (hk) {
    tmp+=(*hk)*(1-coef)*coef;
    tmp+=(*(hk+1))*coef;
    tmp.AXp=hk->AXp*coef;
    tmp.AZp=(hk->AZp*(1-coef)+(hk-hInitData->SizeX+1)->AZp*coef)*coef;
    tmp.AZm=(hk->AZm*(1-coef)+(hk+hInitData->SizeX+1)->AZm*coef)*coef;
  }
  Res.AXp=tmp.RightAXp();

  tmp=(*this)*(1-coef);                                   // --- AXm, wrong AXp!
  tmp+=(*(this-1))*coef;
  tmp.AXm=AXm;
  tmp.AZp=AZp*(1-coef)+(this-hInitData->SizeX-1)->AZp*coef;
  tmp.AZm=AZm*(1-coef)+(this+hInitData->SizeX-1)->AZm*coef;
  if (hk) {
    tmp+=(*hk)*(1-coef)*coef;
    tmp+=(*(hk-hInitData->SizeX))*coef;
    tmp.AXm=hk->AXm*coef;
    tmp.AZp=(hk->AZp*(1-coef)+(hk-hInitData->SizeX-1)->AZp*coef)*coef;
    tmp.AZm=(hk->AZm*(1-coef)+(hk+hInitData->SizeX-1)->AZm*coef)*coef;
  }
  Res.AXm=tmp.RightAXm();
  Res.hInitData=hInitData;
  *hkNew=Res*hInitData->dt;
}

void cOptic2D_Media_Dump::GetData(cOptic2D_Media *hSource) {
  AZp=abs(hSource->AZp); AZm=abs(hSource->AZm);
  AXp=abs(hSource->AXp); AXm=abs(hSource->AXm);
  PZp=abs(hSource->PZp); PZm=abs(hSource->PZm);
  PXp=abs(hSource->PXp); PXm=abs(hSource->PXm);
 /* R0=hSource->R0; RZ=hSource->RZ; RX=hSource->RX;
  Rp=hSource->Rp; Rm=hSource->Rm;*/
}
