//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquOpticTask.h"
#include <valarray>

//---------------------------------------------------------------------------

#pragma package(smart_init)

OPTIC_TYPE qqPZp0,qqPZp1,qqPZp2,qqAZp0,qqAZp1,qqAZp2,qqV1,qqV2,qqGamma;

void cEquOpticBaseTask::Step() {
  //StepNoMedia();
/*
  if(InitData.UseFastP)
    StepNoMediaEuler();
  else
    StepNoMediaRunge2();
*//*
  if(InitData.UseFastP)
    StepNoMediaMatrix();
  else
    StepNoMediaRunge2();
*/
  int ix=30,iz=15;
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  StepNoMediaMatrix();
//  if(!InitData.MicrowaveModeLength) StepNoMediaMatrix();
//  else StepNoMediaMatrixMMode();
  
  if(!InitData.fNoMedia) {
    if(InitData.UseFastP) {
      if(!InitData.fTE) StepMediaRunge2FastP();
      else StepMediaRunge2FastPTE();
    } else {
      if(InitData.UseSimpleEC) {
        StepMediaRunge2FixedNonlinear();
      } else {
        if(!InitData.fTE)
          if(InitData.UseSimpleEC) StepMediaRunge2NoR2();
            else StepMediaRunge2();
        else StepMediaRunge2TE();
      }
    }
  }
 /*
  if(RungeStep==2) {
    qqPZp0=abs(hM[0*Size+IndexXZ(ix,iz)].PZp);
    qqPZp1=abs(hM[1*Size+IndexXZ(ix,iz)].PZp);
    qqPZp2=abs(hM[2*Size+IndexXZ(ix,iz)].PZp);
    qqAZp0=abs(hAZp[IndexXZ(ix,iz)]);
    qqAZp1=abs(hAZp[1*Size+IndexXZ(ix,iz)]);
    qqAZp2=abs(hAZp[2*Size+IndexXZ(ix,iz)]);
    qqV1=qqPZp0-qqGamma*qqAZp0;
    qqV2=InitData.beta*qqAZp0-InitData.PCoef*qqPZp0;
  }*/
}

void cEquOpticBaseTask::FillMatrixes() {
  OPTIC_TYPE dtT=InitData.dt;
  valarray< complex<OPTIC_TYPE> > VA[4];
  complex<OPTIC_TYPE> kk(0,-0.5*dtT*InitData.ACoef);
  kk *= polar(1,InitData.APhase*M_PI);

  for(int i=0;i<4;i++) {                                     //-------- 4x4
    VA[i].resize(8,0);
    VA[i][i]=1;
    VA[i][i+4]=1;
  }
  VA[0][2]=VA[0][3]=VA[1][2]=VA[1][3]=VA[2][0]=VA[2][1]=VA[3][0]=VA[3][1]=kk;

  for(int i=0;i<3;i++) {
    VA[i]=VA[i]/VA[i][i];
    for(int j=i+1;j<4;j++)
      VA[j]=VA[j]-VA[i]*VA[j][i];
  }
  VA[3]=VA[3]/VA[3][3];
  for(int i=1;i<4;i++)
    for(int j=0;j<i;j++)
      VA[j]=VA[j]-VA[i]*VA[j][i];

  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++)
      Matrix4[i*4+j]=VA[i][4+j];
                                                             //-------- 3x3
  for(int i=0;i<3;i++) {
    VA[i]=complex<OPTIC_TYPE>(0);
    VA[i][i]=1;
    VA[i][i+4]=1;
  }
  VA[0][2]=VA[1][2]=VA[2][0]=VA[2][1]=kk;
  for(int i=0;i<2;i++) {
    VA[i]=VA[i]/VA[i][i];
    for(int j=i+1;j<3;j++)
      VA[j]=VA[j]-VA[i]*VA[j][i];
  }
  VA[2]=VA[2]/VA[2][2];
  for(int i=1;i<3;i++)
    for(int j=0;j<i;j++)
      VA[j]=VA[j]-VA[i]*VA[j][i];

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      Matrix3[i*3+j]=VA[i][4+j];
                                                             //-------- 2x2
  complex<OPTIC_TYPE> tmp1mkk=OPTIC_TYPE(1)-kk*kk;
  Matrix2[0]=Matrix2[3]=OPTIC_TYPE(1)/tmp1mkk;
  Matrix2[1]=Matrix2[2]=-Matrix2[0]*kk;
//  complex<OPTIC_TYPE> m1=OPTIC_TYPE(1.0)/tmp1mkk;
//  complex<OPTIC_TYPE> m2=-m1*kk;
  /*
  VA[0]=VA[1]=complex<OPTIC_TYPE>(0);
  VA[0][1]=VA[1][0]=kk;
  VA[0][0]=VA[1][1]=VA[0][4]=VA[1][5]=1;
  for(int i=0;i<1;i++) {
    VA[i]=VA[i]/VA[i][i];
    for(int j=i+1;j<2;j++)
      VA[j]=VA[j]-VA[i]*VA[j][i];
  }
  VA[1]=VA[1]/VA[1][1];
  for(int i=1;i<2;i++)
    for(int j=0;j<i;j++)
      VA[j]=VA[j]-VA[i]*VA[j][i];

  for(int i=0;i<2;i++)
    for(int j=0;j<2;j++)
      Matrix2[i*2+j]=VA[i][4+j];*/
}

/*
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
*/

void cEquOpticBaseTask::StepMediaRunge2NoR2() {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  OPTIC_TYPE dtT2=0.5*dtT;
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1),tmpAXp,tmpAXm,tmpAZp,tmpAZm;
  complex<OPTIC_TYPE> tmpPXp,tmpPXm,tmpPZp,tmpPZm;
  cOptic2D_Media_X tM,tM1;
  complex<OPTIC_TYPE> tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9;
  OPTIC_TYPE rtmp0,rtmp1,rtmp2,rtmp3,rtmp4,rtmp5,rtmp6,rtmp7,rtmp8,rtmp9;
  int istopX=15,istopZ=30;
  int iI;
  cOptic2D_Media_X tmpMedia;
  switch(MediaRungeStep) {
    case 0:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        iI=IndexXZ(ix,iz);
        hM[Size+iI].PXp=((hAXp[iI]*hM[iI].R0*((OPTIC_TYPE)2.0))*InitData.beta -
                        hM[iI].PXp*InitData.PCoef)*dtT2;
        hM[Size+iI].PXm=((hAXm[iI]*hM[iI].R0*((OPTIC_TYPE)2.0))*InitData.beta -
                        hM[iI].PXm*InitData.PCoef)*dtT2;
        hM[Size+iI].PZp=((hAZp[iI]*hM[iI].R0*((OPTIC_TYPE)2.0))*InitData.beta -
                        hM[iI].PZp*InitData.PCoef)*dtT2;
        hM[Size+iI].PZm=((hAZm[iI]*hM[iI].R0*((OPTIC_TYPE)2.0))*InitData.beta -
                        hM[iI].PZm*InitData.PCoef)*dtT2;
        hM[Size+iI].RX=0;
        hM[Size+iI].RZ=0;
        hM[Size+iI].Rp=0;
        hM[Size+iI].Rm=0;
        hM[Size+iI].R0=(InitData.Q-InitData.RCoef*hM[iI].R0
                       -(hAZp[iI]*conj(hM[iI].PZp)
                        +hAXp[iI]*conj(hM[iI].PXp)
                        +hAZm[iI]*conj(hM[iI].PZm)
                        +hAXm[iI]*conj(hM[iI].PXm)).real() )*dtT2;
        hM[Size+iI].R0Generated=(InitData.Q-InitData.RCoef*hM[iI].R0)*dtT2;


        hAXp[Size+iI]+=dtT2*hM[IndexXZ(ix-1,iz)].PXp
                       +dtT2*InitData.AGenCoef*hAXp[IndexXZ(ix-1,iz)]
                       +dtT2*InitData.ABackCoef*hAXm[IndexXZ(ix-1,iz)]
                       +dtT2*InitData.AGenX*hAXp[IndexXZ(ix-1,iz)];
        hAXm[Size+iI]+=dtT2*hM[IndexXZ(ix+1,iz)].PXm
                       +dtT2*InitData.AGenCoef*hAXm[IndexXZ(ix+1,iz)]
                       +dtT2*InitData.ABackCoef*hAXp[IndexXZ(ix+1,iz)]
                       +dtT2*InitData.AGenX*hAXm[IndexXZ(ix+1,iz)];
        hAZp[Size+iI]+=dtT2*hM[IndexXZ(ix,iz-1)].PZp
                       +dtT2*InitData.AGenCoef*hAZp[IndexXZ(ix,iz-1)]
                       +dtT2*InitData.ABackCoef*hAZm[IndexXZ(ix,iz-1)]
                       +dtT2*InitData.AGenZ*hAXp[IndexXZ(ix,iz-1)];
        hAZm[Size+iI]+=dtT2*hM[IndexXZ(ix,iz+1)].PZm
                       +dtT2*InitData.AGenCoef*hAZm[IndexXZ(ix,iz+1)]
                       +dtT2*InitData.ABackCoef*hAZp[IndexXZ(ix,iz+1)]
                       +dtT2*InitData.AGenZ*hAZm[IndexXZ(ix,iz+1)];
      }
    }
    MediaRungeStep++;
    break;
    case 1:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        iI=IndexXZ(ix,iz);
        tmpAXp=hAXp[IndexXZ(ix-1,iz)]+(OPTIC_TYPE)2.0*hAXp[Size+iI];
        tmpAXm=hAXm[IndexXZ(ix+1,iz)]+(OPTIC_TYPE)2.0*hAXm[Size+iI];
        tmpAZp=hAZp[IndexXZ(ix,iz-1)]+(OPTIC_TYPE)2.0*hAZp[Size+iI];
        tmpAZm=hAZm[IndexXZ(ix,iz+1)]+(OPTIC_TYPE)2.0*hAZm[Size+iI];
//        tM=hM[iI]; tM+=hM[Size+iI]*((OPTIC_TYPE)2.0);
//        hM[2*Size+iI]=hM[iI]; hM[2*Size+iI]+=hM[Size+iI];
        tM.M1M2v(hM[iI],hM[Size+iI],2);
        hM[2*Size+iI].M1M2v(hM[iI],hM[Size+iI],1);
//        hM[2*Size+iI]=hM[iI]; hM[2*Size+iI]+=hM[Size+iI];


        hM[2*Size+iI].PXp+=((tmpAXp*tM.R0*((OPTIC_TYPE)2.0))*InitData.beta -
                             tM.PXp*InitData.PCoef)*dtT2;
        hM[2*Size+iI].PXm+=((tmpAXm*tM.R0*((OPTIC_TYPE)2.0))*InitData.beta -
                             tM.PXm*InitData.PCoef)*dtT;
        hM[2*Size+iI].PZp+=((tmpAZp*tM.R0*((OPTIC_TYPE)2.0))*InitData.beta -
                             tM.PZp*InitData.PCoef)*dtT;
        hM[2*Size+iI].PZm+=((tmpAZm*tM.R0*((OPTIC_TYPE)2.0))*InitData.beta -
                             tM.PZm*InitData.PCoef)*dtT;
        hM[2*Size+iI].RX+=0;
        hM[2*Size+iI].RZ+=0;
        hM[2*Size+iI].Rp+=0;
        hM[2*Size+iI].Rm+=0;
        hM[2*Size+iI].R0+=(InitData.Q-InitData.RCoef*tM.R0
                          -(tmpAZp*conj(tM.PZp)
                          +tmpAXp*conj(tM.PXp)
                          +tmpAZm*conj(tM.PZm)
                          +tmpAXm*conj(tM.PXm)).real() )*dtT2;
        hM[2*Size+iI].R0Generated+=(InitData.Q-InitData.RCoef*tM.R0)*dtT2;
        hAXp[2*Size+iI]+=dtT2*tM.PXp+dtT2*InitData.AGenCoef*tmpAXp
                                    +dtT2*InitData.ABackCoef*tmpAXm
                                    +dtT2*InitData.AGenX*tmpAXp;
        hAXm[2*Size+iI]+=dtT2*tM.PXm+dtT2*InitData.AGenCoef*tmpAXm
                                    +dtT2*InitData.ABackCoef*tmpAXp
                                    +dtT2*InitData.AGenX*tmpAXm;
        hAZp[2*Size+iI]+=dtT2*tM.PZp+dtT2*InitData.AGenCoef*tmpAZp
                                    +dtT2*InitData.ABackCoef*tmpAZm
                                    +dtT2*InitData.AGenZ*tmpAZp;
        hAZm[2*Size+iI]+=dtT2*tM.PZm+dtT2*InitData.AGenCoef*tmpAZm
                                    +dtT2*InitData.ABackCoef*tmpAZp
                                    +dtT2*InitData.AGenZ*tmpAZm;
      }
    }
    MediaRungeStep++;
    break;
    case 2:
      memcpy(hM+IndexXZ(ix1st,0),hM+2*Size+IndexXZ(ix1st,0),
        (InitData.nStepsZ+2)*sizeof(cOptic2D_Media_X)*(ixLast-ix1st+1));
    MediaRungeStep=0;
    break;
  }
}

void cEquOpticBaseTask::StepMediaRunge2() {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  OPTIC_TYPE dtT2=0.5*dtT;
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1),tmpAXp,tmpAXm,tmpAZp,tmpAZm;
  complex<OPTIC_TYPE> tmpPXp,tmpPXm,tmpPZp,tmpPZm;
  cOptic2D_Media_X tM,tM1;
  complex<OPTIC_TYPE> tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9;
  OPTIC_TYPE rtmp0,rtmp1,rtmp2,rtmp3,rtmp4,rtmp5,rtmp6,rtmp7,rtmp8,rtmp9;
  int istopX=15,istopZ=30;
  int iI;
  cOptic2D_Media_X tmpMedia;
  complex<OPTIC_TYPE> PCoef=InitData.PCoef;
  PCoef=complex<OPTIC_TYPE>(InitData.PCoef,InitData.Delta);
  switch(MediaRungeStep) {
    case 0:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        iI=IndexXZ(ix,iz);
        hM[Size+iI].PXp=((hAXp[iI]*hM[iI].R0*((OPTIC_TYPE)2.0)+
                        hAXm[iI]*hM[iI].RX+hAZp[iI]*conj(hM[iI].Rm)+
                        hAZm[iI]*hM[iI].Rp)*InitData.beta -
                        hM[iI].PXp*PCoef)*dtT2;
        if((istopX==ix)&&(istopZ==iz)) {
          tmp1=hAXp[iI]*hM[iI].R0*((OPTIC_TYPE)2.0)*InitData.beta;
          tmp2=hM[iI].PXp*PCoef;
          tmpPXp=(hAXp[iI]*hM[iI].R0*((OPTIC_TYPE)2.0)+
                hAXm[iI]*hM[iI].RX+hAZp[iI]*conj(hM[iI].Rm)+
                hAZm[iI]*hM[iI].Rp)*InitData.beta;
        }
        hM[Size+iI].PXm=((hAXm[iI]*hM[iI].R0*((OPTIC_TYPE)2.0)+
                        hAXp[iI]*hM[iI].RZ+hAZp[iI]*conj(hM[iI].Rp)+
                        hAZm[iI]*hM[iI].Rm)*InitData.beta -
                        hM[iI].PXm*PCoef)*dtT2;
        hM[Size+iI].PZp=((hAZp[iI]*hM[iI].R0*((OPTIC_TYPE)2.0)+
                        hAZm[iI]*hM[iI].RZ+hAXp[iI]*hM[iI].Rm+
                        hAXm[iI]*hM[iI].Rp)*InitData.beta -
                        hM[iI].PZp*PCoef)*dtT2;
        hM[Size+iI].PZm=((hAZm[iI]*hM[iI].R0*((OPTIC_TYPE)2.0)+
                        hAZp[iI]*conj(hM[iI].RZ)+hAXp[iI]*conj(hM[iI].Rp)+
                        hAXm[iI]*conj(hM[iI].Rm))*InitData.beta -
                        hM[iI].PZm*PCoef)*dtT2;
        hM[Size+iI].RX=(-hM[iI].RX*InitData.RCoef
                        -(hAXp[iI]*conj(hM[iI].PXm)
                        +hM[iI].PXp*conj(hAXm[iI])) )*dtT2;
        hM[Size+iI].RZ=(-hM[iI].RZ*InitData.RCoef
                        -(hAZp[iI]*conj(hM[iI].PZm)
                        +hM[iI].PZp*conj(hAZm[iI])) )*dtT2;
        hM[Size+iI].Rp=(-hM[iI].Rp*InitData.RCoef
                        -(hM[iI].PZp*conj(hAXm[iI])
                        +hAZp[iI]*conj(hM[iI].PXm)
                        +hM[iI].PXp*conj(hAZm[iI])
                        +hAXp[iI]*conj(hM[iI].PZm)) )*dtT2;
        hM[Size+iI].Rm=(-hM[iI].Rm*InitData.RCoef
                        -(hM[iI].PXm*conj(hAZm[iI])
                        +hAZp[iI]*conj(hM[iI].PXp)
                        +hAXm[iI]*conj(hM[iI].PZm)
                        +hM[iI].PZp*conj(hAXp[iI])) )*dtT2;
        hM[Size+iI].R0=(InitData.Q-InitData.RCoef*hM[iI].R0
                       -(hAZp[iI]*conj(hM[iI].PZp)
                        +hAXp[iI]*conj(hM[iI].PXp)
                        +hAZm[iI]*conj(hM[iI].PZm)
                        +hAXm[iI]*conj(hM[iI].PXm)).real() )*dtT2;
        hM[Size+iI].R0Generated=(InitData.Q-InitData.RCoef*hM[iI].R0)*dtT2;


        hAXp[Size+iI]+=dtT2*hM[IndexXZ(ix-1,iz)].PXp
                       +dtT2*InitData.AGenCoef*hAXp[IndexXZ(ix-1,iz)]
                       +dtT2*InitData.ABackCoef*hAXm[IndexXZ(ix-1,iz)]
                       +dtT2*InitData.AGenX*hAXp[IndexXZ(ix-1,iz)];
        hAXm[Size+iI]+=dtT2*hM[IndexXZ(ix+1,iz)].PXm
                       +dtT2*InitData.AGenCoef*hAXm[IndexXZ(ix+1,iz)]
                       +dtT2*InitData.ABackCoef*hAXp[IndexXZ(ix+1,iz)]
                       +dtT2*InitData.AGenX*hAXm[IndexXZ(ix+1,iz)];
        hAZp[Size+iI]+=dtT2*hM[IndexXZ(ix,iz-1)].PZp
                       +dtT2*InitData.AGenCoef*hAZp[IndexXZ(ix,iz-1)]
                       +dtT2*InitData.ABackCoef*hAZm[IndexXZ(ix,iz-1)]
                       +dtT2*InitData.AGenZ*hAXp[IndexXZ(ix,iz-1)];
        hAZm[Size+iI]+=dtT2*hM[IndexXZ(ix,iz+1)].PZm
                       +dtT2*InitData.AGenCoef*hAZm[IndexXZ(ix,iz+1)]
                       +dtT2*InitData.ABackCoef*hAZp[IndexXZ(ix,iz+1)]
                       +dtT2*InitData.AGenZ*hAZm[IndexXZ(ix,iz+1)];
      }
    }
    MediaRungeStep++;
    break;
    case 1:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        iI=IndexXZ(ix,iz);
        tmpAXp=hAXp[IndexXZ(ix-1,iz)]+(OPTIC_TYPE)2.0*hAXp[Size+iI];
        tmpAXm=hAXm[IndexXZ(ix+1,iz)]+(OPTIC_TYPE)2.0*hAXm[Size+iI];
        tmpAZp=hAZp[IndexXZ(ix,iz-1)]+(OPTIC_TYPE)2.0*hAZp[Size+iI];
        tmpAZm=hAZm[IndexXZ(ix,iz+1)]+(OPTIC_TYPE)2.0*hAZm[Size+iI];
//        tM=hM[iI]; tM+=hM[Size+iI]*((OPTIC_TYPE)2.0);
//        hM[2*Size+iI]=hM[iI]; hM[2*Size+iI]+=hM[Size+iI];
        tM.M1M2v(hM[iI],hM[Size+iI],2);
        hM[2*Size+iI].M1M2v(hM[iI],hM[Size+iI],1);
//        hM[2*Size+iI]=hM[iI]; hM[2*Size+iI]+=hM[Size+iI];


        hM[2*Size+iI].PXp+=((tmpAXp*tM.R0*((OPTIC_TYPE)2.0)+
                             tmpAXm*tM.RX+tmpAZp*conj(tM.Rm)+
                             tmpAZm*tM.Rp)*InitData.beta -
                             tM.PXp*PCoef)*dtT2;
        if((istopX==ix)&&(istopZ==iz)) {
          tmp6=tmpAXp*tM.R0*((OPTIC_TYPE)2.0)*InitData.beta;
          tmp7=tM.PXp*InitData.PCoef;
          tmp3=hM[iI].PXp; tmp4=hM[Size+iI].PXp; tmp5=hM[2*Size+iI].PXp;
          rtmp3=abs(tmp3); rtmp4=abs(tmp4); rtmp5=abs(tmp5);
          rtmp6=abs(tmp6); rtmp7=abs(tmp7);
        }
        hM[2*Size+iI].PXm+=((tmpAXm*tM.R0*((OPTIC_TYPE)2.0)+
                             tmpAXp*tM.RZ+tmpAZp*conj(tM.Rp)+
                             tmpAZm*tM.Rm)*InitData.beta -
                             tM.PXm*PCoef)*dtT;
        hM[2*Size+iI].PZp+=((tmpAZp*tM.R0*((OPTIC_TYPE)2.0)+
                             tmpAZm*tM.RZ+tmpAXp*tM.Rm+
                             tmpAXm*tM.Rp)*InitData.beta -
                             tM.PZp*PCoef)*dtT;
        hM[2*Size+iI].PZm+=((tmpAZm*tM.R0*((OPTIC_TYPE)2.0)+
                             tmpAZp*conj(tM.RZ)+tmpAXp*conj(tM.Rp)+
                             tmpAXm*conj(tM.Rm))*InitData.beta -
                             tM.PZm*PCoef)*dtT;
        hM[2*Size+iI].RX+=(-tM.RX*InitData.RCoef
                           -(tmpAXp*conj(tM.PXm)
                           +tM.PXp*conj(tmpAXm)) )*dtT;
        hM[2*Size+iI].RZ+=(-tM.RZ*InitData.RCoef
                           -(tmpAZp*conj(tM.PZm)
                           +tM.PZp*conj(tmpAZm)) )*dtT;
        hM[2*Size+iI].Rp+=(-tM.Rp*InitData.RCoef
                           -(tM.PZp*conj(tmpAXm)
                           +tmpAZp*conj(tM.PXm)
                           +tM.PXp*conj(tmpAZm)
                           +tmpAXp*conj(tM.PZm)) )*dtT;
        hM[2*Size+iI].Rm+=(-tM.Rm*InitData.RCoef
                           -(tM.PXm*conj(tmpAZm)
                           +tmpAZp*conj(tM.PXp)
                           +tmpAXm*conj(tM.PZm)
                           +tM.PZp*conj(tmpAXp)) )*dtT;
        hM[2*Size+iI].R0+=(InitData.Q-InitData.RCoef*tM.R0
                          -(tmpAZp*conj(tM.PZp)
                          +tmpAXp*conj(tM.PXp)
                          +tmpAZm*conj(tM.PZm)
                          +tmpAXm*conj(tM.PXm)).real() )*dtT2;
        hM[2*Size+iI].R0Generated+=(InitData.Q-InitData.RCoef*tM.R0)*dtT2;
        hAXp[2*Size+iI]+=dtT2*tM.PXp+dtT2*InitData.AGenCoef*tmpAXp
                                    +dtT2*InitData.ABackCoef*tmpAXm
                                    +dtT2*InitData.AGenX*tmpAXp;
        hAXm[2*Size+iI]+=dtT2*tM.PXm+dtT2*InitData.AGenCoef*tmpAXm
                                    +dtT2*InitData.ABackCoef*tmpAXp
                                    +dtT2*InitData.AGenX*tmpAXm;
        hAZp[2*Size+iI]+=dtT2*tM.PZp+dtT2*InitData.AGenCoef*tmpAZp
                                    +dtT2*InitData.ABackCoef*tmpAZm
                                    +dtT2*InitData.AGenZ*tmpAZp;
        hAZm[2*Size+iI]+=dtT2*tM.PZm+dtT2*InitData.AGenCoef*tmpAZm
                                    +dtT2*InitData.ABackCoef*tmpAZp
                                    +dtT2*InitData.AGenZ*tmpAZm;
      }
    }
    MediaRungeStep++;
    break;
    case 2:
      memcpy(hM+IndexXZ(ix1st,0),hM+2*Size+IndexXZ(ix1st,0),
        (InitData.nStepsZ+2)*sizeof(cOptic2D_Media_X)*(ixLast-ix1st+1));
    MediaRungeStep=0;
    break;
  }
}

void cEquOpticBaseTask::StepMediaRunge2TE() {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  OPTIC_TYPE dtT2=0.5*dtT;
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1),tmpAXp,tmpAXm,tmpAZp,tmpAZm;
  complex<OPTIC_TYPE> tmpPXp,tmpPXm,tmpPZp,tmpPZm;
  cOptic2D_Media_X tM,tM1;
  complex<OPTIC_TYPE> tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9;
  OPTIC_TYPE rtmp0,rtmp1,rtmp2,rtmp3,rtmp4,rtmp5,rtmp6,rtmp7,rtmp8,rtmp9;
  int istopX=15,istopZ=30;
  int iI;
  cOptic2D_Media_X tmpMedia;
  switch(MediaRungeStep) {
    case 0:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        iI=IndexXZ(ix,iz);
        hM[Size+iI].PXp=0;
        hM[Size+iI].PXm=0;
        hM[Size+iI].PZp=((hAZp[iI]*hM[iI].R0*((OPTIC_TYPE)2.0)+
                        hAZm[iI]*hM[iI].RZ+hAXp[iI]*hM[iI].Rm+
                        hAXm[iI]*hM[iI].Rp)*InitData.beta -
                        hM[iI].PZp*InitData.PCoef)*dtT2;
        hM[Size+iI].PZm=((hAZm[iI]*hM[iI].R0*((OPTIC_TYPE)2.0)+
                        hAZp[iI]*conj(hM[iI].RZ)+hAXp[iI]*conj(hM[iI].Rp)+
                        hAXm[iI]*conj(hM[iI].Rm))*InitData.beta -
                        hM[iI].PZm*InitData.PCoef)*dtT2;
        hM[Size+iI].RX=(-hM[iI].RX*InitData.RCoef
                        -(hAXp[iI]*conj(hM[iI].PXm)
                        +hM[iI].PXp*conj(hAXm[iI])) )*dtT2;
        hM[Size+iI].RZ=(-hM[iI].RZ*InitData.RCoef
                        -(hAZp[iI]*conj(hM[iI].PZm)
                        +hM[iI].PZp*conj(hAZm[iI])) )*dtT2;
        hM[Size+iI].Rp=(-hM[iI].Rp*InitData.RCoef
                        -(hM[iI].PZp*conj(hAXm[iI])
                        +hAZp[iI]*conj(hM[iI].PXm)
                        +hM[iI].PXp*conj(hAZm[iI])
                        +hAXp[iI]*conj(hM[iI].PZm)) )*dtT2;
        hM[Size+iI].Rm=(-hM[iI].Rm*InitData.RCoef
                        -(hM[iI].PXm*conj(hAZm[iI])
                        +hAZp[iI]*conj(hM[iI].PXp)
                        +hAXm[iI]*conj(hM[iI].PZm)
                        +hM[iI].PZp*conj(hAXp[iI])) )*dtT2;
        hM[Size+iI].R0=(InitData.Q-InitData.RCoef*hM[iI].R0
                       -(hAZp[iI]*conj(hM[iI].PZp)
                        +hAXp[iI]*conj(hM[iI].PXp)
                        +hAZm[iI]*conj(hM[iI].PZm)
                        +hAXm[iI]*conj(hM[iI].PXm)).real() )*dtT2;
        hM[Size+iI].R0Generated=(InitData.Q-InitData.RCoef*hM[iI].R0)*dtT2;


        hAXp[Size+iI]+=dtT2*hM[IndexXZ(ix-1,iz)].PXp
                       +dtT2*InitData.AGenCoef*hAXp[IndexXZ(ix-1,iz)]
                       +dtT2*InitData.ABackCoef*hAXm[IndexXZ(ix-1,iz)]
                       +dtT2*InitData.AGenX*hAXp[IndexXZ(ix-1,iz)];
        hAXm[Size+iI]+=dtT2*hM[IndexXZ(ix+1,iz)].PXm
                       +dtT2*InitData.AGenCoef*hAXm[IndexXZ(ix+1,iz)]
                       +dtT2*InitData.ABackCoef*hAXp[IndexXZ(ix+1,iz)]
                       +dtT2*InitData.AGenX*hAXm[IndexXZ(ix+1,iz)];
        hAZp[Size+iI]+=dtT2*hM[IndexXZ(ix,iz-1)].PZp
                       +dtT2*InitData.AGenCoef*hAZp[IndexXZ(ix,iz-1)]
                       +dtT2*InitData.ABackCoef*hAZm[IndexXZ(ix,iz-1)]
                       +dtT2*InitData.AGenZ*hAXp[IndexXZ(ix,iz-1)];
        hAZm[Size+iI]+=dtT2*hM[IndexXZ(ix,iz+1)].PZm
                       +dtT2*InitData.AGenCoef*hAZm[IndexXZ(ix,iz+1)]
                       +dtT2*InitData.ABackCoef*hAZp[IndexXZ(ix,iz+1)]
                       +dtT2*InitData.AGenZ*hAZm[IndexXZ(ix,iz+1)];
      }
    }
    MediaRungeStep++;
    break;
    case 1:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        iI=IndexXZ(ix,iz);
        tmpAXp=hAXp[IndexXZ(ix-1,iz)]+(OPTIC_TYPE)2.0*hAXp[Size+iI];
        tmpAXm=hAXm[IndexXZ(ix+1,iz)]+(OPTIC_TYPE)2.0*hAXm[Size+iI];
        tmpAZp=hAZp[IndexXZ(ix,iz-1)]+(OPTIC_TYPE)2.0*hAZp[Size+iI];
        tmpAZm=hAZm[IndexXZ(ix,iz+1)]+(OPTIC_TYPE)2.0*hAZm[Size+iI];
//        tM=hM[iI]; tM+=hM[Size+iI]*((OPTIC_TYPE)2.0);
//        hM[2*Size+iI]=hM[iI]; hM[2*Size+iI]+=hM[Size+iI];
        tM.M1M2v(hM[iI],hM[Size+iI],2);
        hM[2*Size+iI].M1M2v(hM[iI],hM[Size+iI],1);
//        hM[2*Size+iI]=hM[iI]; hM[2*Size+iI]+=hM[Size+iI];


        hM[2*Size+iI].PXp=0;
        hM[2*Size+iI].PXm=0;
        hM[2*Size+iI].PZp+=((tmpAZp*tM.R0*((OPTIC_TYPE)2.0)+
                             tmpAZm*tM.RZ+tmpAXp*tM.Rm+
                             tmpAXm*tM.Rp)*InitData.beta -
                             tM.PZp*InitData.PCoef)*dtT;
        hM[2*Size+iI].PZm+=((tmpAZm*tM.R0*((OPTIC_TYPE)2.0)+
                             tmpAZp*conj(tM.RZ)+tmpAXp*conj(tM.Rp)+
                             tmpAXm*conj(tM.Rm))*InitData.beta -
                             tM.PZm*InitData.PCoef)*dtT;
        hM[2*Size+iI].RX+=(-tM.RX*InitData.RCoef
                           -(tmpAXp*conj(tM.PXm)
                           +tM.PXp*conj(tmpAXm)) )*dtT;
        hM[2*Size+iI].RZ+=(-tM.RZ*InitData.RCoef
                           -(tmpAZp*conj(tM.PZm)
                           +tM.PZp*conj(tmpAZm)) )*dtT;
        hM[2*Size+iI].Rp+=(-tM.Rp*InitData.RCoef
                           -(tM.PZp*conj(tmpAXm)
                           +tmpAZp*conj(tM.PXm)
                           +tM.PXp*conj(tmpAZm)
                           +tmpAXp*conj(tM.PZm)) )*dtT;
        hM[2*Size+iI].Rm+=(-tM.Rm*InitData.RCoef
                           -(tM.PXm*conj(tmpAZm)
                           +tmpAZp*conj(tM.PXp)
                           +tmpAXm*conj(tM.PZm)
                           +tM.PZp*conj(tmpAXp)) )*dtT;
        hM[2*Size+iI].R0+=(InitData.Q-InitData.RCoef*tM.R0
                          -(tmpAZp*conj(tM.PZp)
                          +tmpAXp*conj(tM.PXp)
                          +tmpAZm*conj(tM.PZm)
                          +tmpAXm*conj(tM.PXm)).real() )*dtT2;
        hM[2*Size+iI].R0Generated+=(InitData.Q-InitData.RCoef*tM.R0)*dtT2;
        hAXp[2*Size+iI]+=dtT2*tM.PXp+dtT2*InitData.AGenCoef*tmpAXp
                                    +dtT2*InitData.ABackCoef*tmpAXm
                                    +dtT2*InitData.AGenX*tmpAXp;
        hAXm[2*Size+iI]+=dtT2*tM.PXm+dtT2*InitData.AGenCoef*tmpAXm
                                    +dtT2*InitData.ABackCoef*tmpAXp
                                    +dtT2*InitData.AGenX*tmpAXm;
        hAZp[2*Size+iI]+=dtT2*tM.PZp+dtT2*InitData.AGenCoef*tmpAZp
                                    +dtT2*InitData.ABackCoef*tmpAZm
                                    +dtT2*InitData.AGenZ*tmpAZp;
        hAZm[2*Size+iI]+=dtT2*tM.PZm+dtT2*InitData.AGenCoef*tmpAZm
                                    +dtT2*InitData.ABackCoef*tmpAZp
                                    +dtT2*InitData.AGenZ*tmpAZm;
      }
    }
    MediaRungeStep++;
    break;
    case 2:
      memcpy(hM+IndexXZ(ix1st,0),hM+2*Size+IndexXZ(ix1st,0),
        (InitData.nStepsZ+2)*sizeof(cOptic2D_Media_X)*(ixLast-ix1st+1));
    MediaRungeStep=0;
    break;
  }
}

void cEquOpticBaseTask::StepMediaRunge2FixedNonlinear() {
int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  OPTIC_TYPE dtT2=0.5*dtT;
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1);
  OPTIC_TYPE one=1, two=2;
  complex<OPTIC_TYPE> tmpAXp,tmpAXm,tmpAZp,tmpAZm,tmpPXp,tmpPXm,tmpPZp,tmpPZm,tSum;
  cOptic2D_Media_X tM,tM1;
  int iI,iI2;
  cOptic2D_Media_X tmpMedia;
  switch(MediaRungeStep) {
    case 0:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        iI  = IndexXZ(ix,iz);

        iI2 = IndexXZ(ix-1,iz);
        tmpPXp=hAXp[iI2];
        tmpPXp-=hAXp[iI2]*(hAXp[iI2]*conj(hAXp[iI2])+
                           hAXm[iI2]*conj(hAXm[iI2])*two+
                           hAZm[iI2]*conj(hAZm[iI2])*two+
                           hAZp[iI2]*conj(hAZp[iI2])*two
                           )*InitData.beta*InitData.RCoef;
        tmpPXp-=hAZp[iI2]*hAZm[iI2]*conj(hAXm[iI2])*two*
                        InitData.beta*InitData.RCoef;
        tmpPXp*=InitData.beta;

        iI2 = IndexXZ(ix+1,iz);
        tmpPXm=hAXm[iI2];
        tmpPXm-=hAXm[iI2]*(hAXp[iI2]*conj(hAXp[iI2])*two+
                           hAXm[iI2]*conj(hAXm[iI2])+
                           hAZm[iI2]*conj(hAZm[iI2])*two+
                           hAZp[iI2]*conj(hAZp[iI2])*two
                           )*InitData.beta*InitData.RCoef;
        tmpPXm-=hAZp[iI2]*hAZm[iI2]*conj(hAXp[iI2])*two*
                        InitData.beta*InitData.RCoef;
        tmpPXm*=InitData.beta;

        iI2 = IndexXZ(ix,iz-1);
        tmpPZp=hAZp[iI2];
        tmpPZp-=hAZp[iI2]*(hAXp[iI2]*conj(hAXp[iI2])*two+
                          hAXm[iI2]*conj(hAXm[iI2])*two+
                          hAZm[iI2]*conj(hAZm[iI2])*two+
                          hAZp[iI2]*conj(hAZp[iI2])
                          )*InitData.beta*InitData.RCoef;
        tmpPZp-=hAXp[iI2]*hAXm[iI2]*conj(hAZm[iI2])*two*
                        InitData.beta*InitData.RCoef;
        tmpPZp*=InitData.beta;
//        tmpPZp=hAZp[iI2]*InitData.beta;
        if(InitData.fTE) {
          tmpPZp = hAZp[iI2];
          tmpPZp-=hAZp[iI2]*(hAZm[iI2]*conj(hAZm[iI2])*two+
                            hAZp[iI2]*conj(hAZp[iI2])
                            )*InitData.beta*InitData.RCoef;
          tmpPZp*=InitData.beta;
        }

        iI2 = IndexXZ(ix,iz+1);
        tmpPZm=hAZm[iI2];
        tmpPZm-=hAZm[iI2]*(hAXp[iI2]*conj(hAXp[iI2])*two+
                          hAXm[iI2]*conj(hAXm[iI2])*two+
                          hAZm[iI2]*conj(hAZm[iI2])+
                          hAZp[iI2]*conj(hAZp[iI2])*two
                          )*InitData.beta*InitData.RCoef;
        tmpPZm-=hAXp[iI2]*hAXm[iI2]*conj(hAZp[iI2])*two*
                        InitData.beta*InitData.RCoef;
        tmpPZm*=InitData.beta;

        if(InitData.fTE) {
          tmpPZm=hAZm[iI2];
          tmpPZm-=hAZm[iI2]*(hAZm[iI2]*conj(hAZm[iI2])+
                            hAZp[iI2]*conj(hAZp[iI2])*two
                            )*InitData.beta*InitData.RCoef;
          tmpPZm*=InitData.beta;
          tmpPXp=0;
          tmpPXm=0;
        }


        hAXp[Size+iI]+=dtT2*tmpPXp
                      +dtT2*InitData.AGenCoef*hAXp[IndexXZ(ix-1,iz)]
                      +dtT2*InitData.ABackCoef*hAXm[IndexXZ(ix-1,iz)]
                      +dtT2*InitData.AGenX*hAXp[IndexXZ(ix-1,iz)];
        hAXm[Size+iI]+=dtT2*tmpPXm
                      +dtT2*InitData.AGenCoef*hAXm[IndexXZ(ix+1,iz)]
                      +dtT2*InitData.ABackCoef*hAXp[IndexXZ(ix+1,iz)]
                      +dtT2*InitData.AGenX*hAXm[IndexXZ(ix+1,iz)];
        hAZp[Size+iI]+=dtT2*tmpPZp
                      +dtT2*InitData.AGenCoef*hAZp[IndexXZ(ix,iz-1)]
                      +dtT2*InitData.ABackCoef*hAZm[IndexXZ(ix,iz-1)]
                      +dtT2*InitData.AGenZ*hAZp[IndexXZ(ix,iz-1)];
        hAZm[Size+iI]+=dtT2*tmpPZm
                      +dtT2*InitData.AGenCoef*hAZm[IndexXZ(ix,iz+1)]
                      +dtT2*InitData.ABackCoef*hAZm[IndexXZ(ix,iz+1)]
                      +dtT2*InitData.AGenZ*hAZm[IndexXZ(ix,iz+1)];
      }
    }
    MediaRungeStep++;
    break;
    case 1:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        iI=IndexXZ(ix,iz);
        tmpAXp=hAXp[IndexXZ(ix-1,iz)]+(OPTIC_TYPE)2.0*hAXp[Size+iI];
        tmpAXm=hAXm[IndexXZ(ix+1,iz)]+(OPTIC_TYPE)2.0*hAXm[Size+iI];
        tmpAZp=hAZp[IndexXZ(ix,iz-1)]+(OPTIC_TYPE)2.0*hAZp[Size+iI];
        tmpAZm=hAZm[IndexXZ(ix,iz+1)]+(OPTIC_TYPE)2.0*hAZm[Size+iI];
        tSum = (tmpAXp*conj(tmpAXp)+tmpAXm*conj(tmpAXm)+
                tmpAZp*conj(tmpAZp)+tmpAZm*conj(tmpAZm))*two;

        tmpPXp = tmpAXp*InitData.beta*(one-(tSum-tmpAXp*conj(tmpAXp))*InitData.beta*InitData.RCoef);
        tmpPXm = tmpAXm*InitData.beta*(one-(tSum-tmpAXm*conj(tmpAXm))*InitData.beta*InitData.RCoef);
        tmpPZp = tmpAZp*InitData.beta*(one-(tSum-tmpAZp*conj(tmpAZp))*InitData.beta*InitData.RCoef);
        tmpPZm = tmpAZm*InitData.beta*(one-(tSum-tmpAZm*conj(tmpAZm))*InitData.beta*InitData.RCoef);

        tmpPXp-= two*tmpAZp*tmpAZm*conj(tmpAXm)*InitData.beta*InitData.beta*InitData.RCoef;
        tmpPXm-= two*tmpAZp*tmpAZm*conj(tmpAXp)*InitData.beta*InitData.beta*InitData.RCoef;
        tmpPZp-= two*tmpAXp*tmpAXm*conj(tmpAZm)*InitData.beta*InitData.beta*InitData.RCoef;              -
        tmpPZm-= two*tmpAXp*tmpAXm*conj(tmpAZp)*InitData.beta*InitData.beta*InitData.RCoef;

        if(InitData.fTE) {
          tSum = (tmpAZp*conj(tmpAZp)+tmpAZm*conj(tmpAZm))*two;
          tmpPXp = 0;
          tmpPXm = 0;
          tmpPZp = tmpAZp*InitData.beta*(one-(tSum-tmpAZp*conj(tmpAZp))*InitData.beta*InitData.RCoef);
          tmpPZm = tmpAZm*InitData.beta*(one-(tSum-tmpAZm*conj(tmpAZm))*InitData.beta*InitData.RCoef);
        }
/*
        tmpPXp = tmpAXp*InitData.beta;
        tmpPXm = tmpAXm*InitData.beta;
        tmpPZp = tmpAZp*InitData.beta;
        tmpPZm = tmpAZm*InitData.beta;
*/

        hAXp[2*Size+iI]+=dtT2*tmpPXp+dtT2*InitData.AGenCoef*tmpAXp
                                    +dtT2*InitData.ABackCoef*tmpAXm
                                    +dtT2*InitData.AGenX*tmpAXp;
        hAXm[2*Size+iI]+=dtT2*tmpPXm+dtT2*InitData.AGenCoef*tmpAXm
                                    +dtT2*InitData.ABackCoef*tmpAXp
                                    +dtT2*InitData.AGenX*tmpAXm;
        hAZp[2*Size+iI]+=dtT2*tmpPZp+dtT2*InitData.AGenCoef*tmpAZp
                                    +dtT2*InitData.ABackCoef*tmpAZm
                                    +dtT2*InitData.AGenZ*tmpAZp;
        hAZm[2*Size+iI]+=dtT2*tmpPZm+dtT2*InitData.AGenCoef*tmpAZm
                                    +dtT2*InitData.ABackCoef*tmpAZp
                                    +dtT2*InitData.AGenZ*tmpAZm;
      }
    }
    MediaRungeStep++;
    break;
    case 2:
      memcpy(hM+IndexXZ(ix1st,0),hM+2*Size+IndexXZ(ix1st,0),
        (InitData.nStepsZ+2)*sizeof(cOptic2D_Media_X)*(ixLast-ix1st+1));
    MediaRungeStep=0;
    break;
  }
}

void cEquOpticBaseTask::StepMediaRunge2FastP() {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  OPTIC_TYPE dtT2=0.5*dtT;
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1);
  complex<OPTIC_TYPE> tmpAXp,tmpAXm,tmpAZp,tmpAZm,tmpPXp,tmpPXm,tmpPZp,tmpPZm;
  complex<OPTIC_TYPE> RCoef2=InitData.RCoef;
  if(InitData.RCoef2>=0) RCoef2=InitData.RCoef2;
  cOptic2D_Media_X tM,tM1;
  int iI;
  cOptic2D_Media_X tmpMedia;
  complex<OPTIC_TYPE> PCoef = complex<OPTIC_TYPE>(InitData.PCoef,InitData.Delta);
  complex<OPTIC_TYPE> beta_PCoef=InitData.beta/PCoef;
  switch(MediaRungeStep) {
    case 0:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        iI=IndexXZ(ix,iz);
        if(InitData.fForceLinear) {
          hM[iI].R0=InitData.Q;
          hM[iI].RX=hM[iI].RZ=0;
          hM[iI].Rp=hM[iI].Rm=InitData.QGrid;
        }
        hM[Size+iI].PXp=tmpPXp=(hAXp[iI]*hM[iI].R0*((OPTIC_TYPE)2.0)+
                hAXm[iI]*hM[iI].RX+hAZp[iI]*conj(hM[iI].Rm)+
                hAZm[iI]*hM[iI].Rp)*beta_PCoef;
        hM[Size+iI].PXm=tmpPXm=(hAXm[iI]*hM[iI].R0*((OPTIC_TYPE)2.0)+
                hAXp[iI]*hM[iI].RX+hAZp[iI]*conj(hM[iI].Rp)+
                hAZm[iI]*hM[iI].Rm)*beta_PCoef;
        hM[Size+iI].PZp=tmpPZp=(hAZp[iI]*hM[iI].R0*((OPTIC_TYPE)2.0)+
                hAZm[iI]*hM[iI].RZ+hAXp[iI]*hM[iI].Rm+
                hAXm[iI]*hM[iI].Rp)*beta_PCoef;
        hM[Size+iI].PZm=tmpPZm=(hAZm[iI]*hM[iI].R0*((OPTIC_TYPE)2.0)+
                hAZp[iI]*conj(hM[iI].RZ)+hAXp[iI]*conj(hM[iI].Rp)+
                hAXm[iI]*conj(hM[iI].Rm))*beta_PCoef;
        hM[Size+iI].RX=(-hM[iI].RX*RCoef2
                        -(hAXp[iI]*conj(tmpPXm)
                        +tmpPXp*conj(hAXm[iI])) )*dtT2;
        hM[Size+iI].RZ=(-hM[iI].RZ*RCoef2
                        -(hAZp[iI]*conj(tmpPZm)
                        +tmpPZp*conj(hAZm[iI])) )*dtT2;
        hM[Size+iI].Rp=(InitData.QGrid-hM[iI].Rp*RCoef2
                        -(tmpPZp*conj(hAXm[iI])
                        +hAZp[iI]*conj(tmpPXm)
                        +tmpPXp*conj(hAZm[iI])
                        +hAXp[iI]*conj(tmpPZm)) )*dtT2;
        hM[Size+iI].Rm=(InitData.QGrid-hM[iI].Rm*RCoef2
                        -(tmpPXm*conj(hAZm[iI])
                        +hAZp[iI]*conj(tmpPXp)
                        +hAXm[iI]*conj(tmpPZm)
                        +tmpPZp*conj(hAXp[iI])) )*dtT2;
        hM[Size+iI].R0=(InitData.Q-InitData.RCoef*hM[iI].R0
                       -(hAZp[iI]*conj(tmpPZp)
                        +hAXp[iI]*conj(tmpPXp)
                        +hAZm[iI]*conj(tmpPZm)
                        +hAXm[iI]*conj(tmpPXm)).real() )*dtT2;
//        hM[Size+iI].Rp=hM[Size+iI].Rm=0;
//        hM[Size+iI].RX=hM[Size+iI].RZ=0;
        hM[Size+iI].R0Generated=(InitData.Q-InitData.RCoef*hM[iI].R0)*dtT2;


        hAXp[Size+iI]+=dtT2*hM[IndexXZ(ix-1,iz)].PXp
                      +dtT2*InitData.AGenCoef*hAXp[IndexXZ(ix-1,iz)]
                      +dtT2*InitData.ABackCoef*hAXm[IndexXZ(ix-1,iz)]*iii
                      +dtT2*InitData.AGenX*hAXp[IndexXZ(ix-1,iz)];
        hAXm[Size+iI]+=dtT2*hM[IndexXZ(ix+1,iz)].PXm
                      +dtT2*InitData.AGenCoef*hAXm[IndexXZ(ix+1,iz)]
                      +dtT2*InitData.ABackCoef*hAXp[IndexXZ(ix+1,iz)]*iii
                      +dtT2*InitData.AGenX*hAXm[IndexXZ(ix+1,iz)];
        hAZp[Size+iI]+=dtT2*hM[IndexXZ(ix,iz-1)].PZp
                      +dtT2*InitData.AGenCoef*hAZp[IndexXZ(ix,iz-1)]
                      +dtT2*InitData.ABackCoef*hAZm[IndexXZ(ix,iz-1)]*iii
                      +dtT2*InitData.AGenZ*hAZp[IndexXZ(ix,iz-1)];
        hAZm[Size+iI]+=dtT2*hM[IndexXZ(ix,iz+1)].PZm
                      +dtT2*InitData.AGenCoef*hAZm[IndexXZ(ix,iz+1)]
                      +dtT2*InitData.ABackCoef*hAZm[IndexXZ(ix,iz+1)]*iii
                      +dtT2*InitData.AGenZ*hAZm[IndexXZ(ix,iz+1)];
      }
    }
    MediaRungeStep++;
    break;
    case 1:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        iI=IndexXZ(ix,iz);
        tmpAXp=hAXp[IndexXZ(ix-1,iz)]+(OPTIC_TYPE)2.0*hAXp[Size+iI];
        tmpAXm=hAXm[IndexXZ(ix+1,iz)]+(OPTIC_TYPE)2.0*hAXm[Size+iI];
        tmpAZp=hAZp[IndexXZ(ix,iz-1)]+(OPTIC_TYPE)2.0*hAZp[Size+iI];
        tmpAZm=hAZm[IndexXZ(ix,iz+1)]+(OPTIC_TYPE)2.0*hAZm[Size+iI];
//        tM=hM[iI]; tM+=hM[Size+iI]*((OPTIC_TYPE)2.0);
//        hM[2*Size+iI]=hM[iI]; hM[2*Size+iI]+=hM[Size+iI];
        tM.M1M2v(hM[iI],hM[Size+iI],2);
        hM[2*Size+iI].M1M2v(hM[iI],hM[Size+iI],1);
//        hM[2*Size+iI]=hM[iI]; hM[2*Size+iI]+=hM[Size+iI];

/*
        hM[2*Size+iI].PXp=tmpPXp=(hAXp[iI]*tM.R0*((OPTIC_TYPE)2.0)+
                hAXm[iI]*tM.RX+hAZp[iI]*conj(tM.Rm)+
                hAZm[iI]*tM.Rp)*InitData.beta;
        hM[2*Size+iI].PXm=tmpPXm=(hAXm[iI]*tM.R0*((OPTIC_TYPE)2.0)+
                hAXp[iI]*tM.RX+hAZp[iI]*conj(tM.Rp)+
                hAZm[iI]*tM.Rm)*InitData.beta;
        hM[2*Size+iI].PZp=tmpPZp=(hAZp[iI]*tM.R0*((OPTIC_TYPE)2.0)+
                hAZm[iI]*tM.RZ+hAXp[iI]*tM.Rm+
                hAXm[iI]*tM.Rp)*InitData.beta;
        hM[2*Size+iI].PZm=tmpPZm=(hAZm[iI]*tM.R0*((OPTIC_TYPE)2.0)+
                hAZp[iI]*conj(tM.RZ)+hAXp[iI]*conj(tM.Rp)+
                hAXm[iI]*conj(tM.Rm))*InitData.beta;
*/
        if(InitData.fForceLinear) {
          tM.R0=InitData.Q;
          tM.RX=tM.RZ=0;
          tM.Rp=tM.Rm=InitData.QGrid;
        }
        hM[2*Size+iI].PXp=tmpPXp=(tmpAXp*tM.R0*((OPTIC_TYPE)2.0)+
                tmpAXm*tM.RX+tmpAZp*conj(tM.Rm)+
                tmpAZm*tM.Rp)*beta_PCoef;
        hM[2*Size+iI].PXm=tmpPXm=(tmpAXm*tM.R0*((OPTIC_TYPE)2.0)+
                tmpAXp*tM.RX+tmpAZp*conj(tM.Rp)+
                tmpAZm*tM.Rm)*beta_PCoef;
        hM[2*Size+iI].PZp=tmpPZp=(tmpAZp*tM.R0*((OPTIC_TYPE)2.0)+
                tmpAZm*tM.RZ+tmpAXp*tM.Rm+
                tmpAXm*tM.Rp)*beta_PCoef;
        hM[2*Size+iI].PZm=tmpPZm=(tmpAZm*tM.R0*((OPTIC_TYPE)2.0)+
                tmpAZp*conj(tM.RZ)+tmpAXp*conj(tM.Rp)+
                tmpAXm*conj(tM.Rm))*beta_PCoef;
        hM[2*Size+iI].RX+=(-tM.RX*RCoef2
                           -(tmpAXp*conj(tmpPXm)
                           +tmpPXp*conj(tmpAXm)) )*dtT;
        hM[2*Size+iI].RZ+=(-tM.RZ*RCoef2
                           -(tmpAZp*conj(tmpPZm)
                           +tmpPZp*conj(tmpAZm)) )*dtT;
        hM[2*Size+iI].Rp+=(InitData.QGrid-tM.Rp*RCoef2
                           -(tmpPZp*conj(tmpAXm)
                           +tmpAZp*conj(tmpPXm)
                           +tmpPXp*conj(tmpAZm)
                           +tmpAXp*conj(tmpPZm)) )*dtT;
        hM[2*Size+iI].Rm+=(InitData.QGrid-tM.Rm*RCoef2
                           -(tmpPXm*conj(tmpAZm)
                           +tmpAZp*conj(tmpPXp)
                           +tmpAXm*conj(tmpPZm)
                           +tmpPZp*conj(tmpAXp)) )*dtT;
        hM[2*Size+iI].R0+=(InitData.Q-InitData.RCoef*tM.R0
                          -(tmpAZp*conj(tmpPZp)
                          +tmpAXp*conj(tmpPXp)
                          +tmpAZm*conj(tmpPZm)
                          +tmpAXm*conj(tmpPXm)).real() )*dtT2;
        hM[2*Size+iI].R0Generated+=(InitData.Q-InitData.RCoef*tM.R0)*dtT2;
//        hM[2*Size+iI].Rp=hM[2*Size+iI].Rm=0;
//        hM[2*Size+iI].RX=hM[2*Size+iI].RZ=0;
        hAXp[2*Size+iI]+=dtT2*tmpPXp+dtT2*InitData.AGenCoef*tmpAXp
                                    +dtT2*InitData.ABackCoef*tmpAXm*iii
                                    +dtT2*InitData.AGenX*tmpAXp;
        hAXm[2*Size+iI]+=dtT2*tmpPXm+dtT2*InitData.AGenCoef*tmpAXm
                                    +dtT2*InitData.ABackCoef*tmpAXp*iii
                                    +dtT2*InitData.AGenX*tmpAXm;
        hAZp[2*Size+iI]+=dtT2*tmpPZp+dtT2*InitData.AGenCoef*tmpAZp
                                    +dtT2*InitData.ABackCoef*tmpAZm*iii
                                    +dtT2*InitData.AGenZ*tmpAZp;
        hAZm[2*Size+iI]+=dtT2*tmpPZm+dtT2*InitData.AGenCoef*tmpAZm
                                    +dtT2*InitData.ABackCoef*tmpAZp*iii
                                    +dtT2*InitData.AGenZ*tmpAZm;
      }
    }
    MediaRungeStep++;
    break;
    case 2:
      //if(!InitData.fForceLinear)
        memcpy(hM+IndexXZ(ix1st,0),hM+2*Size+IndexXZ(ix1st,0),
             (InitData.nStepsZ+2)*sizeof(cOptic2D_Media_X)*(ixLast-ix1st+1));
    MediaRungeStep=0;
    break;
  }
}

void cEquOpticBaseTask::StepMediaRunge2FastPTE() {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  OPTIC_TYPE dtT2=0.5*dtT;
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1);
  complex<OPTIC_TYPE> tmpAXp,tmpAXm,tmpAZp,tmpAZm,tmpPXp,tmpPXm,tmpPZp,tmpPZm;
  cOptic2D_Media_X tM,tM1;
  int iI;
  cOptic2D_Media_X tmpMedia;
  complex<OPTIC_TYPE> RCoef2=InitData.RCoef;
  if(InitData.RCoef2>=0) RCoef2=InitData.RCoef2;
  switch(MediaRungeStep) {
    case 0:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        iI=IndexXZ(ix,iz);
        hM[Size+iI].PXp=tmpPXp=0;
        hM[Size+iI].PXm=tmpPXm=0;

        hM[Size+iI].PZp=tmpPZp=(hAZp[iI]*hM[iI].R0*((OPTIC_TYPE)2.0)+
                hAZm[iI]*hM[iI].RZ)*InitData.beta;
        hM[Size+iI].PZm=tmpPZm=(hAZm[iI]*hM[iI].R0*((OPTIC_TYPE)2.0)+
                hAZp[iI]*conj(hM[iI].RZ))*InitData.beta;

        hM[Size+iI].RZ=(-hM[iI].RZ*RCoef2
                        -(hAZp[iI]*conj(tmpPZm)
                        +tmpPZp*conj(hAZm[iI])) )*dtT2;

        hM[Size+iI].RX = hM[Size+iI].Rp = hM[Size+iI].Rm = 0;
        hM[Size+iI].R0=(InitData.Q-InitData.RCoef*hM[iI].R0
                       -(hAZp[iI]*conj(tmpPZp)
                        +hAXp[iI]*conj(tmpPXp)
                        +hAZm[iI]*conj(tmpPZm)
                        +hAXm[iI]*conj(tmpPXm)).real() )*dtT2;
//        hM[Size+iI].Rp=hM[Size+iI].Rm=0;
//        hM[Size+iI].RX=hM[Size+iI].RZ=0;
        hM[Size+iI].R0Generated=(InitData.Q-InitData.RCoef*hM[iI].R0)*dtT2;


        hAXp[Size+iI]+=dtT2*InitData.AGenCoef*hAXp[IndexXZ(ix-1,iz)]
                      +dtT2*InitData.ABackCoef*hAXm[IndexXZ(ix-1,iz)]
                      +dtT2*InitData.AGenX*hAXp[IndexXZ(ix-1,iz)];
        hAXm[Size+iI]+=dtT2*InitData.AGenCoef*hAXm[IndexXZ(ix+1,iz)]
                      +dtT2*InitData.ABackCoef*hAXp[IndexXZ(ix+1,iz)]
                      +dtT2*InitData.AGenX*hAXm[IndexXZ(ix+1,iz)];
        hAZp[Size+iI]+=dtT2*hM[IndexXZ(ix,iz-1)].PZp
                      +dtT2*InitData.AGenCoef*hAZp[IndexXZ(ix,iz-1)]
                      +dtT2*InitData.ABackCoef*hAZm[IndexXZ(ix,iz-1)]
                      +dtT2*InitData.AGenZ*hAZp[IndexXZ(ix,iz-1)];
        hAZm[Size+iI]+=dtT2*hM[IndexXZ(ix,iz+1)].PZm
                      +dtT2*InitData.AGenCoef*hAZm[IndexXZ(ix,iz+1)]
                      +dtT2*InitData.ABackCoef*hAZm[IndexXZ(ix,iz+1)]
                      +dtT2*InitData.AGenZ*hAZm[IndexXZ(ix,iz+1)];
      }
    }
    MediaRungeStep++;
    break;
    case 1:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        iI=IndexXZ(ix,iz);
        tmpAXp=hAXp[IndexXZ(ix-1,iz)]+(OPTIC_TYPE)2.0*hAXp[Size+iI];
        tmpAXm=hAXm[IndexXZ(ix+1,iz)]+(OPTIC_TYPE)2.0*hAXm[Size+iI];
        tmpAZp=hAZp[IndexXZ(ix,iz-1)]+(OPTIC_TYPE)2.0*hAZp[Size+iI];
        tmpAZm=hAZm[IndexXZ(ix,iz+1)]+(OPTIC_TYPE)2.0*hAZm[Size+iI];
//        tM=hM[iI]; tM+=hM[Size+iI]*((OPTIC_TYPE)2.0);
//        hM[2*Size+iI]=hM[iI]; hM[2*Size+iI]+=hM[Size+iI];
        tM.M1M2v(hM[iI],hM[Size+iI],2);
        hM[2*Size+iI].M1M2v(hM[iI],hM[Size+iI],1);
//        hM[2*Size+iI]=hM[iI]; hM[2*Size+iI]+=hM[Size+iI];

        hM[2*Size+iI].PXp=tmpPXp=0;
        hM[2*Size+iI].PXm=tmpPXm=0;

        hM[2*Size+iI].PZp=tmpPZp=(tmpAZp*tM.R0*((OPTIC_TYPE)2.0)+
                tmpAZm*tM.RZ)*InitData.beta;
        hM[2*Size+iI].PZm=tmpPZm=(tmpAZm*tM.R0*((OPTIC_TYPE)2.0)+
                tmpAZp*conj(tM.RZ))*InitData.beta;


        hM[2*Size+iI].RX = hM[2*Size+iI].Rp = hM[2*Size+iI].Rm = 0;
        hM[2*Size+iI].RZ+=(-tM.RZ*RCoef2
                           -(tmpAZp*conj(tmpPZm)
                           +tmpPZp*conj(tmpAZm)) )*dtT;
        hM[2*Size+iI].R0+=(InitData.Q-InitData.RCoef*tM.R0
                          -(tmpAZp*conj(tmpPZp)
                          +tmpAZm*conj(tmpPZm)).real() )*dtT2;

        hM[2*Size+iI].R0Generated+=(InitData.Q-InitData.RCoef*tM.R0)*dtT2;
//        hM[2*Size+iI].Rp=hM[2*Size+iI].Rm=0;
//        hM[2*Size+iI].RX=hM[2*Size+iI].RZ=0;

        hAXp[2*Size+iI]+=dtT2*InitData.AGenCoef*tmpAXp
                                    +dtT2*InitData.ABackCoef*tmpAXm
                                    +dtT2*InitData.AGenX*tmpAXp;
        hAXm[2*Size+iI]+=dtT2*InitData.AGenCoef*tmpAXm
                                    +dtT2*InitData.ABackCoef*tmpAXp
                                    +dtT2*InitData.AGenX*tmpAXm;
        hAZp[2*Size+iI]+=dtT2*tmpPZp+dtT2*InitData.AGenCoef*tmpAZp
                                    +dtT2*InitData.ABackCoef*tmpAZm
                                    +dtT2*InitData.AGenZ*tmpAZp;
        hAZm[2*Size+iI]+=dtT2*tmpPZm+dtT2*InitData.AGenCoef*tmpAZm
                                    +dtT2*InitData.ABackCoef*tmpAZp
                                    +dtT2*InitData.AGenZ*tmpAZm;
      }
    }
    MediaRungeStep++;
    break;
    case 2:
      memcpy(hM+IndexXZ(ix1st,0),hM+2*Size+IndexXZ(ix1st,0),
        (InitData.nStepsZ+2)*sizeof(cOptic2D_Media_X)*(ixLast-ix1st+1));
    MediaRungeStep=0;
    break;
  }
}


void cEquOpticBaseTask::StepNoMediaMatrix() {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1),tmpAXp,tmpAXm,tmpAZp,tmpAZm;
//  complex<OPTIC_TYPE> AXpk1,AXpk2,AXmk1,AXmk2,AZpk1,AZpk2,AZmk1,AZmk2;
//  complex<OPTIC_TYPE> AXpk2Arr[100],AXpk1Arr[100],AXpk12Arr[100];
  OPTIC_TYPE dtT2=0.5*dtT;
  int iznb1=1,iznb2=InitData.nStepsZ,ixnb1=ix1st,ixnb2=ixLast;
  switch(RungeStep) {
    case 0:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        hAXp[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*
                  (hAZp[IndexXZ(ix-1,iz)]+hAZm[IndexXZ(ix-1,iz)]);
        hAXm[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*
                  (hAZp[IndexXZ(ix+1,iz)]+hAZm[IndexXZ(ix+1,iz)]);
        hAZp[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*
                  (hAXp[IndexXZ(ix,iz-1)]+hAXm[IndexXZ(ix,iz-1)]);
        hAZm[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*
                  (hAXp[IndexXZ(ix,iz+1)]+hAXm[IndexXZ(ix,iz+1)]);
      }
    }
    RungeStep++;
    break;
    case 1:
    iznb1=2; iznb2=InitData.nStepsZ-1;
    ixnb1=ix1st; ixnb2=ixLast;

    if( (!InitData.LockX)&&(ix1st==1) ) {
      ixnb1=2;
      for(int iz=iznb1;iz<=iznb2;iz++) {
        tmpAXm=hAXm[IndexXZ(2,iz)]+hAXm[Size+IndexXZ(1,iz)];
        tmpAZp=hAZp[IndexXZ(1,iz-1)]+hAZp[Size+IndexXZ(1,iz)];
        tmpAZm=hAZm[IndexXZ(1,iz+1)]+hAZm[Size+IndexXZ(1,iz)];
        hAXp[2*Size+IndexXZ(1,iz)]=hAXp[IndexXZ(0,iz)];
        hAXm[2*Size+IndexXZ(1,iz)]=tmpAZp*Matrix3[6]+tmpAZm*Matrix3[7]+
                                   tmpAXm*Matrix3[8];
        hAZp[2*Size+IndexXZ(1,iz)]=tmpAZp*Matrix3[0]+tmpAZm*Matrix3[1]+
                                   tmpAXm*Matrix3[2];
        hAZm[2*Size+IndexXZ(1,iz)]=tmpAZp*Matrix3[3]+tmpAZm*Matrix3[4]+
                                   tmpAXm*Matrix3[5];
      }
      hAXp[2*Size+IndexXZ(1,1)]=hAXp[IndexXZ(0,1)];
      hAZp[2*Size+IndexXZ(1,1)]=hAZp[IndexXZ(1,0)];
      tmpAXm=hAXm[IndexXZ(2,1)]+hAXm[Size+IndexXZ(1,1)];
      tmpAZm=hAZm[IndexXZ(1,2)]+hAZm[Size+IndexXZ(1,1)];
      hAXm[2*Size+IndexXZ(1,1)]=tmpAZm*Matrix2[1]+tmpAXm*Matrix2[0];
      hAZm[2*Size+IndexXZ(1,1)]=tmpAZm*Matrix2[0]+tmpAXm*Matrix2[1];

      hAXp[2*Size+IndexXZ(1,InitData.nStepsZ)]=hAXp[IndexXZ(0,InitData.nStepsZ)];
      hAZm[2*Size+IndexXZ(1,InitData.nStepsZ)]=hAZm[IndexXZ(1,InitData.nStepsZ+1)];
      tmpAXm=hAXm[IndexXZ(2,InitData.nStepsZ)]+hAXm[Size+IndexXZ(1,InitData.nStepsZ)];
      tmpAZp=hAZp[IndexXZ(1,InitData.nStepsZ-1)]+hAZp[Size+IndexXZ(1,InitData.nStepsZ)];
      hAXm[2*Size+IndexXZ(1,InitData.nStepsZ)]=tmpAZp*Matrix2[1]+tmpAXm*Matrix2[0];
      hAZp[2*Size+IndexXZ(1,InitData.nStepsZ)]=tmpAZp*Matrix2[0]+tmpAXm*Matrix2[1];
    }


    if( (!InitData.LockX)&&(ixLast==InitData.nStepsX) ) {
      ixnb2=InitData.nStepsX-1;
      for(int iz=iznb1;iz<=iznb2;iz++) {
        tmpAXp=hAXp[IndexXZ(InitData.nStepsX-1,iz)]+hAXp[Size+IndexXZ(InitData.nStepsX,iz)];
        tmpAZp=hAZp[IndexXZ(InitData.nStepsX,iz-1)]+hAZp[Size+IndexXZ(InitData.nStepsX,iz)];
        tmpAZm=hAZm[IndexXZ(InitData.nStepsX,iz+1)]+hAZm[Size+IndexXZ(InitData.nStepsX,iz)];
        hAXm[2*Size+IndexXZ(InitData.nStepsX,iz)]=hAXm[IndexXZ(InitData.nStepsX+1,iz)];
        hAXp[2*Size+IndexXZ(InitData.nStepsX,iz)]=tmpAZp*Matrix3[6]+tmpAZm*Matrix3[7]+
                                   tmpAXp*Matrix3[8];
        hAZp[2*Size+IndexXZ(InitData.nStepsX,iz)]=tmpAZp*Matrix3[0]+tmpAZm*Matrix3[1]+
                                   tmpAXp*Matrix3[2];
        hAZm[2*Size+IndexXZ(InitData.nStepsX,iz)]=tmpAZp*Matrix3[3]+tmpAZm*Matrix3[4]+
                                   tmpAXp*Matrix3[5];
      }
      hAXm[2*Size+IndexXZ(InitData.nStepsX,1)]=hAXm[IndexXZ(InitData.nStepsX+1,1)];
      hAZp[2*Size+IndexXZ(InitData.nStepsX,1)]=hAZp[IndexXZ(InitData.nStepsX,0)];
      tmpAXp=hAXp[IndexXZ(InitData.nStepsX-1,1)]+hAXp[Size+IndexXZ(InitData.nStepsX,1)];
      tmpAZm=hAZm[IndexXZ(InitData.nStepsX,2)]+hAZm[Size+IndexXZ(InitData.nStepsX,1)];
      hAXp[2*Size+IndexXZ(InitData.nStepsX,1)]=tmpAZm*Matrix2[1]+tmpAXp*Matrix2[0];
      hAZm[2*Size+IndexXZ(InitData.nStepsX,1)]=tmpAZm*Matrix2[0]+tmpAXp*Matrix2[1];

      hAXm[2*Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)]=hAXm[IndexXZ(InitData.nStepsX+1,InitData.nStepsZ)];
      hAZm[2*Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)]=hAZm[IndexXZ(InitData.nStepsX,InitData.nStepsZ+1)];
      tmpAXp=hAXp[IndexXZ(InitData.nStepsX-1,InitData.nStepsZ)]+hAXp[Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)];
      tmpAZp=hAZp[IndexXZ(InitData.nStepsX,InitData.nStepsZ-1)]+hAZp[Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)];
      hAXp[2*Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)]=tmpAZp*Matrix2[1]+tmpAXp*Matrix2[0];
      hAZp[2*Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)]=tmpAZp*Matrix2[0]+tmpAXp*Matrix2[1];
    }

    for(int ix=ixnb1;ix<=ixnb2;ix++) {
      tmpAXp=hAXp[IndexXZ(ix-1,1)]+hAXp[Size+IndexXZ(ix,1)];
      tmpAXm=hAXm[IndexXZ(ix+1,1)]+hAXm[Size+IndexXZ(ix,1)];
      tmpAZm=hAZm[IndexXZ(ix,2)]+hAZm[Size+IndexXZ(ix,1)];
      hAXp[2*Size+IndexXZ(ix,1)]=tmpAXp*Matrix3[0]+tmpAXm*Matrix3[1]+
                                  tmpAZm*Matrix3[2];
      hAXm[2*Size+IndexXZ(ix,1)]=tmpAXp*Matrix3[3]+tmpAXm*Matrix3[4]+
                                  tmpAZm*Matrix3[5];
      hAZm[2*Size+IndexXZ(ix,1)]=tmpAXp*Matrix3[6]+tmpAXm*Matrix3[7]+
                                  tmpAZm*Matrix3[8];
      hAZp[2*Size+IndexXZ(ix,1)]=hAZp[IndexXZ(ix,0)];
      tmpAXp=hAXp[IndexXZ(ix-1,InitData.nStepsZ)]+hAXp[Size+IndexXZ(ix,InitData.nStepsZ)];
      tmpAXm=hAXm[IndexXZ(ix+1,InitData.nStepsZ)]+hAXm[Size+IndexXZ(ix,InitData.nStepsZ)];
      tmpAZp=hAZp[IndexXZ(ix,InitData.nStepsZ-1)]+hAZp[Size+IndexXZ(ix,InitData.nStepsZ)];
      hAXp[2*Size+IndexXZ(ix,InitData.nStepsZ)]=tmpAXp*Matrix3[0]+tmpAXm*Matrix3[1]+
                                  tmpAZp*Matrix3[2];
      hAXm[2*Size+IndexXZ(ix,InitData.nStepsZ)]=tmpAXp*Matrix3[3]+tmpAXm*Matrix3[4]+
                                  tmpAZp*Matrix3[5];
      hAZp[2*Size+IndexXZ(ix,InitData.nStepsZ)]=tmpAXp*Matrix3[6]+tmpAXm*Matrix3[7]+
                                  tmpAZp*Matrix3[8];
      hAZm[2*Size+IndexXZ(ix,InitData.nStepsZ)]=hAZm[IndexXZ(ix,InitData.nStepsZ+1)];
    }

    if( (!InitData.LockX)&&(InitData.Z0RefCoef) ) { iznb1=1; iznb2=InitData.nStepsZ; ixnb1=ix1st; ixnb2=ixLast; }
              
    for(int ix=ixnb1;ix<=ixnb2;ix++) {
      for(int iz=iznb1;iz<=iznb2;iz++) {
        tmpAXp=hAXp[IndexXZ(ix-1,iz)]+hAXp[Size+IndexXZ(ix,iz)];
        tmpAXm=hAXm[IndexXZ(ix+1,iz)]+hAXm[Size+IndexXZ(ix,iz)];
        tmpAZp=hAZp[IndexXZ(ix,iz-1)]+hAZp[Size+IndexXZ(ix,iz)];
        tmpAZm=hAZm[IndexXZ(ix,iz+1)]+hAZm[Size+IndexXZ(ix,iz)];
        hAXp[2*Size+IndexXZ(ix,iz)]=tmpAXp*Matrix4[0]+tmpAXm*Matrix4[1]+
                                    tmpAZp*Matrix4[2]+tmpAZm*Matrix4[3];
        hAXm[2*Size+IndexXZ(ix,iz)]=tmpAXp*Matrix4[4]+tmpAXm*Matrix4[5]+
                                    tmpAZp*Matrix4[6]+tmpAZm*Matrix4[7];
        hAZp[2*Size+IndexXZ(ix,iz)]=tmpAXp*Matrix4[8]+tmpAXm*Matrix4[9]+
                                    tmpAZp*Matrix4[10]+tmpAZm*Matrix4[11];
        hAZm[2*Size+IndexXZ(ix,iz)]=tmpAXp*Matrix4[12]+tmpAXm*Matrix4[13]+
                                    tmpAZp*Matrix4[14]+tmpAZm*Matrix4[15];
      }
    }
    RungeStep++;
    break;

    case 2:
    memcpy(hAXp+IndexXZ(ix1st,0),hAXp+2*Size+IndexXZ(ix1st,0),
      (InitData.nStepsZ+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    memcpy(hAXm+IndexXZ(ix1st,0),hAXm+2*Size+IndexXZ(ix1st,0),
      (InitData.nStepsZ+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    memcpy(hAZp+IndexXZ(ix1st,0),hAZp+2*Size+IndexXZ(ix1st,0),
      (InitData.nStepsZ+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    memcpy(hAZm+IndexXZ(ix1st,0),hAZm+2*Size+IndexXZ(ix1st,0),
      (InitData.nStepsZ+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    RungeStep=0;
  }
}

void cEquOpticBaseTask::StepNoMediaMatrixMMode() {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1),tmpAXp,tmpAXm,tmpAZp,tmpAZm;
//  complex<OPTIC_TYPE> AXpk1,AXpk2,AXmk1,AXmk2,AZpk1,AZpk2,AZmk1,AZmk2;
//  complex<OPTIC_TYPE> AXpk2Arr[100],AXpk1Arr[100],AXpk12Arr[100];
  OPTIC_TYPE dtT2=0.5*dtT;
  int iznb1=1,iznb2=InitData.nStepsZ,ixnb1=ix1st,ixnb2=ixLast;
  switch(RungeStep) {
    case 0:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        hAXp[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*
                  (hAZp[IndexXZ(ix-1,iz)]+hAZm[IndexXZ(ix-1,iz)]);
        hAXm[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*
                  (hAZp[IndexXZ(ix+1,iz)]+hAZm[IndexXZ(ix+1,iz)]);
        hAZp[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*
                  (hAXp[IndexXZ(ix,iz-1)]+hAXm[IndexXZ(ix,iz-1)]);
        hAZm[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*
                  (hAXp[IndexXZ(ix,iz+1)]+hAXm[IndexXZ(ix,iz+1)]);
      }
    }
    RungeStep++;
    break;
    case 1:
    iznb1=1; iznb2=InitData.nStepsZ-1;
    ixnb1=ix1st; ixnb2=ixLast;

    if( (!InitData.LockX)&&(ix1st==1) ) {
      ixnb1=2;
      for(int iz=iznb1;iz<=iznb2;iz++) {
        tmpAXm=hAXm[IndexXZ(2,iz)]+hAXm[Size+IndexXZ(1,iz)];
        tmpAZp=hAZp[IndexXZ(1,iz-1)]+hAZp[Size+IndexXZ(1,iz)];
        tmpAZm=hAZm[IndexXZ(1,iz+1)]+hAZm[Size+IndexXZ(1,iz)];
        hAXp[2*Size+IndexXZ(1,iz)]=hAXp[IndexXZ(0,iz)];
        hAXm[2*Size+IndexXZ(1,iz)]=tmpAZp*Matrix3[6]+tmpAZm*Matrix3[7]+
                                   tmpAXm*Matrix3[8];
        hAZp[2*Size+IndexXZ(1,iz)]=tmpAZp*Matrix3[0]+tmpAZm*Matrix3[1]+
                                   tmpAXm*Matrix3[2];
        hAZm[2*Size+IndexXZ(1,iz)]=tmpAZp*Matrix3[3]+tmpAZm*Matrix3[4]+
                                   tmpAXm*Matrix3[5];
      }
      hAXp[2*Size+IndexXZ(1,1)]=hAXp[IndexXZ(0,1)];
      hAZp[2*Size+IndexXZ(1,1)]=hAZp[IndexXZ(1,0)];
      tmpAXm=hAXm[IndexXZ(2,1)]+hAXm[Size+IndexXZ(1,1)];
      tmpAZm=hAZm[IndexXZ(1,2)]+hAZm[Size+IndexXZ(1,1)];
      hAXm[2*Size+IndexXZ(1,1)]=tmpAZm*Matrix2[1]+tmpAXm*Matrix2[0];
      hAZm[2*Size+IndexXZ(1,1)]=tmpAZm*Matrix2[0]+tmpAXm*Matrix2[1];

      hAXp[2*Size+IndexXZ(1,InitData.nStepsZ)]=hAXp[IndexXZ(0,InitData.nStepsZ)];
      hAZm[2*Size+IndexXZ(1,InitData.nStepsZ)]=hAZm[IndexXZ(1,InitData.nStepsZ+1)];
      tmpAXm=hAXm[IndexXZ(2,InitData.nStepsZ)]+hAXm[Size+IndexXZ(1,InitData.nStepsZ)];
      tmpAZp=hAZp[IndexXZ(1,InitData.nStepsZ-1)]+hAZp[Size+IndexXZ(1,InitData.nStepsZ)];
      hAXm[2*Size+IndexXZ(1,InitData.nStepsZ)]=tmpAZp*Matrix2[1]+tmpAXm*Matrix2[0];
      hAZp[2*Size+IndexXZ(1,InitData.nStepsZ)]=tmpAZp*Matrix2[0]+tmpAXm*Matrix2[1];
    }


    if( (!InitData.LockX)&&(ixLast==InitData.nStepsX) ) {
      ixnb2=InitData.nStepsX-1;
      for(int iz=iznb1;iz<=iznb2;iz++) {
        tmpAXp=hAXp[IndexXZ(InitData.nStepsX-1,iz)]+hAXp[Size+IndexXZ(InitData.nStepsX,iz)];
        tmpAZp=hAZp[IndexXZ(InitData.nStepsX,iz-1)]+hAZp[Size+IndexXZ(InitData.nStepsX,iz)];
        tmpAZm=hAZm[IndexXZ(InitData.nStepsX,iz+1)]+hAZm[Size+IndexXZ(InitData.nStepsX,iz)];
        hAXm[2*Size+IndexXZ(InitData.nStepsX,iz)]=hAXm[IndexXZ(InitData.nStepsX+1,iz)];
        hAXp[2*Size+IndexXZ(InitData.nStepsX,iz)]=tmpAZp*Matrix3[6]+tmpAZm*Matrix3[7]+
                                   tmpAXp*Matrix3[8];
        hAZp[2*Size+IndexXZ(InitData.nStepsX,iz)]=tmpAZp*Matrix3[0]+tmpAZm*Matrix3[1]+
                                   tmpAXp*Matrix3[2];
        hAZm[2*Size+IndexXZ(InitData.nStepsX,iz)]=tmpAZp*Matrix3[3]+tmpAZm*Matrix3[4]+
                                   tmpAXp*Matrix3[5];
      }
      hAXm[2*Size+IndexXZ(InitData.nStepsX,1)]=hAXm[IndexXZ(InitData.nStepsX+1,1)];
      hAZp[2*Size+IndexXZ(InitData.nStepsX,1)]=hAZp[IndexXZ(InitData.nStepsX,0)];
      tmpAXp=hAXp[IndexXZ(InitData.nStepsX-1,1)]+hAXp[Size+IndexXZ(InitData.nStepsX,1)];
      tmpAZm=hAZm[IndexXZ(InitData.nStepsX,2)]+hAZm[Size+IndexXZ(InitData.nStepsX,1)];
      hAXp[2*Size+IndexXZ(InitData.nStepsX,1)]=tmpAZm*Matrix2[1]+tmpAXp*Matrix2[0];
      hAZm[2*Size+IndexXZ(InitData.nStepsX,1)]=tmpAZm*Matrix2[0]+tmpAXp*Matrix2[1];

      hAXm[2*Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)]=hAXm[IndexXZ(InitData.nStepsX+1,InitData.nStepsZ)];
      hAZm[2*Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)]=hAZm[IndexXZ(InitData.nStepsX,InitData.nStepsZ+1)];
      tmpAXp=hAXp[IndexXZ(InitData.nStepsX-1,InitData.nStepsZ)]+hAXp[Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)];
      tmpAZp=hAZp[IndexXZ(InitData.nStepsX,InitData.nStepsZ-1)]+hAZp[Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)];
      hAXp[2*Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)]=tmpAZp*Matrix2[1]+tmpAXp*Matrix2[0];
      hAZp[2*Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)]=tmpAZp*Matrix2[0]+tmpAXp*Matrix2[1];
    }

    for(int ix=ixnb1;ix<=ixnb2;ix++) {
      tmpAXp=hAXp[IndexXZ(ix-1,1)]+hAXp[Size+IndexXZ(ix,1)];
      tmpAXm=hAXm[IndexXZ(ix+1,1)]+hAXm[Size+IndexXZ(ix,1)];
      tmpAZm=hAZm[IndexXZ(ix,2)]+hAZm[Size+IndexXZ(ix,1)];
      hAXp[2*Size+IndexXZ(ix,1)]=tmpAXp*Matrix3[0]+tmpAXm*Matrix3[1]+
                                  tmpAZm*Matrix3[2];
      hAXm[2*Size+IndexXZ(ix,1)]=tmpAXp*Matrix3[3]+tmpAXm*Matrix3[4]+
                                  tmpAZm*Matrix3[5];
      hAZm[2*Size+IndexXZ(ix,1)]=tmpAXp*Matrix3[6]+tmpAXm*Matrix3[7]+
                                  tmpAZm*Matrix3[8];
      hAZp[2*Size+IndexXZ(ix,1)]=hAZp[IndexXZ(ix,0)];
      tmpAXp=hAXp[IndexXZ(ix-1,InitData.nStepsZ)]+hAXp[Size+IndexXZ(ix,InitData.nStepsZ)];
      tmpAXm=hAXm[IndexXZ(ix+1,InitData.nStepsZ)]+hAXm[Size+IndexXZ(ix,InitData.nStepsZ)];
      tmpAZp=hAZp[IndexXZ(ix,InitData.nStepsZ-1)]+hAZp[Size+IndexXZ(ix,InitData.nStepsZ)];
      hAXp[2*Size+IndexXZ(ix,InitData.nStepsZ)]=tmpAXp*Matrix3[0]+tmpAXm*Matrix3[1]+
                                  tmpAZp*Matrix3[2];
      hAXm[2*Size+IndexXZ(ix,InitData.nStepsZ)]=tmpAXp*Matrix3[3]+tmpAXm*Matrix3[4]+
                                  tmpAZp*Matrix3[5];
      hAZp[2*Size+IndexXZ(ix,InitData.nStepsZ)]=tmpAXp*Matrix3[6]+tmpAXm*Matrix3[7]+
                                  tmpAZp*Matrix3[8];
      hAZm[2*Size+IndexXZ(ix,InitData.nStepsZ)]=hAZm[IndexXZ(ix,InitData.nStepsZ+1)];
    }
              
    for(int ix=ixnb1;ix<=ixnb2;ix++) {
      for(int iz=iznb1;iz<=iznb2;iz++) {
        tmpAXp=hAXp[IndexXZ(ix-1,iz)]+hAXp[Size+IndexXZ(ix,iz)];
        tmpAXm=hAXm[IndexXZ(ix+1,iz)]+hAXm[Size+IndexXZ(ix,iz)];
        tmpAZp=hAZp[IndexXZ(ix,iz-1)]+hAZp[Size+IndexXZ(ix,iz)];
        tmpAZm=hAZm[IndexXZ(ix,iz+1)]+hAZm[Size+IndexXZ(ix,iz)];
        hAXp[2*Size+IndexXZ(ix,iz)]=tmpAXp*Matrix4[0]+tmpAXm*Matrix4[1]+
                                    tmpAZp*Matrix4[2]+tmpAZm*Matrix4[3];
        hAXm[2*Size+IndexXZ(ix,iz)]=tmpAXp*Matrix4[4]+tmpAXm*Matrix4[5]+
                                    tmpAZp*Matrix4[6]+tmpAZm*Matrix4[7];
        hAZp[2*Size+IndexXZ(ix,iz)]=tmpAXp*Matrix4[8]+tmpAXm*Matrix4[9]+
                                    tmpAZp*Matrix4[10]+tmpAZm*Matrix4[11];
        hAZm[2*Size+IndexXZ(ix,iz)]=tmpAXp*Matrix4[12]+tmpAXm*Matrix4[13]+
                                    tmpAZp*Matrix4[14]+tmpAZm*Matrix4[15];
      }
    }
    RungeStep++;
    break;

    case 2:
    memcpy(hAXp+IndexXZ(ix1st,0),hAXp+2*Size+IndexXZ(ix1st,0),
      (InitData.nStepsZ+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    memcpy(hAXm+IndexXZ(ix1st,0),hAXm+2*Size+IndexXZ(ix1st,0),
      (InitData.nStepsZ+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    memcpy(hAZp+IndexXZ(ix1st,0),hAZp+2*Size+IndexXZ(ix1st,0),
      (InitData.nStepsZ+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    memcpy(hAZm+IndexXZ(ix1st,0),hAZm+2*Size+IndexXZ(ix1st,0),
      (InitData.nStepsZ+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    RungeStep=0;
  }
}

void cEquOpticBaseTask::StepNoMediaEuler() {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1),tmpAXp,tmpAXm,tmpAZp,tmpAZm;
  OPTIC_TYPE dtT2=0.5*dtT;
  if(RungeStep==0) {
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {

        tmpAXp=hAXp[IndexXZ(ix-1,iz)];
        tmpAXm=hAXm[IndexXZ(ix+1,iz)];
        tmpAZp=hAZp[IndexXZ(ix,iz-1)];
        tmpAZm=hAZm[IndexXZ(ix,iz+1)];

        hAXp[Size+IndexXZ(ix,iz)]=tmpAXp+dtT*iii*InitData.ACoef*(tmpAZp+tmpAZm);
        hAXm[Size+IndexXZ(ix,iz)]=tmpAXm+dtT*iii*InitData.ACoef*(tmpAZp+tmpAZm);
        hAZp[Size+IndexXZ(ix,iz)]=tmpAZp+dtT*iii*InitData.ACoef*(tmpAXp+tmpAXm);
        hAZm[Size+IndexXZ(ix,iz)]=tmpAZm+dtT*iii*InitData.ACoef*(tmpAXp+tmpAXm);
/*
        hAXp[Size+IndexXZ(ix,iz)]=tmpAXp+dtT*iii*InitData.ACoef*
                  (hAZp[IndexXZ(ix-1,iz-1)]+hAZm[IndexXZ(ix-1,iz+1)]);
        hAXm[Size+IndexXZ(ix,iz)]=tmpAXm+dtT*iii*InitData.ACoef*
                  (hAZp[IndexXZ(ix+1,iz-1)]+hAZm[IndexXZ(ix+1,iz+1)]);
        hAZp[Size+IndexXZ(ix,iz)]=tmpAZp+dtT*iii*InitData.ACoef*
                  (hAXp[IndexXZ(ix-1,iz-1)]+hAXm[IndexXZ(ix+1,iz-1)]);
        hAZm[Size+IndexXZ(ix,iz)]=tmpAZm+dtT*iii*InitData.ACoef*
                  (hAXp[IndexXZ(ix-1,iz+1)]+hAXm[IndexXZ(ix+1,iz+1)]);
*/
      }
    }
    RungeStep++;
  } else {
    memcpy(hAXp+IndexXZ(ix1st,0),hAXp+Size+IndexXZ(ix1st,0),
      (InitData.nStepsZ+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    memcpy(hAXm+IndexXZ(ix1st,0),hAXm+Size+IndexXZ(ix1st,0),
      (InitData.nStepsZ+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    memcpy(hAZp+IndexXZ(ix1st,0),hAZp+Size+IndexXZ(ix1st,0),
      (InitData.nStepsZ+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    memcpy(hAZm+IndexXZ(ix1st,0),hAZm+Size+IndexXZ(ix1st,0),
      (InitData.nStepsZ+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    RungeStep=0;
  }
}

void cEquOpticBaseTask::StepNoMediaRunge2() {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1),tmpAXp,tmpAXm,tmpAZp,tmpAZm;
  complex<OPTIC_TYPE> AXpk1,AXpk2,AXmk1,AXmk2,AZpk1,AZpk2,AZmk1,AZmk2;
  complex<OPTIC_TYPE> AXpk2Arr[100],AXpk1Arr[100],AXpk12Arr[100];
  OPTIC_TYPE dtT2=0.5*dtT;
//  complex<OPTIC_TYPE> AZp10[100],AZp10k1[100],AZp10k2[100],AZp10R[100];
  switch(RungeStep) {
    case 0:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        /*
        tmpAXp=hAXp[IndexXZ(ix-1,iz)];
        tmpAXm=hAXm[IndexXZ(ix+1,iz)];
        tmpAZp=hAZp[IndexXZ(ix,iz-1)];
        tmpAZm=hAZm[IndexXZ(ix,iz+1)];
        hAXp[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*(tmpAZp+tmpAZm);
        hAXm[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*(tmpAZp+tmpAZm);
        hAZp[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*(tmpAXp+tmpAXm);
        hAZm[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*(tmpAXp+tmpAXm);
        */
        hAXp[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*
                  (hAZp[IndexXZ(ix-1,iz)]+hAZm[IndexXZ(ix-1,iz)]);
        hAXm[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*
                  (hAZp[IndexXZ(ix+1,iz)]+hAZm[IndexXZ(ix+1,iz)]);
        hAZp[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*
                  (hAXp[IndexXZ(ix,iz-1)]+hAXm[IndexXZ(ix,iz-1)]);
        hAZm[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*
                  (hAXp[IndexXZ(ix,iz+1)]+hAXm[IndexXZ(ix,iz+1)]);
      }
    }
    RungeStep++;
    break;
    case 1:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        tmpAXp=hAXp[IndexXZ(ix-1,iz)]+(OPTIC_TYPE)2.0*hAXp[Size+IndexXZ(ix,iz)];
        tmpAXm=hAXm[IndexXZ(ix+1,iz)]+(OPTIC_TYPE)2.0*hAXm[Size+IndexXZ(ix,iz)];
        tmpAZp=hAZp[IndexXZ(ix,iz-1)]+(OPTIC_TYPE)2.0*hAZp[Size+IndexXZ(ix,iz)];
        tmpAZm=hAZm[IndexXZ(ix,iz+1)]+(OPTIC_TYPE)2.0*hAZm[Size+IndexXZ(ix,iz)];
        hAXp[2*Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*(tmpAZp+tmpAZm)+
                             hAXp[Size+IndexXZ(ix,iz)]+hAXp[IndexXZ(ix-1,iz)];
        hAXm[2*Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*(tmpAZp+tmpAZm)+
                             hAXm[Size+IndexXZ(ix,iz)]+hAXm[IndexXZ(ix+1,iz)];
        hAZp[2*Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*(tmpAXp+tmpAXm)+
                             hAZp[Size+IndexXZ(ix,iz)]+hAZp[IndexXZ(ix,iz-1)];
        hAZm[2*Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*(tmpAXp+tmpAXm)+
                             hAZm[Size+IndexXZ(ix,iz)]+hAZm[IndexXZ(ix,iz+1)];
        if(iz<100) {
          AXpk1Arr[iz]=hAXp[Size+IndexXZ(ix,iz)];
          AXpk2Arr[iz]=dtT2*iii*InitData.ACoef*(tmpAZp+tmpAZm);
          AXpk12Arr[iz]=AXpk1Arr[iz]+AXpk2Arr[iz];
        }
/*
        AXpk1=hAXp[Size+IndexXZ(ix,iz)];
        AXpk2=dtT2*iii*InitData.ACoef*(tmpAZp+tmpAZm);
        AXmk1=hAXm[Size+IndexXZ(ix,iz)];
        AXmk2=dtT2*iii*InitData.ACoef*(tmpAZp+tmpAZm);
        AZpk1=hAZp[Size+IndexXZ(ix,iz)];
        AZpk2=dtT2*iii*InitData.ACoef*(tmpAXp+tmpAXm);
        AZmk1=hAZm[Size+IndexXZ(ix,iz)];
        AZmk2=dtT2*iii*InitData.ACoef*(tmpAXp+tmpAXm);

        hAXp[2*Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*(tmpAZp+tmpAZm)+
                             hAXp[Size+IndexXZ(ix,iz)]+hAXp[IndexXZ(ix-1,iz)];
        hAXm[2*Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*(tmpAZp+tmpAZm)+
                             hAXm[Size+IndexXZ(ix,iz)]+hAXm[IndexXZ(ix+1,iz)];
        hAZp[2*Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*(tmpAXp+tmpAXm)+
                             hAZp[Size+IndexXZ(ix,iz)]+hAZp[IndexXZ(ix,iz-1)];
        hAZm[2*Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*(tmpAXp+tmpAXm)+
                             hAZm[Size+IndexXZ(ix,iz)]+hAZm[IndexXZ(ix,iz+1)];

        hAXp[2*Size+IndexXZ(ix,iz)]=AXpk1+AXpk2+hAXp[IndexXZ(ix-1,iz)];
        hAXm[2*Size+IndexXZ(ix,iz)]=AXmk1+AXmk2+hAXm[IndexXZ(ix+1,iz)];
        hAZp[2*Size+IndexXZ(ix,iz)]=AZpk1+AZpk2+hAZp[IndexXZ(ix,iz-1)];
        hAZm[2*Size+IndexXZ(ix,iz)]=AZmk1+AZmk2+hAZm[IndexXZ(ix,iz+1)];

        if(iz==10) {
          AZp10[ix]=hAZp[IndexXZ(ix,iz-1)];
          AZp10k1[ix]=AZpk1;
          AZp10k2[ix]=AZpk2;
          AZp10R[ix]=hAZp[2*Size+IndexXZ(ix,iz)];
        }
*/
/*
        hAXp[2*Size+IndexXZ(ix,iz)]=tmpAXp;
        hAXm[2*Size+IndexXZ(ix,iz)]=tmpAXm;
        hAZp[2*Size+IndexXZ(ix,iz)]=tmpAZp;
        hAZm[2*Size+IndexXZ(ix,iz)]=tmpAZm;
*/
      }
    }
    RungeStep++;
    break;

    case 2:
    memcpy(hAXp+IndexXZ(ix1st,0),hAXp+2*Size+IndexXZ(ix1st,0),
      (InitData.nStepsZ+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    memcpy(hAXm+IndexXZ(ix1st,0),hAXm+2*Size+IndexXZ(ix1st,0),
      (InitData.nStepsZ+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    memcpy(hAZp+IndexXZ(ix1st,0),hAZp+2*Size+IndexXZ(ix1st,0),
      (InitData.nStepsZ+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    memcpy(hAZm+IndexXZ(ix1st,0),hAZm+2*Size+IndexXZ(ix1st,0),
      (InitData.nStepsZ+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    RungeStep=0;
  }
}

void cEquOpticBaseTask::StepNoMedia() {
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1);
  /*
  if(RungeStep==0) {
    for()
  } */
}

