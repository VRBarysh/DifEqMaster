//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquOpticTask_MultiMedia.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

OPTIC_TYPE cEquOpticTask_MultiMedia::MediaDistr(int im) {
  int tm=InitData.nMedia/2;
  OPTIC_TYPE tmpC=InitData.DeltaCoef/InitData.MediaGaussCoef;
  return MediaDistNormCoef*exp(-tmpC*tmpC*(im-tm)*(im-tm));
}

OPTIC_TYPE cEquOpticTask_MultiMedia::MediaDelta(int im) {
  int tm=InitData.nMedia/2;
  return (InitData.DeltaCoef*(im-tm));
}

void cEquOpticTask_MultiMedia::MakeMediaDistrNormCoef() {
  OPTIC_TYPE tmp=0, tmpC=InitData.DeltaCoef/InitData.MediaGaussCoef;
  int tm=InitData.nMedia/2;
  for(int im=0; im<InitData.nMedia; im++) tmp+=exp(-tmpC*tmpC*(im-tm)*(im-tm));
  MediaDistNormCoef = 1.0/tmp;
}

void cEquOpticTask_MultiMedia::StepMediaRunge2() {
  int SizeM=(InitData.nStepsX+2)*(InitData.nStepsZ+2)*InitData.nMedia;
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  OPTIC_TYPE dtT2=0.5*dtT;
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1),tmpAXp,tmpAXm,tmpAZp,tmpAZm;
  complex<OPTIC_TYPE> tmpPXp,tmpPXm,tmpPZp,tmpPZm;
  complex<OPTIC_TYPE> PCoef;
  cOptic2D_Media_X tM,tM1;
  complex<OPTIC_TYPE> tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9;
  OPTIC_TYPE rtmp0,rtmp1,rtmp2,rtmp3,rtmp4,rtmp5,rtmp6,rtmp7,rtmp8,rtmp9;
  int iI,iM;
  cOptic2D_Media_X tmpMedia;
  MakeMediaDistrNormCoef();
  switch(MediaRungeStep) {
    case 0:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        tmpPXp=tmpPXm=tmpPZp=tmpPZm=0;
        for(int im=0; im<InitData.nMedia; im++) {
          PCoef=complex<OPTIC_TYPE>(InitData.PCoef,MediaDelta(im));
          iI=IndexXZ(ix,iz);
          iM=IndexXZM(ix,iz,im);
          hM[SizeM+iM].PXp=((hAXp[iI]*hM[iM].R0*((OPTIC_TYPE)2.0)+
                          hAXm[iI]*hM[iM].RX+hAZp[iI]*conj(hM[iM].Rm)+
                          hAZm[iI]*hM[iM].Rp)*InitData.beta -
                          hM[iM].PXp*PCoef)*dtT2;
          hM[SizeM+iM].PXm=((hAXm[iI]*hM[iM].R0*((OPTIC_TYPE)2.0)+
                          hAXp[iI]*hM[iM].RZ+hAZp[iI]*conj(hM[iM].Rp)+
                          hAZm[iI]*hM[iM].Rm)*InitData.beta -
                          hM[iI].PXm*PCoef)*dtT2;
          hM[SizeM+iM].PZp=((hAZp[iI]*hM[iM].R0*((OPTIC_TYPE)2.0)+
                          hAZm[iI]*hM[iM].RZ+hAXp[iI]*hM[iM].Rm+
                          hAXm[iI]*hM[iM].Rp)*InitData.beta -
                          hM[iM].PZp*PCoef)*dtT2;
          hM[SizeM+iM].PZm=((hAZm[iI]*hM[iM].R0*((OPTIC_TYPE)2.0)+
                          hAZp[iI]*conj(hM[iM].RZ)+hAXp[iI]*conj(hM[iM].Rp)+
                          hAXm[iI]*conj(hM[iM].Rm))*InitData.beta -
                          hM[iM].PZm*PCoef)*dtT2;
          hM[SizeM+iM].RX=(-hM[iM].RX*InitData.RCoef
                          -(hAXp[iI]*conj(hM[iM].PXm)
                          +hM[iM].PXp*conj(hAXm[iI])) )*dtT2;
          hM[SizeM+iM].RZ=(-hM[iM].RZ*InitData.RCoef
                          -(hAZp[iI]*conj(hM[iM].PZm)
                          +hM[iM].PZp*conj(hAZm[iI])) )*dtT2;
          hM[SizeM+iM].Rp=(-hM[iM].Rp*InitData.RCoef
                          -(hM[iM].PZp*conj(hAXm[iI])
                          +hAZp[iI]*conj(hM[iM].PXm)
                          +hM[iM].PXp*conj(hAZm[iI])
                          +hAXp[iI]*conj(hM[iM].PZm)) )*dtT2;
          hM[SizeM+iM].Rm=(-hM[iM].Rm*InitData.RCoef
                          -(hM[iM].PXm*conj(hAZm[iI])
                          +hAZp[iI]*conj(hM[iM].PXp)
                          +hAXm[iI]*conj(hM[iM].PZm)
                          +hM[iM].PZp*conj(hAXp[iI])) )*dtT2;
          hM[SizeM+iM].R0=(InitData.Q-InitData.RCoef*hM[iM].R0
                         -(hAZp[iI]*conj(hM[iM].PZp)
                          +hAXp[iI]*conj(hM[iM].PXp)
                          +hAZm[iI]*conj(hM[iM].PZm)
                          +hAXm[iI]*conj(hM[iM].PXm)).real() )*dtT2;
          hM[SizeM+iM].R0Generated=(InitData.Q-InitData.RCoef*hM[iM].R0)*dtT2;
          tmpPXp+=hM[IndexXZM(ix-1,iz,im)].PXp*MediaDistr(im);
          tmpPXm+=hM[IndexXZM(ix+1,iz,im)].PXm*MediaDistr(im);
          tmpPZp+=hM[IndexXZM(ix,iz-1,im)].PZp*MediaDistr(im);
          tmpPZm+=hM[IndexXZM(ix,iz+1,im)].PZm*MediaDistr(im);
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
                       +dtT2*InitData.AGenZ*hAXp[IndexXZ(ix,iz-1)];
        hAZm[Size+iI]+=dtT2*tmpPZm
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
        tmpPXp=tmpPXm=tmpPZp=tmpPZm=0;
        tmpAXp=hAXp[IndexXZ(ix-1,iz)]+(OPTIC_TYPE)2.0*hAXp[Size+iI];
        tmpAXm=hAXm[IndexXZ(ix+1,iz)]+(OPTIC_TYPE)2.0*hAXm[Size+iI];
        tmpAZp=hAZp[IndexXZ(ix,iz-1)]+(OPTIC_TYPE)2.0*hAZp[Size+iI];
        tmpAZm=hAZm[IndexXZ(ix,iz+1)]+(OPTIC_TYPE)2.0*hAZm[Size+iI];
        for(int im=0; im<InitData.nMedia; im++) {
          PCoef=complex<OPTIC_TYPE>(InitData.PCoef,MediaDelta(im));
          iM=IndexXZM(ix,iz,im);
          tM.M1M2v(hM[iM],hM[SizeM+iM],2);
          hM[2*SizeM+iM].M1M2v(hM[iM],hM[SizeM+iM],1);
          hM[2*SizeM+iM].PXp+=((tmpAXp*tM.R0*((OPTIC_TYPE)2.0)+
                             tmpAXm*tM.RX+tmpAZp*conj(tM.Rm)+
                             tmpAZm*tM.Rp)*InitData.beta -
                             tM.PXp*PCoef)*dtT2;
          hM[2*SizeM+iM].PXm+=((tmpAXm*tM.R0*((OPTIC_TYPE)2.0)+
                             tmpAXp*tM.RZ+tmpAZp*conj(tM.Rp)+
                             tmpAZm*tM.Rm)*InitData.beta -
                             tM.PXm*PCoef)*dtT;
          hM[2*SizeM+iM].PZp+=((tmpAZp*tM.R0*((OPTIC_TYPE)2.0)+
                             tmpAZm*tM.RZ+tmpAXp*tM.Rm+
                             tmpAXm*tM.Rp)*InitData.beta -
                             tM.PZp*PCoef)*dtT;
          hM[2*SizeM+iM].PZm+=((tmpAZm*tM.R0*((OPTIC_TYPE)2.0)+
                             tmpAZp*conj(tM.RZ)+tmpAXp*conj(tM.Rp)+
                             tmpAXm*conj(tM.Rm))*InitData.beta -
                             tM.PZm*PCoef)*dtT;
          hM[2*SizeM+iM].RX+=(-tM.RX*InitData.RCoef
                             -(tmpAXp*conj(tM.PXm)
                             +tM.PXp*conj(tmpAXm)) )*dtT;
          hM[2*SizeM+iM].RZ+=(-tM.RZ*InitData.RCoef
                             -(tmpAZp*conj(tM.PZm)
                             +tM.PZp*conj(tmpAZm)) )*dtT;
          hM[2*SizeM+iM].Rp+=(-tM.Rp*InitData.RCoef
                             -(tM.PZp*conj(tmpAXm)
                             +tmpAZp*conj(tM.PXm)
                             +tM.PXp*conj(tmpAZm)
                             +tmpAXp*conj(tM.PZm)) )*dtT;
          hM[2*SizeM+iM].Rm+=(-tM.Rm*InitData.RCoef
                             -(tM.PXm*conj(tmpAZm)
                             +tmpAZp*conj(tM.PXp)
                             +tmpAXm*conj(tM.PZm)
                             +tM.PZp*conj(tmpAXp)) )*dtT;
          hM[2*SizeM+iM].R0+=(InitData.Q-InitData.RCoef*tM.R0
                             -(tmpAZp*conj(tM.PZp)
                             +tmpAXp*conj(tM.PXp)
                            +tmpAZm*conj(tM.PZm)
                            +tmpAXm*conj(tM.PXm)).real() )*dtT2;
          hM[2*SizeM+iM].R0Generated+=(InitData.Q-InitData.RCoef*tM.R0)*dtT2;

          tmpPXp+=tM.PXp*MediaDistr(im);
          tmpPXm+=tM.PXm*MediaDistr(im);
          tmpPZp+=tM.PZp*MediaDistr(im);
          tmpPZm+=tM.PZm*MediaDistr(im);
        }

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

      memcpy(hM+IndexXZM(ix1st,0,0),hM+2*SizeM+IndexXZM(ix1st,0,0),
        (InitData.nStepsZ+2)*sizeof(cOptic2D_Media_X)*(ixLast-ix1st+1)*InitData.nMedia);

    MediaRungeStep=0;
    break;
  }
}

/*
void cEquOpticTask_MultiMedia::StepMediaRunge2FastP() {
  int SizeM=(InitData.nStepsX+2)*(InitData.nStepsZ+2)*InitData.nMedia;
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  OPTIC_TYPE dtT2=0.5*dtT;
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1),tmpAXp,tmpAXm,tmpAZp,tmpAZm;
  complex<OPTIC_TYPE> tmpPXp,tmpPXm,tmpPZp,tmpPZm,tmpPXp2,tmpPXm2,tmpPZp2,tmpPZm2;
  complex<OPTIC_TYPE> PCoef,beta_PCoef;
  cOptic2D_Media_X tM,tM1;
  complex<OPTIC_TYPE> tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9;
  OPTIC_TYPE rtmp0,rtmp1,rtmp2,rtmp3,rtmp4,rtmp5,rtmp6,rtmp7,rtmp8,rtmp9;
  int iI,iM;
  cOptic2D_Media_X tmpMedia;
  MakeMediaDistrNormCoef();
  switch(MediaRungeStep) {
    case 0:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        tmpPXp2=tmpPXm2=tmpPZp2=tmpPZm2=0;
        for(int im=0; im<InitData.nMedia; im++) {
          beta_PCoef=InitData.beta/complex<OPTIC_TYPE>(InitData.PCoef,InitData.Delta+MediaDelta(im));
          iI=IndexXZ(ix,iz);
          iM=IndexXZM(ix,iz,im);
//        hM[Size+iI].PXp=tmpPXp=(hAXp[iI]*hM[iI].R0*((OPTIC_TYPE)2.0)+
//                hAXm[iI]*hM[iI].RX+hAZp[iI]*conj(hM[iI].Rm)+
//                hAZm[iI]*hM[iI].Rp)*beta_PCoef;
          hM[SizeM+iM].PXp=tmpPXp=(hAXp[iI]*hM[iM].R0*((OPTIC_TYPE)2.0)+
                          hAXm[iI]*hM[iM].RX+hAZp[iI]*conj(hM[iM].Rm)+
                          hAZm[iI]*hM[iM].Rp)*beta_PCoef;
          hM[SizeM+iM].PXm=tmpPXm=(hAXm[iI]*hM[iM].R0*((OPTIC_TYPE)2.0)+
                          hAXp[iI]*hM[iM].RZ+hAZp[iI]*conj(hM[iM].Rp)+
                          hAZm[iI]*hM[iM].Rm)*beta_PCoef;
          hM[SizeM+iM].PZp=tmpPZp=(hAZp[iI]*hM[iM].R0*((OPTIC_TYPE)2.0)+
                          hAZm[iI]*hM[iM].RZ+hAXp[iI]*hM[iM].Rm+
                          hAXm[iI]*hM[iM].Rp)*beta_PCoef;
          hM[SizeM+iM].PZm=tmpPZm=(hAZm[iI]*hM[iM].R0*((OPTIC_TYPE)2.0)+
                          hAZp[iI]*conj(hM[iM].RZ)+hAXp[iI]*conj(hM[iM].Rp)+
                          hAXm[iI]*conj(hM[iM].Rm))*beta_PCoef;
          hM[SizeM+iM].RX=(-hM[iM].RX*InitData.RCoef
                          -(hAXp[iI]*conj(tmpPXm)
                          +tmpPXp*conj(hAXm[iI])) )*dtT2;
          hM[SizeM+iM].RZ=(-hM[iM].RZ*InitData.RCoef
                          -(hAZp[iI]*conj(tmpPZm)
                          +tmpPZp*conj(hAZm[iI])) )*dtT2;
          hM[SizeM+iM].Rp=(-hM[iM].Rp*InitData.RCoef
                          -(tmpPZp*conj(hAXm[iI])
                          +hAZp[iI]*conj(tmpPXm)
                          +tmpPXp*conj(hAZm[iI])
                          +hAXp[iI]*conj(tmpPZm)) )*dtT2;
          hM[SizeM+iM].Rm=(-hM[iM].Rm*InitData.RCoef
                          -(hM[iM].PXm*conj(hAZm[iI])
                          +hAZp[iI]*conj(tmpPXp)
                          +hAXm[iI]*conj(tmpPZm)
                          +tmpPZp*conj(hAXp[iI])) )*dtT2;
          hM[SizeM+iM].R0=(InitData.Q-InitData.RCoef*hM[iM].R0
                         -(hAZp[iI]*conj(tmpPZp)
                          +hAXp[iI]*conj(tmpPXp)
                          +hAZm[iI]*conj(tmpPZm)
                          +hAXm[iI]*conj(tmpPXm)).real() )*dtT2;
          hM[SizeM+iM].R0Generated=(InitData.Q-InitData.RCoef*hM[iM].R0)*dtT2;
          iI=IndexXZ(ix-1,iz); iM=IndexXZM(ix-1,iz,im);
          tmpPXp2+=(hAXp[iI]*hM[iM].R0*((OPTIC_TYPE)2.0)+
                          hAXm[iI]*hM[iM].RX+hAZp[iI]*conj(hM[iM].Rm)+
                          hAZm[iI]*hM[iM].Rp)*beta_PCoef*MediaDistr(im);
          iI=IndexXZ(ix+1,iz); iM=IndexXZM(ix+1,iz,im);
          tmpPXm2+=(hAXm[iI]*hM[iM].R0*((OPTIC_TYPE)2.0)+
                          hAXp[iI]*hM[iM].RZ+hAZp[iI]*conj(hM[iM].Rp)+
                          hAZm[iI]*hM[iM].Rm)*beta_PCoef*MediaDistr(im);
          iI=IndexXZ(ix,iz-1); iM=IndexXZM(ix,iz-1,im);
          tmpPZp2+=(hAZp[iI]*hM[iM].R0*((OPTIC_TYPE)2.0)+
                          hAZm[iI]*hM[iM].RZ+hAXp[iI]*hM[iM].Rm+
                          hAXm[iI]*hM[iM].Rp)*beta_PCoef*MediaDistr(im);
          iI=IndexXZ(ix,iz+1); iM=IndexXZM(ix,iz+1,im);
          tmpPZm2+=(hAZm[iI]*hM[iM].R0*((OPTIC_TYPE)2.0)+
                          hAZp[iI]*conj(hM[iM].RZ)+hAXp[iI]*conj(hM[iM].Rp)+
                          hAXm[iI]*conj(hM[iM].Rm))*beta_PCoef*MediaDistr(im);
        }

        hAXp[Size+iI]+=dtT2*tmpPXp2
                         +dtT2*InitData.AGenCoef*hAXp[IndexXZ(ix-1,iz)]
                         +dtT2*InitData.ABackCoef*hAXm[IndexXZ(ix-1,iz)]
                         +dtT2*InitData.AGenX*hAXp[IndexXZ(ix-1,iz)];
        hAXm[Size+iI]+=dtT2*tmpPXm2
                       +dtT2*InitData.AGenCoef*hAXm[IndexXZ(ix+1,iz)]
                       +dtT2*InitData.ABackCoef*hAXp[IndexXZ(ix+1,iz)]
                       +dtT2*InitData.AGenX*hAXm[IndexXZ(ix+1,iz)];
        hAZp[Size+iI]+=dtT2*tmpPZp2
                       +dtT2*InitData.AGenCoef*hAZp[IndexXZ(ix,iz-1)]
                       +dtT2*InitData.ABackCoef*hAZm[IndexXZ(ix,iz-1)]
                       +dtT2*InitData.AGenZ*hAXp[IndexXZ(ix,iz-1)];
        hAZm[Size+iI]+=dtT2*tmpPZm2
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
        tmpPXp=tmpPXm=tmpPZp=tmpPZm=0;
        tmpAXp=hAXp[IndexXZ(ix-1,iz)]+(OPTIC_TYPE)2.0*hAXp[Size+iI];
        tmpAXm=hAXm[IndexXZ(ix+1,iz)]+(OPTIC_TYPE)2.0*hAXm[Size+iI];
        tmpAZp=hAZp[IndexXZ(ix,iz-1)]+(OPTIC_TYPE)2.0*hAZp[Size+iI];
        tmpAZm=hAZm[IndexXZ(ix,iz+1)]+(OPTIC_TYPE)2.0*hAZm[Size+iI];
        for(int im=0; im<InitData.nMedia; im++) {
          PCoef=complex<OPTIC_TYPE>(InitData.PCoef,MediaDelta(im));
          beta_PCoef=InitData.beta/complex<OPTIC_TYPE>(InitData.PCoef,InitData.Delta+MediaDelta(im));
          iM=IndexXZM(ix,iz,im);
          tM.M1M2v(hM[iM],hM[SizeM+iM],2);
          hM[2*SizeM+iM].M1M2v(hM[iM],hM[SizeM+iM],1);
          hM[2*SizeM+iM].PXp+=((tmpAXp*tM.R0*((OPTIC_TYPE)2.0)+
                             tmpAXm*tM.RX+tmpAZp*conj(tM.Rm)+
                             tmpAZm*tM.Rp)*InitData.beta -
                             tM.PXp*PCoef)*dtT2;
          hM[2*SizeM+iM].PXm+=((tmpAXm*tM.R0*((OPTIC_TYPE)2.0)+
                             tmpAXp*tM.RZ+tmpAZp*conj(tM.Rp)+
                             tmpAZm*tM.Rm)*InitData.beta -
                             tM.PXm*PCoef)*dtT;
          hM[2*SizeM+iM].PZp+=((tmpAZp*tM.R0*((OPTIC_TYPE)2.0)+
                             tmpAZm*tM.RZ+tmpAXp*tM.Rm+
                             tmpAXm*tM.Rp)*InitData.beta -
                             tM.PZp*PCoef)*dtT;
          hM[2*SizeM+iM].PZm+=((tmpAZm*tM.R0*((OPTIC_TYPE)2.0)+
                             tmpAZp*conj(tM.RZ)+tmpAXp*conj(tM.Rp)+
                             tmpAXm*conj(tM.Rm))*InitData.beta -
                             tM.PZm*PCoef)*dtT;
          hM[2*SizeM+iM].RX+=(-tM.RX*InitData.RCoef
                             -(tmpAXp*conj(tM.PXm)
                             +tM.PXp*conj(tmpAXm)) )*dtT;
          hM[2*SizeM+iM].RZ+=(-tM.RZ*InitData.RCoef
                             -(tmpAZp*conj(tM.PZm)
                             +tM.PZp*conj(tmpAZm)) )*dtT;
          hM[2*SizeM+iM].Rp+=(-tM.Rp*InitData.RCoef
                             -(tM.PZp*conj(tmpAXm)
                             +tmpAZp*conj(tM.PXm)
                             +tM.PXp*conj(tmpAZm)
                             +tmpAXp*conj(tM.PZm)) )*dtT;
          hM[2*SizeM+iM].Rm+=(-tM.Rm*InitData.RCoef
                             -(tM.PXm*conj(tmpAZm)
                             +tmpAZp*conj(tM.PXp)
                             +tmpAXm*conj(tM.PZm)
                             +tM.PZp*conj(tmpAXp)) )*dtT;
          hM[2*SizeM+iM].R0+=(InitData.Q-InitData.RCoef*tM.R0
                             -(tmpAZp*conj(tM.PZp)
                             +tmpAXp*conj(tM.PXp)
                            +tmpAZm*conj(tM.PZm)
                            +tmpAXm*conj(tM.PXm)).real() )*dtT2;
          hM[2*SizeM+iM].R0Generated+=(InitData.Q-InitData.RCoef*tM.R0)*dtT2;

          tmpPXp+=tM.PXp*MediaDistr(im);
          tmpPXm+=tM.PXm*MediaDistr(im);
          tmpPZp+=tM.PZp*MediaDistr(im);
          tmpPZm+=tM.PZm*MediaDistr(im);
        }

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

      memcpy(hM+IndexXZM(ix1st,0,0),hM+2*SizeM+IndexXZM(ix1st,0,0),
        (InitData.nStepsZ+2)*sizeof(cOptic2D_Media_X)*(ixLast-ix1st+1)*InitData.nMedia);

    MediaRungeStep=0;
    break;
  }
}
*/
