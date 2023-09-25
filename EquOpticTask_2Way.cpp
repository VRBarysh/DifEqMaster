//---------------------------------------------------------------------------                                            Z0RefCoef


#pragma hdrstop

#include "EquOpticTask_2Way.h"
#include <fstream.h>

//---------------------------------------------------------------------------

#pragma package(smart_init)

void cEquOpticBaseTask_2Way::Step() {
  StepNoMediaMatrix();
  if(!InitData.fNoMedia) {
    if(InitData.UseFastP) {
      StepMediaRunge2FastP();
    } else {
      StepMediaRunge2();
    }
  }
}


OPTIC_TYPE cEquOpticBaseTask_2Way::BetaF(int ix){
int s2=(InitData.nStepsX-InitData.nMediaRowSize)/InitData.nMediaRows;
int ss=s2-InitData.nMediaRowSize;
int sr=(InitData.nStepsX-InitData.nMediaRowSize) % InitData.nMediaRows;
int dx=ix-1-(s2+1)*sr;
int res = dx>0 ? (dx % s2 <InitData.nMediaRowSize) : ( (ix-1) % (s2+1) <InitData.nMediaRowSize );
return res ? InitData.beta : 0;
}

void cEquOpticBaseTask_2Way::StepMediaRunge2() {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZB+2);
  int SizeBragg=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  OPTIC_TYPE dtT2=0.5*dtT;
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1);
  complex<OPTIC_TYPE> tmpAXp,tmpAXm,tmpAZp,tmpAZm,tmpPXp,tmpPXm,tmpPZp,tmpPZm;
  complex<OPTIC_TYPE> ZRef = polar(InitData.Z0RefCoef,InitData.Z0RefPhase);
  cOptic2D_Media_X_2Way tM,tM1;
  int iI;
  cOptic2D_Media_X tmpMedia;
  OPTIC_TYPE beta=InitData.beta;
  complex<OPTIC_TYPE> PCoef=InitData.PCoef;
  switch(MediaRungeStep) {
    case 0:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      //if(isBragg(ix)) {
        PCoef=complex<OPTIC_TYPE>(InitData.PCoef,InitData.Delta);
        for(int iz=1;iz<=InitData.nStepsZB;iz++) {
          iI=IndexXZ_M(ix,iz);

          hM2[Size+iI].PZp=tmpPZp=(hAZpMedia[iI]*hM2[iI].R0*((OPTIC_TYPE)2.0)+
                  hAZmMedia[iI]*hM2[iI].RZ)*beta;
          hM2[Size+iI].PZm=tmpPZm=(hAZmMedia[iI]*hM2[iI].R0*((OPTIC_TYPE)2.0)+
                  hAZpMedia[iI]*conj(hM2[iI].RZ))*beta;

          hM2[Size+iI].PZp=((hAZpMedia[iI]*hM2[iI].R0*((OPTIC_TYPE)2.0)+
                          hAZmMedia[iI]*hM2[iI].RZ)*beta -
                          hM2[iI].PZp*PCoef)*dtT2;
          hM2[Size+iI].PZm=((hAZmMedia[iI]*hM2[iI].R0*((OPTIC_TYPE)2.0)+
                          hAZpMedia[iI]*conj(hM2[iI].RZ))*beta -
                          hM2[iI].PZm*PCoef)*dtT2;

          hM2[Size+iI].RZ=(-hM2[iI].RZ*InitData.RCoef
                          -(hAZpMedia[iI]*conj(hM2[iI].PZm)
                          +hM2[iI].PZp*conj(hAZmMedia[iI])) )*dtT2;
          hM2[Size+iI].R0=(InitData.Q-InitData.RCoef*hM2[iI].R0
                          -(hAZpMedia[iI]*conj(hM2[iI].PZp)
                          +hAZmMedia[iI]*conj(hM2[iI].PZm)).real() )*dtT2;
          hM2[Size+iI].R0Generated=(InitData.Q-InitData.RCoef*hM2[iI].R0)*dtT2;

          hAZpMedia[Size+iI]+=dtT2*hM2[IndexXZ_M(ix,iz-1)].PZp
                        +dtT2*InitData.AGenCoef*hAZpMedia[IndexXZ_M(ix,iz-1)]
                        +dtT2*InitData.AGenZ*hAZpMedia[IndexXZ_M(ix,iz-1)];
          hAZmMedia[Size+iI]+=dtT2*hM2[IndexXZ_M(ix,iz+1)].PZm
                        +dtT2*InitData.AGenCoef*hAZmMedia[IndexXZ_M(ix,iz+1)]
                        +dtT2*InitData.AGenZ*hAZmMedia[IndexXZ_M(ix,iz+1)];
        }
//      }
    }
    MediaRungeStep++;
    break;
    case 1:
    for(int ix=ix1st;ix<=ixLast;ix++) {
//      if(isBragg(ix)) {
        PCoef=complex<OPTIC_TYPE>(InitData.PCoef,InitData.Delta);
        for(int iz=1;iz<=InitData.nStepsZB;iz++) {
          iI=IndexXZ_M(ix,iz);
          tmpAZp=hAZpMedia[IndexXZ_M(ix,iz-1)]+(OPTIC_TYPE)2.0*hAZpMedia[Size+iI];
          tmpAZm=hAZmMedia[IndexXZ_M(ix,iz+1)]+(OPTIC_TYPE)2.0*hAZmMedia[Size+iI];
          tM.M1M2v(hM2[iI],hM2[Size+iI],2);
          hM2[2*Size+iI].M1M2v(hM2[iI],hM2[Size+iI],1);

          tmpPZp=(tmpAZp*tM.R0*((OPTIC_TYPE)2.0)+
                  tmpAZm*tM.RZ)*beta;
          tmpPZm=(tmpAZm*tM.R0*((OPTIC_TYPE)2.0)+
                  tmpAZp*conj(tM.RZ))*beta;

          hM2[2*Size+iI].PZp+=((tmpAZp*tM.R0*((OPTIC_TYPE)2.0)+
                               tmpAZm*tM.RZ)*beta -
                               tM.PZp*PCoef)*dtT;
          hM2[2*Size+iI].PZm+=((tmpAZm*tM.R0*((OPTIC_TYPE)2.0)+
                               tmpAZp*conj(tM.RZ))*beta -
                               tM.PZm*PCoef)*dtT;

          tmpPZm=tM.PZm; tmpPZp=tM.PZp;

          hM2[2*Size+iI].RZ+=(-tM.RZ*InitData.RCoef
                             -(tmpAZp*conj(tmpPZm)
                             +tmpPZp*conj(tmpAZm)) )*dtT;
          hM2[2*Size+iI].R0+=(InitData.Q-InitData.RCoef*tM.R0
                            -(tmpAZp*conj(tmpPZp)
                            +tmpAZm*conj(tmpPZm)).real() )*dtT2;
          hM2[2*Size+iI].R0Generated+=(InitData.Q-InitData.RCoef*tM.R0)*dtT2;
          hAZpMedia[2*Size+iI]+=dtT2*tmpPZp+dtT2*InitData.AGenCoef*tmpAZp
                                      +dtT2*InitData.AGenZ*tmpAZp;
          hAZmMedia[2*Size+iI]+=dtT2*tmpPZm+dtT2*InitData.AGenCoef*tmpAZm
                                      +dtT2*InitData.AGenZ*tmpAZm;
        }
/*
      } else {
        for(int iz=1;iz<=InitData.nStepsZB;iz++) {
          iI=IndexXZ_M(ix,iz);
          hM2[2*Size+iI]=hM2[iI];
        }
      }*/
    }
    MediaRungeStep++;
    break;
    case 2:
      memcpy(hM2+IndexXZ_M(ix1st,0),hM2+2*Size+IndexXZ_M(ix1st,0),
        (InitData.nStepsZB+2)*sizeof(cOptic2D_Media_X_2Way)*(ixLast-ix1st+1));
    MediaRungeStep=0;
    break;
  }
}

/*
void cEquOpticBaseTask_2Way::StepMediaRunge2() {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZB+2);
  int SizeBragg=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  complex<OPTIC_TYPE> ZRef = polar(InitData.Z0RefCoef,InitData.Z0RefPhase);
  OPTIC_TYPE dtT2=0.5*dtT;
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1),tmpAXp,tmpAXm,tmpAZp,tmpAZm;
  complex<OPTIC_TYPE> tmpPXp,tmpPXm,tmpPZp,tmpPZm;
  cOptic2D_Media_X_2Way tM,tM1;
  complex<OPTIC_TYPE> tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9;
  OPTIC_TYPE rtmp0,rtmp1,rtmp2,rtmp3,rtmp4,rtmp5,rtmp6,rtmp7,rtmp8,rtmp9;
  //int istopX=15,istopZ=30;
  int iI;
  OPTIC_TYPE beta=InitData.beta;
  cOptic2D_Media_X_2Way tmpMedia;
  complex<OPTIC_TYPE> PCoef = complex<OPTIC_TYPE>(InitData.PCoef,InitData.Delta);
  switch(MediaRungeStep) {
    case 0:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      if(InitData.nMediaRows) beta=BetaF(ix);
      for(int iz=1;iz<=InitData.nStepsZB;iz++) {
        iI=IndexXZ_M(ix,iz);
        hM2[Size+iI].PZp=((hAZpMedia[iI]*hM2[iI].R0*((OPTIC_TYPE)2.0)+
                        hAZmMedia[iI]*hM2[iI].RZ)*beta -
                        hM2[iI].PZp*PCoef)*dtT2;
        hM2[Size+iI].PZm=((hAZmMedia[iI]*hM2[iI].R0*((OPTIC_TYPE)2.0)+
                        hAZpMedia[iI]*conj(hM2[iI].RZ))*beta -
                        hM2[iI].PZm*PCoef)*dtT2;
        hM2[Size+iI].RZ=(-hM2[iI].RZ*InitData.RCoef
                        -(hAZp[iI]*conj(hM2[iI].PZm)
                        +hM2[iI].PZp*conj(hAZm[iI])) )*dtT2;
        hM2[Size+iI].R0=(InitData.Q-InitData.RCoef*hM2[iI].R0
                       -(hAZpMedia[iI]*conj(hM2[iI].PZp)
                        +hAZmMedia[iI]*conj(hM2[iI].PZm)).real() )*dtT2;
        hM2[Size+iI].R0Generated=(InitData.Q-InitData.RCoef*hM2[iI].R0)*dtT2;

        hAZpMedia[Size+iI]+=dtT2*hM2[IndexXZ_M(ix,iz-1)].PZp
                       +dtT2*InitData.AGenCoef*hAZpMedia[IndexXZ_M(ix,iz-1)]
                       +dtT2*InitData.AGenZ*hAZpMedia[IndexXZ_M(ix,iz-1)];
        hAZmMedia[Size+iI]+=dtT2*hM2[IndexXZ_M(ix,iz+1)].PZm
                       +dtT2*InitData.AGenCoef*hAZmMedia[IndexXZ_M(ix,iz+1)]
                       +dtT2*InitData.AGenZ*hAZmMedia[IndexXZ_M(ix,iz+1)];
      }
      hAZpMedia[Size+IndexXZ_M(ix,1)]+=ZRef*(dtT2*hM2[Size+iI].PZm
          +dtT2*InitData.AGenCoef*hAZpMedia[IndexXZ_M(ix,InitData.nStepsZB)]
          +dtT2*InitData.AGenZ*hAZpMedia[IndexXZ_M(ix,InitData.nStepsZB)]);
      hAZp[SizeBragg+IndexXZ_M(ix,1)]+=dtT2*hM2[Size+IndexXZ_M(ix,InitData.nStepsZB)].PZp
          +dtT2*InitData.AGenCoef*hAZpMedia[IndexXZ_M(ix,InitData.nStepsZB)]
          +dtT2*InitData.AGenZ*hAZpMedia[IndexXZ_M(ix,InitData.nStepsZB)];
    }
    MediaRungeStep++;
    break;
    case 1:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      if(InitData.nMediaRows) beta=BetaF(ix);
      for(int iz=1;iz<=InitData.nStepsZB;iz++) {
        iI=IndexXZ_M(ix,iz);
        tmpAZp=hAZpMedia[IndexXZ_M(ix,iz-1)]+(OPTIC_TYPE)2.0*hAZpMedia[Size+iI];
        tmpAZm=hAZmMedia[IndexXZ_M(ix,iz+1)]+(OPTIC_TYPE)2.0*hAZmMedia[Size+iI];
//        tM=hM2[iI]; tM+=hM2[Size+iI]*((OPTIC_TYPE)2.0);
//        hM2[2*Size+iI]=hM2[iI]; hM2[2*Size+iI]+=hM2[Size+iI];
        tM.M1M2v(hM2[iI],hM2[Size+iI],2);
        hM2[2*Size+iI].M1M2v(hM2[iI],hM2[Size+iI],1);
//        hM2[2*Size+iI]=hM2[iI]; hM2[2*Size+iI]+=hM2[Size+iI];


        hM2[2*Size+iI].PZp+=((tmpAZp*tM.R0*((OPTIC_TYPE)2.0)+
                             tmpAZm*tM.RZ)*beta -
                             tM.PZp*PCoef)*dtT;
        hM2[2*Size+iI].PZm+=((tmpAZm*tM.R0*((OPTIC_TYPE)2.0)+
                             tmpAZp*conj(tM.RZ))*beta -
                             tM.PZm*PCoef)*dtT;
        hM2[2*Size+iI].RZ+=(-tM.RZ*InitData.RCoef
                           -(tmpAZp*conj(tM.PZm)
                           +tM.PZp*conj(tmpAZm)) )*dtT;
        hM2[2*Size+iI].R0+=(InitData.Q-InitData.RCoef*tM.R0
                          -(tmpAZp*conj(tM.PZp)
                          +tmpAZm*conj(tM.PZm)).real() )*dtT2;
        hM2[2*Size+iI].R0Generated+=(InitData.Q-InitData.RCoef*tM.R0)*dtT2;
        hAZpMedia[2*Size+iI]+=dtT2*tM.PZp+dtT2*InitData.AGenCoef*tmpAZp
                                    +dtT2*InitData.AGenZ*tmpAZp;
        hAZmMedia[2*Size+iI]+=dtT2*tM.PZm+dtT2*InitData.AGenCoef*tmpAZm
                                    +dtT2*InitData.AGenZ*tmpAZm;
      }
    }
    MediaRungeStep++;
    break;
    case 2:
      memcpy(hM2+IndexXZ_M(ix1st,0),hM2+2*Size+IndexXZ_M(ix1st,0),
        (InitData.nStepsZ+2)*sizeof(cOptic2D_Media_X)*(ixLast-ix1st+1));
    MediaRungeStep=0;
    break;
  }
}
*/


void cEquOpticBaseTask_2Way::StepMediaRunge2FastP() {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZB+2);
  int SizeBragg=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  OPTIC_TYPE dtT2=0.5*dtT;
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1);
  complex<OPTIC_TYPE> tmpAXp,tmpAXm,tmpAZp,tmpAZm,tmpPXp,tmpPXm,tmpPZp,tmpPZm;
  complex<OPTIC_TYPE> ZRef = polar(InitData.Z0RefCoef,InitData.Z0RefPhase);
  cOptic2D_Media_X_2Way tM,tM1;
  int iI;
  cOptic2D_Media_X tmpMedia;
  OPTIC_TYPE beta=InitData.beta;
  switch(MediaRungeStep) {
    case 0:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      if(InitData.nMediaRows) beta=BetaF(ix);
      for(int iz=1;iz<=InitData.nStepsZB;iz++) {
        iI=IndexXZ_M(ix,iz);
        hM2[Size+iI].PZp=tmpPZp=(hAZpMedia[iI]*hM2[iI].R0*((OPTIC_TYPE)2.0)+
                hAZmMedia[iI]*hM2[iI].RZ)*beta;
        hM2[Size+iI].PZm=tmpPZm=(hAZmMedia[iI]*hM2[iI].R0*((OPTIC_TYPE)2.0)+
                hAZpMedia[iI]*conj(hM2[iI].RZ))*beta;
        hM2[Size+iI].RZ=(-hM2[iI].RZ*InitData.RCoef
                        -(hAZpMedia[iI]*conj(tmpPZm)
                        +tmpPZp*conj(hAZmMedia[iI])) )*dtT2;
        hM2[Size+iI].R0=(InitData.Q-InitData.RCoef*hM2[iI].R0
                       -(hAZpMedia[iI]*conj(tmpPZp)
                        +hAZmMedia[iI]*conj(tmpPZm)).real() )*dtT2;
        hM2[Size+iI].R0Generated=(InitData.Q-InitData.RCoef*hM2[iI].R0)*dtT2;

        hAZpMedia[Size+iI]+=dtT2*hM2[IndexXZ_M(ix,iz-1)].PZp
                      +dtT2*InitData.AGenCoef*hAZpMedia[IndexXZ_M(ix,iz-1)]
                      +dtT2*InitData.AGenZ*hAZpMedia[IndexXZ_M(ix,iz-1)];
        hAZmMedia[Size+iI]+=dtT2*hM2[IndexXZ_M(ix,iz+1)].PZm
                      +dtT2*InitData.AGenCoef*hAZmMedia[IndexXZ_M(ix,iz+1)]
                      +dtT2*InitData.AGenZ*hAZmMedia[IndexXZ_M(ix,iz+1)];
      }
      /*
      hAZpMedia[Size+IndexXZ_M(ix,1)]+=ZRef*(dtT2*hM2[Size+iI].PZm
          +dtT2*InitData.AGenCoef*hAZpMedia[IndexXZ_M(ix,InitData.nStepsZB)]
          +dtT2*InitData.AGenZ*hAZpMedia[IndexXZ_M(ix,InitData.nStepsZB)]);
      hAZp[SizeBragg+IndexXZ(ix,1)]+=dtT2*hM2[IndexXZ_M(ix,InitData.nStepsZ)].PZp
          +dtT2*InitData.AGenCoef*hAZpMedia[IndexXZ_M(ix,InitData.nStepsZ)]
          +dtT2*InitData.AGenZ*hAZpMedia[IndexXZ_M(ix,InitData.nStepsZ)];
      */
    }
    MediaRungeStep++;
    break;
    case 1:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      if(InitData.nMediaRows) beta=BetaF(ix);
      for(int iz=1;iz<=InitData.nStepsZB;iz++) {
        iI=IndexXZ_M(ix,iz);
        tmpAZp=hAZpMedia[IndexXZ_M(ix,iz-1)]+(OPTIC_TYPE)2.0*hAZpMedia[Size+iI];
        tmpAZm=hAZmMedia[IndexXZ_M(ix,iz+1)]+(OPTIC_TYPE)2.0*hAZmMedia[Size+iI];
//        tM=hM2[iI]; tM+=hM2[Size+iI]*((OPTIC_TYPE)2.0);
//        hM2[2*Size+iI]=hM2[iI]; hM2[2*Size+iI]+=hM2[Size+iI];
        tM.M1M2v(hM2[iI],hM2[Size+iI],2);
        hM2[2*Size+iI].M1M2v(hM2[iI],hM2[Size+iI],1);
//        hM2[2*Size+iI]=hM2[iI]; hM2[2*Size+iI]+=hM2[Size+iI];

        hM2[2*Size+iI].PZp=tmpPZp=(tmpAZp*tM.R0*((OPTIC_TYPE)2.0)+
                tmpAZm*tM.RZ)*beta;
        hM2[2*Size+iI].PZm=tmpPZm=(tmpAZm*tM.R0*((OPTIC_TYPE)2.0)+
                tmpAZp*conj(tM.RZ))*beta;
        hM2[2*Size+iI].RZ+=(-tM.RZ*InitData.RCoef
                           -(tmpAZp*conj(tmpPZm)
                           +tmpPZp*conj(tmpAZm)) )*dtT;
        hM2[2*Size+iI].R0+=(InitData.Q-InitData.RCoef*tM.R0
                          -(tmpAZp*conj(tmpPZp)
                          +tmpAZm*conj(tmpPZm)).real() )*dtT2;
        hM2[2*Size+iI].R0Generated+=(InitData.Q-InitData.RCoef*tM.R0)*dtT2;
        hAZpMedia[2*Size+iI]+=dtT2*tmpPZp+dtT2*InitData.AGenCoef*tmpAZp
                                    +dtT2*InitData.AGenZ*tmpAZp;
        hAZmMedia[2*Size+iI]+=dtT2*tmpPZm+dtT2*InitData.AGenCoef*tmpAZm
                                    +dtT2*InitData.AGenZ*tmpAZm;
      }
    }
    MediaRungeStep++;
    break;
    case 2:
      memcpy(hM2+IndexXZ_M(ix1st,0),hM2+2*Size+IndexXZ_M(ix1st,0),
        (InitData.nStepsZB+2)*sizeof(cOptic2D_Media_X_2Way)*(ixLast-ix1st+1));
    MediaRungeStep=0;
    break;
  }
}

void cEquOpticBaseTask_2Way::StepNoMediaMatrix() {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  int SizeMedia=(InitData.nStepsX+2)*(InitData.nStepsZB+2);
  complex<OPTIC_TYPE> ZRef = polar(InitData.Z0RefCoef,InitData.Z0RefPhase);
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1),tmpAXp,tmpAXm,tmpAZp,tmpAZm;
  complex<OPTIC_TYPE> iiiA=iii*InitData.ACoef*polar(1,InitData.APhase);
//  complex<OPTIC_TYPE> AXpk1,AXpk2,AXmk1,AXmk2,AZpk1,AZpk2,AZmk1,AZmk2;
//  complex<OPTIC_TYPE> AXpk2Arr[100],AXpk1Arr[100],AXpk12Arr[100];
  OPTIC_TYPE dtT2=0.5*dtT;
  int iznb1=1,iznb2=InitData.nStepsZ,ixnb1=ix1st,ixnb2=ixLast;

  int id1,id2,id3,id4,id5,id6,id7,id8,id9;
  complex<OPTIC_TYPE> dd1,dd2,dd3,dd4,dd5,dd6,dd7,dd8,dd9;
  OPTIC_TYPE ddr1,ddr2,ddr3,ddr4;

  id1=IndexXZ(1,InitData.nStepsZ/3+1);
  id2=IndexXZ(2,InitData.nStepsZ/3+1);
  id3=IndexXZ(3,InitData.nStepsZ/3+1);
  id4=IndexXZ(1,InitData.nStepsZ/3);
  id5=IndexXZ(2,InitData.nStepsZ/3);
  id6=IndexXZ(3,InitData.nStepsZ/3);
  id7=IndexXZ(1,InitData.nStepsZ/3-1);
  id8=IndexXZ(2,InitData.nStepsZ/3-1);
  id9=IndexXZ(3,InitData.nStepsZ/3-1);

  switch(RungeStep) {
    case 0:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        hAXp[Size+IndexXZ(ix,iz)]=dtT2*iiiA*
                  (hAZp[IndexXZ(ix-1,iz)]+hAZm[IndexXZ(ix-1,iz)]);
        hAXm[Size+IndexXZ(ix,iz)]=dtT2*iiiA*
                  (hAZp[IndexXZ(ix+1,iz)]+hAZm[IndexXZ(ix+1,iz)]);
        hAZp[Size+IndexXZ(ix,iz)]=dtT2*iiiA*
                  (hAXp[IndexXZ(ix,iz-1)]+hAXm[IndexXZ(ix,iz-1)]);
        hAZm[Size+IndexXZ(ix,iz)]=dtT2*iiiA*
                  (hAXp[IndexXZ(ix,iz+1)]+hAXm[IndexXZ(ix,iz+1)]);
      }
    
      for(int iz=1; iz<=InitData.nStepsZB; iz++){
        hAZpMedia[SizeMedia+IndexXZ_M(ix,iz)]=0;//hAZpMedia[2*Size+IndexXZ(ix,iz-1)];
        hAZmMedia[SizeMedia+IndexXZ_M(ix,iz)]=0;//hAZmMedia[2*Size+IndexXZ(ix,iz+1)];
      }
      hAZpMedia[SizeMedia+IndexXZ_M(ix,1)]=0;//hAZmMedia[2*Size+IndexXZ(ix,1)];
      hAZmMedia[SizeMedia+IndexXZ_M(ix,InitData.nStepsZB)]=dtT2*iiiA*
                  (hAXp[IndexXZ(ix,1)]+hAXm[IndexXZ(ix,1)]);//+hAZm[IndexXZ(ix,1);
    }

    id1=IndexXZ(1,InitData.nStepsZ/3+1);
    id2=IndexXZ(2,InitData.nStepsZ/3+1);
    id3=IndexXZ(3,InitData.nStepsZ/3+1);
    id4=IndexXZ(1,InitData.nStepsZ/3);
    id5=IndexXZ(2,InitData.nStepsZ/3);
    id6=IndexXZ(3,InitData.nStepsZ/3);
    id7=IndexXZ(1,InitData.nStepsZ/3-1);
    id8=IndexXZ(2,InitData.nStepsZ/3-1);
    id9=IndexXZ(3,InitData.nStepsZ/3-1);

    if(fDebug) {
      dd1=hAZp[Size+id4];
      dd2=hAZp[Size+id5];
      dd3=hAZp[id4];
      dd4=hAZp[id5];
      dd5=hAXp[id4];
      dd6=hAXm[id4];
      dd7=hAXp[id5];
      dd8=hAXm[id5];
      //dd4=hAZp[IndexXZ(2,InitData.nStepsZ/3+1)];
      ddr1=abs(dd2)-abs(dd1); ddr2=abs(dd4)-abs(dd3);
      ddr1/=InitData.dt; ddr2/=InitData.dt;
    }

    RungeStep++;
    break;
    case 1:
    iznb1=1; iznb2=InitData.nStepsZ;
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
      /*
      hAXp[2*Size+IndexXZ(1,1)]=hAXp[IndexXZ(0,1)];
      hAZp[2*Size+IndexXZ(1,1)]=hAZp[IndexXZ(1,0)];
      tmpAXm=hAXm[IndexXZ(2,1)]+hAXm[Size+IndexXZ(1,1)];
      tmpAZm=hAZm[IndexXZ(1,2)]+hAZm[Size+IndexXZ(1,1)];
      hAXm[2*Size+IndexXZ(1,1)]=tmpAZm*Matrix2[1]+tmpAXm*Matrix2[0];
      hAZm[2*Size+IndexXZ(1,1)]=tmpAZm*Matrix2[0]+tmpAXm*Matrix2[1];
      */

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
      /*
      hAXm[2*Size+IndexXZ(InitData.nStepsX,1)]=hAXm[IndexXZ(InitData.nStepsX+1,1)];
      hAZp[2*Size+IndexXZ(InitData.nStepsX,1)]=hAZp[IndexXZ(InitData.nStepsX,0)];
      tmpAXp=hAXp[IndexXZ(InitData.nStepsX-1,1)]+hAXp[Size+IndexXZ(InitData.nStepsX,1)];
      tmpAZm=hAZm[IndexXZ(InitData.nStepsX,2)]+hAZm[Size+IndexXZ(InitData.nStepsX,1)];
      hAXp[2*Size+IndexXZ(InitData.nStepsX,1)]=tmpAZm*Matrix2[1]+tmpAXp*Matrix2[0];
      hAZm[2*Size+IndexXZ(InitData.nStepsX,1)]=tmpAZm*Matrix2[0]+tmpAXp*Matrix2[1];
      */

      hAXm[2*Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)]=hAXm[IndexXZ(InitData.nStepsX+1,InitData.nStepsZ)];
      hAZm[2*Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)]=hAZm[IndexXZ(InitData.nStepsX,InitData.nStepsZ+1)];
      tmpAXp=hAXp[IndexXZ(InitData.nStepsX-1,InitData.nStepsZ)]+hAXp[Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)];
      tmpAZp=hAZp[IndexXZ(InitData.nStepsX,InitData.nStepsZ-1)]+hAZp[Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)];
      hAXp[2*Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)]=tmpAZp*Matrix2[1]+tmpAXp*Matrix2[0];
      hAZp[2*Size+IndexXZ(InitData.nStepsX,InitData.nStepsZ)]=tmpAZp*Matrix2[0]+tmpAXp*Matrix2[1];
    }

    /*
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
    */

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
        if(fDebug) {
          if(IndexXZ(ix,iz)==id4) {
            dd1=tmpAZp;
            dd2=tmpAXp*Matrix4[8]+tmpAXm*Matrix4[9];
            dd3=tmpAZm*Matrix4[11];
          }
          if(IndexXZ(ix,iz)==id5) {
            dd4=tmpAZp;
            dd5=tmpAXp*Matrix4[8]+tmpAXm*Matrix4[9];
            dd6=tmpAZm*Matrix4[11];
          }
          if(IndexXZ(ix,iz)==id6) {
            dd7=tmpAZp;
            dd8=tmpAXp*Matrix4[8]+tmpAXm*Matrix4[9];
            dd9=tmpAZm*Matrix4[11];
          }
        }
      }
    }
    //  int tmpi1,tmpi2,tmpi3;
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1; iz<=InitData.nStepsZB; iz++){
        hAZpMedia[2*SizeMedia+IndexXZ_M(ix,iz)]=hAZpMedia[SizeMedia+IndexXZ_M(ix,iz)]+
                                         hAZpMedia[IndexXZ_M(ix,iz-1)];
/*        tmpi1=2*SizeMedia+IndexXZ_M(ix,iz);
        tmpi2=IndexXZ_M(ix,iz+1);
        tmpi3= SizeMedia+IndexXZ_M(ix,iz);
*/
        hAZmMedia[2*SizeMedia+IndexXZ_M(ix,iz)]=hAZmMedia[SizeMedia+IndexXZ_M(ix,iz)]+
                                         hAZmMedia[IndexXZ_M(ix,iz+1)];
//        hAZmMedia[tmpi1]=hAZmMedia[tmpi3]+hAZmMedia[tmpi2];
      }

//      hAZpMedia[2*SizeMedia+IndexXZ_M(ix,1)]=hAZpMedia[SizeMedia+IndexXZ_M(ix,1)]+
//                                      ZRef*hAZmMedia[IndexXZ_M(ix,1)];
/*
      tmpi1=2*SizeMedia+IndexXZ_M(ix,InitData.nStepsZB);
      tmpi2=IndexXZ(ix,1);
      tmpi3= SizeMedia+IndexXZ_M(ix,InitData.nStepsZB);
*/
//      hAZmMedia[2*SizeMedia+IndexXZ_M(ix,InitData.nStepsZB)]=
//          hAZmMedia[SizeMedia+IndexXZ_M(ix,InitData.nStepsZB)]+hAZm[IndexXZ(ix,1)];
//      hAZmMedia[tmpi1]=hAZmMedia[tmpi3]+hAZm[tmpi2];
    }

    id1=IndexXZ(1,InitData.nStepsZ/3+1);
    id2=IndexXZ(2,InitData.nStepsZ/3+1);
    id3=IndexXZ(3,InitData.nStepsZ/3+1);
    id4=IndexXZ(1,InitData.nStepsZ/3);
    id5=IndexXZ(2,InitData.nStepsZ/3);
    id6=IndexXZ(3,InitData.nStepsZ/3);
    id7=IndexXZ(1,InitData.nStepsZ/3-1);
    id8=IndexXZ(2,InitData.nStepsZ/3-1);
    id9=IndexXZ(3,InitData.nStepsZ/3-1);

    if(fDebug) {
      dd1=hAZp[id4];
      dd2=hAZp[id1];
      dd3=hAZp[id5];
      dd4=hAZp[id2];
      dd5=hAXp[id4];
      dd6=hAXm[id4];
      dd7=hAXp[id5];
      dd8=hAXm[id5];
      //dd4=hAZp[IndexXZ(2,InitData.nStepsZ/3+1)];
      ddr1=abs(dd2)-abs(dd1); ddr2=abs(dd4)-abs(dd3); 
      ddr1/=InitData.dt; ddr2/=InitData.dt;
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

    memcpy(hAZpMedia+IndexXZ_M(ix1st,0),hAZpMedia+2*SizeMedia+IndexXZ_M(ix1st,0),
      (InitData.nStepsZB+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    memcpy(hAZmMedia+IndexXZ_M(ix1st,0),hAZmMedia+2*SizeMedia+IndexXZ_M(ix1st,0),
      (InitData.nStepsZB+2)*sizeof(complex<OPTIC_TYPE>)*(ixLast-ix1st+1));
    RungeStep=0;
  }
}

void cEquOpticBaseTask_2Way::StepNoMediaEuler() {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1),tmpAXp,tmpAXm,tmpAZp,tmpAZm;
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

void cEquOpticBaseTask_2Way::StepNoMediaRunge2() {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1),tmpAXp,tmpAXm,tmpAZp,tmpAZm;
  complex<OPTIC_TYPE> iiiA=iii*InitData.ACoef*polar(1,InitData.APhase);
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
        hAXp[Size+IndexXZ(ix,iz)]=dtT2*iiiA*
                  (hAZp[IndexXZ(ix-1,iz)]+hAZm[IndexXZ(ix-1,iz)]);
        hAXm[Size+IndexXZ(ix,iz)]=dtT2*iiiA*
                  (hAZp[IndexXZ(ix+1,iz)]+hAZm[IndexXZ(ix+1,iz)]);
        hAZp[Size+IndexXZ(ix,iz)]=dtT2*iiiA*
                  (hAXp[IndexXZ(ix,iz-1)]+hAXm[IndexXZ(ix,iz-1)]);
        hAZm[Size+IndexXZ(ix,iz)]=dtT2*iiiA*
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
        hAXp[2*Size+IndexXZ(ix,iz)]=dtT2*iiiA*(tmpAZp+tmpAZm)+
                             hAXp[Size+IndexXZ(ix,iz)]+hAXp[IndexXZ(ix-1,iz)];
        hAXm[2*Size+IndexXZ(ix,iz)]=dtT2*iiiA*(tmpAZp+tmpAZm)+
                             hAXm[Size+IndexXZ(ix,iz)]+hAXm[IndexXZ(ix+1,iz)];
        hAZp[2*Size+IndexXZ(ix,iz)]=dtT2*iiiA*(tmpAXp+tmpAXm)+
                             hAZp[Size+IndexXZ(ix,iz)]+hAZp[IndexXZ(ix,iz-1)];
        hAZm[2*Size+IndexXZ(ix,iz)]=dtT2*iiiA*(tmpAXp+tmpAXm)+
                             hAZm[Size+IndexXZ(ix,iz)]+hAZm[IndexXZ(ix,iz+1)];
        if(iz<100) {
          AXpk1Arr[iz]=hAXp[Size+IndexXZ(ix,iz)];
          AXpk2Arr[iz]=dtT2*iiiA*(tmpAZp+tmpAZm);
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
