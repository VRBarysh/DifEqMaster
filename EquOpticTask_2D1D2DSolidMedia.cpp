//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquOpticTask_2D1D2DSolidMedia.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

void cEquOpticTask_2D1D2DSolidMedia::StepNoMediaMatrix() {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1),tmpAXp,tmpAXm,tmpAZp,tmpAZm;
//  complex<OPTIC_TYPE> AXpk1,AXpk2,AXmk1,AXmk2,AZpk1,AZpk2,AZmk1,AZmk2;
//  complex<OPTIC_TYPE> AXpk2Arr[100],AXpk1Arr[100],AXpk12Arr[100];
  OPTIC_TYPE dtT2=0.5*dtT;
  int iznb1=1,iznb2=InitData.nStepsZ,ixnb1=ix1st,ixnb2=ixLast;
  complex<OPTIC_TYPE> ACoef2 = InitData.ACoef;
  if(InitData.ACoef2) ACoef2 = polar(InitData.ACoef2,InitData.APhase);
  complex<OPTIC_TYPE> ACoef1D = polar(InitData.ACoef1D,InitData.APhase1D);
  complex<OPTIC_TYPE> ACoef1Dih2 = dtT2*iii*ACoef1D;
  complex<OPTIC_TYPE> ACoef1Dih2Construct = complex<OPTIC_TYPE>(1,0)/(complex<OPTIC_TYPE>(1,0)-dtT2*iii*dtT2*iii*InitData.ACoef1D*InitData.ACoef1D);

  ACoef2 = polar(InitData.ACoef2,InitData.APhase);  //  ------ DEBUG
  switch(RungeStep) {
    case 0:
    for(int ix=ix1st;ix<=ixLast;ix++) {
      for(int iz=1;iz<=InitData.nStepsZB;iz++) {
        hAXp[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*
                  (hAZp[IndexXZ(ix-1,iz)]+hAZm[IndexXZ(ix-1,iz)]);
        hAXm[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*
                  (hAZp[IndexXZ(ix+1,iz)]+hAZm[IndexXZ(ix+1,iz)]);
        hAZp[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*
                  (hAXp[IndexXZ(ix,iz-1)]+hAXm[IndexXZ(ix,iz-1)]);
        hAZm[Size+IndexXZ(ix,iz)]=dtT2*iii*InitData.ACoef*
                  (hAXp[IndexXZ(ix,iz+1)]+hAXm[IndexXZ(ix,iz+1)]);
      }

      for(int iz=InitData.nStepsZB+1; iz<=InitData.nStepsZ-InitData.nStepsZB2;iz++) {
//      for(int iz=InitData.nStepsZB+1; iz<=InitData.nStepsZ;iz++) {
        hAXp[Size+IndexXZ(ix,iz)]=0;
        hAXm[Size+IndexXZ(ix,iz)]=0;
        hAZp[Size+IndexXZ(ix,iz)]=dtT2*iii*ACoef1D*(hAZm[IndexXZ(ix,iz-1)])+dtT2*InitData.AGenCoef1D*hAZp[Size+IndexXZ(ix,iz-1)];
        hAZm[Size+IndexXZ(ix,iz)]=dtT2*iii*ACoef1D*(hAZp[IndexXZ(ix,iz+1)])+dtT2*InitData.AGenCoef1D*hAZm[Size+IndexXZ(ix,iz+1)];;
      }

      for(int iz=InitData.nStepsZ-InitData.nStepsZB2+1;iz<=InitData.nStepsZ;iz++) {
        hAXp[Size+IndexXZ(ix,iz)]=dtT2*iii*ACoef2*
                  (hAZp[IndexXZ(ix-1,iz)]+hAZm[IndexXZ(ix-1,iz)]);
        hAXm[Size+IndexXZ(ix,iz)]=dtT2*iii*ACoef2*
                  (hAZp[IndexXZ(ix+1,iz)]+hAZm[IndexXZ(ix+1,iz)]);
        hAZp[Size+IndexXZ(ix,iz)]=dtT2*iii*ACoef2*
                  (hAXp[IndexXZ(ix,iz-1)]+hAXm[IndexXZ(ix,iz-1)]);
        hAZm[Size+IndexXZ(ix,iz)]=dtT2*iii*ACoef2*
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
      for(int iz=iznb1; iz<=InitData.nStepsZB;iz++) {
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
      for(int iz=InitData.nStepsZB+1; iz<=InitData.nStepsZ-InitData.nStepsZB2;iz++) {
//      for(int iz=InitData.nStepsZB+1; iz<=InitData.nStepsZ;iz++) {
        tmpAZp=hAZp[IndexXZ(1,iz-1)]+hAZp[Size+IndexXZ(1,iz)];
        tmpAZm=hAZm[IndexXZ(1,iz+1)]+hAZm[Size+IndexXZ(1,iz)];
        hAXp[2*Size+IndexXZ(1,iz)]=0;//hAXp[IndexXZ(0,iz)];
        hAXm[2*Size+IndexXZ(1,iz)]=0;//hAXm[IndexXZ(2,iz)];
        hAZp[2*Size+IndexXZ(1,iz)]=(tmpAZp+tmpAZm*ACoef1Dih2)*ACoef1Dih2Construct*(1+InitData.AGenCoef1D*dtT2);
        hAZm[2*Size+IndexXZ(1,iz)]=(tmpAZm+tmpAZp*ACoef1Dih2)*ACoef1Dih2Construct*(1+InitData.AGenCoef1D*dtT2);
      }
      for(int iz=InitData.nStepsZ-InitData.nStepsZB2+1;iz<=iznb2;iz++) {
        tmpAXm=hAXm[IndexXZ(2,iz)]+hAXm[Size+IndexXZ(1,iz)];
        tmpAZp=hAZp[IndexXZ(1,iz-1)]+hAZp[Size+IndexXZ(1,iz)];
        tmpAZm=hAZm[IndexXZ(1,iz+1)]+hAZm[Size+IndexXZ(1,iz)];
        hAXp[2*Size+IndexXZ(1,iz)]=hAXp[IndexXZ(0,iz)];
        hAXm[2*Size+IndexXZ(1,iz)]=tmpAZp*Matrix3V2[6]+tmpAZm*Matrix3V2[7]+
                                   tmpAXm*Matrix3V2[8];
        hAZp[2*Size+IndexXZ(1,iz)]=tmpAZp*Matrix3V2[0]+tmpAZm*Matrix3V2[1]+
                                   tmpAXm*Matrix3V2[2];
        hAZm[2*Size+IndexXZ(1,iz)]=tmpAZp*Matrix3V2[3]+tmpAZm*Matrix3V2[4]+
                                   tmpAXm*Matrix3V2[5];
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
      hAXm[2*Size+IndexXZ(1,InitData.nStepsZ)]=tmpAZp*Matrix2V2[1]+tmpAXm*Matrix2V2[0];
      hAZp[2*Size+IndexXZ(1,InitData.nStepsZ)]=tmpAZp*Matrix2V2[0]+tmpAXm*Matrix2V2[1];
    }


    if( (!InitData.LockX)&&(ixLast==InitData.nStepsX) ) {
      ixnb2=InitData.nStepsX-1;
      for(int iz=iznb1; iz<=InitData.nStepsZB;iz++) {
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
      for(int iz=InitData.nStepsZB+1; iz<=InitData.nStepsZ-InitData.nStepsZB2;iz++) {
//      for(int iz=InitData.nStepsZB+1; iz<=InitData.nStepsZ;iz++) {
        tmpAZp=hAZp[IndexXZ(InitData.nStepsX,iz-1)]+hAZp[Size+IndexXZ(InitData.nStepsX,iz)];
        tmpAZm=hAZm[IndexXZ(InitData.nStepsX,iz+1)]+hAZm[Size+IndexXZ(InitData.nStepsX,iz)];
        hAXp[2*Size+IndexXZ(InitData.nStepsX,iz)]=0;//hAXp[IndexXZ(0,iz)];
        hAXm[2*Size+IndexXZ(InitData.nStepsX,iz)]=0;//hAXm[IndexXZ(2,iz)];
        hAZp[2*Size+IndexXZ(InitData.nStepsX,iz)]=(tmpAZp+tmpAZm*ACoef1Dih2)*ACoef1Dih2Construct*(1+InitData.AGenCoef1D*dtT2);
        hAZm[2*Size+IndexXZ(InitData.nStepsX,iz)]=(tmpAZm+tmpAZp*ACoef1Dih2)*ACoef1Dih2Construct*(1+InitData.AGenCoef1D*dtT2);
      }

      for(int iz=InitData.nStepsZ-InitData.nStepsZB2+1;iz<=iznb2;iz++) {
        tmpAXp=hAXp[IndexXZ(InitData.nStepsX-1,iz)]+hAXp[Size+IndexXZ(InitData.nStepsX,iz)];
        tmpAZp=hAZp[IndexXZ(InitData.nStepsX,iz-1)]+hAZp[Size+IndexXZ(InitData.nStepsX,iz)];
        tmpAZm=hAZm[IndexXZ(InitData.nStepsX,iz+1)]+hAZm[Size+IndexXZ(InitData.nStepsX,iz)];
        hAXm[2*Size+IndexXZ(InitData.nStepsX,iz)]=hAXm[IndexXZ(InitData.nStepsX+1,iz)];
        hAXp[2*Size+IndexXZ(InitData.nStepsX,iz)]=tmpAZp*Matrix3V2[6]+tmpAZm*Matrix3V2[7]+
                                   tmpAXp*Matrix3V2[8];
        hAZp[2*Size+IndexXZ(InitData.nStepsX,iz)]=tmpAZp*Matrix3V2[0]+tmpAZm*Matrix3V2[1]+
                                   tmpAXp*Matrix3V2[2];
        hAZm[2*Size+IndexXZ(InitData.nStepsX,iz)]=tmpAZp*Matrix3V2[3]+tmpAZm*Matrix3V2[4]+
                                   tmpAXp*Matrix3V2[5];
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
      hAXp[2*Size+IndexXZ(ix,InitData.nStepsZ)]=tmpAXp*Matrix3V2[0]+tmpAXm*Matrix3V2[1]+
                                  tmpAZp*Matrix3V2[2];
      hAXm[2*Size+IndexXZ(ix,InitData.nStepsZ)]=tmpAXp*Matrix3V2[3]+tmpAXm*Matrix3V2[4]+
                                  tmpAZp*Matrix3V2[5];
      hAZp[2*Size+IndexXZ(ix,InitData.nStepsZ)]=tmpAXp*Matrix3V2[6]+tmpAXm*Matrix3V2[7]+
                                  tmpAZp*Matrix3V2[8];
      hAZm[2*Size+IndexXZ(ix,InitData.nStepsZ)]=hAZm[IndexXZ(ix,InitData.nStepsZ+1)];
    }

    if( (!InitData.LockX)&&(InitData.Z0RefCoef) ) { iznb1=1; iznb2=InitData.nStepsZ; ixnb1=ix1st; ixnb2=ixLast; }

    for(int ix=ixnb1;ix<=ixnb2;ix++) {
      for(int iz=iznb1; iz<=InitData.nStepsZB;iz++) {
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
      for(int iz=InitData.nStepsZB+1; iz<=InitData.nStepsZ-InitData.nStepsZB2;iz++) {
//      for(int iz=InitData.nStepsZB+1; iz<=InitData.nStepsZ;iz++) {
        tmpAZp=hAZp[IndexXZ(ix,iz-1)]+hAZp[Size+IndexXZ(ix,iz)];
        tmpAZm=hAZm[IndexXZ(ix,iz+1)]+hAZm[Size+IndexXZ(ix,iz)];
        hAXp[2*Size+IndexXZ(ix,iz)]=0;//hAXp[IndexXZ(0,iz)];
        hAXm[2*Size+IndexXZ(ix,iz)]=0;//hAXm[IndexXZ(2,iz)];
        hAZp[2*Size+IndexXZ(ix,iz)]=(tmpAZp+tmpAZm*ACoef1Dih2)*ACoef1Dih2Construct*(1+InitData.AGenCoef1D*dtT2);
        hAZm[2*Size+IndexXZ(ix,iz)]=(tmpAZm+tmpAZp*ACoef1Dih2)*ACoef1Dih2Construct*(1+InitData.AGenCoef1D*dtT2);
      }

      for(int iz=InitData.nStepsZ-InitData.nStepsZB2+1;iz<=iznb2;iz++) {
        tmpAXp=hAXp[IndexXZ(ix-1,iz)]+hAXp[Size+IndexXZ(ix,iz)];
        tmpAXm=hAXm[IndexXZ(ix+1,iz)]+hAXm[Size+IndexXZ(ix,iz)];
        tmpAZp=hAZp[IndexXZ(ix,iz-1)]+hAZp[Size+IndexXZ(ix,iz)];
        tmpAZm=hAZm[IndexXZ(ix,iz+1)]+hAZm[Size+IndexXZ(ix,iz)];
        hAXp[2*Size+IndexXZ(ix,iz)]=tmpAXp*Matrix4V2[0]+tmpAXm*Matrix4V2[1]+
                                    tmpAZp*Matrix4V2[2]+tmpAZm*Matrix4V2[3];
        hAXm[2*Size+IndexXZ(ix,iz)]=tmpAXp*Matrix4V2[4]+tmpAXm*Matrix4V2[5]+
                                    tmpAZp*Matrix4V2[6]+tmpAZm*Matrix4V2[7];
        hAZp[2*Size+IndexXZ(ix,iz)]=tmpAXp*Matrix4V2[8]+tmpAXm*Matrix4V2[9]+
                                    tmpAZp*Matrix4V2[10]+tmpAZm*Matrix4V2[11];
        hAZm[2*Size+IndexXZ(ix,iz)]=tmpAXp*Matrix4V2[12]+tmpAXm*Matrix4V2[13]+
                                    tmpAZp*Matrix4V2[14]+tmpAZm*Matrix4V2[15];
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