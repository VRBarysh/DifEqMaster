//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquOpticTask_2BraggsSolidMedia.h"
#include <valarray>

//---------------------------------------------------------------------------

#pragma package(smart_init)

void cEquOpticTask_2BraggsSolidMedia::FillMatrixes() {
  OPTIC_TYPE dtT=InitData.dt;
  valarray< complex<OPTIC_TYPE> > VA[4];
  complex<OPTIC_TYPE> kk(0,-0.5*dtT*InitData.ACoef);
//  kk *= polar(1,InitData.APhase*M_PI);
  complex<OPTIC_TYPE> kk2(0,-0.5*dtT*InitData.ACoef2);
  kk2 *= polar(1,InitData.APhase*M_PI);

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

  for(int i=0;i<4;i++) {                                     //-------- 4x4 V2
    VA[i].resize(8,0);
    VA[i][i]=1;
    VA[i][i+4]=1;
  }
  VA[0][2]=VA[0][3]=VA[1][2]=VA[1][3]=VA[2][0]=VA[2][1]=VA[3][0]=VA[3][1]=kk2;

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
      Matrix4V2[i*4+j]=VA[i][4+j];
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
                                                               //-------- 3x3 V2
  for(int i=0;i<3;i++) {
    VA[i]=complex<OPTIC_TYPE>(0);
    VA[i][i]=1;
    VA[i][i+4]=1;
  }
  VA[0][2]=VA[1][2]=VA[2][0]=VA[2][1]=kk2;
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
      Matrix3V2[i*3+j]=VA[i][4+j];
                                                             //-------- 2x2
  complex<OPTIC_TYPE> tmp1mkk=OPTIC_TYPE(1)-kk*kk;
  Matrix2[0]=Matrix2[3]=OPTIC_TYPE(1)/tmp1mkk;
  Matrix2[1]=Matrix2[2]=-Matrix2[0]*kk;

  tmp1mkk=OPTIC_TYPE(1)-kk2*kk2;
  Matrix2V2[0]=Matrix2V2[3]=OPTIC_TYPE(1)/tmp1mkk;
  Matrix2V2[1]=Matrix2V2[2]=-Matrix2V2[0]*kk2;



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


void cEquOpticTask_2BraggsSolidMedia::StepNoMediaMatrix() {
  int Size=(InitData.nStepsX+2)*(InitData.nStepsZ+2);
  complex<OPTIC_TYPE> iii=complex<OPTIC_TYPE>(0,1),tmpAXp,tmpAXm,tmpAZp,tmpAZm;
//  complex<OPTIC_TYPE> AXpk1,AXpk2,AXmk1,AXmk2,AZpk1,AZpk2,AZmk1,AZmk2;
//  complex<OPTIC_TYPE> AXpk2Arr[100],AXpk1Arr[100],AXpk12Arr[100];
  OPTIC_TYPE dtT2=0.5*dtT;
  int iznb1=1,iznb2=InitData.nStepsZ,ixnb1=ix1st,ixnb2=ixLast;
  complex<OPTIC_TYPE> ACoef2 = InitData.ACoef;
  if(InitData.ACoef2) ACoef2 = polar(InitData.ACoef2,InitData.APhase);

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
        hAZp[Size+IndexXZ(ix,iz)]=0;
        hAZm[Size+IndexXZ(ix,iz)]=0;
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
        hAXp[2*Size+IndexXZ(1,iz)]=0;//hAXp[IndexXZ(0,iz)];
        hAXm[2*Size+IndexXZ(1,iz)]=0;//hAXm[IndexXZ(2,iz)];
        hAZp[2*Size+IndexXZ(1,iz)]=hAZp[IndexXZ(1,iz-1)];
        hAZm[2*Size+IndexXZ(1,iz)]=hAZm[IndexXZ(1,iz+1)];
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
        hAXp[2*Size+IndexXZ(InitData.nStepsX,iz)]=0;//hAXp[IndexXZ(0,iz)];
        hAXm[2*Size+IndexXZ(InitData.nStepsX,iz)]=0;//hAXm[IndexXZ(2,iz)];
        hAZp[2*Size+IndexXZ(InitData.nStepsX,iz)]=hAZp[IndexXZ(1,iz-1)];
        hAZm[2*Size+IndexXZ(InitData.nStepsX,iz)]=hAZm[IndexXZ(1,iz+1)];
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
        hAZp[2*Size+IndexXZ(ix,iz)]=hAZp[IndexXZ(ix,iz-1)];
        hAZm[2*Size+IndexXZ(ix,iz)]=hAZm[IndexXZ(ix,iz+1)];
        hAXp[2*Size+IndexXZ(ix,iz)]=0;
        hAXm[2*Size+IndexXZ(ix,iz)]=0;
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