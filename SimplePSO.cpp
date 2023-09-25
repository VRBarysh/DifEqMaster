//---------------------------------------------------------------------------


#pragma hdrstop

#include <stdlib.h>
#include "SimplePSO.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)
/*
inline PSOBaseType XAjust(PSOBaseType Val) {
  int iVal=floor(Val);
  double tmp1=Val-iVal;
  RevertSpeed=(iVal)&1;
  return (RevertSpeed)? 1.0-tmp1 : tmp1;
}
*/

void CPSOPoint::Init(CPSOInitData *hInit){
  PSOBaseType VC= hInit->VCoef > 0.000000000001 ? hInit->VCoef : 0.3;
  for(int i=0;i<Size;i++) {
    X[i]=hInit->hPointData[i].vMin +
         (hInit->hPointData[i].vMax-hInit->hPointData[i].vMin)*rand()/RAND_MAX;
    V[i]=VC*(hInit->hPointData[i].vMax-hInit->hPointData[i].vMin)*(2.0*(double)rand()/RAND_MAX-1);
    XMax[i]=X[i];
  }
  MaxVal=NO_VALUE;
  LifeTime=0;
}

void CPSOPoint::Park(TVecPSO *hMax,TVecPSO *hSpot) {
  for(int i=0;i<Size;i++) {
    X[i]=(*hMax)[i] - (*hSpot)[i]+ ((*hSpot)[i])*2.0*rand()/RAND_MAX;
    XMax[i]=X[i];
    V[i]=0;
  }
  MaxVal=NO_VALUE;
  LifeTime=0;
}

//---------------------------------------------------------------------------

void CPSOSimple::LoadInitData(CPSOInitData *hInit){
  if( (hInit->nPoints!=InitData.nPoints) ) {
    delete [] hPoints; hPoints= new CPSOPoint[hInit->nPoints];
  }
  InitData=*hInit;
  for(int i=0;i<InitData.nPoints;i++) {
   hPoints[i].SetSize(InitData.nComponents);
   hPoints[i].Init(hInit);
  }
  tmpArray1.resize(InitData.nComponents);
  tmpArray2.resize(InitData.nComponents);
  tmpArray3.resize(InitData.nComponents);
  BorderXMin.resize(InitData.nComponents);
  BorderXMax.resize(InitData.nComponents);
  AreaSize.resize(InitData.nComponents);
  GlobalMaxX.resize(InitData.nComponents);
  SpotSize.resize(InitData.nComponents);
  SpotSize2.resize(InitData.nComponents);
  P1Array.resize(InitData.nComponents,1);
  M1Array.resize(InitData.nComponents,-1);
  for (int i=0;i<InitData.nComponents;i++) {
    BorderXMin[i]=InitData.hPointData[i].vMin;
    BorderXMax[i]=InitData.hPointData[i].vMax;
    SpotSize[i]=InitData.SpotSize;
    SpotSize2[i]=InitData.SpotSize*0.5;
  }
  AreaSize=BorderXMax-BorderXMin;
  GlobalMaxX=hPoints[0].X;
  GlobalMaxVal=NO_VALUE;
  nStepCurrent=0; nPointCurrent=0;
  StoreData();
}

void CPSOSimple::LoadCoefData(CPSOInitData *hInit){
  InitData.Omega=hInit->Omega;
  InitData.c1=hInit->c1;
  InitData.c2=hInit->c2;
}

void CPSOSimple::StoreData(){
  for(int i=0;i<InitData.nComponents;i++)
    (*InitData.hPointData[i].hV)=hPoints[nPointCurrent].X[i];
}

void CPSOSimple::MoveX(int PointIndex) {  // valarray
  tmpArray1=hPoints[PointIndex].X + hPoints[PointIndex].V;
  tmpArray2=BorderXMax-tmpArray1;

  bool f1=tmpArray2.min()>0;
  tmpArray2=tmpArray1-BorderXMin;
  if( (f1) && (tmpArray2.min()>0) ) {
    hPoints[PointIndex].X=tmpArray1;
    return;
  }

  tmpArray2/=AreaSize;
  //tmpArray2.apply(XAjust);
  for(int i=0;i<InitData.nComponents;i++) {
    int iVal=floor(tmpArray2[i]);
    double tmp1=tmpArray2[i]-iVal;
    bool RevertSpeed=(iVal)&1;
    tmpArray2[i]= (RevertSpeed)? 1.0-tmp1 : tmp1;
    if (RevertSpeed) hPoints[PointIndex].V[i]=-hPoints[PointIndex].V[i];
  }
  hPoints[PointIndex].X=BorderXMin+AreaSize*tmpArray2;
}

void CPSOSimple::MovePoints() {
  for(int i=0;i<InitData.nPoints;i++) {
    for(int j=0;j<InitData.nComponents;j++) {
      tmpArray1[j]=InitData.c1*rand()/RAND_MAX;
      tmpArray2[j]=InitData.c2*rand()/RAND_MAX;
    }
    tmpArray3 = tmpArray2*(hPoints[i].XMax-hPoints[i].X) +
                tmpArray1*(GlobalMaxX-hPoints[i].X);
    MoveX(i);
    hPoints[i].V=hPoints[i].V*InitData.Omega + tmpArray3;
  }
}

void CPSOSimple::LoadFData(PSOBaseType NewValue) {
  int fArea; 
  if(NewValue > hPoints[nPointCurrent].MaxVal ) {
    hPoints[nPointCurrent].MaxVal=NewValue;
    hPoints[nPointCurrent].XMax=hPoints[nPointCurrent].X;
    hPoints[nPointCurrent].LifeTime=0;
    if ( NewValue > GlobalMaxVal ) {
      GlobalMaxVal=NewValue; GlobalMaxX=hPoints[nPointCurrent].X;
    }
  } else hPoints[nPointCurrent].LifeTime++;
  if( (InitData.MaxLifeTime) &&
      (InitData.MaxLifeTime<hPoints[nPointCurrent].LifeTime) ) {
    fArea=1;
    tmpArray1=hPoints[nPointCurrent].X-(GlobalMaxX-SpotSize);
    tmpArray2=hPoints[nPointCurrent].X-(GlobalMaxX+SpotSize);
    for(int j=0;j<InitData.nComponents;j++)
      if (!((tmpArray1[j]>0)&&(tmpArray2[j]<0))) fArea = 0;
    if (!fArea) hPoints[nPointCurrent].Init(&InitData);
      else hPoints[nPointCurrent].LifeTime=0;
  }
  nPointCurrent= (nPointCurrent + 1) % InitData.nPoints;
  if(!nPointCurrent) {
    nStepCurrent++;
    MovePoints();
  }
  StoreData();
}

double CPSOSimple::GetDoneFrac() {
  return( ((double)InitData.nPoints*nStepCurrent+nPointCurrent)/InitData.nSteps/InitData.nPoints);
}

PSOBaseType CPSOSimple::GetAvgSwarmSize() {
  PSOBaseType tmp=0;
  for(int i=0;i<InitData.nPoints;i++) {
    tmpArray1=hPoints[i].X-GlobalMaxX;
    tmpArray1*=tmpArray1;
    tmp+=tmpArray1.sum();
  }
  tmp/=InitData.nPoints*InitData.nComponents;
  return tmp;
}

