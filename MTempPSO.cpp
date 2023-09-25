//---------------------------------------------------------------------------


#pragma hdrstop

#include "MTempPSO.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

void CPSOMultiMaxTmp::LoadInitData(CPSOInitData *hInit){
  if( (hInit->nPoints!=InitData.nPoints) ) {
    delete [] hPoints; hPoints= new CPSOPoint[hInit->nPoints];
  }
  if ((hInit->nGroups!=InitData.nGroups)) {
    delete [] hMaxPoint; hMaxPoint= new CPSOMaxPointTmp[hInit->nGroups];
    delete [] hMaxIndex; hMaxIndex= new int[hInit->nGroups];
  }
  InitData=*hInit;
  for(int i=0;i<InitData.nPoints;i++) {
   hPoints[i].SetSize(InitData.nComponents);
   hPoints[i].Init(hInit);
  }
  tmpArray1.resize(InitData.nComponents);
  tmpArray2.resize(InitData.nComponents);
  tmpArray3.resize(InitData.nComponents);
  tmpArray4.resize(InitData.nComponents);
  BorderXMin.resize(InitData.nComponents);
  BorderXMax.resize(InitData.nComponents);
  AreaSize.resize(InitData.nComponents);
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
  for(int i=0;i<InitData.nGroups;i++) {
    hMaxPoint[i].X.resize(InitData.nComponents);
    hMaxPoint[i].Value=NO_VALUE;
    hMaxPoint[i].X=hPoints[i*InitData.nPoints/InitData.nGroups].X;
    hMaxIndex[i]=i;
  }
  nStepCurrent=0; nPointCurrent=0;
  StoreData();
}

void CPSOMultiMaxTmp::LoadCoefData(CPSOInitData *hInit){
  InitData.Omega=hInit->Omega;
  InitData.c1=hInit->c1;
  InitData.c2=hInit->c2;
}

void CPSOMultiMaxTmp::StoreData(){
  for(int i=0;i<InitData.nComponents;i++)
    (*InitData.hPointData[i].hV)=hPoints[nPointCurrent].X[i];
}



double CPSOMultiMaxTmp::GetDoneFrac() {
  return( ((double)InitData.nPoints*nStepCurrent+nPointCurrent)/InitData.nSteps/InitData.nPoints);
}

PSOBaseType CPSOMultiMaxTmp::GetAvgSwarmSize() {
  PSOBaseType tmp=0;
  for(int i=0;i<InitData.nPoints;i++) {
    tmpArray1=hPoints[i].X-hMaxPoint[Group(i)].Value;
    tmpArray1*=tmpArray1;
    tmp+=tmpArray1.sum();
  }
  tmp/=InitData.nPoints*InitData.nComponents;
  return tmp;           
}

 void CPSOMultiMaxTmp::MoveX(int PointIndex) {  // valarray
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

void CPSOMultiMaxTmp::MovePoints() {
  for(int i=0;i<InitData.nPoints;i++) {
    tmpArray4=0;
    if(InitData.c3m>0) {
      for (int j=0; j<InitData.nGroups; j++) {
        tmpArray1=hPoints[i].X-(hMaxPoint[j].X-SpotSize2);
        tmpArray2=hPoints[i].X-(hMaxPoint[j].X+SpotSize2);
        hMaxPoint[j].fArea=1;
        for(int k=0;k<InitData.nComponents;k++)
          if (!((tmpArray1[k]>0)&&(tmpArray2[k]<0))) hMaxPoint[j].fArea = 0;
        if((hMaxPoint[j].fArea)&&(Group(i)!=j)) {
          tmpArray1=(hPoints[i].X-hMaxPoint[j].X);
          for(int k=0;k<InitData.nComponents;k++)
            tmpArray2[k]=tmpArray1[k]>0 ? SpotSize2[k]-tmpArray1[k] : -SpotSize2[k]-tmpArray1[k];
          tmpArray4+=InitData.c3m*tmpArray2;
        }
      }
    }
    for(int j=0;j<InitData.nComponents;j++) {
      tmpArray1[j]=InitData.c1*rand()/RAND_MAX;
      tmpArray2[j]=InitData.c2*rand()/RAND_MAX;
    }
    tmpArray3 = tmpArray2*(hPoints[i].XMax-hPoints[i].X) +
                tmpArray1*(hMaxPoint[Group(i)].X-hPoints[i].X)+tmpArray4;
    MoveX(i);
    hPoints[i].V=hPoints[i].V*InitData.Omega + tmpArray3;
  }
}

void CPSOMultiMaxTmp::LoadFData(PSOBaseType NewValue) {
  int fArea;
  if(NewValue > hPoints[nPointCurrent].MaxVal ) {
    hPoints[nPointCurrent].MaxVal=NewValue;
    hPoints[nPointCurrent].XMax=hPoints[nPointCurrent].X;
    hPoints[nPointCurrent].LifeTime=0;
    /*
    if ( NewValue > GlobalMaxVal ) {
      GlobalMaxVal=NewValue; GlobalMaxX=hPoints[nPointCurrent].X;
    } */
  } else hPoints[nPointCurrent].LifeTime++;
  if( (InitData.MaxLifeTime) &&
      (InitData.MaxLifeTime<hPoints[nPointCurrent].LifeTime) ) {
    fArea=1;
    tmpArray1=hPoints[nPointCurrent].X;//-(GlobalMaxX-SpotSize);
    tmpArray2=hPoints[nPointCurrent].X;//-(GlobalMaxX+SpotSize);
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