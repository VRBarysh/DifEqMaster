//---------------------------------------------------------------------------


#pragma hdrstop

#include <stdlib.h>
#include "MultiMaxPSO.h"

//---------------------------------------------------------------------------

void CPSOMultiMax::LoadInitData(CPSOInitData *hInit){
  if( (hInit->nPoints!=InitData.nPoints) ) {
    delete [] hPoints; hPoints= new CPSOPoint[hInit->nPoints];
  }
  if ((hInit->nGroups!=InitData.nGroups)) {
    delete [] hMaxPoint; hMaxPoint= new CPSOMaxPoint[hInit->nGroups];
    delete [] hMaxIndex; hMaxIndex= new int[hInit->nGroups];
  }
  delete [] hBlackList; hBlackList=new TVecPSO[MAXBLACKLISTSIZE]; nBlack=0;
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

void CPSOMultiMax::LoadCoefData(CPSOInitData *hInit){
  InitData.Omega=hInit->Omega;
  InitData.c1=hInit->c1;
  InitData.c2=hInit->c2;
}

void CPSOMultiMax::StoreData(){
  for(int i=0;i<InitData.nComponents;i++)
    (*InitData.hPointData[i].hV)=hPoints[nPointCurrent].X[i];
}

void CPSOMultiMax::BlackListUpdate() {
  if (!InitData.BlackTimeInterval) return;
  if ( ((nStepCurrent % InitData.BlackTimeInterval) == InitData.BlackTimeInterval - 1) &&
       (nPointCurrent==0) ) {
    for(int i=0;i<InitData.nGroups;i++) {
      if ( (hMaxPoint[i].Value<InitData.BlackListValue) && (nBlack<MAXBLACKLISTSIZE)) {
        hBlackList[nBlack].resize(InitData.nComponents);
        hBlackList[nBlack]=hMaxPoint[i].X;
        nBlack++;
        hMaxPoint[i].Value=NO_VALUE;
        for(int j=G1st(i);j<=GLast(i);j++) hPoints[j].Init(&InitData);
      }
    }
  }
}

int CPSOMultiMax::BlackListed(int nPoint) {
  int fSpot;
  if (!InitData.BlackTimeInterval) return 0;
  for (int i=0; i<nBlack; i++) {
    tmpArray1=hPoints[nPoint].X-(hBlackList[i]-InitData.BlackSpotSize);
    tmpArray2=hPoints[nPoint].X-(hBlackList[i]+InitData.BlackSpotSize);
    fSpot=1;
    for(int k=0;k<InitData.nComponents;k++)
      if (!((tmpArray1[k]>0)&&(tmpArray2[k]<0))) fSpot = 0;
    if(fSpot) return 1;
  }
  return(0);
}

void CPSOMultiMax::MoveX(int PointIndex) {  // valarray
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

void CPSOMultiMax::MovePoints() {
  for(int i=0;i<InitData.nPoints;i++) {
    if(InitData.c3m>0) {
      tmpArray4=0;
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

void CPSOMultiMax::LoadFData(PSOBaseType NewValue) {
  int fMax=-1,fArea,fValid=1,i1stDead=-1,nAliveGroups=InitData.nGroups,itmp1,fMaxAreaCalculated=0;
  if(NewValue > hPoints[nPointCurrent].MaxVal ) {
    if(InitData.LocalsInSpots) {
      for (int i=0; i<InitData.nGroups; i++) {
        tmpArray1=hPoints[nPointCurrent].X-(hMaxPoint[i].X-SpotSize);
        tmpArray2=hPoints[nPointCurrent].X-(hMaxPoint[i].X+SpotSize);
        hMaxPoint[i].fArea=1;
        for(int j=0;j<InitData.nComponents;j++)
          if (!((tmpArray1[j]>0)&&(tmpArray2[j]<0))) hMaxPoint[i].fArea = 0;
        if ( (hMaxPoint[i].fArea)&&(i!=Group(nPointCurrent)) ) fValid=0;
      }
      fMaxAreaCalculated=1;
    }
    if ( (fValid)&&(!BlackListed(nPointCurrent)) ) {
      hPoints[nPointCurrent].MaxVal=NewValue;
      hPoints[nPointCurrent].XMax=hPoints[nPointCurrent].X;
      hPoints[nPointCurrent].LifeTime=0;
    } else hPoints[nPointCurrent].LifeTime++;
  } else hPoints[nPointCurrent].LifeTime++;
  fValid=1;
  for (int i=0; (i<InitData.nGroups)&&(fMax==-1); i++)
    if (NewValue > hMaxPoint[hMaxIndex[i]].Value)
      fMax=i;
  if (fMax!=-1) {
    for (int i=0; i<InitData.nGroups; i++) {
      if(!fMaxAreaCalculated) {
        tmpArray1=hPoints[nPointCurrent].X-(hMaxPoint[hMaxIndex[i]].X-SpotSize);
        tmpArray2=hPoints[nPointCurrent].X-(hMaxPoint[hMaxIndex[i]].X+SpotSize);
        hMaxPoint[hMaxIndex[i]].fArea=1;
        for(int j=0;j<InitData.nComponents;j++)
          if (!((tmpArray1[j]>0)&&(tmpArray2[j]<0))) hMaxPoint[hMaxIndex[i]].fArea = 0;
      }
      if ( (hMaxPoint[hMaxIndex[i]].fArea)&&(i<fMax) ) fValid=0;
      if ( (hMaxPoint[hMaxIndex[i]].fArea)&&(i>=fMax)&&(i1stDead==-1) ) i1stDead=i;
    }
    if((fValid)&&(!BlackListed(nPointCurrent))) {
      i1stDead = i1stDead>-1 ? i1stDead : InitData.nGroups-1;
      itmp1=hMaxIndex[i1stDead];
      for (int i=i1stDead; i>fMax; i--) hMaxIndex[i]=hMaxIndex[i-1];
      hMaxIndex[fMax]=itmp1;
      hMaxPoint[hMaxIndex[fMax]].X=hPoints[nPointCurrent].X;
      hMaxPoint[hMaxIndex[fMax]].Value=NewValue;
      for(int i=0;i<InitData.nParked;i++) {
        hPoints[hMaxIndex[fMax]*InitData.nPoints/InitData.nGroups+i].Park(&hMaxPoint[hMaxIndex[fMax]].X,&SpotSize);
      }
      for (int i=i1stDead+1; i<nAliveGroups; i++)
        while ((hMaxPoint[hMaxIndex[i]].fArea)&&(i<nAliveGroups)) {
          itmp1=hMaxIndex[i];
          for(int j=i; j<nAliveGroups-1; j++) hMaxIndex[j]=hMaxIndex[j+1];
          hMaxIndex[nAliveGroups-1]=itmp1;
          nAliveGroups--;
        }
      for(int i=nAliveGroups; i<InitData.nGroups; i++) {
        hMaxPoint[hMaxIndex[i]].Value=NO_VALUE;
      }
    }
  }
  if( (InitData.MaxLifeTime) &&
      (InitData.MaxLifeTime<hPoints[nPointCurrent].LifeTime) ) {
    fArea=1;
    tmpArray1=hPoints[nPointCurrent].X-(hMaxPoint[Group(nPointCurrent)].X-SpotSize);
    tmpArray2=hPoints[nPointCurrent].X-(hMaxPoint[Group(nPointCurrent)].X+SpotSize);
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
  BlackListUpdate();
  StoreData();
}

double CPSOMultiMax::GetDoneFrac() {
  return( ((double)InitData.nPoints*nStepCurrent+nPointCurrent)/InitData.nSteps/InitData.nPoints);
}

PSOBaseType CPSOMultiMax::GetAvgSwarmSize() {
  PSOBaseType tmp=0;
  for(int i=0;i<InitData.nPoints;i++) {
    tmpArray1=hPoints[i].X-hMaxPoint[Group(i)].Value;
    tmpArray1*=tmpArray1;
    tmp+=tmpArray1.sum();
  }
  tmp/=InitData.nPoints*InitData.nComponents;
  return tmp;           
}

