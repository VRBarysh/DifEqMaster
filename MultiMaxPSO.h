//---------------------------------------------------------------------------

#ifndef MultiMaxPSOH
#define MultiMaxPSOH
//---------------------------------------------------------------------------
#define MAXBLACKLISTSIZE    1000

#include "SimplePSO.h"

class CPSOMaxPoint {
public:
  TVecPSO X;
  PSOBaseType Value;
  int fArea;
};

class CPSOMultiMax {
public:
  CPSOMultiMax(){ hPoints=NULL; nStepCurrent=0; nPointCurrent=0; nBlack=0;
                InitData.nSteps=0; InitData.nPoints=0; InitData.nComponents=0;
                InitData.nGroups=0; hMaxPoint=NULL; hMaxIndex=NULL; hBlackList=NULL; }
  ~CPSOMultiMax(){ delete [] hPoints; hPoints=NULL;
                   delete [] hMaxPoint; hMaxPoint=NULL; hMaxIndex=NULL;
                   delete [] hBlackList; hBlackList=NULL;}

  void LoadInitData(CPSOInitData *hInit);
  void LoadCoefData(CPSOInitData *hInit);
  void LoadFData(PSOBaseType NewValue);
  void MovePoints();
  void MoveX(int PointIndex);
  void BlackListUpdate();
  int BlackListed(int nPoint);
  PSOBaseType GetAvgSwarmSize();
  inline void StoreData();
  inline int Group(int nPoint) { return nPoint*InitData.nGroups/InitData.nPoints; }
  inline int G1st(int nGroup) { return nGroup*InitData.nPoints/InitData.nGroups; }
  inline int GLast(int nGroup) { return (nGroup+1)*InitData.nPoints/InitData.nGroups-1; }
  double GetDoneFrac();
  PSOBaseType GetGlobalMaxVal(int nMax) { return hMaxPoint[hMaxIndex[nMax]].Value; }
  void SetGlobalMaxVal(int nMax,PSOBaseType newMaxVal) { hMaxPoint[hMaxIndex[nMax]].Value=newMaxVal; }
  TVecPSO *GetGlobalMaxX(int nMax) { return &(hMaxPoint[hMaxIndex[nMax]].X); }
  TVecPSO *GetSpotSize() { return &(SpotSize); }

//private:
public:   // ------------------------------ Debug
CPSOInitData InitData;
CPSOPoint *hPoints;
CPSOMaxPoint *hMaxPoint;
TVecPSO *hBlackList;
TVecPSO BorderXMin,BorderXMax,AreaSize,SpotSize,SpotSize2;
TVecPSO tmpArray1,tmpArray2,tmpArray3,tmpArray4,P1Array,M1Array;
int nStepCurrent,nPointCurrent,nGroupCurrent,*hMaxIndex,nBlack;
};

#endif
