//---------------------------------------------------------------------------

#ifndef MTempPSOH
#define MTempPSOH
//---------------------------------------------------------------------------


#include "SimplePSO.h"

class CPSOMaxPointTmp {
public:
  TVecPSO X;
  PSOBaseType Value;
  int fArea;
};

class CPSOMultiMaxTmp {
public:
  CPSOMultiMaxTmp(){ hPoints=NULL; nStepCurrent=0; nPointCurrent=0;
                InitData.nSteps=0; InitData.nPoints=0; InitData.nComponents=0;
                InitData.nGroups=0; hMaxPoint=NULL; hMaxIndex=NULL; }
  ~CPSOMultiMaxTmp(){ delete [] hPoints; hPoints=NULL;
                   delete [] hMaxPoint; hMaxPoint=NULL; hMaxIndex=NULL; }

  void LoadInitData(CPSOInitData *hInit);
  void LoadCoefData(CPSOInitData *hInit);
  void LoadFData(PSOBaseType NewValue);
  void MovePoints();
  void MoveX(int PointIndex);
  PSOBaseType GetAvgSwarmSize();
  inline void StoreData();
  inline int Group(int nPoint) { return nPoint*InitData.nGroups/InitData.nPoints; }
  double GetDoneFrac();
  PSOBaseType GetGlobalMaxVal(int nMax) { return hMaxPoint[hMaxIndex[nMax]].Value; }
  void SetGlobalMaxVal(int nMax,PSOBaseType newMaxVal) { hMaxPoint[hMaxIndex[nMax]].Value=newMaxVal; }
  TVecPSO *GetGlobalMaxX(int nMax) { return &(hMaxPoint[hMaxIndex[nMax]].X); }
  TVecPSO *GetSpotSize() { return &(SpotSize); }

//private:
public:   // ------------------------------ Debug
CPSOInitData InitData;
CPSOPoint *hPoints;
CPSOMaxPointTmp *hMaxPoint;
TVecPSO BorderXMin,BorderXMax,AreaSize,SpotSize,SpotSize2;
TVecPSO tmpArray1,tmpArray2,tmpArray3,tmpArray4,P1Array,M1Array;
int nStepCurrent,nPointCurrent,nGroupCurrent,*hMaxIndex;
};

#endif
