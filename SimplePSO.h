//---------------------------------------------------------------------------

#ifndef SimplePSOH
#define SimplePSOH
//---------------------------------------------------------------------------

#include <valarray>
#define NO_VALUE  (-1.0e+305) 

using namespace std;

typedef long double PSOBaseType;
typedef valarray<PSOBaseType> TVecPSO;

class CPSOInitComp {
public:
PSOBaseType *hV,vMin,vMax;
};

class CPSOInitData {
public:
  CPSOInitData() { hPointData=NULL; }
int nSteps,nPoints,nComponents,nGroups,nParked,MaxLifeTime,LocalsInSpots,BlackTimeInterval;
PSOBaseType Omega,c1,c2,VCoef,SpotSize,c3m,BlackSpotSize,BlackListValue;
CPSOInitComp *hPointData;
};

class CPSOPoint {
public:
  void SetSize(int size) { X.resize(size); V.resize(size); XMax.resize(size); Size=size; }
  void Init(CPSOInitData *hInit);
  void Park(TVecPSO *hMax,TVecPSO *hSpot);

  TVecPSO X,V,XMax;
  PSOBaseType MaxVal;
  int Size,LifeTime;
};

//valarray<double> V
//V.resize(10)

class CPSOSimple {
public:
  CPSOSimple(){ hPoints=NULL; nStepCurrent=0; nPointCurrent=0;
                InitData.nSteps=0; InitData.nPoints=0; InitData.nComponents=0;
                GlobalMaxVal=NO_VALUE; InitData.nGroups=0; }
  ~CPSOSimple(){ delete [] hPoints; hPoints=NULL; }
  
  void LoadInitData(CPSOInitData *hInit);
  void LoadCoefData(CPSOInitData *hInit);
  void LoadFData(PSOBaseType NewValue);
  void MovePoints();
  void MoveX(int PointIndex);
  PSOBaseType GetAvgSwarmSize();
  inline void StoreData();
  double GetDoneFrac();
  PSOBaseType GetGlobalMaxVal() { return GlobalMaxVal; }
  void SetGlobalMaxVal(PSOBaseType newMaxVal) { GlobalMaxVal=newMaxVal; }
  TVecPSO *GetGlobalMaxX() { return &GlobalMaxX; }

//private:
public:   // ------------------------------ Debug
CPSOInitData InitData;
CPSOPoint *hPoints;
TVecPSO GlobalMaxX,BorderXMin,BorderXMax,AreaSize,SpotSize,SpotSize2;
TVecPSO tmpArray1,tmpArray2,tmpArray3,P1Array,M1Array;
int nStepCurrent,nPointCurrent,nGroupCurrent;
PSOBaseType GlobalMaxVal;
};

#endif
