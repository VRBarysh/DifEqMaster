//---------------------------------------------------------------------------

#ifndef EquOptic2D01H
#define EquOptic2D01H

#include <complex.h>
#include <vcl\math.hpp>
#include "EquBaseMaster.h"

typedef double OPTIC_TYPE;
//---------------------------------------------------------------------------

class cDifEqu_InitDataRec_Optic2D : public cDifEqu_InitDataRec
{
public:
int nStepsX,nStepsZ,SizeX,SizeZ,nStepsT,MediaSize,FFTSize,FFTScale,DumpInt,DumpReq,nAHarm,nStepsXMedia,nStepsXBragg;
int nStepsZB,nStepsZB2,SizeZB,MediaSizeB,nStepsDelay,DelayPos,fNoMedia,fTE,nStepsX2,fForceLinear,f2BraggsSolidMedia,f2D1D2DSolidMedia;
int nTSubSteps,UseSimpleEC,UseFastP,LockX,UseMirrorDelay,Coaxial,nMediaRows,nMediaRowSize,nMedia;
int nThreads;
int PurgeNoiseDelay;
long double rPolNoise,dz,PCoef,RCoef,beta,ACoef,A0,A1,A2,AGenCoef,PCoefFast,Q,AGenZ,AGenX,ABackCoef,RCoef2;
long double ACoef1D,APhase1D,AGenCoef1D;
long double Z0RefCoef,Z0RefPhase,Z1RefCoef,DiffractionCoef,EWrong,StartingInv,QGrid,ACoef2,APhase;
long double MaxV[5],SaveBitmapInterval,BitmapNorm,MicrowaveModeLength,Delta,DeltaCoef,MediaGaussCoef;
};

template<class fT>
class cDifEqu_InitDataRec_Optic2D_T : public cDifEqu_InitDataRec
{
public:
int nStepsX,nStepsZ,SizeX,SizeZ,nStepsT,MediaSize,FFTSize,FFTScale,DumpInt,DumpReq,nAHarm,nStepsXMedia,nStepsXBragg;
int nStepsZB,nStepsZB2,SizeZB,MediaSizeB,nStepsDelay,DelayPos,fNoMedia,fTE,nStepsX2,fForceLinear,f2BraggsSolidMedia,f2D1D2DSolidMedia;
int nTSubSteps,UseSimpleEC,UseFastP,LockX,UseMirrorDelay,Coaxial,nMediaRows,nMediaRowSize,nMedia;
int nThreads;
int PurgeNoiseDelay;
fT rPolNoise,dz,PCoef,RCoef,beta,ACoef,A0,A1,A2,AGenCoef,PCoefFast,Q,AGenZ,AGenX,ABackCoef,RCoef2;
fT ACoef1D,APhase1D,AGenCoef1D;
fT Z0RefCoef,Z0RefPhase,Z1RefCoef,DiffractionCoef,EWrong,StartingInv,QGrid,ACoef2,APhase;
fT MaxV[5],SaveBitmapInterval,BitmapNorm,MicrowaveModeLength,Delta,DeltaCoef,MediaGaussCoef;
};

class cOptic2D_Media
{
public:
  cOptic2D_Media() {}
  ~cOptic2D_Media() {}

//  cOptic2D_Media operator+(cOptic2D_Media &m1);
  cOptic2D_Media operator+(const cOptic2D_Media &m1);
  void operator+=(const cOptic2D_Media &m1);
  void operator+=(const cOptic2D_Media m1);
  cOptic2D_Media operator*(long double m);

  void MakeK(cOptic2D_Media *hk, cOptic2D_Media *hkNew, long double coef);
  void MakeNew(cOptic2D_Media *hBase, cOptic2D_Media *hkArray, cOptic2D_Media *hNew);

  void MakeNewSeveralSteps(cOptic2D_Media *hNew,const int i,const int j);
  void MakeNewSimple(cOptic2D_Media *hNew,const int i,const int j);
  void MakeNewSimpleEC(cOptic2D_Media *hNew,const int i,const int j);
  void MakeNewSimpleEC1(cOptic2D_Media *hNew,const int i,const int j);

cOptic2D_Media RightPart();
cOptic2D_Media RightPartMedia();
inline TCplxLong RightAZp();
inline TCplxLong RightAZm();
inline TCplxLong RightAXp();
inline TCplxLong RightAXm();
inline TCplxLong RightPZp();
inline TCplxLong RightPZm();
inline TCplxLong RightPXp();
inline TCplxLong RightPXm();
inline TCplxLong RightPZpShort();
inline TCplxLong RightPZmShort();
inline TCplxLong RightPXpShort();
inline TCplxLong RightPXmShort();
inline long double RightR0();
inline TCplxLong RightRZ();
inline TCplxLong RightRX();
inline TCplxLong RightRp();
inline TCplxLong RightRm();

inline TCplxLong RightPZpSimple();
inline TCplxLong RightPZmSimple();
inline TCplxLong RightPXpSimple();
inline TCplxLong RightPXmSimple();
inline TCplxLong RightPZpShortSimple();
inline TCplxLong RightPZmShortSimple();
inline TCplxLong RightPXpShortSimple();
inline TCplxLong RightPXmShortSimple();
inline long double RightR0Simple();
inline TCplxLong RightRZSimple();
inline TCplxLong RightRXSimple();
inline TCplxLong RightRpSimple();
inline TCplxLong RightRmSimple();


  TCplxLong AZp,AZm,AXp,AXm,PZp,PZm,PXp,PXm,RZ,RX,Rp,Rm; // r2z, r2x, rz+x, rz-x
  long double R0;   // r0
  cDifEqu_InitDataRec_Optic2D *hInitData;
};

class cOptic2D_Media_Dump
{
public:
  cOptic2D_Media_Dump () {}
  ~cOptic2D_Media_Dump () {}
  void GetData(cOptic2D_Media *hSource);
  
  float AZp,AZm,AXp,AXm,PZp,PZm,PXp,PXm;
  float R0,RZ,RX,Rp,Rm;   // r0, r2z, r2x, rz+x, rz-x
};

class cOptic2D_AbsStore {
public:
  void GetData(cOptic2D_Media *hSource) {
    AZp=abs(hSource->AZp); AZm=abs(hSource->AZm);
    AXp=abs(hSource->AXp); AXm=abs(hSource->AXm); }
  float AZp,AZm,AXp,AXm,PZp,PZm,PXp,PXm;
};

#endif
