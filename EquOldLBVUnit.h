//---------------------------------------------------------------------------

#ifndef EquationUnitH
#define EquationUnitH
//---------------------------------------------------------------------------
#include <math.h>
#include <complex.h>
#include <iostream.h>
#include <fstream.h>

#include "Dump25DElModel.h"
#include "Dump25DElForm.h"

#define  CLASS_SIMPLE_LOV          0
#define  CLASS_GAMMA_LOV           1
#define  CLASS_SIMPLE_LBV          2
#define  CLASS_GAMMA_LBV           3
#define  CLASS_25D_LOV             8
#define  CLASS_25D_LOV_GAMMA       9
#define  CLASS_25D_BASIC          10
#define  CLASS_25D_GAMMA          11
#define  CLASS_25D_FREE_LOV       (16+8+1)
#define  CLASS_25D_FREE_LBV       (16+8+2+1)

#define  REAL_TYPE       (long double)


class cEqu_OldLBVBase : public cDifEqu_DataBankIO
{
public:
  cEqu_OldLBVBase(){ hIdata=hAdata=hIdata1=hAdata1=hInterpolated=NULL; OutPoint=0;
                 ClassType=CLASS_SIMPLE_LOV; hDumpStorage=NULL; hElPresent=NULL;
                 hEldata=hElk1=hElk2=hElk3=hElk4=hElNew=NULL; nEls=0; nSteps=0;
                 hSplineCoefB=hSplineCoefC=hSplineCoefD=NULL; G=+1; SplineRequest=0; GammaCoef=1;}
  virtual ~cEqu_OldLBVBase()
        { delete [] hIdata; delete [] hIdata1;
          delete [] hAdata; delete [] hAdata1;
          delete [] hEldata; delete [] hElNew;
          delete [] hElk1; delete [] hElk2;
          delete [] hElk3; delete [] hElk4;
          delete [] hInterpolated; delete [] hSplineCoefB;
          delete [] hSplineCoefC;  delete [] hSplineCoefD;
          delete [] hElPresent; }

        virtual void LoadInitData(cDifEqu_InitDataRec *hInitData);
        void InitResData(int pResRange,REAL_TYPE *hA)
          { ResRange=pResRange; hResDataA=hA; }

        void PercentStep(ofstream *hOutFile=NULL);
        REAL_TYPE CalcMaxA();
        REAL_TYPE CalcMaxA_FirstPeak();
        virtual void AfterLastStep() {}
        virtual void Step();
        virtual void CheckIfElStopped(REAL_TYPE *hEl,int i) {}
        REAL_TYPE InputField(REAL_TYPE t)
          { return( t>TInput ? 0 : aInput*sin(M_PI*t*TInput_inv)*sin(M_PI*t*TInput_inv) ); }
        void cDiffurMain::Step_Euler();
        virtual complex<REAL_TYPE> Calc_I(REAL_TYPE * hEls);
        void StoreData(REAL_TYPE * hData,int type);
        void FillSplineArray(int ofs);
        REAL_TYPE kappa(REAL_TYPE z);

        void SPLINE_TEST(int nPoints,REAL_TYPE *hTestData);
        REAL_TYPE SPLINE_TEST_Y(REAL_TYPE x);

        virtual REAL_TYPE RightTH0(REAL_TYPE *hEl);
        virtual REAL_TYPE RightTH1(REAL_TYPE *hEl, complex<REAL_TYPE> A);
        virtual REAL_TYPE Gamma()
          {return((t>=delta) ? gamma0+RTH1Coef*alpha*(t-delta) : gamma0);}
        virtual void InitEls();

        virtual void SecondaryInit(){}
//        virtual void SecondaryInit_Gamma(REAL_TYPE pCPierce=0.001,REAL_TYPE pgamma0=3.0) {}
        virtual void SecondaryInit_Gamma(REAL_TYPE pCPierce=0.001,REAL_TYPE pgamma0=3.0);

        int GetElDroppedSum() { return(ElDroppedSum); }
        double GetElDroppedMaxFrac() { return(ElDroppedMaxFrac); }
        double GetElDroppedFrac() { return( (double)ElDroppedSum/(double)(nEls*(nCStep+1))); }
        double GetAvgGamma() {return(AvgGamma);}
        double GetAvgAvgGamma() {return(nGammaSteps ? AvgAvgGamma/nGammaSteps : 0);}

        void   SetGammaFileHandlers(ofstream *hAvg,ofstream *hMax)
          {hGammaFileAvg=hAvg; hGammaFileMax=hMax; }

        REAL_TYPE GetDebugValue(int nVal) { return(rDebug[nVal]); }

        void NewRequest(cDumpRequest *hRequest)
          { DumpRequest.NewRequest(hRequest); }

        int GetClassType() { return(ClassType); }

        REAL_TYPE r() { return( t<rTime ? rValue : 0 ); }
        REAL_TYPE GetAInput()  { if (aInput>0) return( t<TInput ? aInput : 0 );
          else return( t<TInput ? -aInput*sin(M_PI*t/TInput)*sin(M_PI*t/TInput) : 0 ); }
protected:
cDifEqu_InitDataRec_25DElModel InitData;

REAL_TYPE             dt,t,L,tMax,tMax_Res,rTime,rValue,delta,alpha,G,GammaCoef;
REAL_TYPE             TInput,TInput_inv,aInput,p,Z0,Z0_inv,Ckappa;
int                   nSteps,nEls,nMovingEls,nCStep,SplineRequest,OutPoint;
complex<REAL_TYPE>    *hIdata,*hAdata,*hIdata1,*hAdata1,*hInterpolated;
REAL_TYPE             *hEldata,*hElk1,*hElk2,*hElk3,*hElk4,*hElNew;
int                   *hElPresent;
REAL_TYPE             *hResDataA;
REAL_TYPE             *hSplineCoefB,*hSplineCoefC,*hSplineCoefD;
REAL_TYPE             ElDroppedMaxFrac;
double                AvgGamma,MaxGamma,AvgAvgGamma;
int                   nGammaSteps;
int                   ElDroppedSum;
int                   ResIndex,ResRange;
int                   PercentageCounter;
int                   j;
REAL_TYPE             backOutData;
cDumpRequest          DumpRequest;
int                   ClassType;
bool                  CalcPower;
bool                  bKappaLinear;

REAL_TYPE             rDebug[10];

REAL_TYPE gamma0,sqgamma0,CPierce,CPierce_inv,beta0,beta0_inv,RTH1Coef;
ofstream     *hGammaFileAvg,*hGammaFileMax;
};

class cDiffurMain_Gamma : public cDiffurMain
{
public:
  cDiffurMain_Gamma(){ hIdata=hAdata=hIdata1=hAdata1=hInterpolated=NULL; OutPoint=0;
                 ClassType=CLASS_GAMMA_LOV; hDumpStorage=NULL; hElPresent=NULL;
                 hEldata=hElk1=hElk2=hElk3=hElk4=hElNew=NULL; nEls=0; nSteps=0;
                 hSplineCoefB=hSplineCoefC=hSplineCoefD=NULL; G=+1; SplineRequest=0; GammaCoef=1;}
  cDiffurMain_Gamma(cDumpStorage *phDumpStorage)
               { hIdata=hAdata=hIdata1=hAdata1=hInterpolated=NULL; OutPoint=0;
                 ClassType=CLASS_GAMMA_LOV; hDumpStorage=phDumpStorage; hElPresent=NULL;
                 hEldata=hElk1=hElk2=hElk3=hElk4=hElNew=NULL; nEls=0; nSteps=0;
                 hSplineCoefB=hSplineCoefC=hSplineCoefD=NULL; G=+1; SplineRequest=0; GammaCoef=1;}
  ~cDiffurMain_Gamma()
        { delete [] hIdata; delete [] hIdata1;
          delete [] hAdata; delete [] hAdata1;
          delete [] hEldata; delete [] hElNew;
          delete [] hElk1; delete [] hElk2;
          delete [] hElk3; delete [] hElk4;
          delete [] hInterpolated; delete [] hSplineCoefB;
          delete [] hSplineCoefC;  delete [] hSplineCoefD;
          delete [] hElPresent; }
  virtual REAL_TYPE RightTH0(REAL_TYPE *hEl);
  virtual REAL_TYPE RightTH1(REAL_TYPE *hEl, complex<REAL_TYPE> A);
  virtual void InitEls();
  virtual void CheckIfElStopped(REAL_TYPE *hEl,int i);
  virtual void SecondaryInit_Gamma(REAL_TYPE pCPierce=0.001,REAL_TYPE pgamma0=3.0);

protected:
};

class cDiffurMain_LBV : public cDiffurMain
{
public:
  cDiffurMain_LBV(){ hIdata=hAdata=hIdata1=hAdata1=hInterpolated=NULL; OutPoint=0;
                 ClassType=CLASS_SIMPLE_LBV; hDumpStorage=NULL; hElPresent=NULL;
                 hEldata=hElk1=hElk2=hElk3=hElk4=hElNew=NULL; nEls=0; nSteps=0;
                 hSplineCoefB=hSplineCoefC=hSplineCoefD=NULL; G=-1; SplineRequest=0; GammaCoef=1;}
  cDiffurMain_LBV(cDumpStorage *phDumpStorage)
               { hIdata=hAdata=hIdata1=hAdata1=hInterpolated=NULL; OutPoint=0;
                 ClassType=CLASS_SIMPLE_LBV; hDumpStorage=phDumpStorage; hElPresent=NULL;
                 hEldata=hElk1=hElk2=hElk3=hElk4=hElNew=NULL; nEls=0; nSteps=0;
                 hSplineCoefB=hSplineCoefC=hSplineCoefD=NULL; G=-1; SplineRequest=0; GammaCoef=1;}
  virtual void Step();
  virtual void AfterLastStep();
  virtual void SecondaryInit()
          { G=1; OutPoint=nSteps; }

protected:
};

class cDiffurMain_LBV_Gamma : public cDiffurMain_LBV
{
public:

  cDiffurMain_LBV_Gamma(){ hIdata=hAdata=hIdata1=hAdata1=hInterpolated=NULL; OutPoint=0;
                 ClassType=CLASS_GAMMA_LBV; hDumpStorage=NULL; hElPresent=NULL;
                 hEldata=hElk1=hElk2=hElk3=hElk4=hElNew=NULL; nEls=0; nSteps=0;
                 hSplineCoefB=hSplineCoefC=hSplineCoefD=NULL; G=-1; SplineRequest=0; GammaCoef=1;}
  cDiffurMain_LBV_Gamma(cDumpStorage *phDumpStorage)
               { hIdata=hAdata=hIdata1=hAdata1=hInterpolated=NULL; OutPoint=0;
                 ClassType=CLASS_GAMMA_LBV; hDumpStorage=phDumpStorage; hElPresent=NULL;
                 hEldata=hElk1=hElk2=hElk3=hElk4=hElNew=NULL; nEls=0; nSteps=0;
                 hSplineCoefB=hSplineCoefC=hSplineCoefD=NULL; G=-1; SplineRequest=0; GammaCoef=1;}
                 
  virtual REAL_TYPE RightTH0(REAL_TYPE *hEl);
  virtual REAL_TYPE RightTH1(REAL_TYPE *hEl, complex<REAL_TYPE> A);
  virtual void InitEls();
  virtual void CheckIfElStopped(REAL_TYPE *hEl,int i);
  virtual void SecondaryInit_Gamma(REAL_TYPE pCPierce=0.001,REAL_TYPE pgamma0=3.0);
protected:
};

#endif
