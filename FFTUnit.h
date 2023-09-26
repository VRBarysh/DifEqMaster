//---------------------------------------------------------------------------

#ifndef FFTUnitH
#define FFTUnitH

#include <complex>
//---------------------------------------------------------------------------
using namespace std;

template <class FFT_TYPE>
class CFFTrans
{
public:
  CFFTrans() { hRevTable=NULL; }
  CFFTrans(int piSMax) { hRevTable=NULL; MakeRevTable(piSMax); }
  ~CFFTrans() { delete hRevTable; }

  void MakeRevTable(int piSMax);
  void BitReverseCopy(complex<FFT_TYPE> *dest, complex<FFT_TYPE> *src);
  void TransForward(complex<FFT_TYPE> *dest, complex<FFT_TYPE> *src);
  void TransBackward(complex<FFT_TYPE> *dest, complex<FFT_TYPE> *src);

protected:
int *hRevTable;
int iSMax,iSize;
};

template <class FFT_TYPE>
void CFFTrans<FFT_TYPE>::MakeRevTable(int piSMax) {
  iSMax=piSMax; iSize=1 << iSMax;
  delete hRevTable;
  hRevTable = new int[ iSize ];
  int tmpBit,tmpVal;
  for(int i=0;i<iSize;i++) {
    hRevTable[i]=0; tmpVal=i;
    for(int j=0;j<iSMax;j++) {
      tmpBit=tmpVal&1;
      tmpVal>>=1;
      hRevTable[i]+= tmpBit << (iSMax-j-1);
    }
  }
}

template <class FFT_TYPE>
void CFFTrans<FFT_TYPE>::BitReverseCopy(complex<FFT_TYPE> *hDest, complex<FFT_TYPE> *hSrc){
  FFT_TYPE mul=1.0/sqrt((FFT_TYPE)iSize);
  for(int i=0;i<iSize;i++) hDest[hRevTable[i]]=hSrc[i]*mul;
}

template <class FFT_TYPE>
void CFFTrans<FFT_TYPE>::TransForward(complex<FFT_TYPE> *hDest, complex<FFT_TYPE> *hSrc) {
  BitReverseCopy(hDest,hSrc);
  int m;
  complex<FFT_TYPE> omega,omegaM,tmp1,tmp2;
  for(int i=1;i<=iSMax;i++) {
    m = 1 << i;
    omegaM=complex<FFT_TYPE>(cos(2*M_PI/m),sin(2*M_PI/m));
    omega=1;
    for(int j=0; j<(m>>1);j++) {
      for(int k=j;k<iSize-(m>>1);k+=m ) {
        tmp1=omega*hDest[k+(m>>1)];
        tmp2=hDest[k];
        hDest[k]=tmp1+tmp2;
        hDest[k+(m>>1)]=tmp2-tmp1;
      }
      omega*=omegaM;
    }
  }
}

template <class FFT_TYPE>
void CFFTrans<FFT_TYPE>::TransBackward(complex<FFT_TYPE> *hDest, complex<FFT_TYPE> *hSrc) {
  BitReverseCopy(hDest,hSrc);
  int m;
  complex<FFT_TYPE> omega,omegaM,tmp1,tmp2;
  for(int i=1;i<=iSMax;i++) {
    m = 1 << i;
    omegaM=complex<FFT_TYPE>(cos(2*M_PI/m),-sin(2*M_PI/m));
    omega=1;
    for(int j=0; j<(m>>1);j++) {
      for(int k=j;k<iSize-(m>>1);k+=m ) {
        tmp1=omega*hDest[k+(m>>1)];
        tmp2=hDest[k];
        hDest[k]=tmp1+tmp2;
        hDest[k+(m>>1)]=tmp2-tmp1;
      }
      omega*=omegaM;
    }
  }
}

#endif
