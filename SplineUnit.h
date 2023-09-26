//---------------------------------------------------------------------------

#ifndef SplineUnitH
#define SplineUnitH

#include <complex.h>
//---------------------------------------------------------------------------

using namespace std;

template <class T>
class CSpline {
public:
  CSpline() { hB=hC=hD=NULL; hcY=NULL; hX=hY=NULL; }
  void MakeSpline(int n,T *Y, T *B, T *C, T *D, T dx, T *X = NULL);
  void MakeSpline(int n,complex<T> *Y, T *B, T *C, T *D, T dx, T *X = NULL);
  inline T GetYIter(T x);
  inline complex<T> GetYIterC(T x);
  void IterReset() { index=0; }
protected:
int Size,index;
T *hB,*hC,*hD,*hX,*hY,d;
complex<T> *hcY;
};

template <class T>
T CSpline<T>::GetYIter(T x) {
  T tmpX;
  if(hX) {
    while((x>hX[index+1])&&(index<Size-1)) index++;
    tmpX=x-hX[index];
    return hY[index]+hB[index]*tmpX+hC[index]*tmpX*tmpX+hD[index]*tmpX*tmpX*tmpX;
  } else {
    while((x>d*(index+1))&&(index<Size-1)) index++;
    tmpX=x-d*index;
    return hY[index]+hB[index]*tmpX+hC[index]*tmpX*tmpX+hD[index]*tmpX*tmpX*tmpX;
  }
}

template <class T>
complex<T> CSpline<T>::GetYIterC(T x) {
  T tmpX;
  if(hX) {
    while((x>hX[index+1])&&(index<Size-1)) index++;
    tmpX=x-hX[index];
    return hcY[index] +
           complex<T>( hB[index]*tmpX+
                       hC[index]*tmpX*tmpX+
                       hD[index]*tmpX*tmpX*tmpX,
                       hB[Size+index]*tmpX+
                       hC[Size+index]*tmpX*tmpX+
                       hD[Size+index]*tmpX*tmpX*tmpX);
  } else {
    while((x>d*(index+1))&&(index<Size-1)) index++;
    tmpX=x-d*index;
    return hcY[index]+
           complex<T>( hB[index]*tmpX+
                       hC[index]*tmpX*tmpX+
                       hD[index]*tmpX*tmpX*tmpX,
                       hB[Size+index]*tmpX+
                       hC[Size+index]*tmpX*tmpX+
                       hD[Size+index]*tmpX*tmpX*tmpX);
  }
}

template <class T>
void CSpline<T>::MakeSpline(int n,T *Y, T *B, T *C, T *D, T dx, T *X) {
  Size=n; hY=Y; hX=X; hB=B; hC=C; hD=D; d=dx;
  T tmpV,d2x=dx*2,d3x=dx*3;
  int nm1=n-1,i1;
  index=0;
  if(Size<2) return;
  if(X) {
    if(Size<3) {
      B[0]=(Y[1]-Y[0])/(X[1]-X[0]);
      C[0]=0; D[0]=0; B[1]=B[0]; C[1]=0; D[0]=0;
      return;
    }
    D[0]=X[1]-X[0];
    C[1]=(Y[1]-Y[0])/D[0];
    for(int i=1;i<nm1;i++) {
      D[i]=X[i+1]-X[i];
      B[i]=2*(D[i-1]+D[i]);
      C[i+1]=(Y[i+1]-Y[i])/D[i];
      C[i]=C[i+1]-C[i];
    }
    B[0]=-D[0];
    B[nm1]=-D[n-2];
    C[0]=0;
    C[nm1]=0;
    if(n>3) {
      C[0]=C[2]/(X[3]-X[1])-C[1]/(X[2]-X[0]);
      C[nm1]=C[n-2]/(X[nm1]-X[n-3])-C[n-3]/(X[n-2]-X[n-4]);
      C[0]=C[0]*D[0]*D[0]/(X[3]-X[0]);
      C[nm1]=-C[nm1]*D[n-2]*D[n-2]/(X[nm1]-X[n-4]);
    }
    for(int i=1;i<n;i++) {
      tmpV=D[i-1]/B[i-1];
      B[i]=B[i]-tmpV*D[i-1];
      C[i]=C[i]-tmpV*C[i-1];
    }
    C[nm1]=C[nm1]/B[nm1];
    for(int i=n-2;i>=0;i--) {
      C[i]=(C[i]-D[i]*C[i+1])/B[i];
    }
    B[nm1]=(Y[nm1]-Y[n-2])/D[n-2]+D[n-2]*(C[n-2]+2*C[nm1]);
    for(int i=0;i<nm1;i++) {
      B[i]=(Y[i+1]-Y[i])/D[i]-D[i]*(C[i+1]+2*C[i]);
      D[i]=(C[i+1]-C[i])/D[i];
      C[i]=3*C[i];
    }
    C[nm1]=3*C[nm1];
    D[nm1]=D[n-2];
  } else {
    if(Size<3) {
      B[0]=(Y[1]-Y[0])/dx;
      C[0]=0; D[0]=0; B[1]=B[0]; C[1]=0; D[0]=0;
      return;
    }
    D[0]=dx;
    C[1]=(Y[1]-Y[0])/D[0];
    for(int i=1;i<nm1;i++) {
      D[i]=dx;
      B[i]=2*(D[i-1]+D[i]);
      C[i+1]=(Y[i+1]-Y[i])/D[i];
      C[i]=C[i+1]-C[i];
    }
    B[0]=-D[0];
    B[nm1]=-D[n-2];
    C[0]=0;
    C[nm1]=0;
    if(n>3) {
      C[0]=(C[2]-C[1])/d2x;
      C[nm1]=(C[n-2]-C[n-3])/d2x;
      C[0]=C[0]*D[0]*D[0]/d3x;
      C[nm1]=-C[nm1]*D[n-2]*D[n-2]/d3x;
    }
    for(int i=1;i<n;i++) {
      tmpV=D[i-1]/B[i-1];
      B[i]=B[i]-tmpV*D[i-1];
      C[i]=C[i]-tmpV*C[i-1];
    }
    C[nm1]=C[nm1]/B[nm1];
    for(int i=n-2;i>=0;i--) {
      C[i]=(C[i]-D[i]*C[i+1])/B[i];
    }
    B[nm1]=(Y[nm1]-Y[n-2])/D[n-2]+D[n-2]*(C[n-2]+2*C[nm1]);
    for(int i=0;i<nm1;i++) {
      B[i]=(Y[i+1]-Y[i])/D[i]-D[i]*(C[i+1]+2*C[i]);
      D[i]=(C[i+1]-C[i])/D[i];
      C[i]=3*C[i];
    }
    C[nm1]=3*C[nm1];
    D[nm1]=D[n-2];
  }
}

template <class T>
void CSpline<T>::MakeSpline(int n,complex<T> *Y,
                     T *B, T *C, T *D, T dx, T *X) {
  Size=n; hcY=Y; hX=X; hB=B; hC=C; hD=D; d=dx;
  T tmpV,d2x=dx*2,d3x=dx*3;
  int nm1=n-1,i1;
  index=0;
  if(Size<2) return;
  if(X) {                                              //  ------- X[i]

                                                       // ------ re
    D[0]=X[1]-X[0];
    C[1]=(Y[1].real()-Y[0].real())/D[0];
    for(int i=1;i<nm1;i++) {
      D[i]=X[i+1]-X[i];
      B[i]=2*(D[i-1]+D[i]);
      C[i+1]=(Y[i+1].real()-Y[i].real())/D[i];
      C[i]=C[i+1]-C[i];
    }
    B[0]=-D[0];
    B[nm1]=-D[n-2];
    C[0]=0;
    C[nm1]=0;
    if(n>3) {
      C[0]=C[2]/(X[3]-X[1])-C[1]/(X[2]-X[0]);
      C[nm1]=C[n-2]/(X[nm1]-X[n-3])-C[n-3]/(X[n-2]-X[n-4]);
      C[0]=C[0]*D[0]*D[0]/(X[3]-X[0]);
      C[nm1]=-C[nm1]*D[n-2]*D[n-2]/(X[nm1]-X[n-4]);
    }
    for(int i=1;i<n;i++) {
      tmpV=D[i-1]/B[i-1];
      B[i]=B[i]-tmpV*D[i-1];
      C[i]=C[i]-tmpV*C[i-1];
    }
    C[nm1]=C[nm1]/B[nm1];
    for(int i=n-2;i>=0;i--) {
      C[i]=(C[i]-D[i]*C[i+1])/B[i];
    }
    B[nm1]=(Y[nm1].real()-Y[n-2].real())/D[n-2]+D[n-2]*(C[n-2]+2*C[nm1]);
    for(int i=0;i<nm1;i++) {
      B[i]=(Y[i+1].real()-Y[i].real())/D[i]-D[i]*(C[i+1]+2*C[i]);
      D[i]=(C[i+1]-C[i])/D[i];
      C[i]=3*C[i];
    }
    C[nm1]=3*C[nm1];
    D[nm1]=D[n-2];

                                                       // ------ im
    D[Size]=X[1]-X[0];
    C[Size+1]=(Y[1].imag()-Y[0].imag())/D[Size];
    for(int i=1;i<nm1;i++) {
      D[Size+i]=X[i+1]-X[i];
      B[Size+i]=2*(D[Size+i-1]+D[Size+i]);
      C[Size+i+1]=(Y[i+1].imag()-Y[i].imag())/D[Size+i];
      C[Size+i]=C[Size+i+1]-C[Size+i];
    }
    B[Size]=-D[Size];
    B[Size+nm1]=-D[Size+n-2];
    C[Size]=0;
    C[Size+nm1]=0;
    if(n>3) {
      C[Size]=C[Size+2]/(X[3]-X[1])-C[Size+1]/(X[2]-X[0]);
      C[Size+nm1]=C[Size+n-2]/(X[nm1]-X[n-3])-C[Size+n-3]/(X[n-2]-X[n-4]);
      C[Size]=C[Size]*D[Size]*D[Size]/(X[3]-X[0]);
      C[Size+nm1]=-C[Size+nm1]*D[Size+n-2]*D[Size+n-2]/(X[nm1]-X[n-4]);
    }
    for(int i=Size+1;i<Size+n;i++) {
      tmpV=D[i-1]/B[i-1];
      B[i]=B[i]-tmpV*D[i-1];
      C[i]=C[i]-tmpV*C[i-1];
    }
    C[Size+nm1]=C[Size+nm1]/B[Size+nm1];
    for(int i=Size+n-2;i>=Size;i--) {
      C[i]=(C[i]-D[i]*C[i+1])/B[i];
    }
    B[Size+nm1]=(Y[nm1].imag()-Y[n-2].imag())/D[Size+n-2]+D[Size+n-2]*(C[Size+n-2]+2*C[Size+nm1]);
    for(int i=Size;i<Size+nm1;i++) {
      B[i]=(Y[i+1-Size].imag()-Y[i-Size].imag())/D[i]-D[i]*(C[i+1]+2*C[i]);
      D[i]=(C[i+1]-C[i])/D[i];
      C[i]=3*C[i];
    }
    C[Size+nm1]=3*C[Size+nm1];
    D[Size+nm1]=D[Size+n-2];

  } else {                                                  //  ----------- dx


    D[0]=dx;                                                // ------------ re
    C[1]=(Y[1].real()-Y[0].real())/D[0];
    for(int i=1;i<nm1;i++) {
      D[i]=dx;
      B[i]=2*(D[i-1]+D[i]);
      C[i+1]=(Y[i+1].real()-Y[i].real())/D[i];
      C[i]=C[i+1]-C[i];
    }
    B[0]=-D[0];
    B[nm1]=-D[n-2];
    C[0]=0;
    C[nm1]=0;
    if(n>3) {
      C[0]=(C[2]-C[1])/d2x;
      C[nm1]=(C[n-2]-C[n-3])/d2x;
      C[0]=C[0]*D[0]*D[0]/d3x;
      C[nm1]=-C[nm1]*D[n-2]*D[n-2]/d3x;
    }
    for(int i=1;i<n;i++) {
      tmpV=D[i-1]/B[i-1];
      B[i]=B[i]-tmpV*D[i-1];
      C[i]=C[i]-tmpV*C[i-1];
    }
    C[nm1]=C[nm1]/B[nm1];
    for(int i=n-2;i>=0;i--) {
      C[i]=(C[i]-D[i]*C[i+1])/B[i];
    }
    B[nm1]=(Y[nm1].real()-Y[n-2].real())/D[n-2]+D[n-2]*(C[n-2]+2*C[nm1]);
    for(int i=0;i<nm1;i++) {
      B[i]=(Y[i+1].real()-Y[i].real())/D[i]-D[i]*(C[i+1]+2*C[i]);
      D[i]=(C[i+1]-C[i])/D[i];
      C[i]=3*C[i];
    }
    C[nm1]=3*C[nm1];
    D[nm1]=D[n-2];

    D[Size]=dx;                                               // ---------- im
    C[Size+1]=(Y[1].imag()-Y[0].imag())/D[Size+0];
    for(int i=Size+1;i<Size+nm1;i++) {
      D[i]=dx;
      B[i]=2*(D[i-1]+D[i]);
      C[i+1]=(Y[i+1-Size].imag()-Y[i-Size].imag())/D[i];
      C[i]=C[i+1]-C[i];
    }
    B[Size]=-D[Size];
    B[Size+nm1]=-D[Size+n-2];
    C[Size]=0;
    C[Size+nm1]=0;
    if(n>3) {
      C[Size]=(C[Size+2]-C[Size+1])/d2x;
      C[Size+nm1]=(C[Size+n-2]-C[Size+n-3])/d2x;
      C[Size]=C[Size]*D[Size]*D[Size]/d3x;
      C[Size+nm1]=-C[Size+nm1]*D[Size+n-2]*D[Size+n-2]/d3x;
    }
    for(int i=Size+1;i<Size+n;i++) {
      tmpV=D[i-1]/B[i-1];
      B[i]=B[i]-tmpV*D[i-1];
      C[i]=C[i]-tmpV*C[i-1];
    }
    C[Size+nm1]=C[Size+nm1]/B[Size+nm1];
    for(int i=Size+n-2;i>=Size;i--) {
      C[i]=(C[i]-D[i]*C[i+1])/B[i];
    }
    B[Size+nm1]=(Y[nm1].imag()-Y[n-2].imag())/D[Size+n-2]+D[Size+n-2]*(C[Size+n-2]+2*C[Size+nm1]);
    for(int i=Size;i<Size+nm1;i++) {
      B[i]=(Y[i+1-Size].imag()-Y[i-Size].imag())/D[i]-D[i]*(C[i+1]+2*C[i]);
      D[i]=(C[i+1]-C[i])/D[i];
      C[i]=3*C[i];
    }
    C[Size+nm1]=3*C[Size+nm1];
    D[Size+nm1]=D[Size+n-2];
  }
}

#endif
