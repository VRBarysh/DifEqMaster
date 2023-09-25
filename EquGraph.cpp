//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquGraph.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

void cEquIOGraph::GetData(int Line,int dataSize, double *hData) {
  double c1,c2;
  int    i1,i2;
  Min=Max=Avg=0; if(!dataSize) { Size=0; return; }
  if ( (dataSize<=MemSize)&&(ShortMode==EQU_GRAPH_MODE_SHORT) ) {
    for(int i=0; i<dataSize; i++) Data[Line*MemSize+i]=hData[i]; Size=dataSize; return; }
  Size=MemSize;
  for(int i=0;i<MemSize-1;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    Data[Line*MemSize+i]=c1*hData[i1]+c2*hData[i2];
    if(i==0) Min=Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]>Max) Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]<Min) Min=Data[Line*MemSize+i];
    Avg+=Data[Line*MemSize+i];
  }
  Avg/=MemSize-1;
  Data[Line*MemSize+MemSize-1]=hData[dataSize-1];
}

void cEquIOGraph::GetData(int Line,int dataSize, long double *hData) {
  double c1,c2;
  int    i1,i2;
  Min=Max=Avg=0; if(!dataSize) { Size=0; return; }
  if ( (dataSize<=MemSize)&&(ShortMode==EQU_GRAPH_MODE_SHORT) ) {
    for(int i=0; i<dataSize; i++) Data[Line*MemSize+i]=hData[i]; Size=dataSize; return; }
  Size=MemSize;
  for(int i=0;i<MemSize-1;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    Data[Line*MemSize+i]=c1*hData[i1]+c2*hData[i2];
    if(i==0) Min=Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]>Max) Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]<Min) Min=Data[Line*MemSize+i];
    Avg+=Data[Line*MemSize+i];
  }
  Avg/=MemSize-1;
  Data[MemSize-1]=hData[dataSize-1];
}


void cEquIOGraph::GetDataAbs(int Line,int dataSize, complex<float> *hData) {
  double c1,c2;
  int    i1,i2;
  Min=Max=Avg=0; if(!dataSize) { Size=0; return; }
  if ( (dataSize<=MemSize)&&(ShortMode==EQU_GRAPH_MODE_SHORT) ) {
    for(int i=0; i<dataSize; i++) Data[Line*MemSize+i]=abs(hData[i]); Size=dataSize; return; }
  Size=MemSize;
  for(int i=0;i<MemSize-1;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    Data[Line*MemSize+i]=c1*abs(hData[i1])+c2*abs(hData[i2]);
    if(i==0) Min=Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]>Max) Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]<Min) Min=Data[Line*MemSize+i];
    Avg+=Data[Line*MemSize+i];
  }
  Avg/=MemSize-1;
  Data[Line*MemSize+MemSize-1]=abs(hData[dataSize-1]);
}

void cEquIOGraph::GetDataAbs(int Line,int dataSize, complex<double> *hData) {
  double c1,c2;
  int    i1,i2;
  Min=Max=Avg=0; if(!dataSize) { Size=0; return; }
  if ( (dataSize<=MemSize)&&(ShortMode==EQU_GRAPH_MODE_SHORT) ) {
    for(int i=0; i<dataSize; i++) Data[Line*MemSize+i]=abs(hData[i]); Size=dataSize; return; }
  Size=MemSize;
  for(int i=0;i<MemSize-1;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    Data[Line*MemSize+i]=c1*abs(hData[i1])+c2*abs(hData[i2]);
    if(i==0) Min=Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]>Max) Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]<Min) Min=Data[Line*MemSize+i];
    Avg+=Data[Line*MemSize+i];
  }
  Avg/=MemSize-1;
  Data[Line*MemSize+MemSize-1]=abs(hData[dataSize-1]);
}

void cEquIOGraph::GetDataAbs(int Line,int dataSize, complex<long double> *hData) {
  double c1,c2;
  int    i1,i2;
  long double tmp1,tmp2;
  Min=Max=Avg=0; if(!dataSize) { Size=0; return; }
  if ( (dataSize<=MemSize)&&(ShortMode==EQU_GRAPH_MODE_SHORT) ) {
    for(int i=0; i<dataSize; i++) Data[Line*MemSize+i]=abs(hData[i]); Size=dataSize; return; }
  Size=MemSize;
  for(int i=0;i<MemSize-1;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    tmp1=abs(hData[i1]);
    tmp2=abs(hData[i2]);
    Data[Line*MemSize+i]=c1*abs(hData[i1])+c2*abs(hData[i2]);
    if(i==0) Min=Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]>Max) Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]<Min) Min=Data[Line*MemSize+i];
    Avg+=Data[Line*MemSize+i];
  }
  Avg/=MemSize-1;
  Data[Line*MemSize+MemSize-1]=abs(hData[dataSize-1]);
}

void cEquIOGraph::GetDataArg(int Line,int dataSize, complex<float> *hData) {
  double c1,c2;
  double tmpA1,tmpA2;
  int    i1,i2;
  Min=Max=Avg=0; if(!dataSize) { Size=0; return; }
  if ( (dataSize<=MemSize)&&(ShortMode==EQU_GRAPH_MODE_SHORT) ) {
    for(int i=0; i<dataSize; i++) Data[Line*MemSize+i]=arg1(hData[i]); Size=dataSize; return; }
  Size=MemSize;
  for(int i=0;i<MemSize-1;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    tmpA1 = 0; if(i1>0) tmpA1 = arg1(hData[i1]*conj(hData[i1-1]));
    tmpA2 = 0; if(i2>0) tmpA2 = arg1(hData[i2]*conj(hData[i2-1]));
    Data[Line*MemSize+i]=c1*tmpA1+c2*tmpA2;
    if(i==0) Min=Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]>Max) Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]<Min) Min=Data[Line*MemSize+i];
    Avg+=Data[Line*MemSize+i];
  }
  Avg/=MemSize-1;
  Data[Line*MemSize+MemSize-1]=arg1(hData[dataSize-1]*conj(hData[dataSize-2]));
}

void cEquIOGraph::GetDataArg(int Line,int dataSize, complex<double> *hData) {
  double tmpA1,tmpA2;
  double c1,c2;
  int    i1,i2;
  Min=Max=Avg=0; if(!dataSize) { Size=0; return; }
  if ( (dataSize<=MemSize)&&(ShortMode==EQU_GRAPH_MODE_SHORT) ) {
    for(int i=0; i<dataSize; i++) Data[Line*MemSize+i]=arg1(hData[i]); Size=dataSize; return; }
  Size=MemSize;
  for(int i=0;i<MemSize-1;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    tmpA1 = 0; if(i1>0) tmpA1 = arg1(hData[i1]*conj(hData[i1-1]));
    tmpA2 = 0; if(i2>0) tmpA2 = arg1(hData[i2]*conj(hData[i2-1]));
    Data[Line*MemSize+i]=c1*tmpA1+c2*tmpA2;
    if(i==0) Min=Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]>Max) Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]<Min) Min=Data[Line*MemSize+i];
    Avg+=Data[Line*MemSize+i];
  }
  Avg/=MemSize-1;
  Data[Line*MemSize+MemSize-1]=arg1(hData[dataSize-1]*conj(hData[dataSize-2]));
}

void cEquIOGraph::GetDataArg(int Line,int dataSize, complex<long double> *hData) {
  double tmpA1,tmpA2;
  double c1,c2;
  int    i1,i2;
  long double tmp1,tmp2;
  Min=Max=Avg=0; if(!dataSize) { Size=0; return; }
  if ( (dataSize<=MemSize)&&(ShortMode==EQU_GRAPH_MODE_SHORT) ) {
    for(int i=0; i<dataSize; i++) Data[Line*MemSize+i]=arg1(hData[i]); Size=dataSize; return; }
  Size=MemSize;
  for(int i=0;i<MemSize-1;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    tmp1=arg1(hData[i1]);
    tmp2=arg1(hData[i2]);
    tmpA1 = 0; if(i1>0) tmpA1 = arg1(hData[i1]*conj(hData[i1-1]));
    tmpA2 = 0; if(i2>0) tmpA2 = arg1(hData[i2]*conj(hData[i2-1]));
    Data[Line*MemSize+i]=c1*tmpA1+c2*tmpA2;
    if(i==0) Min=Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]>Max) Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]<Min) Min=Data[Line*MemSize+i];
    Avg+=Data[Line*MemSize+i];
  }
  Avg/=MemSize-1;
  Data[Line*MemSize+MemSize-1]=arg1(hData[dataSize-1]*conj(hData[dataSize-2]));
}

void cEquIOGraph::GetDataRe(int Line,int dataSize, complex<float> *hData) {
  double c1,c2;
  int    i1,i2;
  Min=Max=Avg=0; if(!dataSize) { Size=0; return; }
  if ( (dataSize<=MemSize)&&(ShortMode==EQU_GRAPH_MODE_SHORT) ) {
    for(int i=0; i<dataSize; i++) Data[Line*MemSize+i]=hData[i].real(); Size=dataSize; return; }
  Size=MemSize;
  for(int i=0;i<MemSize-1;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    Data[Line*MemSize+i]=c1*hData[i1].real()+c2*hData[i2].real();
    if(i==0) Min=Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]>Max) Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]<Min) Min=Data[Line*MemSize+i];
    Avg+=Data[Line*MemSize+i];
  }
  Avg/=MemSize-1;
  Data[Line*MemSize+MemSize-1]=hData[dataSize-1].real();
}

void cEquIOGraph::GetDataRe(int Line,int dataSize, complex<double> *hData) {
  double c1,c2;
  int    i1,i2;
  Min=Max=Avg=0; if(!dataSize) { Size=0; return; }
  if ( (dataSize<=MemSize)&&(ShortMode==EQU_GRAPH_MODE_SHORT) ) {
    for(int i=0; i<dataSize; i++) Data[Line*MemSize+i]=hData[i].real(); Size=dataSize; return; }
  Size=MemSize;
  for(int i=0;i<MemSize-1;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    Data[Line*MemSize+i]=c1*hData[i1].real()+c2*hData[i2].real();
    if(i==0) Min=Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]>Max) Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]<Min) Min=Data[Line*MemSize+i];
    Avg+=Data[Line*MemSize+i];
  }
  Avg/=MemSize-1;
  Data[Line*MemSize+MemSize-1]=hData[dataSize-1].real();
}

void cEquIOGraph::GetDataRe(int Line,int dataSize, complex<long double> *hData) {
  double c1,c2;
  int    i1,i2;
  Min=Max=Avg=0; if(!dataSize) { Size=0; return; }
  if ( (dataSize<=MemSize)&&(ShortMode==EQU_GRAPH_MODE_SHORT) ) {
    for(int i=0; i<dataSize; i++) Data[Line*MemSize+i]=hData[i].real(); Size=dataSize; return; }
  Size=MemSize;
  for(int i=0;i<MemSize-1;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    Data[Line*MemSize+i]=c1*hData[i1].real()+c2*hData[i2].real();
    if(i==0) Min=Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]>Max) Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]<Min) Min=Data[Line*MemSize+i];
    Avg+=Data[Line*MemSize+i];
  }
  Avg/=MemSize-1;
  Data[Line*MemSize+MemSize-1]=hData[dataSize-1].real();
}

void cEquIOGraph::GetDataIm(int Line,int dataSize, complex<float> *hData) {
  double c1,c2;
  int    i1,i2;
  Min=Max=Avg=0; if(!dataSize) { Size=0; return; }
  if ( (dataSize<=MemSize)&&(ShortMode==EQU_GRAPH_MODE_SHORT) ) {
    for(int i=0; i<dataSize; i++) Data[Line*MemSize+i]=hData[i].imag(); Size=dataSize; return; }
  Size=MemSize;
  for(int i=0;i<MemSize-1;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    Data[Line*MemSize+i]=c1*hData[i1].imag()+c2*hData[i2].imag();
    if(i==0) Min=Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]>Max) Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]<Min) Min=Data[Line*MemSize+i];
    Avg+=Data[Line*MemSize+i];
  }
  Avg/=MemSize-1;
  Data[Line*MemSize+MemSize-1]=hData[dataSize-1].imag();
}

void cEquIOGraph::GetDataIm(int Line,int dataSize, complex<double> *hData) {
  double c1,c2;
  int    i1,i2;
  Min=Max=Avg=0; if(!dataSize) { Size=0; return; }
  if ( (dataSize<=MemSize)&&(ShortMode==EQU_GRAPH_MODE_SHORT) ) {
    for(int i=0; i<dataSize; i++) Data[Line*MemSize+i]=hData[i].imag(); Size=dataSize; return; }
  Size=MemSize;
  for(int i=0;i<MemSize-1;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    Data[Line*MemSize+i]=c1*hData[i1].imag()+c2*hData[i2].imag();
    if(i==0) Min=Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]>Max) Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]<Min) Min=Data[Line*MemSize+i];
    Avg+=Data[Line*MemSize+i];
  }
  Avg/=MemSize-1;
  Data[MemSize-1]=hData[dataSize-1].imag();
}

void cEquIOGraph::GetDataIm(int Line,int dataSize, complex<long double> *hData) {
  double c1,c2;
  int    i1,i2;
  Min=Max=Avg=0; if(!dataSize) { Size=0; return; }
  if ( (dataSize<=MemSize)&&(ShortMode==EQU_GRAPH_MODE_SHORT) ) {
    for(int i=0; i<dataSize; i++) Data[Line*MemSize+i]=hData[i].imag(); Size=dataSize; return; }
  Size=MemSize;
  for(int i=0;i<MemSize-1;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    Data[Line*MemSize+i]=c1*hData[i1].imag()+c2*hData[i2].imag();
    if(i==0) Min=Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]>Max) Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]<Min) Min=Data[Line*MemSize+i];
    Avg+=Data[Line*MemSize+i];
  }
  Avg/=MemSize-1;
  Data[Line*MemSize+MemSize-1]=hData[dataSize-1].imag();
}

void cEquIOGraph::GetDataSq(int Line,int dataSize, complex<float> *hData) {
  double c1,c2;
  int    i1,i2;
  Min=Max=Avg=0; if(!dataSize) { Size=0; return; }
  if ( (dataSize<=MemSize)&&(ShortMode==EQU_GRAPH_MODE_SHORT) ) {
    for(int i=0; i<dataSize; i++) Data[Line*MemSize+i]=hData[i].real()*hData[i].real()
                                  +hData[i].imag()*hData[i].imag();
  Size=dataSize; return; }
  Size=MemSize;
  for(int i=0;i<MemSize-1;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    Data[Line*MemSize+i]=c1*(hData[i1].real()*hData[i1].real()+hData[i1].imag()*hData[i1].imag())
           +c2*(hData[i2].real()*hData[i2].real()+hData[i2].imag()*hData[i2].imag());
    if(i==0) Min=Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]>Max) Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]<Min) Min=Data[Line*MemSize+i];
    Avg+=Data[Line*MemSize+i];
  }
  Avg/=MemSize-1;
  Data[Line*MemSize+MemSize-1]=hData[dataSize-1].real()*hData[dataSize-1].real()
                 +hData[dataSize-1].imag()*hData[dataSize-1].imag();
}

void cEquIOGraph::GetDataSq(int Line,int dataSize, complex<double> *hData) {
  double c1,c2;
  int    i1,i2;
  Min=Max=Avg=0; if(!dataSize) { Size=0; return; }
  if ( (dataSize<=MemSize)&&(ShortMode==EQU_GRAPH_MODE_SHORT) ) {
    for(int i=0; i<dataSize; i++) Data[Line*MemSize+i]=hData[i].real()*hData[i].real()
                                  +hData[i].imag()*hData[i].imag();
  Size=dataSize; return; }
  Size=MemSize;
  for(int i=0;i<MemSize-1;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    Data[Line*MemSize+i]=c1*(hData[i1].real()*hData[i1].real()+hData[i1].imag()*hData[i1].imag())
           +c2*(hData[i2].real()*hData[i2].real()+hData[i2].imag()*hData[i2].imag());
    if(i==0) Min=Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]>Max) Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]<Min) Min=Data[Line*MemSize+i];
    Avg+=Data[Line*MemSize+i];
  }
  Avg/=MemSize-1;
  Data[Line*MemSize+MemSize-1]=hData[dataSize-1].real()*hData[dataSize-1].real()
                 +hData[dataSize-1].imag()*hData[dataSize-1].imag();
}

void cEquIOGraph::GetDataSq(int Line,int dataSize, complex<long double> *hData) {
  double c1,c2;
  int    i1,i2;
  Min=Max=Avg=0; if(!dataSize) { Size=0; return; }
  if ( (dataSize<=MemSize)&&(ShortMode==EQU_GRAPH_MODE_SHORT) ) {
    for(int i=0; i<dataSize; i++) Data[Line*MemSize+i]=hData[i].real()*hData[i].real()
                                  +hData[i].imag()*hData[i].imag();
  Size=dataSize; return; }
  Size=MemSize;
  for(int i=0;i<MemSize-1;i++) {
    c2=(double)i*(dataSize-1)/(Size-1);
    i1=c2; i2=i1+1;
    c2-=i1; c1=1-c2;
    Data[Line*MemSize+i]=c1*(hData[i1].real()*hData[i1].real()+hData[i1].imag()*hData[i1].imag())
           +c2*(hData[i2].real()*hData[i2].real()+hData[i2].imag()*hData[i2].imag());
    if(i==0) Min=Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]>Max) Max=Data[Line*MemSize+i];
    if(Data[Line*MemSize+i]<Min) Min=Data[Line*MemSize+i];
    Avg+=Data[Line*MemSize+i];
  }
  Avg/=MemSize-1;
  Data[Line*MemSize+MemSize-1]=hData[dataSize-1].real()*hData[dataSize-1].real()
                 +hData[dataSize-1].imag()*hData[dataSize-1].imag();
}
