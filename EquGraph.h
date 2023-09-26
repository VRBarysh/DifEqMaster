//---------------------------------------------------------------------------

#ifndef EquGraphH
#define EquGraphH

#define EQU_GRAPH_SIZE                  600
#define EQU_GRAPH_MODE_SCALE            0
#define EQU_GRAPH_MODE_SHORT            1

#include <complex.h>
//---------------------------------------------------------------------------

using namespace std;

class cEquIOGraph
{
public:
  cEquIOGraph() {Size=0; ShortMode=0; Min=Max=0; MemSize=0; Number2D=-1;
                 Data = NULL; }
  cEquIOGraph(int pMemSize ) {Size=0; ShortMode=0; Min=Max=0; Number2D=-1;
                    MemSize=pMemSize; Data = new double[MemSize];}
  ~cEquIOGraph() { delete [] Data; }

  void SetShortMode(int NewMode) { ShortMode=NewMode; }
  void Init(int pnLines,int pMemSize)
    { if(!MemSize) { nLines=pnLines; MemSize=pMemSize;
         Data = new double[nLines*MemSize]; FillZero();} Number2D=-1; }
  void FillZero() {if(Data) for(int i=0;i<MemSize;i++) Data[i]=0;}

  inline double arg1(complex<float> v)
    { double q=abs(v); return (q<0.000001)? 0 : arg(v); }
  inline double arg1(complex<double> v)
    { double q=abs(v); return (q<0.00000000001)? 0 : arg(v); }
  inline double arg1(complex<long double> v)
    { double q=abs(v); return (q<0.00000000001)? 0 : arg(v); }

  void GetData(int Line,int dataSize, double *hData);
  void GetData(int Line,int dataSize, long double *hData);

  void GetDataAbs(int Line,int dataSize, complex<float> *hData);
  void GetDataAbs(int Line,int dataSize, complex<double> *hData);
  void GetDataAbs(int Line,int dataSize, complex<long double> *hData);

  void GetDataArg(int Line,int dataSize, complex<float> *hData);
  void GetDataArg(int Line,int dataSize, complex<double> *hData);
  void GetDataArg(int Line,int dataSize, complex<long double> *hData);

  void GetDataRe(int Line,int dataSize, complex<float> *hData);
  void GetDataRe(int Line,int dataSize, complex<double> *hData);
  void GetDataRe(int Line,int dataSize, complex<long double> *hData);

  void GetDataIm(int Line,int dataSize, complex<float> *hData);
  void GetDataIm(int Line,int dataSize, complex<double> *hData);
  void GetDataIm(int Line,int dataSize, complex<long double> *hData);

  void GetDataSq(int Line,int dataSize, complex<float> *hData);
  void GetDataSq(int Line,int dataSize, complex<double> *hData);
  void GetDataSq(int Line,int dataSize, complex<long double> *hData);

  double Value(int Line,int index) { return(Data[Line*MemSize+index]); }
  double *Data;
  double Min,Max,Avg;
  int Size, ShortMode, MemSize, nLines, Number2D;
protected:
};

#endif
