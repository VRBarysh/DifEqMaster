//---------------------------------------------------------------------------

#ifndef DumpStorageH
#define DumpStorageH
//---------------------------------------------------------------------------
#include <vcl.h>

class cDumpPlane {
public:
  cDumpPlane() { hData=NULL; SizeZ=SizeT=nStepT=0;}
  ~cDumpPlane() { delete [] hData;}
  MakePlane(int pSizeZ,pSizeT)
    { if( SizeZ*Size<pSizeZ*pSizeT ) {
        delete [] hData; hData=new long double[pSizeZ*pSizeT];
        SizeZ=pSizeZ; SizeT=pSizeT; } }

AnsiString Name;
long double *hData;
int SizeZ,SizeT,nStepT;
protected:
};

class cDumpStorage {
public:
  cDumpStorage() {}
  ~cDumpStorage() {}
protected:

};

#endif
