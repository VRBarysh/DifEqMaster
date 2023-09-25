//---------------------------------------------------------------------------


#pragma hdrstop

#include <io.h>
#include <fcntl.h>
#include <sys\stat.h>
#include "Dump2DOpticModel.h"
#include "CephesLib.h"


//---------------------------------------------------------------------------

#pragma package(smart_init)

cDump2DOpticModel::cDump2DOpticModel()// : cEquBaseDump()
{
  DumpEls=0; Init.nStepsX=Init.nStepsZ=Init.nStepsT=0;
  hMediaDump=NULL;
  DumpChanged=1;
}

void cDump2DOpticModel::LoadInitData(cDifEqu_InitDataRec_Optic2D *hInit) {
  hInit->nStepsT=hInit->tMax/hInit->dt;
  DumpSizeX=hInit->nStepsX/hInit->DumpInt;
  DumpSizeZ=hInit->nStepsX/hInit->DumpInt;
  DumpSizeT=hInit->nStepsT/hInit->DumpInt;
  delete [] hMediaDump;
  hMediaDump = new cOptic2D_Media_Dump[DumpSizeX*DumpSizeZ*DumpSizeT];
  Init=*hInit;
  Step=0;
  DumpEnabled=DumpReady=1;
}

void cDump2DOpticModel::Dump(cOptic2D_Media *hMedia) {
  if(!DumpEnabled) return;
  if(Step/Init.DumpInt>Init.nStepsT-1) return;
  if(Step%Init.DumpInt==0)
    for(int i=0;i<DumpSizeZ;i++)
      for(int j=0;j<DumpSizeX;j++)
        hMediaDump[i*DumpSizeX+j].GetData(hMedia+Init.SizeX*i*Init.DumpInt+j*Init.DumpInt);
  Step++;
}
/*
TCplxLong cDump2DOpticModel::GetAVal(int ZStep,int XStep,int TStep) {
  TCplxLong dI,Res=0;
  TCplxLong alpha,beta;
  TCplxLong ACoef=TCplxLong(0,1);
  ACoef=-sqrt(ACoef/(long double)M_PI);
  for(int i=0;i<ZStep;i++) {
    dI=hI[(Init.nSteps+1)*TStep+ZStep-i]-hI[(Init.nSteps+1)*TStep+ZStep-i-1];
    alpha=(hI[(Init.nSteps+1)*TStep+ZStep-i]+dI*(long double)i);
    beta=dI/Init.dt;
    Res+=hXProfileA[i*nXSteps+XStep]*alpha-hXProfileB[i*nXSteps+XStep]*beta;
  }
  Res*=ACoef;
  return(Res);
}

void cDump2DOpticModel::MakeXProfileCoef(long double pdx, int pnXSteps) {
  if( (dx==pdx) && (nXSteps==pnXSteps) ) return;
  if(pnXSteps!=nXSteps) {
    delete [] hXProfileA; delete [] hXProfileB;
    hXProfileA=new TCplxLong[pnXSteps*Init.nSteps];
    hXProfileB=new TCplxLong[pnXSteps*Init.nSteps];
  }
  dx=pdx; nXSteps=pnXSteps;
  long double tmpRe=0,tmpIm=0;
  long double sqPI=sqrt(M_PI);
  long double sqdt=sqrt(Init.dt),sq2=sqrt((long double)2);
  long double x,x3,a,cs,sn,d1,d2,sqd1,sqd2,cs1,cs2,sn1,sn2;
  double tmpFC,tmpFS,tmpFC2,tmpFS2;
  for(int j=0;j<nXSteps;j++) {
    x=dx*j; x3=x*x*x; a=-0.25*x*x/Init.dt; sn=sin(a); cs=cos(a);
    fresnl(sq2*x/(2*sqPI*sqdt),&tmpFS,&tmpFC);
    tmpRe=-0.5*x*sq2*sqPI+2*sqdt*cs+x*sq2*sqPI*tmpFS;
    tmpIm=-0.5*x*sq2*sqPI+2*sqdt*sn+x*sq2*sqPI*tmpFC;
    hXProfileA[j]=TCplxLong(tmpRe,tmpIm);
    tmpRe=-0.5*x3*sq2*sqPI+4*sqdt*sqdt*sqdt*cs+2*x*x*sqdt*sn+x3*sq2*sqPI*tmpFC;
    tmpIm=0.5*x3*sq2*sqPI+4*sqdt*sqdt*sqdt*sn-2*x*x*sqdt*cs-x3*sq2*sqPI*tmpFS;
    tmpIm/=6;  tmpRe/=6;
    hXProfileB[j]=TCplxLong(tmpRe,tmpIm);
    for(int i=1;i<Init.nSteps;i++) {
      d1=(i+1)*Init.dt; sqd1=sqrt(d1); d2=i*Init.dt; sqd2=sqrt(d2);
      a=-0.25*x*x/d1; cs1=cos(a); sn1=sin(a);
      a=0.25*x*x/d2; cs2=cos(a); sn2=sin(a);
      fresnl(sq2*x/(2*sqPI*sqd1),&tmpFS,&tmpFC);
      fresnl(sq2*x/(2*sqPI*sqd2),&tmpFS2,&tmpFC2);
      tmpRe=2*sqd1*cs1+x*sq2*sqPI*tmpFS-2*sqd2*cs2-x*sq2*sqPI*tmpFS2;
      tmpIm=2*sqd1*sn1+x*sq2*sqPI*tmpFC+2*sqd2*sn2-x*sq2*sqPI*tmpFC2;
      hXProfileA[i*nXSteps+j]=TCplxLong(tmpRe,tmpIm);
      tmpRe=4*sqd1*sqd1*sqd1*cs1+2*x*x*sqd1*sn1+x3*sq2*sqPI*tmpFC
           -4*sqd2*sqd2*sqd2*cs2+2*x*x*sqd2*sn2-x3*sq2*sqPI*tmpFC2;
      tmpIm=4*sqd2*sqd2*sqd2*sn2+2*x*x*sqd2*cs2+x3*sq2*sqPI*tmpFS2
           +4*sqd1*sqd1*sqd1*sn1-2*x*x*sqd1*cs1-x3*sq2*sqPI*tmpFS;
      tmpIm/=6;  tmpRe/=6;
      hXProfileB[i*nXSteps+j]=TCplxLong(tmpRe,tmpIm);
    }
  }
}

void cDump2DOpticModel::DrawZTSlice(Graphics::TBitmap *hBitmap,cGr3DInitData Gr3DInit,DrawDumpInitData DrawInit){
  long double x1,y1,x2,y2,x3,y3,x4,y4,q1,q2;
  TPoint Poly[4];
  cGr3DMaster Gr3D=cGr3DMaster(Gr3DInit);
  hBitmap->Canvas->Brush->Color=clWhite;
  hBitmap->Canvas->FillRect(TRect(0,0,640,480));
  hBitmap->Canvas->Brush->Color=clBlue;
  for(int i=Init.nSteps;i>0;i--)
    for(int j=Step-1;j>0;j--) {
      x1=Gr3D.GetX(i*Init.dt,j*Init.dt,abs(hA[j*(Init.nSteps+1)+i]));
      y1=Gr3D.GetY(i*Init.dt,j*Init.dt,abs(hA[j*(Init.nSteps+1)+i]));
      x2=Gr3D.GetX((i-1)*Init.dt,j*Init.dt,abs(hA[j*(Init.nSteps+1)+(i-1)]));
      y2=Gr3D.GetY((i-1)*Init.dt,j*Init.dt,abs(hA[j*(Init.nSteps+1)+(i-1)]));
      x3=Gr3D.GetX((i-1)*Init.dt,(j-1)*Init.dt,abs(hA[(j-1)*(Init.nSteps+1)+(i-1)]));
      y3=Gr3D.GetY((i-1)*Init.dt,(j-1)*Init.dt,abs(hA[(j-1)*(Init.nSteps+1)+(i-1)]));
      x4=Gr3D.GetX(i*Init.dt,(j-1)*Init.dt,abs(hA[(j-1)*(Init.nSteps+1)+i]));
      y4=Gr3D.GetY(i*Init.dt,(j-1)*Init.dt,abs(hA[(j-1)*(Init.nSteps+1)+i]));
      Poly[0]=TPoint(x1,y1);
      Poly[1]=TPoint(x2,y2);
      Poly[2]=TPoint(x3,y3);
      Poly[3]=TPoint(x4,y4);
      q1=(x1-x2)*(y4-y2)-(x4-x2)*(y1-y2);
      q2=(x3-x2)*(y4-y2)-(x4-x2)*(y3-y2);
      hBitmap->Canvas->Brush->Color= (q1-q2)>0 ? clLime : clBlue;
      if(DrawInit.DrawLines) {
        hBitmap->Canvas->Pen->Color=clBlack;
        hBitmap->Canvas->MoveTo(x1,y1);
        hBitmap->Canvas->LineTo(x2,y2);
        hBitmap->Canvas->LineTo(x3,y3);
      } else hBitmap->Canvas->Polygon(Poly,4);
    }
}

void cDump2DOpticModel::DrawXTSlice(Graphics::TBitmap *hBitmap,cGr3DInitData Gr3DInit,DrawDumpInitData DrawInit){
  long double x1,y1,x2,y2,x3,y3,x4,y4,q1,q2;
  TPoint Poly[5];
  int zStep=DrawInit.t/Init.dt;
  cGr3DMaster Gr3D=cGr3DMaster(Gr3DInit);
  hBitmap->Canvas->Brush->Color=clWhite;
  hBitmap->Canvas->FillRect(TRect(0,0,640,480));
  if(zStep>Init.nSteps) return;
  hBitmap->Canvas->Brush->Color=clBlue;
  if ( (DrawInit.dx!=dx) || (DrawInit.nXSteps!=nXSteps))
                         MakeXProfileCoef(DrawInit.dx,DrawInit.nXSteps);
  int iBack=nXSteps-1,jBack;
  for(int i=nXSteps-2;i>=0;i=i-DrawInit.Skip > 0 ? i-DrawInit.Skip : (i>0?0:-1) ) {
    jBack=Step-1;
    for(int j=Step-2;j>=0;j=j-DrawInit.Skip > 0 ? j-DrawInit.Skip : (j>0?0:-1)) {
      x1=Gr3D.GetX(iBack*dx,jBack*Init.dt,abs(GetAVal(zStep,iBack,jBack)));
      y1=Gr3D.GetY(iBack*dx,jBack*Init.dt,abs(GetAVal(zStep,iBack,jBack)));
      x2=Gr3D.GetX(i*dx,jBack*Init.dt,abs(GetAVal(zStep,i,jBack)));
      y2=Gr3D.GetY(i*dx,jBack*Init.dt,abs(GetAVal(zStep,i,jBack)));
      x3=Gr3D.GetX(i*dx,j*Init.dt,abs(GetAVal(zStep,i,j)));
      y3=Gr3D.GetY(i*dx,j*Init.dt,abs(GetAVal(zStep,i,j)));
      x4=Gr3D.GetX(iBack*dx,j*Init.dt,abs(GetAVal(zStep,iBack,j)));
      y4=Gr3D.GetY(iBack*dx,j*Init.dt,abs(GetAVal(zStep,iBack,j)));
      Poly[0]=TPoint(x1,y1);
      Poly[1]=TPoint(x2,y2);
      Poly[2]=TPoint(x3,y3);
      Poly[3]=TPoint(x4,y4);
      q1=(x1-x2)*(y4-y2)-(x4-x2)*(y1-y2);
      q2=(x3-x2)*(y4-y2)-(x4-x2)*(y3-y2);
      hBitmap->Canvas->Brush->Color= (q1-q2)>0 ? clLime : clBlue;
      if(DrawInit.DrawLines) {
        hBitmap->Canvas->Pen->Color=clBlack;
        hBitmap->Canvas->MoveTo(x1,y1);
        hBitmap->Canvas->LineTo(x2,y2);
        hBitmap->Canvas->LineTo(x3,y3);
      } else {
        hBitmap->Canvas->Polygon(Poly,3);
      }
      jBack=j;
    }
    iBack=i;
  }
  if(DrawInit.DoubleSided) {
    iBack=0;
    for(int i=1;i<=nXSteps-1;i=i+DrawInit.Skip < nXSteps-1 ? i+DrawInit.Skip : (i<nXSteps-1?nXSteps-1:nXSteps) ) {
      jBack=0;
      for(int j=1;j<=Step-1;j=j+DrawInit.Skip < Step-1 ? j+DrawInit.Skip : (j<Step-1?Step-1:Step)) {
        x1=Gr3D.GetX(-iBack*dx,jBack*Init.dt,abs(GetAVal(zStep,iBack,jBack)));
        y1=Gr3D.GetY(-iBack*dx,jBack*Init.dt,abs(GetAVal(zStep,iBack,jBack)));
        x2=Gr3D.GetX(-i*dx,jBack*Init.dt,abs(GetAVal(zStep,i,jBack)));
        y2=Gr3D.GetY(-i*dx,jBack*Init.dt,abs(GetAVal(zStep,i,jBack)));
        x3=Gr3D.GetX(-i*dx,j*Init.dt,abs(GetAVal(zStep,i,j)));
        y3=Gr3D.GetY(-i*dx,j*Init.dt,abs(GetAVal(zStep,i,j)));
        x4=Gr3D.GetX(-iBack*dx,j*Init.dt,abs(GetAVal(zStep,iBack,j)));
        y4=Gr3D.GetY(-iBack*dx,j*Init.dt,abs(GetAVal(zStep,iBack,j)));
        Poly[0]=TPoint(x1,y1);
        Poly[1]=TPoint(x2,y2);
        Poly[2]=TPoint(x3,y3);
        Poly[3]=TPoint(x4,y4);
        q1=(x1-x2)*(y4-y2)-(x4-x2)*(y1-y2);
        q2=(x3-x2)*(y4-y2)-(x4-x2)*(y3-y2);
        hBitmap->Canvas->Brush->Color= (q1-q2)>0 ? clLime : clBlue;
        if(DrawInit.DrawLines) {
          hBitmap->Canvas->Pen->Color=clBlack;
          hBitmap->Canvas->MoveTo(x1,y1);
          hBitmap->Canvas->LineTo(x2,y2);
          hBitmap->Canvas->LineTo(x3,y3);
        } else hBitmap->Canvas->Polygon(Poly,3);
        jBack=j;
      }
      iBack=i;
    }
  }
}

void cDump2DOpticModel::DrawXZSlice(Graphics::TBitmap *hBitmap,cGr3DInitData Gr3DInit,DrawDumpInitData DrawInit){
  long double x1,y1,x2,y2,x3,y3,x4,y4,q1,q2;
  TPoint Poly[5];
  int tStep=DrawInit.t/Init.dt;
  cGr3DMaster Gr3D=cGr3DMaster(Gr3DInit);
  hBitmap->Canvas->Brush->Color=clWhite;
  hBitmap->Canvas->FillRect(TRect(0,0,640,480));
  if(tStep>Step-1) return;
  hBitmap->Canvas->Brush->Color=clBlue;
  if ( (DrawInit.dx!=dx) || (DrawInit.nXSteps!=nXSteps))
                         MakeXProfileCoef(DrawInit.dx,DrawInit.nXSteps);
  int iBack=nXSteps-1,jBack;
  for(int i=nXSteps-2;i>=0;i=i-DrawInit.Skip > 0 ? i-DrawInit.Skip : (i>0?0:-1) ) {
    jBack=Init.nSteps-1;
    for(int j=Init.nSteps-2;j>=0;j=j-DrawInit.Skip > 0 ? j-DrawInit.Skip : (j>0?0:-1)) {
      x1=Gr3D.GetX(iBack*dx,jBack*Init.dt,abs(GetAVal(jBack,iBack,tStep)));
      y1=Gr3D.GetY(iBack*dx,jBack*Init.dt,abs(GetAVal(jBack,iBack,tStep)));
      x2=Gr3D.GetX(i*dx,jBack*Init.dt,abs(GetAVal(jBack,i,tStep)));
      y2=Gr3D.GetY(i*dx,jBack*Init.dt,abs(GetAVal(jBack,i,tStep)));
      x3=Gr3D.GetX(i*dx,j*Init.dt,abs(GetAVal(j,i,tStep)));
      y3=Gr3D.GetY(i*dx,j*Init.dt,abs(GetAVal(j,i,tStep)));
      x4=Gr3D.GetX(iBack*dx,j*Init.dt,abs(GetAVal(j,iBack,tStep)));
      y4=Gr3D.GetY(iBack*dx,j*Init.dt,abs(GetAVal(j,iBack,tStep)));
      Poly[0]=TPoint(x1,y1);
      Poly[1]=TPoint(x2,y2);
      Poly[2]=TPoint(x3,y3);
      Poly[3]=TPoint(x4,y4);
      q1=(x1-x2)*(y4-y2)-(x4-x2)*(y1-y2);
      q2=(x3-x2)*(y4-y2)-(x4-x2)*(y3-y2);
      hBitmap->Canvas->Brush->Color= (q1-q2)>0 ? clLime : clBlue;
      if(DrawInit.DrawLines) {
        hBitmap->Canvas->Pen->Color=clBlack;
        hBitmap->Canvas->MoveTo(x1,y1);
        hBitmap->Canvas->LineTo(x2,y2);
        hBitmap->Canvas->LineTo(x3,y3);
      } else {
        hBitmap->Canvas->Polygon(Poly,3);
      }
      jBack=j;
    }
    iBack=i;
  }
  if(DrawInit.DoubleSided) {
    iBack=0;
    for(int i=1;i<=nXSteps-1;i=i+DrawInit.Skip < nXSteps-1 ? i+DrawInit.Skip : (i<nXSteps-1?nXSteps-1:nXSteps) ) {
      jBack=0;
      for(int j=1;j<=Init.nSteps-1;j=j+DrawInit.Skip < Init.nSteps-1 ? j+DrawInit.Skip : (j<Init.nSteps-1?Init.nSteps-1:Init.nSteps)) {
        x1=Gr3D.GetX(-iBack*dx,jBack*Init.dt,abs(GetAVal(jBack,iBack,tStep)));
        y1=Gr3D.GetY(-iBack*dx,jBack*Init.dt,abs(GetAVal(jBack,iBack,tStep)));
        x2=Gr3D.GetX(-i*dx,jBack*Init.dt,abs(GetAVal(jBack,i,tStep)));
        y2=Gr3D.GetY(-i*dx,jBack*Init.dt,abs(GetAVal(jBack,i,tStep)));
        x3=Gr3D.GetX(-i*dx,j*Init.dt,abs(GetAVal(j,i,tStep)));
        y3=Gr3D.GetY(-i*dx,j*Init.dt,abs(GetAVal(j,i,tStep)));
        x4=Gr3D.GetX(-iBack*dx,j*Init.dt,abs(GetAVal(j,iBack,tStep)));
        y4=Gr3D.GetY(-iBack*dx,j*Init.dt,abs(GetAVal(j,iBack,tStep)));
        Poly[0]=TPoint(x1,y1);
        Poly[1]=TPoint(x2,y2);
        Poly[2]=TPoint(x3,y3);
        Poly[3]=TPoint(x4,y4);
        q1=(x1-x2)*(y4-y2)-(x4-x2)*(y1-y2);
        q2=(x3-x2)*(y4-y2)-(x4-x2)*(y3-y2);
        hBitmap->Canvas->Brush->Color= (q1-q2)>0 ? clLime : clBlue;
        if(DrawInit.DrawLines) {
          hBitmap->Canvas->Pen->Color=clBlack;
          hBitmap->Canvas->MoveTo(x1,y1);
          hBitmap->Canvas->LineTo(x2,y2);
          hBitmap->Canvas->LineTo(x3,y3);
        } else hBitmap->Canvas->Polygon(Poly,3);
        jBack=j;
      }
      iBack=i;
    }
  }
}

void cDump2DOpticModel::SaveCURR(char *hFileName){
double t,dt=Init.L/Init.nSteps;
  _fmode=O_BINARY;
  int fCURR=creat("CURR.dat", S_IWRITE);
  char fEoLine[3]; fEoLine[0]=13; fEoLine[1]=10;  fEoLine[2]=' ';
  AnsiString tmpStr;
  int nZLines;
  int tmp2=0;
  tmpStr.sprintf("%G %d %d %G",(double)dt,1,Init.nSteps,(double)0);
  write(fCURR,tmpStr.c_str(),tmpStr.Length());
  write(fCURR,fEoLine,2);
  for(int j=0; (j<Step) ;j++) {
    for(int i=0;i<=Init.nSteps;i++){
      tmp2++;
      tmpStr.sprintf("%G",(double)hI[ (Init.nSteps+1)*j+i ].real());
      write(fCURR,tmpStr.c_str(),tmpStr.Length());
      if(tmp2==1) {
        write(fCURR,fEoLine,2);
        tmp2=0;
      }
        else write(fCURR,fEoLine+2,1);
    }
    write(fCURR,fEoLine,2);
    for(int i=0;i<=Init.nSteps;i++){
      tmp2++;
      tmpStr.sprintf("%G",(double)hI[ (Init.nSteps+1)*j+i ].imag());
      write(fCURR,tmpStr.c_str(),tmpStr.Length());
      if(tmp2==1) {
        write(fCURR,fEoLine,2);
        tmp2=0;
      }
        else write(fCURR,fEoLine+2,1);
    }
  }
  close(fCURR);
}

void cDump2DOpticModel::SaveZTSlice(DrawDumpInitData DrawInit) {
  double z,t;
  _fmode=O_BINARY;
  int fSurface=creat("ZTSlice.dat", S_IWRITE);
  char fEoLine[2]; fEoLine[0]=13; fEoLine[1]=10;
  AnsiString tmpStr;
  int nZLines;
  int tmp1=0;
  for(int i=0;i<=Init.nSteps;i++){
    z=i*Init.dt;
    for(int j=0; (j<Step) ;j++) {
      t=j*Init.dt;
      tmpStr.sprintf("%G,%G,%G",(double)z,(double)t,(double)abs(hA[ (Init.nSteps+1)*j+i ]));
      tmp1=write(fSurface,tmpStr.c_str(),tmpStr.Length());
      tmp1=write(fSurface,fEoLine,2);
      tmp1=tmp1;
    }
  }
  close(fSurface);
}

void cDump2DOpticModel::SaveXZSlice(DrawDumpInitData DrawInit) {
  double z,x;
  _fmode=O_BINARY;
  int tStep=DrawInit.t/Init.dt;
  if(tStep>Step-1) return;
  int fSurface=creat("ZTSlice.dat", S_IWRITE);
  char fEoLine[2]; fEoLine[0]=13; fEoLine[1]=10;
  if ( (DrawInit.dx!=dx) || (DrawInit.nXSteps!=nXSteps))
                         MakeXProfileCoef(DrawInit.dx,DrawInit.nXSteps);
  AnsiString tmpStr;
  int nZLines;
  int tmp1=0;
  for(int i=0;i<=Init.nSteps;i++){
    z=i*Init.dt;
    for(int j=0; (j<nXSteps) ;j++) {
      x=j*dx;
      tmpStr.sprintf("%G,%G,%G",(double)x,(double)z,(double)abs(GetAVal(i,j,tStep)));
      tmp1=write(fSurface,tmpStr.c_str(),tmpStr.Length());
      tmp1=write(fSurface,fEoLine,2);
      tmpStr.sprintf("%G,%G,%G",-(double)x,(double)z,(double)abs(GetAVal(i,j,tStep)));
      tmp1=write(fSurface,tmpStr.c_str(),tmpStr.Length());
      tmp1=write(fSurface,fEoLine,2);
      tmp1=tmp1;
    }
  }
  close(fSurface);
}

void cDump2DOpticModel::SaveXTSlice(DrawDumpInitData DrawInit) {
  double t,x;
  _fmode=O_BINARY;
  int zStep=DrawInit.t/Init.dt;
  if(zStep>Init.nSteps) return;
  int fSurface=creat("ZTSlice.dat", S_IWRITE);
  char fEoLine[2]; fEoLine[0]=13; fEoLine[1]=10;
  if ( (DrawInit.dx!=dx) || (DrawInit.nXSteps!=nXSteps))
                         MakeXProfileCoef(DrawInit.dx,DrawInit.nXSteps);
  AnsiString tmpStr;
  int nZLines;
  int tmp1=0;
  for(int i=0;i<Step;i++){
    t=i*Init.dt;
    for(int j=0; (j<nXSteps) ;j++) {
      x=j*dx;
      tmpStr.sprintf("%G,%G,%G",(double)x,(double)t,(double)abs(GetAVal(zStep,j,i)));
      tmp1=write(fSurface,tmpStr.c_str(),tmpStr.Length());
      tmp1=write(fSurface,fEoLine,2);
      tmpStr.sprintf("%G,%G,%G",-(double)x,(double)t,(double)abs(GetAVal(zStep,j,i)));
      tmp1=write(fSurface,tmpStr.c_str(),tmpStr.Length());
      tmp1=write(fSurface,fEoLine,2);
      tmp1=tmp1;
    }
  }
  close(fSurface);
}
*/
/*
void cDumpStorage::Init_Diffur(int pnSteps, int pnEls,
                               REAL_TYPE pL, REAL_TYPE ptMax)
{
int pnStepsT=pnSteps*ptMax/pL;
  if(!DumpRequired) return;
  if((nSteps+1)*nStepsT!=(pnSteps+1)*pnStepsT)
  {
    delete [] hADump; delete [] hIDump;
    hADump = new complex<float>[(pnSteps+1)*pnStepsT];
    hIDump = new complex<float>[(pnSteps+1)*pnStepsT];
  }
  if( (nEls*(nSteps+1)*nStepsT!=pnEls*(pnSteps+1)*pnSteps) & (DumpElectrons) )
  {
    delete [] hElDump; hElDump = new float[ 2*pnEls*(pnSteps+1)*pnStepsT ];
  }
  nSteps=pnSteps; nStepsT=pnStepsT; nEls=pnEls; L=pL; EndTime=ptMax;
  Dumping=DUMPING_DIFFUR; Step=StepT=0;
  DumpReady=false;
}

void cDumpStorage::Init_Opti(int Type,
            int pnOptiSteps1,int pnOptiSteps2,int pnOptiSteps3)
{
  if(!DumpRequired) return;
  if( nOptiSteps1*nOptiSteps2*nOptiSteps3!=
      pnOptiSteps1*pnOptiSteps2*pnOptiSteps3 )
  {
    delete [] hMaxADump;
    hMaxADump = new float[pnOptiSteps1*pnOptiSteps2*pnOptiSteps3];
  }
  nOptiSteps1=pnOptiSteps1; nOptiSteps2=pnOptiSteps2; nOptiSteps3=pnOptiSteps3;
  Dumping=Type; OptiStep1=OptiStep2=OptiStep3=0;
  DumpReady=false;
}

void cDumpStorage::DumpOpti(REAL_TYPE MaxA)
{
  if(!DumpRequired) return;
  if( (Dumping>DUMPING_DIFFUR) && (!DumpReady) )
  {
    hMaxADump[ (nOptiSteps2*OptiStep3+OptiStep2)*nOptiSteps1+OptiStep1 ]=MaxA;
    if(++OptiStep1==nOptiSteps1)
    {
      OptiStep1=0;
      if(++OptiStep2==nOptiSteps2)
      {
        OptiStep2=0;
        if(++OptiStep3==nOptiSteps3) DumpReady=true;
      }
    }
  }
}

void cDumpStorage::DumpDiffur( complex<REAL_TYPE> *hA, complex<REAL_TYPE> *hI,
                               REAL_TYPE *hEls)
{
  if(!DumpRequired) return;
  if( (Dumping==DUMPING_DIFFUR) && (!DumpReady) )
  {
    hADump[ (nSteps+1)*StepT+Step ]=*hA;
    hIDump[ (nSteps+1)*StepT+Step ]=*hI;
    if (DumpElectrons)
      for(int i=0;i<nEls;i++)
      {
        hElDump[ 2*(((nSteps+1)*StepT+Step)*nEls+i)]=hEls[2*i];
        hElDump[ 2*(((nSteps+1)*StepT+Step)*nEls+i)+1]=hEls[2*i+1];
      }
    if(++Step==nSteps+1)
    {
      Step=0;
      if(++StepT==nStepsT) DumpReady=true;
    }
  }
}

void cDumpStorage::Save_AUDump(char *hFileName,float time)
{
 if(!DumpReady) return;
 ofstream outFile = ofstream(hFileName,ios::out);
 StepT=time*(float)(nStepsT-1)/EndTime;
 MaxA=0;
 for(int j=0;j<nSteps+1;j++)
  if(MaxA<abs(hADump[ (nSteps+1)*StepT+j ])) MaxA=abs(hADump[ (nSteps+1)*StepT+j ]);
 if (MaxA<0.0000000001) MaxA=0.0000000001;
 MaxU=0;
 for(int j=0;j<nSteps+1;j++)
   for(int i=0;i<nEls;i++)
     if (MaxU<hElDump[ 2*(((nSteps+1)*StepT+j)*nEls+i)+1])
       MaxU=hElDump[ 2*(((nSteps+1)*StepT+j)*nEls+i)+1];
 if (MaxU<0.0001) MaxU=0.0001;
 for(int j=0;j<nSteps+1;j++)
 {
   outFile << L*(float)j/(float)nSteps << "  ";
   outFile << MaxU*abs(hADump[ (nSteps+1)*StepT+j ])/MaxA + MaxU*1.2<< "  ";
   for(int i=0;i<nEls;i++)
     outFile << hElDump[ 2*(((nSteps+1)*StepT+j)*nEls+i)+1] << "  ";
   outFile << endl;
 }
}

void cDumpStorage::Save_AIDump(char *hFileName,float time)
{
 if(!DumpReady) return;
 ofstream outFile = ofstream(hFileName,ios::out);
 StepT=time*(float)(nStepsT-1)/EndTime;
 MaxA=0;
 for(int j=0;j<nSteps+1;j++)
  if(MaxA<abs(hADump[ (nSteps+1)*StepT+j ])) MaxA=abs(hADump[ (nSteps+1)*StepT+j ]);
 if (MaxA<0.0000000001) MaxA=0.0000000001;
 MaxI=0;
 for(int j=0;j<nSteps+1;j++)
  if(MaxI<abs(hIDump[ (nSteps+1)*StepT+j ])) MaxI=abs(hIDump[ (nSteps+1)*StepT+j ]);
 if (MaxI<0.0000000001) MaxI=0.0000000001;
 for(int j=0;j<nSteps+1;j++)
 {
   outFile << L*(float)j/(float)nSteps << "  ";
   outFile << abs(hADump[ (nSteps+1)*StepT+j ])<< "  ";
   outFile << abs(hIDump[ (nSteps+1)*StepT+j ])*MaxA/MaxI - MaxA*1.2<< "  ";
   outFile << endl;
 }
}

void cDumpStorage::Save_AIUDump(char *hFileName,float time)
{
 if(!(DumpReady & DumpElectrons)) return;
 ofstream outFile = ofstream(hFileName,ios::out);
 StepT=time*(float)nStepsT/EndTime;
 MaxA=0;
 for(int j=0;j<nSteps+1;j++)
  if(MaxA<abs(hADump[ (nSteps+1)*StepT+j ])) MaxA=abs(hADump[ (nSteps+1)*StepT+j ]);
 if (MaxA<0.0000000001) MaxA=0.0000000001;
 MaxI=0;
 for(int j=0;j<nSteps+1;j++)
  if(MaxI<abs(hIDump[ (nSteps+1)*StepT+j ])) MaxI=abs(hIDump[ (nSteps+1)*StepT+j ]);
 if (MaxI<0.0000000001) MaxI=0.0000000001;
 MaxU=0;
 for(int j=0;j<nSteps+1;j++)
   for(int i=0;i<nEls;i++)
     if (MaxU<hElDump[ 2*(((nSteps+1)*StepT+j)*nEls+i)+1])
       MaxU=hElDump[ 2*(((nSteps+1)*StepT+j)*nEls+i)+1];
 if (MaxU<0.0001) MaxU=0.0001;
 for(int j=0;j<nSteps+1;j++)
 {
   outFile << L*(float)j/(float)nSteps << "  ";
   outFile << MaxU*abs(hADump[ (nSteps+1)*StepT+j ])/MaxA + MaxU*2.2<< "  ";
   outFile << MaxU*abs(hIDump[ (nSteps+1)*StepT+j ])/MaxI + MaxU*1.2<< "  ";
   for(int i=0;i<nEls;i++)
     outFile << hElDump[ 2*(((nSteps+1)*StepT+j)*nEls+i)+1] << "  ";
   outFile << endl;
 }
}

void cDumpStorage::Draw_AUDump(TImage *hImage,float time)
{
 if(!DumpReady) return;
 StepT=time*(float)(nStepsT-1)/EndTime;
 StepT = StepT<nStepsT ? StepT : nStepsT-1;
 double realMaxA,realMaxI;
 MaxA=realMaxA=0;
 for(int j=0;j<nSteps+1;j++)
  if(MaxA<abs(hADump[ (nSteps+1)*StepT+j ])) realMaxA=MaxA=abs(hADump[ (nSteps+1)*StepT+j ]);
 if (MaxA<0.0000000001) MaxA=0.0000000001;
 realMaxI=MaxI=0;
 for(int j=0;j<nSteps+1;j++)
  if(MaxI<abs(hIDump[ (nSteps+1)*StepT+j ])) realMaxI=MaxI=abs(hIDump[ (nSteps+1)*StepT+j ]);
 if (MaxI<0.0000000001) MaxI=0.0000000001;
 MaxU=0;
 MinU=0;
 if (DumpElectrons) {
   for(int j=0;j<nSteps+1;j++)
     for(int i=0;i<nEls;i++)
     {
       if (MaxU<hElDump[ 2*(((nSteps+1)*StepT+j)*nEls+i)+1])
         MaxU=hElDump[ 2*(((nSteps+1)*StepT+j)*nEls+i)+1];
       if (MinU>hElDump[ 2*(((nSteps+1)*StepT+j)*nEls+i)+1])
         MinU=hElDump[ 2*(((nSteps+1)*StepT+j)*nEls+i)+1];
     }
   if (MaxU<0.0001) MaxU=0.0001; if (MinU>-0.0001) MinU=-0.0001;
 }
 hImage->Canvas->FillRect(TRect(0,0,640,480));
 hImage->Canvas->MoveTo(640,120);
 hImage->Canvas->LineTo(0,120);
 hImage->Canvas->MoveTo(0,120-100*abs(hADump[ (nSteps+1)*StepT ])/MaxA);
 for(int j=0;j<nSteps+1;j++)
   hImage->Canvas->LineTo(640*j/nSteps,120-100*abs(hADump[ (nSteps+1)*StepT+j ])/MaxA);
 hImage->Canvas->MoveTo(640,350);
 hImage->Canvas->LineTo(0,350);
 if (DumpElectrons) {
   for(int i=0;i<nEls;i++)
   {
     hImage->Canvas->MoveTo(0,350-120*hElDump[ 2*(((nSteps+1)*StepT)*nEls+i)+1]/(MaxU-MinU));
     for(int j=0;j<nSteps+1;j++)
       hImage->Canvas->LineTo(640*j/nSteps,350-120*hElDump[ 2*(((nSteps+1)*StepT+j)*nEls+i)+1]/(MaxU-MinU));
   }
 }
 hImage->Canvas->TextOutA(0,0,"Maximum A: "+AnsiString(realMaxA));
// hImage->Canvas->TextOutA(10,25,"Maximum I: "+AnsiString(realMaxI));
}

void cDumpStorage::Draw_AIUDump(TImage *hImage,float time)
{
 if(!DumpReady) return;
 StepT=time*(float)(nStepsT-1)/EndTime;
 double realMaxA,realMaxI;
 MaxA=realMaxA=0;
 for(int j=0;j<nSteps+1;j++)
  if(MaxA<abs(hADump[ (nSteps+1)*StepT+j ])) realMaxA=MaxA=abs(hADump[ (nSteps+1)*StepT+j ]);
 if (MaxA<0.0000000001) MaxA=0.0000000001;
 realMaxI=MaxI=0;
 for(int j=0;j<nSteps+1;j++)
  if(MaxI<abs(hIDump[ (nSteps+1)*StepT+j ])) realMaxI=MaxI=abs(hIDump[ (nSteps+1)*StepT+j ]);
 if (MaxI<0.0000000001) MaxI=0.0000000001;
 MaxU=0;
 MinU=0;
 if (DumpElectrons) {
   for(int j=0;j<nSteps+1;j++)
     for(int i=0;i<nEls;i++)
     {
       if (MaxU<hElDump[ 2*(((nSteps+1)*StepT+j)*nEls+i)+1])
         MaxU=hElDump[ 2*(((nSteps+1)*StepT+j)*nEls+i)+1];
       if (MinU>hElDump[ 2*(((nSteps+1)*StepT+j)*nEls+i)+1])
         MinU=hElDump[ 2*(((nSteps+1)*StepT+j)*nEls+i)+1];
     }
   if (MaxU<0.0001) MaxU=0.0001; if (MinU>-0.0001) MinU=-0.0001;
 }
 hImage->Canvas->FillRect(TRect(0,0,640,480));
 hImage->Canvas->MoveTo(640,110);
 hImage->Canvas->LineTo(0,110);
 hImage->Canvas->MoveTo(0,110-100*abs(hADump[ (nSteps+1)*StepT ])/MaxA);
 for(int j=0;j<nSteps+1;j++)
   hImage->Canvas->LineTo(640*j/nSteps,110-100*abs(hADump[ (nSteps+1)*StepT+j ])/MaxA);

 hImage->Canvas->MoveTo(640,230);
 hImage->Canvas->LineTo(0,230);
 hImage->Canvas->MoveTo(0,230-100*abs(hIDump[ (nSteps+1)*StepT ])/MaxI);
 for(int j=0;j<nSteps+1;j++)
   hImage->Canvas->LineTo(640*j/nSteps,230-100*abs(hIDump[ (nSteps+1)*StepT+j ])/MaxI);
   
 hImage->Canvas->MoveTo(640,360);
 hImage->Canvas->LineTo(0,360);
 if (DumpElectrons) {
   for(int i=0;i<nEls;i++)
   {
     hImage->Canvas->MoveTo(0,350-120*hElDump[ 2*(((nSteps+1)*StepT)*nEls+i)+1]/(MaxU-MinU));
     for(int j=0;j<nSteps+1;j++)
       hImage->Canvas->LineTo(640*j/nSteps,350-120*hElDump[ 2*(((nSteps+1)*StepT+j)*nEls+i)+1]/(MaxU-MinU));
   }
 }
 hImage->Canvas->TextOutA(0,0,"Maximum A: "+AnsiString(realMaxA));
// hImage->Canvas->TextOutA(10,25,"Maximum I: "+AnsiString(realMaxI));
}

complex<float> cDumpStorage::GetAValue(REAL_TYPE z,REAL_TYPE t)
{
  if (!DumpReady) return( complex<float>(0,0) );
  StepT=t*(double)nStepsT/EndTime;
  int StepZ=z*(double)(nSteps+1)/L;
  return( hADump[ (nSteps+1)*StepT+StepZ ] );
}
complex<float> cDumpStorage::GetIValue(REAL_TYPE z,REAL_TYPE t)
{
  if (!DumpReady) return( complex<float>(0,0) );
  StepT=t*(double)nStepsT/EndTime;
  int StepZ=z*(double)(nSteps+1)/L;
  return( hIDump[ (nSteps+1)*StepT+StepZ ] );
}

void cDumpStorage::SaveSurfaceZT(char *hFileName) {
  double z,t,dt=L/nSteps;
  _fmode=O_BINARY;
  int fSurface=creat(hFileName, S_IWRITE);
  char fEoLine[2]; fEoLine[0]=13; fEoLine[1]=10;
  AnsiString tmpStr;
  int nZLines;
  int tmp1=0;
  for(int i=0;i<=nSteps;i++){
    z=i*dt;
    for(int j=0; (j<nStepsT) ;j++) {
      t=j*dt;
      tmpStr.sprintf("%G,%G,%G",z,t,abs(hADump[ (nSteps+1)*j+i ]));
      tmp1=write(fSurface,tmpStr.c_str(),tmpStr.Length());
      tmp1=write(fSurface,fEoLine,2);
    }
  }
  close(fSurface);
}

void cDumpStorage::SaveCURR(char *hFileName){
double t,dt=L/nSteps;
  _fmode=O_BINARY;
  int fCURR=creat("CURR.dat", S_IWRITE);
  char fEoLine[3]; fEoLine[0]=13; fEoLine[1]=10;  fEoLine[2]=' ';
  AnsiString tmpStr;
  int nZLines;
  int tmp1=0,tmp2=0;
  tmpStr.sprintf("%G %d %d %G",dt,1,nSteps,0);
  tmp1=write(fCURR,tmpStr.c_str(),tmpStr.Length());
  tmp1=write(fCURR,fEoLine,2);
  for(int j=0; (j<nStepsT) ;j++) {
    for(int i=0;i<=nSteps;i++){
      tmp2++;
      tmpStr.sprintf("%G",hIDump[ (nSteps+1)*j+i ].real());
      tmp1=write(fCURR,tmpStr.c_str(),tmpStr.Length());
      if(tmp2==5) {
        tmp1=write(fCURR,fEoLine,2);
        tmp2=0;
      }
        else tmp1=write(fCURR,fEoLine+2,1);
    }
    tmp1=write(fCURR,fEoLine,2);
    for(int i=0;i<=nSteps;i++){
      tmp2++;
      tmpStr.sprintf("%G",hIDump[ (nSteps+1)*j+i ].imag());
      tmp1=write(fCURR,tmpStr.c_str(),tmpStr.Length());
      if(tmp2==5) {
        tmp1=write(fCURR,fEoLine,2);
        tmp2=0;
      }
        else tmp1=write(fCURR,fEoLine+2,1);
    }
  }
  close(fCURR);
}
*/
