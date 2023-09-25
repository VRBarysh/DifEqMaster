//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquPSO2DBraggSolver.h"
#include <fstream>

//---------------------------------------------------------------------------

#pragma package(smart_init)

cEquPSO2DBraggSolver::cEquPSO2DBraggSolver() : cDifEqu_DataBankIO()
{
  hVariable=NULL; hPointData=NULL; hResult=NULL; hAnalResult=NULL;
}

cEquPSO2DBraggSolver::~cEquPSO2DBraggSolver()// : ~cDifEqu_DataBankIO()
{
  delete [] hPointData; hPointData=NULL;
  delete [] hResult; delete [] hAnalResult;
  hResult=NULL; hAnalResult=NULL;
}

void cEquPSO2DBraggSolver::InitReport(cDifEquReport *hReport) {
  int Lines[EQU_MAX_GRAPHS];
  for(int i=0;i<EQU_MAX_GRAPHS;i++) Lines[i]=1;
  Lines[0]=1; Lines[1]=1; Lines[2]=1; Lines[3]=1;
  hReport->nVals=0; hReport->nGraphs=0;
//  hReport->LoadText(hReportText);
  hReport->SetGraphSize(EQU_GRAPH_SIZE,Lines);
  hReport->SetMaxTo1();
}

void cEquPSO2DBraggSolver::Report(cDifEquReport *hReport) {
  long double MaxFVal[50];
  cDifEqu_InitDataRec_PSOSolver tmpInitData=InitData;
  hReport->nGraphs=2;  hReport->nVals=0;
  hReport->Graph[0].GetData(0,InitData.nScaleSteps,hAnalResult);
  hReport->Graph[1].GetData(0,InitData.nScaleSteps,hResult);
  hReport->ValText[hReport->nVals]="Anal Max   ";
  hReport->Val[hReport->nVals]=hReport->Graph[0].Max;
  hReport->nVals++;
  hReport->ValText[hReport->nVals]="Delta Max   ";
  hReport->Val[hReport->nVals]=hReport->Graph[1].Max;
  hReport->nVals++;

  if( (InitData.PSOData.nComponents!=2) ) {
    for(int i=0;i<InitData.PSOData.nGroups;i++) {
      curDeltaRe=(*PSOOptimizer.GetGlobalMaxX(i))[0];
      curDeltaIm=(*PSOOptimizer.GetGlobalMaxX(i))[1];
      curGammaRe=(*PSOOptimizer.GetGlobalMaxX(i))[2];
      curGammaIm=(*PSOOptimizer.GetGlobalMaxX(i))[3];
      MaxFVal[i]=MainFunction();
    }
    InitData=tmpInitData;
    long double TargetFVal=0;
//  if(InitData.mode==1) {
      TCplxLong gamma=GammaAnal(InitData.mode,InitData.tnx,InitData.tnz);
      TCplxLong delta=DeltaAnal(InitData.mode,InitData.tnx,InitData.tnz);
      curDeltaRe=delta.real();
      curDeltaIm=delta.imag();
      curGammaRe=gamma.real();
      curGammaIm=gamma.imag();
      TargetFVal=MainFunction();
      InitData=tmpInitData;
//    }
    for(int i=0;i<InitData.PSOData.nGroups;i++) {
      hReport->ValText[hReport->nVals]=AnsiString("Delta(") +
               AnsiString(double((*PSOOptimizer.GetGlobalMaxX(i))[0])) + AnsiString(",") +
               AnsiString(double((*PSOOptimizer.GetGlobalMaxX(i))[1])) + AnsiString(")  ") +
               AnsiString("Gamma(") +
               AnsiString(double((*PSOOptimizer.GetGlobalMaxX(i))[2])) + AnsiString(",") +
               AnsiString(double((*PSOOptimizer.GetGlobalMaxX(i))[3])) + AnsiString(")    ");
      hReport->Val[hReport->nVals]=0; hReport->nVals++;
      hReport->ValText[hReport->nVals]="Function Max   ";
      hReport->Val[hReport->nVals]=MaxFVal[i];
      hReport->nVals++;
    }
  } else {
    for(int i=0;i<InitData.PSOData.nGroups;i++) {
      curDeltaRe=0;
      curDeltaIm=(*PSOOptimizer.GetGlobalMaxX(i))[0];
      curGammaRe=0;
      curGammaIm=(*PSOOptimizer.GetGlobalMaxX(i))[1];
      MaxFVal[i]=MainFunction();
    }
    for(int i=0;i<InitData.PSOData.nGroups;i++) {
      hReport->ValText[hReport->nVals]=AnsiString("Delta(") +
               AnsiString(0) + AnsiString(",") +
               AnsiString(double((*PSOOptimizer.GetGlobalMaxX(i))[0])) + AnsiString(")  ") +
               AnsiString("Gamma(") +
               AnsiString(0) + AnsiString(",") +
               AnsiString(double((*PSOOptimizer.GetGlobalMaxX(i))[1])) + AnsiString(")    ");
      hReport->Val[hReport->nVals]=0; hReport->nVals++;
      hReport->ValText[hReport->nVals]="Function Max   ";
      hReport->Val[hReport->nVals]=MaxFVal[i];
      hReport->nVals++;
    }
    InitData=tmpInitData;
  }

//  hReport->nVals=6;
//  hReport->ValText[0]="Delta Re       "; hReport->Val[0]=(*PSOOptimizer.GetGlobalMaxX())[0];
//  hReport->ValText[1]="Delta Im       "; hReport->Val[1]=(*PSOOptimizer.GetGlobalMaxX())[1];
//  hReport->ValText[2]="Gamma Re       "; hReport->Val[2]=(*PSOOptimizer.GetGlobalMaxX())[2];
//  hReport->ValText[3]="Gamma Im       "; hReport->Val[3]=(*PSOOptimizer.GetGlobalMaxX())[3];
//  hReport->ValText[4]="Function Max   "; hReport->Val[4]=MaxFVal;
/*  hReport->ValText[hReport->nVals]="Function Target ";
  hReport->Val[hReport->nVals]=TargetFVal;   hReport->nVals++;
  hReport->ValText[hReport->nVals]="Delta Re   ";
  hReport->Val[hReport->nVals]=delta.real(); hReport->nVals++;
  hReport->ValText[hReport->nVals]="Delta Im   ";
  hReport->Val[hReport->nVals]=delta.imag(); hReport->nVals++;
  hReport->ValText[hReport->nVals]="Gamma Re   ";
  hReport->Val[hReport->nVals]=gamma.real(); hReport->nVals++;
  hReport->ValText[hReport->nVals]="Gamma Im   ";
  hReport->Val[hReport->nVals]=gamma.imag(); hReport->nVals++;
  hReport->ValText[hReport->nVals]="fError     "; hReport->Val[hReport->nVals]=fError;
  hReport->nVals++;*/
}

void cEquPSO2DBraggSolver::FinalReport(cDifEquReport *hReport) {
  Report(hReport);
/*  cDifEqu_InitDataRec_PSOSolver tmpInitData=InitData;
  hReport->nGraphs=0;
  curDeltaRe=(*PSOOptimizer.GetGlobalMaxX())[0];
  curDeltaIm=(*PSOOptimizer.GetGlobalMaxX())[1];
  curGammaRe=(*PSOOptimizer.GetGlobalMaxX())[2];
  curGammaIm=(*PSOOptimizer.GetGlobalMaxX())[3];
  long double tmp=MainFunction();
  InitData=tmpInitData;
  hReport->nVals=5;
  hReport->ValText[0]="Delta Re "; hReport->Val[0]=(*PSOOptimizer.GetGlobalMaxX())[0];
  hReport->ValText[1]="Delta Im "; hReport->Val[1]=(*PSOOptimizer.GetGlobalMaxX())[1];
  hReport->ValText[2]="Gamma Re "; hReport->Val[2]=(*PSOOptimizer.GetGlobalMaxX())[2];
  hReport->ValText[3]="Gamma Im "; hReport->Val[3]=(*PSOOptimizer.GetGlobalMaxX())[3];
  hReport->ValText[4]="Function "; hReport->Val[4]=tmp;*/
/*
  hReport->nGraphs=2;
  hReport->GraphText[0]="OutGraph";
  hReport->Graph[0].GetData(0,Var.nSteps,MaxOut);
  hReport->GraphText[1]="MaxTotalGraph";
  hReport->Graph[1].GetData(0,Var.nSteps,MaxTotal);
  hReport->nVals=4;
  hReport->ValText[0]="Max Out A = ";  hReport->Val[0]=MaxOutMax;
  hReport->ValText[1]="Optimum Value = ";  hReport->Val[1]=BestOutValue;
  hReport->ValText[2]="Max A = ";  hReport->Val[2]=MaxTotalMax;
  hReport->ValText[3]="Optimum Value = ";  hReport->Val[3]=BestTotalValue;

  ofstream FileOut = ofstream("MaxOut.dat",ios::out);
  ofstream FileTotal = ofstream("MaxTotal.dat",ios::out);
  long double value;
  for(int i=0;i<Var.nSteps;i++) {
    value=Var.V1+(Var.V2-Var.V1)*i/(Var.nSteps-1);
    FileOut << value << " " << MaxOut[i] << endl;
    FileTotal << value << " " << MaxTotal[i] << endl;
  }*/
  if(InitData.mode==4) {
    SaveResults();
  }
  if(InitData.mode==3) {
    // SaveDeltaImGammaImPlane();
    SaveDeltaReIm();
    SaveGammaReIm();
  }
}

void cEquPSO2DBraggSolver::SaveResults() {
  ofstream FileOut = ofstream("OptResults.dat",ios::out);
  double Lz;
  for(int i=0;i<InitData.nScaleSteps;i++) {
    Lz=InitData.Lz0-(InitData.Lz0-InitData.LzMin)*i/(InitData.nScaleSteps-1);
    FileOut << Lz << " " << hResult[i] << " " << hAnalResult[i] << endl;
  }
}

void cEquPSO2DBraggSolver::SaveDeltaImGammaImPlane() {
  ofstream FileOut = ofstream("DeltaPlane.dat",ios::out);
  curDeltaRe=0;
  curGammaRe=0;
  for(curDeltaIm=InitData.deltaImMin;curDeltaIm<InitData.deltaImMax;curDeltaIm+=InitData.epsilon)
    for(curGammaIm=InitData.gammaImMin;curGammaIm<InitData.gammaImMax;curGammaIm+=InitData.epsilon) {
      FileOut << curDeltaIm << " " << curGammaIm << " " << MainFunction() << endl;
    }
}

void cEquPSO2DBraggSolver::SaveDeltaReIm() {
  ofstream FileOut = ofstream("PSODeltaReIm.dat",ios::out);
  TCplxLong gamma=GammaAnal(1,InitData.tnx,InitData.tnz);
  TCplxLong delta=DeltaAnal(1,InitData.tnx,InitData.tnz);
  curGammaRe=gamma.real(); curGammaIm=gamma.imag();
  for(int iRe=0;iRe<InitData.nScaleSteps;iRe++)
    for(int iIm=0;iIm<InitData.nScaleSteps;iIm++) {
      curDeltaRe = delta.real() - InitData.epsilon + 2.0*InitData.epsilon*iRe/(InitData.nScaleSteps-1);
      curDeltaIm = delta.imag() - InitData.epsilon + 2.0*InitData.epsilon*iIm/(InitData.nScaleSteps-1);
      FileOut << curDeltaRe << " " << curDeltaIm << " " << MainFunction() << endl;
    }
}

void cEquPSO2DBraggSolver::SaveGammaReIm() {
  ofstream FileOut = ofstream("PSOGammaReIm.dat",ios::out);
  TCplxLong gamma=GammaAnal(1,InitData.tnx,InitData.tnz);
  TCplxLong delta=DeltaAnal(1,InitData.tnx,InitData.tnz);
  curDeltaRe=delta.real(); curDeltaIm=delta.imag();
  for(int iRe=0;iRe<InitData.nScaleSteps;iRe++)
    for(int iIm=0;iIm<InitData.nScaleSteps;iIm++) {
      curGammaRe = gamma.real() - InitData.epsilon + 2.0*InitData.epsilon*iRe/(InitData.nScaleSteps-1);
      curGammaIm = gamma.imag() - InitData.epsilon + 2.0*InitData.epsilon*iIm/(InitData.nScaleSteps-1);
      FileOut << curGammaRe << " " << curGammaIm << " " << MainFunction() << endl;
    }
}

void cEquPSO2DBraggSolver::LoadInitData(cDifEqu_InitDataRec *hInitData)
{
  cDifEqu_InitDataRec_PSOSolver *hInit=(cDifEqu_InitDataRec_PSOSolver *)hInitData;
  InitData=*hInit;
  InitData.Lz0=InitData.Lz; InitData.Lx0=InitData.Lx; 
  cDifEqu_InitDataRec_PSOSolver tmpInitData=InitData;
  if(hPointData) delete [] hPointData;
  if(hResult) delete [] hResult;
  if(hAnalResult) delete [] hAnalResult;
  hResult = new long double[InitData.nScaleSteps];
  hAnalResult = new long double[InitData.nScaleSteps];
  for(int i=0;i<InitData.nScaleSteps;i++) {
    InitData.Lz=tmpInitData.Lz-(tmpInitData.Lz-InitData.LzMin)*i/(InitData.nScaleSteps-1);
    InitData.Lx=tmpInitData.Lx*InitData.Lz/tmpInitData.Lz;
    hResult[i]=0; hAnalResult[i]=DeltaAnal(1,InitData.tnx,InitData.tnz).imag();
  }
  InitData=tmpInitData;
  InitData.PSOData.nComponents=4;
  hPointData = new CPSOInitComp[InitData.PSOData.nComponents];
  hPointData[0].vMin=InitData.deltaReMin;
  hPointData[0].vMax=InitData.deltaReMax;
  hPointData[0].hV=&curDeltaRe;
  hPointData[1].vMin=InitData.deltaImMin;
  hPointData[1].vMax=InitData.deltaImMax;
  hPointData[1].hV=&curDeltaIm;
  hPointData[2].vMin=InitData.gammaReMin;
  hPointData[2].vMax=InitData.gammaReMax;
  hPointData[2].hV=&curGammaRe;
  hPointData[3].vMin=InitData.gammaImMin;
  hPointData[3].vMax=InitData.gammaImMax;
  hPointData[3].hV=&curGammaIm;
  if(InitData.mode==1) {
    TCplxLong gamma=GammaAnal(InitData.mode,InitData.tnx,InitData.tnz);
    TCplxLong delta=DeltaAnal(InitData.mode,InitData.tnx,InitData.tnz);
    hPointData[0].vMin=delta.real()-InitData.epsilon;
    hPointData[0].vMax=delta.real()+InitData.epsilon;
    hPointData[1].vMin=delta.imag()-InitData.epsilon;
    hPointData[1].vMax=delta.imag()+InitData.epsilon;
    hPointData[2].vMin=gamma.real()-InitData.epsilon;
    hPointData[2].vMax=gamma.real()+InitData.epsilon;
    hPointData[3].vMin=gamma.imag()-InitData.epsilon;
    hPointData[3].vMax=gamma.imag()+InitData.epsilon;
  }
  if(InitData.mode==2) {
    TCplxLong gamma=GammaAnal(1,InitData.tnx,InitData.tnz);
    TCplxLong delta=DeltaAnal(1,InitData.tnx,InitData.tnz);
    curDeltaRe=0; curGammaRe=0;
    InitData.PSOData.nComponents=2;
    hPointData[0].hV=&curDeltaIm;
    hPointData[0].vMin=delta.imag()-InitData.epsilon;
    hPointData[0].vMax=delta.imag()+InitData.epsilon;
    hPointData[1].hV=&curGammaIm;
    hPointData[1].vMin=gamma.imag()-InitData.epsilon;
    hPointData[1].vMax=gamma.imag()+InitData.epsilon;
  }
  if(InitData.mode==4) {
    TCplxLong gamma=GammaAnal(1,InitData.tnx,InitData.tnz);
    TCplxLong delta=DeltaAnal(1,InitData.tnx,InitData.tnz);
    curDeltaRe=0; curGammaRe=0;
    InitData.PSOData.nComponents=2;
    hPointData[0].hV=&curDeltaIm;
    hPointData[0].vMin=delta.imag()-InitData.epsilon;
    hPointData[0].vMax=delta.imag()+InitData.epsilon;
    hPointData[1].hV=&curGammaIm;
    hPointData[1].vMin=gamma.imag()-InitData.epsilon;
    hPointData[1].vMax=gamma.imag()+InitData.epsilon;
  }
  InitData.PSOData.hPointData=hPointData;
  PSOOptimizer.LoadInitData(&InitData.PSOData);

  if(InitData.SaveBitmapsInt) PrepareBitmapFolders();

  InitData.dt=1; InitData.tMax=-0.00000001+InitData.PSOData.nSteps*InitData.PSOData.nPoints;
  nTSteps=(InitData.tMax-0.00000001)/InitData.dt;
  t=0; tMax=InitData.PSOData.nPoints*InitData.PSOData.nSteps-0.5; dt=1;
  nTSteps=(InitData.tMax-0.00000001)/InitData.dt;
  if(InitData.mode==4) {
    tMax=InitData.PSOData.nPoints*InitData.PSOData.nSteps*InitData.nScaleSteps-0.5;
    nTSteps=(InitData.tMax-0.00000001)/InitData.dt;
  }
  tStep=0; fError=0; iScaleStep=0; iPSOStep=0; nBitmap=0; LastBitmapSaveTime=-1;
  SaveBitmapInterval=InitData.PSOData.nPoints;
}

TCplxLong cEquPSO2DBraggSolver::GammaAnal(int mode,int nx,int nz) {
  long double pi=M_PI,nxl=nx,nzl=nz,lx=InitData.Lx,lz=InitData.Lz,a=InitData.alpha;
  TCplxLong ii(0,1);
  if(mode==1) {
    if(nz==0) return ii/lz/a;
    if(nx==0) return -ii*lx*a;
    return ((TCplxLong(1,0)+ii*(nxl/nzl/lx-nzl/nxl/lz)/a)*nzl*lx/nxl/lz);
  }
  if (mode==2) {
    return(0);
  }
  return 0;
}

TCplxLong cEquPSO2DBraggSolver::DeltaAnal(int mode,int nx,int nz) {
  long double pi=M_PI,nxl=nx,nzl=nz,lx=InitData.Lx,lz=InitData.Lz,a=InitData.alpha;
  TCplxLong ii(0,1);
  if(mode==1) {
    return TCplxLong(pi*pi*nxl*nzl/2.0/a/lz/lx,pi*pi*(nxl*nxl/lx+nzl*nzl/lz)/2.0/a/a/lz/lx);
  }
  if (mode==2) {
    return(0);
  }
  return 0;
}

long double cEquPSO2DBraggSolver::MainFunction() {
  delta=TCplxLong(curDeltaRe,curDeltaIm); gamma=TCplxLong(curGammaRe,curGammaIm);
  if(abs(gamma)<0.0000000000001) { fError=1; return(NO_VALUE); }
  TCplxLong lambdaX=LambdaX(),lambdaZ=LambdaZ();
  if(abs(delta+lambdaX)<0.0000000000001) { fError=2; return(NO_VALUE); }
  if(abs(delta+lambdaZ)<0.0000000000001) { fError=3; return(NO_VALUE); }
  TCplxLong tmpx1=(delta-lambdaX)/(delta+lambdaX);
  TCplxLong tmpz1=(delta-lambdaZ)/(delta+lambdaZ);
  TCplxLong tmpx2=tmpx1*tmpx1;
  TCplxLong tmpx31=TCplxLong(0,2)*lambdaX*InitData.Lx;
  TCplxLong tmpx32=exp(tmpx31);
  TCplxLong tmpx3= tmpx32 - tmpx2;
  long double tmpx=abs(tmpx3);
//  long double tmpx=abs( exp(TCplxLong(0,2)*lambdaX*InitData.Lx)- tmpx1*tmpx1);
  long double tmpz=abs( exp(TCplxLong(0,2)*lambdaZ*InitData.Lz)- tmpz1*tmpz1);
  return -(tmpx+tmpz);
}

TCplxLong cEquPSO2DBraggSolver::LambdaX() {
  return(sqrt(delta*delta+2.0*InitData.alpha*delta/gamma));
}
TCplxLong cEquPSO2DBraggSolver::LambdaZ() {
  return(sqrt(delta*delta+2.0*InitData.alpha*delta*gamma));
}

void cEquPSO2DBraggSolver::StepRoutine() {
  long double FuncValue;
  SaveBitmaps();
  if (InitData.mode==4) { StepRoutineRepeat(); return; }
  if(PSOOptimizer.GetDoneFrac()<1.0) {
    PSOOptimizer.LoadFData(MainFunction());
  }
}

void cEquPSO2DBraggSolver::StepRoutineRepeat() {
  PSOOptimizer.LoadFData(MainFunction());
  if( (++iPSOStep) == InitData.PSOData.nPoints*InitData.PSOData.nSteps) {
    hResult[iScaleStep]=(*PSOOptimizer.GetGlobalMaxX(0))[0];
    iScaleStep++; iPSOStep=0;
    InitData.Lz=InitData.Lz0-(InitData.Lz0-InitData.LzMin)*iScaleStep/(InitData.nScaleSteps-1);
    InitData.Lx=InitData.Lx0*InitData.Lz/InitData.Lz0;
    hPointData[0].vMin=(*PSOOptimizer.GetGlobalMaxX(0))[0]*(1.0-InitData.epsilon);
    hPointData[0].vMax=(*PSOOptimizer.GetGlobalMaxX(0))[0]*(1.0+InitData.epsilon);
    hPointData[1].vMin=(*PSOOptimizer.GetGlobalMaxX(0))[1]-InitData.epsilon;
    hPointData[1].vMax=(*PSOOptimizer.GetGlobalMaxX(0))[1]+InitData.epsilon;
    PSOOptimizer.LoadInitData(&InitData.PSOData);
  }
}

void cEquPSO2DBraggSolver::SaveBitmaps() {
  int tmpInt;
  int FrameWidth=200;
  int FrameHeight=200;
  int FrameDist=20;
  int ix,iy;
  if(!InitData.SaveBitmapsInt) return;

  if( (LastBitmapSaveTime!=-1)&&
    (t-LastBitmapSaveTime < SaveBitmapInterval) ) return;
  LastBitmapSaveTime=t;

  AnsiString FileName=AnsiString(nBitmap);
  while(FileName.Length()<8) FileName="0"+FileName;
  FileName+=".bmp";

  Graphics::TBitmap *hBitmap = new Graphics::TBitmap;
  hBitmap->Width=FrameWidth*2+FrameDist; hBitmap->Height=FrameHeight;
  hBitmap->Canvas->MoveTo(0,0);
  hBitmap->Canvas->LineTo(FrameWidth,0);
  hBitmap->Canvas->LineTo(FrameWidth,FrameHeight);
  hBitmap->Canvas->LineTo(0,FrameHeight);
  hBitmap->Canvas->LineTo(0,0);
  hBitmap->Canvas->MoveTo(FrameWidth+FrameDist,0);
  hBitmap->Canvas->LineTo(2*FrameWidth+FrameDist,0);
  hBitmap->Canvas->LineTo(2*FrameWidth+FrameDist,FrameHeight);
  hBitmap->Canvas->LineTo(FrameWidth+FrameDist,FrameHeight);
  hBitmap->Canvas->LineTo(FrameWidth+FrameDist,0);
  for(int i=0;i<InitData.PSOData.nPoints;i++) {
    ix=FrameWidth*(PSOOptimizer.hPoints[i].X[0]-hPointData[0].vMin)
                  /(hPointData[0].vMax-hPointData[0].vMin);
    iy=FrameHeight*(PSOOptimizer.hPoints[i].X[1]-hPointData[1].vMin)
                  /(hPointData[1].vMax-hPointData[1].vMin);
    hBitmap->Canvas->Pixels[ix][iy]=clRed;

    ix=FrameWidth+FrameDist+
              FrameWidth*(PSOOptimizer.hPoints[i].X[2]-hPointData[2].vMin)
                  /(hPointData[2].vMax-hPointData[2].vMin);
    iy=FrameHeight*(PSOOptimizer.hPoints[i].X[3]-hPointData[3].vMin)
                  /(hPointData[3].vMax-hPointData[3].vMin);
    hBitmap->Canvas->Pixels[ix][iy]=clRed;
  }
  hBitmap->SaveToFile("Bitmaps\\GammaDelta\\"+FileName);
  delete hBitmap;
  nBitmap++;
}

void cEquPSO2DBraggSolver::PrepareBitmapFolders() {
  TSearchRec SearchRec;
  int findV;
  CreateDir("Bitmaps");
  CreateDir("Bitmaps\\GammaDelta");
  findV = FindFirst("Bitmaps\\GammaDelta\\*.*",0,SearchRec);
  while(!findV) {
    DeleteFile("Bitmaps\\GammaDelta\\"+SearchRec.Name);
    findV=FindNext(SearchRec);
  }
  FindClose(SearchRec);
}


