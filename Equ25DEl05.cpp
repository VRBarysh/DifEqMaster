//---------------------------------------------------------------------------


#pragma hdrstop

#include "Equ25DEl05.h"
#include <fstream.h>

//---------------------------------------------------------------------------

#pragma package(smart_init)

cEqu25DEl_PSOOptimizer::cEqu25DEl_PSOOptimizer() : cDifEqu_DataBankIO()
{
  hEquation=NULL; hVariable=NULL; hPointData=NULL; hKappaData=NULL; hMaxKappaData=NULL;
  hDeltaZData=NULL;
}

cEqu25DEl_PSOOptimizer::~cEqu25DEl_PSOOptimizer()// : ~cDifEqu_DataBankIO()
{
  delete hEquation;  hEquation=NULL; delete [] hPointData; hPointData=NULL;
  delete [] hKappaData; hKappaData=NULL;  delete [] hMaxKappaData; hMaxKappaData=NULL;
  delete [] hDeltaZData; hDeltaZData=NULL; 
}

void cEqu25DEl_PSOOptimizer::InitReport(cDifEquReport *hReport) {
  int Lines[EQU_MAX_GRAPHS];
  for(int i=0;i<EQU_MAX_GRAPHS;i++) Lines[i]=1;
  Lines[0]=1; Lines[1]=1; Lines[2]=1; Lines[3]=1;
  hReport->nVals=0; hReport->nGraphs=0;
//  hReport->LoadText(hReportText);
  hReport->SetGraphSize(EQU_GRAPH_SIZE,Lines);
  hReport->SetMaxTo1();
}

void cEqu25DEl_PSOOptimizer::Report(cDifEquReport *hReport) {
  cDifEqu_InitDataRec_25DElModel tmpInitData=InitData;
  hReport->nGraphs=1;
  if(InitData.PSOMode==2) {                   // Delta(Z)
    memcpy(hMaxKappaData,InitData.hDeltaZData,2*InitData.nDeltaZSteps*sizeof(long double));
    for(int i=0;i<InitData.nDeltaZSteps;i++)
       hMaxKappaData[2*i+1]=(*PSOOptimizer.GetGlobalMaxX())[i];
    for(int i=1;i<InitData.nDeltaZSteps-1;i++)
       hMaxKappaData[2*i]=(*PSOOptimizer.GetGlobalMaxX())[i+InitData.nDeltaZSteps-1];
    NormalizeDeltaZData(hMaxKappaData);
    InitData.hDeltaZData=hMaxKappaData;
    for(int i=0;i<InitData.nSteps;i++)
      MaxKappa[i]=hEquation->DeltaZEXT(InitData.dt*i,0,&InitData);
    hReport->Graph[0].GetData(0,InitData.nSteps,MaxKappa);
    InitData=tmpInitData;
  } else {                                    // Kappa
    for(int i=0;i<InitData.PSOData.nComponents;i++)
      hMaxKappaData[i]=(*PSOOptimizer.GetGlobalMaxX())[i];
    NormalizeKappaData(hMaxKappaData);
    if ((InitData.KappaMode==3)||(InitData.KappaMode==4))
      InitData.hKappaData=hMaxKappaData;
    if (InitData.KappaMode==1) {
      InitData.KappaV1=hMaxKappaData[0];
      InitData.KappaV2=hMaxKappaData[1];
      if(!(InitData.KappaV5<0)) InitData.KappaV4=hMaxKappaData[2];
    }
    for(int i=0;i<InitData.nSteps;i++)
      MaxKappa[i]=hEquation->KappaEXT(InitData.dt*i,0,&InitData);
    hReport->Graph[0].GetData(0,InitData.nSteps,MaxKappa);
    InitData=tmpInitData;
    InitData.hKappaData=hKappaData;
  }
  hReport->nVals=4;
  hReport->ValText[0]="Max A  "; hReport->Val[0]=PSOOptimizer.GetGlobalMaxVal()+FuckFunction(hMaxKappaData);
  hReport->ValText[1]="Max Kappa value ";
    hReport->Val[1]=hReport->Graph[0].Max;
  hReport->ValText[2]="Swarm Size "; hReport->Val[2]=PSOOptimizer.GetAvgSwarmSize();
  hReport->ValText[3]="Fuck Function "; hReport->Val[3]=FuckFunction(hMaxKappaData);
  if (InitData.KappaMode==1) {
    hReport->ValText[hReport->nVals]="Z0 "; hReport->Val[hReport->nVals++]=hMaxKappaData[0];
    hReport->ValText[hReport->nVals]="P  "; hReport->Val[hReport->nVals++]=hMaxKappaData[1];
  }
  if(InitData.PSOMode==1) {
    hReport->nVals=3;
    hReport->ValText[0]="Max A  "; hReport->Val[0]=PSOOptimizer.GetGlobalMaxVal();
    hReport->ValText[1]="V1  "; hReport->Val[1]=(*PSOOptimizer.GetGlobalMaxX())[0];
    hReport->ValText[2]="V2  "; hReport->Val[2]=(*PSOOptimizer.GetGlobalMaxX())[1];
  }
}

void cEqu25DEl_PSOOptimizer::FinalReport(cDifEquReport *hReport) {
  cDifEqu_InitDataRec_25DElModel tmpInitData=InitData;
  hReport->nGraphs=1;
  if(InitData.PSOMode==2) {                   // Delta(Z)
    memcpy(hMaxKappaData,InitData.hDeltaZData,2*sizeof(long double));
    for(int i=0;i<InitData.nDeltaZSteps;i++)
       hMaxKappaData[2*i+1]=(*PSOOptimizer.GetGlobalMaxX())[i];
    for(int i=1;i<InitData.nDeltaZSteps-1;i++)
       hMaxKappaData[2*i]=(*PSOOptimizer.GetGlobalMaxX())[i+InitData.nDeltaZSteps-1];
    NormalizeDeltaZData(hMaxKappaData);
    InitData.hDeltaZData=hMaxKappaData;
    for(int i=0;i<InitData.nSteps;i++)
      MaxKappa[i]=hEquation->DeltaZEXT(InitData.dt*i,0,&InitData);
    hReport->Graph[0].GetData(0,InitData.nSteps,MaxKappa);
    InitData=tmpInitData;
  } else {                                    // Kappa
    for(int i=0;i<InitData.PSOData.nComponents;i++)
      hMaxKappaData[i]=(*PSOOptimizer.GetGlobalMaxX())[i];
    NormalizeKappaData(hMaxKappaData);
    if ((InitData.KappaMode==3)||(InitData.KappaMode==4))
      InitData.hKappaData=hMaxKappaData;
    if (InitData.KappaMode==1) {
      InitData.KappaV1=hMaxKappaData[0];
      InitData.KappaV2=hMaxKappaData[1];
      if(!(InitData.KappaV5<0)) InitData.KappaV4=hMaxKappaData[2];
    }
    for(int i=0;i<InitData.nSteps;i++)
      MaxKappa[i]=hEquation->KappaEXT(InitData.dt*i,0,&InitData);
    hReport->Graph[0].GetData(0,InitData.nSteps,MaxKappa);
    InitData=tmpInitData;
    InitData.hKappaData=hKappaData;
  }
  hReport->nVals=4;
  hReport->ValText[0]="Max A  "; hReport->Val[0]=PSOOptimizer.GetGlobalMaxVal()+FuckFunction(hMaxKappaData);
  hReport->ValText[1]="Max Kappa value ";
  hReport->ValText[2]="Swarm Size "; hReport->Val[2]=PSOOptimizer.GetAvgSwarmSize();
  hReport->ValText[3]="Fuck Function "; hReport->Val[3]=FuckFunction(hMaxKappaData);
  hReport->Val[1]=hReport->Graph[0].Max;
  if (InitData.KappaMode==1) {
    hReport->ValText[hReport->nVals]="Z0 "; hReport->Val[hReport->nVals++]=hMaxKappaData[0];
    hReport->ValText[hReport->nVals]="P  "; hReport->Val[hReport->nVals++]=hMaxKappaData[1];
  }
  ofstream FileOutKappa = ofstream("Kappa.dat",ios::out);
  ofstream FileOutDeltaZ = ofstream("DeltaZ.dat",ios::out);
  long double value;
  if(InitData.PSOMode==2) {
    for(int i=0;i<InitData.nDeltaZSteps;i++) {
      value=hMaxKappaData[2*i];
      FileOutDeltaZ << value << " " << hMaxKappaData[2*i+1] << endl;
    }
  } else {
    for(int i=0;i<InitData.nSteps;i++) {
      value=InitData.dt*i;
      FileOutKappa << value << " " << MaxKappa[i] << endl;
    }
  }

/*
  hReport->nGraphs=0;
  hReport->nVals=3;
  hReport->ValText[0]="Max A  "; hReport->Val[0]=PSOOptimizer.GetGlobalMaxVal();
  hReport->ValText[1]="Max Z0 "; hReport->Val[1]=(*PSOOptimizer.GetGlobalMaxX())[0];
  hReport->ValText[2]="Max P  "; hReport->Val[2]=(*PSOOptimizer.GetGlobalMaxX())[1];
  */
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
  if(InitData.PSOMode==1) {
    hReport->nVals=3;
    hReport->ValText[0]="Max A  "; hReport->Val[0]=PSOOptimizer.GetGlobalMaxVal();
    hReport->ValText[1]="V1  "; hReport->Val[1]=(*PSOOptimizer.GetGlobalMaxX())[0];
    hReport->ValText[2]="V2  "; hReport->Val[2]=(*PSOOptimizer.GetGlobalMaxX())[1];
  }
}

void cEqu25DEl_PSOOptimizer::LoadInitData(cDifEqu_InitDataRec *hInitData)
{
  cDifEqu_InitDataRec_25DElModel *hInit=(cDifEqu_InitDataRec_25DElModel *)hInitData;
  InitData=*hInit;
  cDifEqu_InitDataRec_25DElModel tmpInitData=InitData;
  if(hPointData) delete [] hPointData;
  if(hKappaData) delete [] hKappaData;
  if(hMaxKappaData) delete [] hMaxKappaData;
  InitData.PSOData.nComponents=InitData.KappaV1;
  if (InitData.KappaMode==1) 
    if (InitData.KappaV5<0) {
      InitData.PSOData.nComponents=2;
    } else {
      InitData.PSOData.nComponents=3;
    }
  if(InitData.PSOMode==2) {
    InitData.PSOData.nComponents=2*(InitData.nDeltaZSteps);
    hDeltaZData = new long double[2*InitData.nDeltaZSteps];
  }
  InitData.nKappaSteps=InitData.PSOData.nComponents;
  hPointData = new CPSOInitComp[InitData.PSOData.nComponents];
  hKappaData = new long double[InitData.PSOData.nComponents];
  hMaxKappaData = new long double[InitData.PSOData.nComponents];
  InitData.hKappaData=hKappaData;
  InitData.hDeltaZData=hDeltaZData;
  if( (InitData.KappaMode==3)||(InitData.KappaMode==4) ) {
    for(int i=0;i<InitData.PSOData.nComponents;i++) {
      hPointData[i].vMin=0;
      hPointData[i].vMax=1;
      hPointData[i].hV=hKappaData+i;
    }
  } else if (InitData.KappaMode==1) {
    if (InitData.KappaV5<0) {
      InitData.PSOData.nComponents=2;
      hPointData[0].vMin=0;
      hPointData[0].vMax=InitData.dt*(InitData.nSteps-1);
      hPointData[0].hV=&InitData.KappaV1;
      hPointData[1].vMin=1.1;
      hPointData[1].vMax=100;
      hPointData[1].hV=&InitData.KappaV2;
    } else {
      InitData.PSOData.nComponents=3;
      hPointData[0].vMin=0;
      hPointData[0].vMax=InitData.dt*(InitData.nSteps-1);
      hPointData[0].hV=&InitData.KappaV1;
      hPointData[1].vMin=1.1;
      hPointData[1].vMax=InitData.KappaV6;
      hPointData[1].hV=&InitData.KappaV2;
      hPointData[2].vMin=0;
      hPointData[2].vMax=InitData.KappaV5;
      hPointData[2].hV=&InitData.KappaV4;
    }
  }
//  hPointData[0].hV=&(InitData.KappaV1);
//  hPointData[0].vMin=0;
//  hPointData[0].vMax=InitData.dt*(InitData.nSteps-1);
//  hPointData[1].hV=&(InitData.KappaV2);
//  hPointData[1].vMin=1;
//  hPointData[1].vMax=InitData.KappaV3;
  if(InitData.PSOMode==2) {               // Optimize by delta(Z)
    InitData.PSOData.nComponents=2*(InitData.nDeltaZSteps-1);
    for(int i=0;i<InitData.nDeltaZSteps;i++) {
      hDeltaZData[2*i+1]=0;
      hDeltaZData[2*i]=i*(InitData.nSteps+1)*InitData.dt/(InitData.nDeltaZSteps-1);
      hPointData[i].vMin=-InitData.DeltaZMax;
      hPointData[i].vMax=InitData.DeltaZMax;
      hPointData[i].hV=(hDeltaZData+2*i+1);
    }
    for(int i=1;i<InitData.nDeltaZSteps-1;i++) {
      hPointData[i+InitData.nDeltaZSteps-1].vMin=0;
      hPointData[i+InitData.nDeltaZSteps-1].vMax=1;
      hPointData[i+InitData.nDeltaZSteps-1].hV=(hDeltaZData+2*i);
    }
  }
  if(InitData.PSOMode==1) {               // Optimize by delta
    InitData.PSOData.nComponents=2;
    hPointData[0].vMin=-InitData.deltaV1;
    hPointData[0].vMax=InitData.deltaV1;
    hPointData[0].hV=&InitData.deltaV1;
    hPointData[1].vMin=-InitData.deltaV2;
    hPointData[1].vMax=InitData.deltaV2;
    hPointData[1].hV=&InitData.deltaV2;
  }
  InitData.PSOData.hPointData=hPointData;
  PSOOptimizer.LoadInitData(&InitData.PSOData);
//  if(InitData.KappaV6>0) {               //  ----------------- !
//    InitData.Nu=-(InitData.KappaV6*expl(-0.66*logl(InitData.KappaV4)));
//  }
  if(InitData.KappaV3>1) {
    TVecPSO *hMax=PSOOptimizer.GetGlobalMaxX();
    PSOOptimizer.SetGlobalMaxVal(InitData.KappaV3);
    if((InitData.KappaV8>0)&&(InitData.KappaV9>0)) {
      tmpInitData.KappaMode=1;
      tmpInitData.KappaV1=InitData.KappaV8;
      tmpInitData.KappaV2=InitData.KappaV9;
      tmpInitData.KappaV4=1.0;
      for(int i=0;i<InitData.PSOData.nComponents;i++)
        (*hMax)[i]=hEquation->KappaEXT(
          InitData.dt*i*(InitData.nSteps-1)/(InitData.PSOData.nComponents-1),
          0,&tmpInitData);
    } else
      for(int i=0;i<InitData.PSOData.nComponents;i++) (*hMax)[i]=1;
  }
  nTSteps=(InitData.tMax-0.00000001)/InitData.dt;
  t=0; tMax=nTSteps*InitData.PSOData.nPoints*InitData.PSOData.nSteps-0.5; dt=1;
  tStep=0;
  BestOutTime=BestTotalTime=0;
  MaxOut=0;
  hEquation->LoadInitData(&InitData);
}

void cEqu25DEl_PSOOptimizer::StepRoutine() {
  int omg,imax=0;
  long double tmp;
  hEquation->Step();
  tStep++;
  if(hEquation->GetOutA()>MaxOut) {
    MaxOut=hEquation->GetOutA();
    imax=tStep;
    if(((MaxOut-FuckFunction(hKappaData))>InitData.KappaV7)&&(InitData.KappaV7>0)) {
      omg=3;
    }
  }
  if( (hEquation->GetOutA()<MaxOut)&&(MaxOut>1)&&(InitData.fNeed1stPeak)) tStep=nTSteps;
//  long double tmpMax=hEquation->GetMaxA();
  if(tStep==nTSteps) {
    if(((MaxOut-FuckFunction(hKappaData))>InitData.KappaV7)&&(InitData.KappaV7>0)) {
      omg=3;
    }
//    if(InitData.KappaV6>0)
    PSOOptimizer.LoadFData(MaxOut-FuckFunction(hKappaData));
//    if(InitData.KappaV6>0) {
//      InitData.Nu=-(InitData.KappaV6*expl(-0.66*logl(InitData.KappaV4)));
//    }
    NormalizeKappaData(hKappaData);
    NormalizeDeltaZData(hDeltaZData);
    tStep=0; MaxOut=0;
    hEquation->LoadInitData(&InitData);
  }
}

void cEqu25DEl_PSOOptimizer::NormalizeDeltaZData(long double *hData) {
  long double L=(InitData.nSteps)*InitData.dt;
  long double tmpSum=0,L2;
  if(!(InitData.PSOMode==2)) return;
  for(int i=1;i<InitData.nDeltaZSteps-1;i++) tmpSum+=hData[2*i];
  L2=L*tmpSum/(InitData.nDeltaZSteps-2);
  if(tmpSum==0) {
    for(int i=1;i<InitData.nDeltaZSteps-1;i++) hData[2*i]=i*InitData.dt;
    return;
  }
  for(int i=1;i<InitData.nDeltaZSteps-1;i++)
    hData[2*i]=hData[2*i-2]+L2*hData[2*i]/tmpSum;
  /*
  for(int i=1;i<InitData.nDeltaZSteps-1;i++) {
    hData[2*i]=hData[2*i-2]+hData[2*i]*(L-hData[2*i-2]);
  }
  */
}

void cEqu25DEl_PSOOptimizer::NormalizeKappaData(long double *hData) {
  long double tmpV=0;
  if( (InitData.KappaV2==1)&&((InitData.KappaMode==3)||(InitData.KappaMode==4)) ) {
    for(int i=0;i<InitData.PSOData.nComponents;i++)
      if(tmpV<hData[i]) tmpV=hData[i];
    tmpV=1.0/tmpV;
    for(int i=0;i<InitData.PSOData.nComponents;i++)
      hData[i]*=tmpV;
  }
  if(InitData.KappaV2==2) {
    if(InitData.KappaMode==3) {
      tmpV=-0.5*(hData[0]+hData[InitData.PSOData.nComponents-1]);
      for(int i=0;i<InitData.PSOData.nComponents;i++) tmpV+=hData[i];
      tmpV=1.0/tmpV;
      tmpV*=(InitData.PSOData.nComponents-1);
      for(int i=0;i<InitData.PSOData.nComponents;i++) hData[i]*=tmpV;
    }
    if(InitData.KappaMode==4) {
      for(int i=0;i<InitData.PSOData.nComponents;i++) tmpV+=hData[i];
      tmpV=1.0/tmpV;
      tmpV*=(InitData.PSOData.nComponents);
      for(int i=0;i<InitData.PSOData.nComponents;i++) hData[i]*=tmpV;
    }
    if(InitData.KappaMode==6) {
      for(int i=0;i<InitData.nSteps;i++) {
        MaxKappa[i]=hEquation->KappaEXT(InitData.dt*i,0,&InitData);
        if(tmpV<MaxKappa[i]) tmpV=MaxKappa[i];
      }
      if (tmpV<0) {
        tmpV=1.0/(1.0-tmpV);
        for(int i=0;i<InitData.PSOData.nComponents;i++) hData[i]*=tmpV;
      }
    }
  }
  if(InitData.KappaV2==3) {
    if(InitData.KappaMode==3) {
      for(int i=InitData.PSOData.nComponents-1;i>0;i--) hData[i-1]*=hData[i];
    }
  }
}

long double cEqu25DEl_PSOOptimizer::FuckFunction(long double *hData){
  long double tmpFuck=0,tmpFuck1=0;
    if(InitData.KappaMode==3) {
      for(int i=0;i<InitData.PSOData.nComponents-2;i++) {
        tmpFuck+=(hData[i+2]-2.0*hData[i+1]+hData[i])*
                  (hData[i+2]-2.0*hData[i+1]+hData[i]);
        tmpFuck1+=(hData[i+2]-2.0*hData[i+1]+hData[i])*
                  (hData[i+2]-2.0*hData[i+1]+hData[i])*
                  (hData[i+2]-2.0*hData[i+1]+hData[i])*
                  (hData[i+2]-2.0*hData[i+1]+hData[i]);
      }
    }
    if(InitData.KappaMode==4) {
      for(int i=0;i<InitData.PSOData.nComponents-1;i++) {
        tmpFuck+=(hData[i+1]-hData[i])*(hData[i+1]-hData[i]);
        tmpFuck1+=(hData[i+1]-hData[i])*(hData[i+1]-hData[i])*(hData[i+1]-hData[i])*(hData[i+1]-hData[i]);
      }
    }
    if(InitData.KappaMode==6) {
    }
    return tmpFuck*InitData.KappaV4+tmpFuck1*InitData.KappaV5;
}


