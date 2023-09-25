//---------------------------------------------------------------------------


#pragma hdrstop

#include <fstream> 
#include "Equ2WaveElModel.h"

#define MAX_DUMP_SIZE 268435457

//---------------------------------------------------------------------------

#pragma package(smart_init)

void cEqu2WaveElModel::InitReport(cDifEquReport *hReport) {
  int Lines[EQU_MAX_GRAPHS];
  for(int i=0;i<EQU_MAX_GRAPHS;i++) Lines[i]=1;
  Lines[0]=1; Lines[1]=1; Lines[2]=1; Lines[3]=1; Lines[4]=1;
  hReport->LoadText(hReportText);
  hReport->SetGraphSize(EQU_GRAPH_SIZE,Lines);
  hReport->SetMaxTo1();
}
void cEqu2WaveElModel::Report(cDifEquReport *hReport) {

  hReport->Graph[0].GetDataAbs(0,InitData.nSteps,hAp);
  hReport->Graph[1].GetDataAbs(0,InitData.nSteps,hAm);
  hReport->Graph[2].GetDataAbs(0,InitData.nStepsT-1,hOutA);
  hReport->Graph[3].GetDataAbs(0,InitData.nSteps,hI);
  hReport->Graph[4].GetDataArg(0,InitData.nSteps,hAp);
  hReport->Graph[5].GetDataArg(0,InitData.nStepsT-1,hOutA);
  hReport->Graph[6].GetDataAbs(0,InitData.nStepsT-1,hOutA+InitData.nStepsT);
  hReport->Graph[4].Max=M_PI;   hReport->Graph[5].Max=M_PI;

  hReport->nVals=0;

  hReport->ValText[hReport->nVals]="Max A+(z) ";
  hReport->Val[hReport->nVals++]=hReport->Graph[0].Max;
  hReport->ValText[hReport->nVals]="Max A-(z) ";
  hReport->Val[hReport->nVals++]=hReport->Graph[1].Max;
  hReport->ValText[hReport->nVals]="Max A-(t,z=L) ";
  hReport->Val[hReport->nVals++]=hReport->Graph[2].Max;
  hReport->ValText[hReport->nVals]="Max I(z) ";
  hReport->Val[hReport->nVals++]=hReport->Graph[3].Max;

  CalcEnergy();

  hReport->ValText[hReport->nVals]="EÀ+ ";
  hReport->Val[hReport->nVals++]=EAp;
  hReport->ValText[hReport->nVals]="EÀ- ";
  hReport->Val[hReport->nVals++]=EAm;
  hReport->ValText[hReport->nVals]="EOut ";
  hReport->Val[hReport->nVals++]=EOut;
  hReport->ValText[hReport->nVals]="EE ";
  hReport->Val[hReport->nVals++]=EE;
  hReport->ValText[hReport->nVals]="ESumA ";
  hReport->Val[hReport->nVals++]=EOut+EAp+EAm;
  hReport->ValText[hReport->nVals]="ESum ";
  hReport->Val[hReport->nVals++]=EOut+EAp+EAm+EE;
  hReport->ValText[hReport->nVals]="t = ";
  hReport->Val[hReport->nVals++]=t;
  /*
  hReport->Graph[4].GetDataRe(0,10,hApInter);
  hReport->Graph[5].GetDataIm(0,10,hApInter);
  hReport->Graph[6].GetDataRe(0,10,hAmInter);
  hReport->Graph[7].GetDataIm(0,10,hAmInter);
  */
}
void cEqu2WaveElModel::FinalReport(cDifEquReport *hReport) {
  hReport->nVals=0; hReport->nGraphs=3;
  hReport->Graph[0].GetDataAbs(0,InitData.nSteps,hAp);
  hReport->Graph[1].GetDataAbs(0,InitData.nSteps,hAm);
  hReport->Graph[2].GetDataAbs(0,InitData.nStepsT,hOutA);

  hReport->ValText[hReport->nVals]="Max A+(z) ";
  hReport->Val[hReport->nVals++]=hReport->Graph[0].Max;
  hReport->ValText[hReport->nVals]="Max A-(z) ";
  hReport->Val[hReport->nVals++]=hReport->Graph[1].Max;
  hReport->ValText[hReport->nVals]="Max A-(t,z=L) ";
  hReport->Val[hReport->nVals++]=hReport->Graph[2].Max;

  CalcEnergy();
  hReport->ValText[hReport->nVals]="EÀ+ ";
  hReport->Val[hReport->nVals++]=EAp;
  hReport->ValText[hReport->nVals]="EÀ- ";
  hReport->Val[hReport->nVals++]=EAm;
  hReport->ValText[hReport->nVals]="EOut ";
  hReport->Val[hReport->nVals++]=EOut;
  hReport->ValText[hReport->nVals]="EE ";
  hReport->Val[hReport->nVals++]=EE;
  hReport->ValText[hReport->nVals]="ESumA ";
  hReport->Val[hReport->nVals++]=EOut+EAp+EAm;
  hReport->ValText[hReport->nVals]="ESum ";
  hReport->Val[hReport->nVals++]=EOut+EAp+EAm+EE;

  FinalSave();
}

void cEqu2WaveElModel::FinalSave() {
  long double tmp1;
  int tmpi;

  ofstream FileOutA = ofstream("OutA.dat",ios::out);
  ofstream FileAp = ofstream("Ap.dat",ios::out);
  ofstream FileAm = ofstream("Am.dat",ios::out);
  ofstream FileApzt = ofstream("Apzt.dat",ios::out);
  ofstream FileAmzt = ofstream("Amzt.dat",ios::out);
  ofstream FileApt = ofstream("Apt.dat",ios::out);
  ofstream FileApArg = ofstream("ApArg.dat",ios::out);
  ofstream FileE = ofstream("E.dat",ios::out);
  REAL_TYPE tmpArg;
  complex<float> backA;
  for(int i=0;i<InitData.nSteps;i++) {
    tmpArg=0;
    if(abs(hAp[i])>0.0000000001) tmpArg=arg(hAp[i]);
    FileAp << InitData.dt*i << " " << abs(hAp[i]) << endl;
    FileApArg << InitData.dt*i << " " << tmpArg << endl;
    FileAm << InitData.dt*i << " " << abs(hAm[i]) << endl;
  }
  for(int i=0;i<nStepNumber;i++)
    FileOutA << InitData.dt*i << " " << abs(hOutA[i])
                              << " " << abs(hOutA[i+InitData.nStepsT]) << endl;

  if( (hApDump)&&(hAmDump) ) {
    for(int it=0;it<nStepNumber;it++) {
      if(InitData.TSliceSaveZ>=0) {
        tmpi = InitData.TSliceSaveZ < InitData.nSteps ? InitData.TSliceSaveZ : InitData.nSteps - 1;
        FileApt << InitData.dt*it << " " << hApDump[InitData.nSteps*it+tmpi] << endl;
      }
      FileE << InitData.dt*it << " " << hEnergyData[it].Ap
                              << " " << hEnergyData[it].Am
                              << " " << hEnergyData[it].Out
                              << " " << hEnergyData[it].E
                              << " " << hEnergyData[it].Ap+
                                        hEnergyData[it].Am+
                                        hEnergyData[it].Out+
                                        hEnergyData[it].E << endl;
      backA=hApDump[InitData.nSteps*it];
      for(int iz=0;iz<InitData.nSteps;iz++) {
        tmpArg = 0;
        if (abs(hApDump[InitData.nSteps*it+iz]*conj(backA))>0.000001)
          tmpArg=arg(hApDump[InitData.nSteps*it+iz]*conj(backA));
        FileApzt << InitData.dt*iz << " "<< InitData.dt*it << " " <<
          abs(hApDump[InitData.nSteps*it+iz]) << " " << tmpArg << endl;
        tmpArg = 0;
        if (abs(hAmDump[InitData.nSteps*it+iz])>0.000001)
          tmpArg=arg(hAmDump[InitData.nSteps*it+iz]);
        FileAmzt << InitData.dt*iz << " "<< InitData.dt*it << " " <<
          abs(hAmDump[InitData.nSteps*it+iz]) << " " << tmpArg << endl;
        backA=hApDump[InitData.nSteps*it+iz];
      }
    }
  }
}

void cEqu2WaveElModel::LoadInitData(cDifEqu_InitDataRec *hInitData) {
  cDifEqu_InitDataRec_2WaveElModel *hInit=
             (cDifEqu_InitDataRec_2WaveElModel *)hInitData;
  hInit->nStepsT=(hInit->tMax/hInit->dt)+1;
  if(!hTaskCritSection) hTaskCritSection = new TCriticalSection;
  if(!hTaskEvent) hTaskEvent = new TEvent(NULL,true,false,"");
  if((hInit->nSteps+hInit->nStepsT)!=(InitData.nSteps+InitData.nStepsT)) {
    delete [] hReportBuf; delete [] htmpAp; delete [] htmpAm;
    delete [] hAp; delete [] hAm; delete [] hAhalf; delete [] hEl; delete [] hI;
    for(int i=0;i<6;i++) delete [] hSCoef[i]; delete [] hZ;
    delete [] hApInter; delete [] hAmInter;
    delete [] hOutA; delete [] hAmDump; delete [] hApDump;
    delete [] hEnergyData;
    hAp = new complex<REAL_TYPE>[5*(hInit->nSteps+1)];
    hAm = new complex<REAL_TYPE>[5*(hInit->nSteps+1)];
    hAhalf = new complex<REAL_TYPE>[hInit->nSteps+1];
    hI = new complex<REAL_TYPE>[hInit->nSteps+1];
    htmpAp = new complex<REAL_TYPE>[hInit->nSteps+1];
    htmpAm = new complex<REAL_TYPE>[hInit->nSteps+1];
    hApInter = new complex<REAL_TYPE>[hInit->nSteps+1];
    hAmInter = new complex<REAL_TYPE>[hInit->nSteps+1];
    hEnergyData = new cEnergyData_2WaveEl[hInit->nStepsT];
    if(4*(hInit->nSteps*hInit->nStepsT)<268435457) {
      hAmDump = new complex<float>[hInit->nSteps*hInit->nStepsT];
      hApDump = new complex<float>[hInit->nSteps*hInit->nStepsT];
    }
    for(int i=0;i<6;i++)
      hSCoef[i]=new REAL_TYPE[ 2*(hInit->nSteps) ];
    hZ = new REAL_TYPE[hInit->nSteps];
    hReportBuf = new REAL_TYPE[hInit->nSteps+hInit->nStepsT+1];
    hOutA = new complex<REAL_TYPE>[hInit->nStepsT*2];
  }
  if((hInit->nSteps*hInit->nEls)!=(InitData.nSteps*InitData.nEls)) {
    hEl   = new REAL_TYPE[hInit->nSteps*hInit->nEls*8];
  }

  for(int i=0;i<hInit->nThreads+1-nEqus;i++) {
    hEqus[i] = new cEqu2WaveBaseTask;
  }
  for(int i=0;i<hInit->nThreads-InitData.nThreads;i++) {
    hThreads[i] = new TTaskThread(true);
    hThreads[i]->Priority = tpLowest;
    /*
        tpIdle	The thread executes only when the system is idle. The system will not interrupt other threads to execute a thread with tpIdle priority.
        tpLowest	The thread's priority is two points below normal.
        tpLower	The thread's priority is one point below normal.
        tpNormal	The thread has normal priority.
        tpHigher	The thread's priority is one point above normal.
        tpHighest	The thread's priority is two points above normal.
        tpTimeCritical	The thread gets highest priority.
    */
  }
  InitData=*hInit; nEqus=InitData.nThreads+1;
  int zPart=InitData.nSteps/nEqus;
  int zRest=InitData.nSteps - nEqus*zPart;
  for(int i=0;i<nEqus;i++) {
    hEqus[i]->InitData=InitData;
    hEqus[i]->dtT=InitData.dt;
    hEqus[i]->hAp=hAp;  hEqus[i]->hAm=hAm; hEqus[i]->htmpAp=htmpAp; 
    hEqus[i]->hEl=hEl;  hEqus[i]->hI=hI;
    hEqus[i]->iz1st=i*zPart + (i<zRest?i:zRest);
    hEqus[i]->izLast=(i+1)*zPart + (i<zRest?(i+1):zRest)-1;
  }
  for(int i=0;i<InitData.nThreads;i++) {
    hThreads[i]->Init(hEqus[i+1],hTaskEvent,hTaskCritSection,&ThreadCounter);
  }
  t=0; dt=InitData.dt; tMax=InitData.tMax;
  PrepareMedia();

  nStepNumber=0; t=0; EOut=0;
}


void cEqu2WaveElModel::PrepareMedia() {
  REAL_TYPE randVal;
  if(InitData.RandSeed) srand(InitData.RandSeed);
  for(int iz=0;iz<InitData.nSteps;iz++) {
    hAp[iz]=sin(M_PI*iz/InitData.nSteps); hAm[iz]=0;  hI[iz]=0; hAp[iz]=0;
    randVal=InitData.rndCoef*M_PI*rand()/RAND_MAX;
    for(int ie=0;ie<InitData.nEls;ie++) {
      long double Theta0=2.0*M_PI*ie/InitData.nEls;
      hEl[iz*InitData.nEls*2+ie*2]=Theta0+InitData.r*cos(randVal+Theta0);
      hEl[iz*InitData.nEls*2+ie*2+1]=InitData.deltaV1;
    }
  }
  for(int it=0;it<InitData.nStepsT;it++) hOutA[it]=0;
}

void cEqu2WaveElModel::CalcEnergy() {
  EAp=EAm=EE=0;
  for(int iz=0;iz<InitData.nSteps;iz++) {
    EAm+=abs(hAm[iz])*abs(hAm[iz])*InitData.dt;
    EAp+=abs(hAp[iz])*abs(hAp[iz])*InitData.dt;
    for(int ie=0;ie<InitData.nEls;ie++) {
      EE+=hEl[iz*InitData.nEls*2+2*ie+1]-InitData.deltaV1;
    }
  }
  EE*=-2.0*InitData.dt/(InitData.nEls);
}

void cEqu2WaveElModel::StepRoutine() {
  int ASize = (InitData.nSteps);
  REAL_TYPE dt2=InitData.dt*0.5,dt1=InitData.dt;
  REAL_TYPE onesix = 1.0/6.0;
  complex<REAL_TYPE> tcA[10];
  complex<REAL_TYPE> ai=complex<REAL_TYPE>(0,InitData.ACoef),tmpA,tmpA1,nextAm;

//  hOutA[nStepNumber]=hAm[InitData.nSteps];
//  EOut+=abs(hAm[InitData.nSteps])*abs(hAm[InitData.nSteps])*dt1;
  CalcEnergy();
  hEnergyData[nStepNumber].E=EE;
  hEnergyData[nStepNumber].Ap=EAp;
  hEnergyData[nStepNumber].Am=EAm;
  hEnergyData[nStepNumber].Out=EOut;

  hOutA[InitData.nStepsT+nStepNumber]=hAp[InitData.nSteps-1];

  DoTasks();
  for(int iz=0;iz<InitData.nSteps;iz++) {
    //hI[iz]=0;
    if( (hApDump)&&(hAmDump) ) {
      hAmDump[InitData.nSteps*nStepNumber+iz]=hAm[iz];
      hApDump[InitData.nSteps*nStepNumber+iz]=hAp[iz];
    }
    hAp[ASize+iz]=(hAm[iz]*ai+Fz(iz*InitData.dt)*hI[iz])*dt2;
    hAm[ASize+iz]=(hAp[iz]*ai)*dt2;
    htmpAp[iz] = hAp[iz]+hAp[ASize+iz];
    htmpAm[iz] = hAm[iz]+hAm[ASize+iz];
  }
  /*                           //  ----------------------- Euler Part
  for(int iz=0;iz<InitData.nSteps;iz++) {
    hAp[iz]=hAp[iz]+hAp[ASize+iz]*2.0;
    hAm[iz+1]=htmpAm[iz]+hAm[ASize+iz];
  }
  hAm[0]=0;
  hOutA[nStepNumber]=hAm[InitData.nSteps];
  EOut+=abs(hOutA[nStepNumber])*abs(hOutA[nStepNumber])*dt1;
  return;
  */

  /*                                               // ----------- Runge 2
  hOutA[nStepNumber]=hAm[InitData.nSteps-1];
  EOut+=abs(hOutA[nStepNumber])*abs(hOutA[nStepNumber])*dt1;
  for(int iz=0;iz<InitData.nSteps;iz++) {
    tmpA=0; if(iz) tmpA=(htmpAm[iz-1]+hAm[ASize+iz-1]);//SplineAm.GetYIterC(dt1*(-0.5+iz));
    hAp[2*ASize+iz]=(tmpA*ai+Fz(iz*dt1)*hI[iz])*dt2;
    tmpA=0; if(iz!=InitData.nSteps-1) tmpA=(htmpAp[iz+1]+hAp[ASize+iz+1]);
    hAm[2*ASize+iz]=(tmpA)*ai*dt2;
    htmpAp[iz] = hAp[iz]+hAp[ASize+iz]+hAp[2*ASize+iz];
    htmpAm[iz] = hAm[iz]+hAm[ASize+iz]+hAm[2*ASize+iz];
  }
  hAm[0]=0;
  for(int iz=0;iz<InitData.nSteps;iz++) {
    hAp[iz]=htmpAp[iz];
    hAm[iz+1]=htmpAm[iz];
  }
  return;
  */
/*
  for(int iz=0;iz<InitData.nSteps;iz++) {
    tmpA=(htmpAm[iz])*0.5; if(iz) tmpA+=(htmpAm[iz-1])*0.5;
    hAp[2*ASize+iz]=(tmpA*ai+Fz(iz*dt1)*hI[iz])*dt2;
    tmpA=(htmpAp[iz])*0.5; if(iz!=InitData.nSteps-1) tmpA+=(htmpAp[iz+1])*0.5;
    hAm[2*ASize+iz]=(tmpA)*ai*dt2;
    htmpAp[iz] = hAp[iz]+hAp[2*ASize+iz];
    htmpAm[iz] = hAm[iz]+hAm[2*ASize+iz];
  }
  for(int iz=0;iz<InitData.nSteps;iz++) {
    tmpA=(htmpAm[iz])*0.5; if(iz) tmpA+=(htmpAm[iz-1])*0.5;
    hAp[3*ASize+iz]=(tmpA*ai+Fz(iz*dt1)*hI[iz])*dt1;
    tmpA=(htmpAp[iz])*0.5; if(iz!=InitData.nSteps-1) tmpA+=(htmpAp[iz+1])*0.5;
    hAm[3*ASize+iz]=(tmpA)*ai*dt1;
    htmpAp[iz] = hAp[iz]+hAp[3*ASize+iz];
    htmpAm[iz] = hAm[iz]+hAm[3*ASize+iz];
  }
*/
  DoTasks();
//  SplineAm.MakeSpline(InitData.nSteps,htmpAm,hSCoef[0],hSCoef[1],hSCoef[2],dt,NULL);
//  SplineAp.MakeSpline(InitData.nSteps,htmpAp,hSCoef[3],hSCoef[4],hSCoef[5],dt,NULL);
  for(int iz=0;iz<InitData.nSteps;iz++) {
//    tmpA=(htmpAp[iz])*0.5; if(iz!=InitData.nSteps-1) tmpA=SplineAp.GetYIterC(dt1*(0.5+iz));
    tmpA1=(htmpAp[iz])*0.5; if(iz!=InitData.nSteps-1) tmpA1+=(htmpAp[iz+1])*0.5;
    hApInter[iz]=tmpA1;
//    tmpA=(htmpAm[iz])*0.5; if(iz) tmpA=SplineAm.GetYIterC(dt1*(-0.5+iz));
    tmpA1=(htmpAm[iz])*0.5; if(iz) tmpA1+=(htmpAm[iz-1])*0.5;
    hAmInter[iz]=tmpA1;
  }
  //nextAm = (htmpAm[0])*0.5;
  for(int iz=0;iz<InitData.nSteps;iz++) {
//    tmpA=(htmpAm[iz])*0.5; if(iz) tmpA=SplineAm.GetYIterC(dt1*(-0.5+iz));
//    tmpA1=(htmpAm[iz])*0.5; if(iz) tmpA1+=(htmpAm[iz-1])*0.5;
    hAp[2*ASize+iz]=(hAmInter[iz]*ai+Fz(iz*dt1)*hI[iz])*dt2;
//    nextAm=SplineAm.GetYIterC(dt1*(0.5+iz));

//    tmpA=(htmpAp[iz])*0.5; if(iz!=InitData.nSteps-1) tmpA=SplineAp.GetYIterC(dt1*(0.5+iz));
//    tmpA1=(htmpAp[iz])*0.5; if(iz!=InitData.nSteps-1) tmpA1+=(htmpAp[iz+1])*0.5;
    hAm[2*ASize+iz]=(hApInter[iz])*ai*dt2;

    htmpAp[iz] = hAp[iz]+hAp[2*ASize+iz];
    htmpAm[iz] = hAm[iz]+hAm[2*ASize+iz];
  }

  DoTasks();
//  SplineAm.MakeSpline(InitData.nSteps,htmpAm,hSCoef[0],hSCoef[1],hSCoef[2],dt,NULL);
//  SplineAp.MakeSpline(InitData.nSteps,htmpAp,hSCoef[3],hSCoef[4],hSCoef[5],dt,NULL);
  for(int iz=0;iz<InitData.nSteps;iz++) {
//    tmpA=(htmpAp[iz])*0.5; if(iz!=InitData.nSteps-1) tmpA=SplineAp.GetYIterC(dt1*(0.5+iz));
    tmpA1=(htmpAp[iz])*0.5; if(iz!=InitData.nSteps-1) tmpA1+=(htmpAp[iz+1])*0.5;
    hApInter[iz]=tmpA1;
//    tmpA=(htmpAm[iz])*0.5; if(iz) tmpA=SplineAm.GetYIterC(dt1*(-0.5+iz));
    tmpA1=(htmpAm[iz])*0.5; if(iz) tmpA1+=(htmpAm[iz-1])*0.5;
    hAmInter[iz]=tmpA1;
  }
  //nextAm = (htmpAm[0])*0.5;
  for(int iz=0;iz<InitData.nSteps;iz++) {
//    tmpA=(htmpAm[iz])*0.5; if(iz) tmpA=SplineAm.GetYIterC(dt1*(-0.5+iz));
//    tmpA1=(htmpAm[iz])*0.5; if(iz) tmpA1+=(htmpAm[iz-1])*0.5;
    hAp[3*ASize+iz]=(hAmInter[iz]*ai+Fz(iz*dt1)*hI[iz])*dt1;
//    nextAm=SplineAm.GetYIterC(dt1*(0.5+iz));
/*
    if(iz<InitData.nSteps/2) {
      hApInter[2*iz]=tmpA;
      hApInter[2*iz+1]=htmpAm[iz];
      hAmInter[2*iz]=tmpA1;
      hAmInter[2*iz+1]=htmpAm[iz];
    }
*/
//    tmpA=(htmpAp[iz])*0.5; if(iz!=InitData.nSteps-1) tmpA=SplineAp.GetYIterC(dt1*(0.5+iz));
//    tmpA1=(htmpAp[iz])*0.5; if(iz!=InitData.nSteps-1) tmpA1+=(htmpAp[iz+1])*0.5;
    hAm[3*ASize+iz]=(hApInter[iz])*ai*dt1;

    htmpAp[iz] = hAp[iz]+hAp[3*ASize+iz];
    htmpAm[iz] = hAm[iz]+hAm[3*ASize+iz];
  }


  DoTasks();
  hOutA[nStepNumber]=hAm[InitData.nSteps-1]+onesix*( hAm[ASize+InitData.nSteps-1]*2.0+
                 hAm[2*ASize+InitData.nSteps-1]*4.0+hAm[3*ASize+InitData.nSteps-1]*2.0);
//  hOutA[nStepNumber]=hAp[0];

  for(int iz=0;iz<InitData.nSteps;iz++) {
    tmpA=0; if(iz) tmpA=htmpAm[iz-1];
    hAp[4*ASize+iz]=hAp[iz]+onesix*( hAp[ASize+iz]*2.0+hAp[2*ASize+iz]*4.0+hAp[3*ASize+iz]*2.0+
                        (tmpA*ai+Fz(iz*InitData.dt)*hI[iz])*dt1);
    if(iz<InitData.nSteps-1) hAm[4*ASize+iz+1]=hAm[iz]+onesix*( hAm[ASize+iz]*2.0+
                                hAm[2*ASize+iz]*4.0+hAm[3*ASize+iz]*2.0+
                          (htmpAp[iz+1]*ai)*dt1);
  }
  for(int iz=0;iz<InitData.nSteps;iz++) {
    hAm[iz]=hAm[4*ASize+iz];
    hAp[iz]=hAp[4*ASize+iz];
  }
  hAm[0]=0;
  EOut+=abs(hOutA[nStepNumber])*abs(hOutA[nStepNumber])*dt1;
  hAm[0]=0;
}

void cEqu2WaveElModel::DoTasks() {
  hTaskCritSection->Acquire();
  ThreadCounter=InitData.nThreads;
  hTaskEvent->ResetEvent();
  hTaskCritSection->Release();
  for(int i=0;i<InitData.nThreads;i++) {
    hThreads[i]->Resume();
  }
  hEqus[0]->Step();
  if(InitData.nThreads) {
    if (hTaskEvent->WaitFor(THREAD_WAIT_TIME) != wrSignaled) {
      Application->MessageBoxA("Thread wait time is off", "OMG",MB_OK);
    }
  }
}


/*
void cEquOptic2D_X_base::StepRoutineEulerNoMedia() {
  FillBorders(0);
  DoTasks();
  DoTasks();
}

void cEquOptic2D_X_base::StepRoutineRunge2NoMedia() {
  FillBorders(0);
  DoTasks();
  DoTasks();
  DoTasks();
  CalcOutEnergy();
  //CalcEnergy();
  //hEnergyData[nStepNumber].E=EnergyXP+EnergyXM+EnergyZP+EnergyZM;
  //SaveBitmaps();
  for(int i=0;i<FFTSizeFull-1;i++) hFFTSrc[i]=hFFTSrc[i+1];
  hFFTSrc[FFTSizeFull-1]=hAZp[InitData.SizeZ*InitData.SizeX/2+InitData.SizeX/2];
//  hFFTSrc[FFTSizeFull-1]=hMedia[(InitData.SizeZ-1)*InitData.SizeX+InitData.SizeX/2].AZm;
}

void cEquOptic2D_X_base::CalcEnergy() {
  int iI;
  EnergyXP=EnergyXM=EnergyZP=EnergyZM=EnergyR=EnergyQ=0;
  for(int ix=1;ix<=InitData.nStepsX;ix++)
    for(int iz=1;iz<=InitData.nStepsZ;iz++) {
      iI=IndexXZ(ix,iz);
      EnergyXP+=hAXp[iI].real()*hAXp[iI].real()+hAXp[iI].imag()*hAXp[iI].imag();
      EnergyXM+=hAXm[iI].real()*hAXm[iI].real()+hAXm[iI].imag()*hAXm[iI].imag();
      EnergyZP+=hAZp[iI].real()*hAZp[iI].real()+hAZp[iI].imag()*hAZp[iI].imag();
      EnergyZM+=hAZm[iI].real()*hAZm[iI].real()+hAZm[iI].imag()*hAZm[iI].imag();;
      EnergyR+=hM[iI].R0;
      EnergyQ+=hM[iI].R0Generated;
    }
  EnergyXP*=0.5*InitData.dt*InitData.dt;
  EnergyXM*=0.5*InitData.dt*InitData.dt;
  EnergyZP*=0.5*InitData.dt*InitData.dt;
  EnergyZM*=0.5*InitData.dt*InitData.dt;
  EnergyR*=InitData.dt*InitData.dt;
  EnergyQ*=InitData.dt*InitData.dt;
}

void cEquOptic2D_X_base::CalcOutEnergy() {
  int iI;
  OPTIC_TYPE tmpXP=0,tmpXM=0,tmpZP=0,tmpZM=0,tmpGenZ=0,tmpGenX=0;
  for(int ix=1;ix<=InitData.nStepsX;ix++) {
    tmpZP+= hAZp[IndexXZ(ix,InitData.nStepsZ)].real()
           *hAZp[IndexXZ(ix,InitData.nStepsZ)].real()
           +hAZp[IndexXZ(ix,InitData.nStepsZ)].imag()
           *hAZp[IndexXZ(ix,InitData.nStepsZ)].imag();
    tmpZM+=  hAZm[IndexXZ(ix,1)].real()
            *hAZm[IndexXZ(ix,1)].real()
            +hAZm[IndexXZ(ix,1)].imag()
            *hAZm[IndexXZ(ix,1)].imag();
  }
  for(int iz=1;iz<=InitData.nStepsZ;iz++) {
    tmpXP+= hAXp[IndexXZ(InitData.nStepsX,iz)].real()
           *hAXp[IndexXZ(InitData.nStepsX,iz)].real()
           +hAXp[IndexXZ(InitData.nStepsX,iz)].imag()
           *hAXp[IndexXZ(InitData.nStepsX,iz)].imag();
    tmpXM+= hAXm[IndexXZ(1,iz)].real()
           *hAXm[IndexXZ(1,iz)].real()
           +hAXm[IndexXZ(1,iz)].imag()
           *hAXm[IndexXZ(1,iz)].imag();
  }
  if((InitData.AGenZ)||(InitData.AGenX)) {
    for(int ix=1;ix<=InitData.nStepsX;ix++) {
      for(int iz=1;iz<=InitData.nStepsZ;iz++) {
        tmpGenZ+=hAZp[IndexXZ(ix,iz)].real()
                *hAZp[IndexXZ(ix,iz)].real()
                +hAZp[IndexXZ(ix,iz)].imag()
                *hAZp[IndexXZ(ix,iz)].imag()
                +hAZm[IndexXZ(ix,iz)].real()
                *hAZm[IndexXZ(ix,iz)].real()
                +hAZm[IndexXZ(ix,iz)].imag()
                *hAZm[IndexXZ(ix,iz)].imag();
        tmpGenX+=hAXp[IndexXZ(ix,iz)].real()
                *hAXp[IndexXZ(ix,iz)].real()
                +hAXp[IndexXZ(ix,iz)].imag()
                *hAXp[IndexXZ(ix,iz)].imag()
                +hAXm[IndexXZ(ix,iz)].real()
                *hAXm[IndexXZ(ix,iz)].real()
                +hAXm[IndexXZ(ix,iz)].imag()
                *hAXm[IndexXZ(ix,iz)].imag();
      }
    }
  }
  tmpXP*=0.5*InitData.dt*InitData.dt;  EnergyXPOut+=tmpXP;
  tmpXM*=0.5*InitData.dt*InitData.dt;  EnergyXMOut+=tmpXM;
  tmpZP*=0.5*InitData.dt*InitData.dt;  EnergyZPOut+=tmpZP;
  tmpZM*=0.5*InitData.dt*InitData.dt;  EnergyZMOut+=tmpZM;
  tmpGenZ*=-InitData.dt*InitData.dt*InitData.dt*InitData.AGenZ;
  tmpGenX*=-InitData.dt*InitData.dt*InitData.dt*InitData.AGenX;
  EnergyUPOut+=tmpGenZ+tmpGenX;
  hEnergyData[nStepNumber].OutXp=tmpXP;
  hEnergyData[nStepNumber].OutXm=tmpXM;
  hEnergyData[nStepNumber].OutZp=tmpZP;
  hEnergyData[nStepNumber].OutZm=tmpZM;
  hEnergyData[nStepNumber].OutGenZ=tmpGenZ;
  hEnergyData[nStepNumber].Out=tmpXP+tmpXM+tmpZP+tmpZM;
  if( (hEnergyData[nStepNumber].Out>InitData.MaxV[0])&&
      (LastMaxOut>hEnergyData[nStepNumber].Out)&&
      (!MaxESaved))
              SaveADataMaxE();
  LastMaxOut=hEnergyData[nStepNumber].Out;
}

void cEquOptic2D_X_base::SaveADataMaxE() {
  ofstream FileOutData = ofstream("DataMaxE.txt",ios::out);
  FileOutData << "t = " << t << endl;
  FileOutData << "SizeZ = " << InitData.dt*InitData.nStepsZ << endl;
  FileOutData << "SizeX = " << InitData.dt*InitData.nStepsX << endl;
  FileOutData << "dt = " << InitData.nStepsZ << endl;
  FileOutData << "a = " << InitData.ACoef << endl;
  FileOutData << "beta = " << InitData.beta << endl;
  FileOutData << "Q = " << InitData.Q << endl;
  FileOutData << "PCoef = " << InitData.PCoef << endl;
  FileOutData << "RCoef = " << InitData.RCoef << endl;
  ofstream FileOutAZp = ofstream("AZpMaxE.dat",ios::out);
  ofstream FileOutAXp = ofstream("AXpMaxE.dat",ios::out);
  ofstream FileOutAZm = ofstream("AZmMaxE.dat",ios::out);
  ofstream FileOutAXm = ofstream("AXmMaxE.dat",ios::out);

  ofstream FileOutAZpF = ofstream("AZpFMaxE.dat",ios::out);
  ofstream FileOutAXpF = ofstream("AXpFMaxE.dat",ios::out);
  ofstream FileOutAZmF = ofstream("AZmFMaxE.dat",ios::out);
  ofstream FileOutAXmF = ofstream("AXmFMaxE.dat",ios::out);
  for(int i=0;i<InitData.nStepsZ;i++)
    for(int j=0;j<InitData.nStepsX;j++) {
      FileOutAZm << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hAZm[Index_ZX(i+1,j+1)]) << endl;
      FileOutAXm << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hAXm[Index_ZX(i+1,j+1)]) << endl;
      FileOutAZp << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hAZp[Index_ZX(i+1,j+1)]) << endl;
      FileOutAXp << InitData.dt*i << " " << InitData.dt*j
                 << " " << abs(hAXp[Index_ZX(i+1,j+1)]) << endl;
    }
  for(int iz=1;iz<InitData.nStepsZ;iz++) {
    FileOutAXpF << InitData.dt*iz << " " << arg(hAXp[IndexXZ(InitData.nStepsX,iz)]) << endl;
    FileOutAXmF << InitData.dt*iz << " " <<arg(hAXm[IndexXZ(1,iz)]) << endl;
  }
  for(int ix=1;ix<InitData.nStepsZ;ix++) {
    FileOutAZpF << InitData.dt*ix << " " <<arg(hAZp[IndexXZ(ix,InitData.nStepsZ)]) << endl;
    FileOutAZmF << InitData.dt*ix << " " <<arg(hAZm[IndexXZ(ix,1)]) << endl;
  }
  MaxESaved=1;
}

*/
