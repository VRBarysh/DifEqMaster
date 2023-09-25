//---------------------------------------------------------------------------


#pragma hdrstop

#include "Equ2WaveTask.h"
#include "Dump2WaveElModel.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

REAL_TYPE cEqu2WaveBaseTask::RightPart0(REAL_TYPE El1) {
  return El1;
}

REAL_TYPE cEqu2WaveBaseTask::RightPart1(complex<REAL_TYPE> Ap, REAL_TYPE El0) {
  return (Ap*polar(1,El0)).real();
}

void cEqu2WaveBaseTask::Step() { Step01(); }
void cEqu2WaveBaseTask::Step01() {
  ElSize = InitData.nEls*InitData.nSteps*2;
  REAL_TYPE dt2 = 0.5*InitData.dt, dt1=InitData.dt;
  REAL_TYPE onesix = 1.0/6.0;
  REAL_TYPE tmpE,EE1=0,EE2=0;
  //                             -----------         Euler part
  /*
  if(!RungeStep) {
    for(int iz=iz1st;iz<=izLast;iz++) {
      hI[iz]=0;
      for(int ie=0;ie<InitData.nEls;ie++) {
        hEl[ElSize+iE0(iz,ie)]= dt1*RightPart0(hEl[iE1(iz,ie)]);
        hEl[ElSize+iE1(iz,ie)]= dt1*RightPart1( hAp[iz],hEl[iE0(iz,ie)] );
        EE1+=hEl[iE1(iz,ie)];
        hI[iz]+=polar(1,-hEl[iE0(iz,ie)]);
        hEl[iE0(iz,ie)]+=hEl[ElSize+iE0(iz,ie)];
        hEl[iE1(iz,ie)]+=hEl[ElSize+iE1(iz,ie)];
        EE2+=hEl[iE1(iz,ie)];
      }
      hI[iz]/=InitData.nEls;
      EE1*=-2.0/(InitData.nEls*InitData.nSteps);
      EE2*=-2.0/(InitData.nEls*InitData.nSteps);
    }
  }
  RungeStep=(RungeStep+1)%4;
  return;
                                                                      */
  switch(RungeStep) {
    case 0:
      for(int iz=iz1st;iz<=izLast;iz++) {
        hI[iz]=0;
        for(int ie=0;ie<InitData.nEls;ie++) {
          hEl[ElSize+iE0(iz,ie)] = dt2*RightPart0( hEl[iE1(iz,ie)] );
          hEl[ElSize+iE1(iz,ie)] = dt2*RightPart1( hAp[iz],hEl[iE0(iz,ie)] );
          hI[iz]+=polar(1,-hEl[iE0(iz,ie)]);
        }
        hI[iz]/=InitData.nEls;
      }
      RungeStep++;
    break;
    case 1:
      for(int iz=iz1st;iz<=izLast;iz++) {
        hI[iz]=0;
        for(int ie=0;ie<InitData.nEls;ie++) {
          hEl[2*ElSize+iE0(iz,ie)]= dt2*RightPart0(
             (hEl[iE1(iz,ie)]+hEl[ElSize+iE1(iz,ie)]));
          hEl[2*ElSize+iE1(iz,ie)]=
             dt2*RightPart1( htmpAp[iz],hEl[iE0(iz,ie)]+hEl[ElSize+iE0(iz,ie)] );
          hI[iz]+=polar(1,-hEl[iE0(iz,ie)]-hEl[ElSize+iE0(iz,ie)]);
        }
        hI[iz]/=InitData.nEls;
      }
      RungeStep++;
    break;
    case 2:
      for(int iz=iz1st;iz<=izLast;iz++) {
        hI[iz]=0;
        for(int ie=0;ie<InitData.nEls;ie++) {
          hEl[3*ElSize+iE0(iz,ie)]= dt1*RightPart0(
             (hEl[iE1(iz,ie)]+hEl[2*ElSize+iE1(iz,ie)]));
          hEl[3*ElSize+iE1(iz,ie)]= dt1*
             RightPart1( htmpAp[iz],hEl[iE0(iz,ie)]+hEl[2*ElSize+iE0(iz,ie)] );
          hI[iz]+=polar(1,-hEl[iE0(iz,ie)]-hEl[2*ElSize+iE0(iz,ie)]);
        }
        hI[iz]/=InitData.nEls;
      }
      RungeStep++;
    break;
    case 3:
      for(int iz=iz1st;iz<=izLast;iz++) {
        hI[iz]=0;
        for(int ie=0;ie<InitData.nEls;ie++) {
          hI[iz]+=polar(1,-hEl[iE0(iz,ie)]-hEl[3*ElSize+iE0(iz,ie)]);

          tmpE = hEl[iE0(iz,ie)] +
             onesix*( hEl[ElSize+iE0(iz,ie)]*2.0+
                      hEl[2*ElSize+iE0(iz,ie)]*4.0+
                      hEl[3*ElSize+iE0(iz,ie)]*2.0+
                      dt1 * RightPart0(
              hEl[iE1(iz,ie)]+hEl[3*ElSize+iE1(iz,ie)]));

          hEl[iE1(iz,ie)]= hEl[iE1(iz,ie)]+
             onesix*( hEl[ElSize+iE1(iz,ie)]*2.0+
                      hEl[2*ElSize+iE1(iz,ie)]*4.0+
                      hEl[3*ElSize+iE1(iz,ie)]*2.0+
                      dt1 * RightPart1( htmpAp[iz],
                         hEl[iE0(iz,ie)]+hEl[3*ElSize+iE0(iz,ie)] ));
          hEl[iE0(iz,ie)]=tmpE;
        }
        hI[iz]/=InitData.nEls;
      }
      RungeStep=0;
    break;
  }
}
