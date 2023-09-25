//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquOpticTask_Wide1DBragg.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

void cEquOpticTask_Wide1DBragg::Step() {
}

void cEquOpticTask_Wide1DBragg::StepRunge4FastP() {
    complex<OPTIC_TYPE> iAlpha=complex<OPTIC_TYPE>(0,InitData.ACoef);
//  MakeAField
    for(int iz=0; iz<InitData.nStepsZ; iz++) {
      FFT.TransBackward(hApReal,hAp);  hApReal[InitData.nAHarm]=hApReal[0];
      FFT.TransBackward(hAmReal,hAm);  hAmReal[InitData.nAHarm]=hApReal[0];
//  Make K1 for Ap Am
      for(int ix=0; ix<InitData.nAHarm; ix++) {
        if(IsBragg(ix)) {
          htmpApReal[ix]=0.5*dtT*iAlpha*hAmReal[ix];
          htmpAmReal[ix]=0.5*dtT*iAlpha*hApReal[ix];
        } else {
          htmpApReal[ix]=htmpAmReal[ix]=0;
        }
      }
    }
//  --------------- need to add media effect on htmpApReal, htmpAmReal here
    FFT.TransForward(hAmReal,hAm);
}
