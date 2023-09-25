//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquOpticTask_BMB.h"
#include <fstream.h>

//---------------------------------------------------------------------------

#pragma package(smart_init)

void cEquOpticBaseTask_BMB::Step() {
  int RungeStepBack=RungeStep;
  OPTIC_TYPE backACoef, backAPhase;
  hAXp=hAXp1; hAXm=hAXm1; hAZp=hAZp1; hAZm=hAZm1;
  StepNoMediaMatrix();
  if(InitData.ACoef2) {
    backACoef=InitData.ACoef;
    backAPhase=InitData.APhase;
    InitData.ACoef=InitData.ACoef2;
    InitData.APhase=0;
    FillMatrixes();
  }
  RungeStep=RungeStepBack;
  hAXp=hAXp2; hAXm=hAXm2; hAZp=hAZp2; hAZm=hAZm2;
  StepNoMediaMatrix();
  if(InitData.ACoef2) {
    InitData.ACoef=backACoef;
    InitData.APhase=backAPhase;
    FillMatrixes();
  }
  if(!InitData.fNoMedia) {
    if(InitData.UseFastP) {
      StepMediaRunge2FastP();
    } else {
      StepMediaRunge2();
    }
  }
}

