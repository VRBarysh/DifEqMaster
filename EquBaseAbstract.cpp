//---------------------------------------------------------------------------


#pragma hdrstop

#include "EquBaseAbstract.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

void cDifEquReport::LoadText(cDifEquReportText *hReportText) {
  nGraphs=hReportText->nGraphs;
  for(int i=0;i<nGraphs;i++) GraphText[i]=hReportText->GraphText[i];
  nVals=hReportText->nVals;
  for(int i=0;i<nVals;i++) ValText[i]=hReportText->ValText[i];
}

void cDifEquReport::SetGraphSize(int GraphSize,int *hLines) {
  for(int i=0;i<EQU_MAX_GRAPHS;i++)
    if(hLines[i]) Graph[i].Init(hLines[i],GraphSize);
}


