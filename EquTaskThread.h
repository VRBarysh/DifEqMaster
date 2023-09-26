//---------------------------------------------------------------------------

#ifndef EquTaskThreadH
#define EquTaskThreadH

#include "EquBaseTask.h"
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <vcl\syncobjs.hpp>
//---------------------------------------------------------------------------

#define MAX_THREADS      20
#define THREAD_WAIT_TIME 100000


class TTaskThread : public TThread
{            
private:
  cEquBaseTask        *hEquation;
  TEvent              *hEvent;
  TCriticalSection    *hCritSection;
  int                 *hCounter;
protected:
        void __fastcall Execute();
public:
        __fastcall TTaskThread(bool CreateSuspended);
  void Init(cEquBaseTask *hEq, TEvent *hEv, TCriticalSection *hCr, int *hCounter);
};
//---------------------------------------------------------------------------
#endif
