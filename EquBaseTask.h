//---------------------------------------------------------------------------

#ifndef EquBaseTaskH
#define EquBaseTaskH
//---------------------------------------------------------------------------

class cEquBaseTask {
public:
  virtual void Step()=0;
  virtual int Finished()=0;
protected:
};

#endif
