//---------------------------------------------------------------------------

#ifndef Gr3DPlaneH
#define Gr3DPlaneH
//---------------------------------------------------------------------------

class cGr3DInitData
{
public:
long double PosX,PosY,ScaleX,ScaleY,ScaleZ,Angle1,Angle2;
protected:
};


class cGr3DMaster
{
public:
  cGr3DMaster() {}
  ~cGr3DMaster() {}
  cGr3DMaster(cGr3DInitData pInit) {Init=pInit;}
  long double GetX(long double x,long double y, long double z);
  long double GetY(long double x,long double y, long double z);
  inline int FrontPlane(long double x1, long double y1,
                        long double x2, long double y2,
                        long double x3, long double y3,
                        long double x4, long double y4);

cGr3DInitData Init;
protected:
};

#endif
