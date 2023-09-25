//---------------------------------------------------------------------------


#pragma hdrstop

#include <math.h>
#include "Gr3DPlane.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)


long double cGr3DMaster::GetX(long double x,long double y, long double z) {
  long double X=Init.PosX - Init.ScaleY*y*sin(Init.Angle1)+Init.ScaleX*x*cos(Init.Angle1);
  return(X);
}

long double cGr3DMaster::GetY(long double x,long double y, long double z) {
  long double Y=Init.PosY - Init.ScaleX*x*sin(Init.Angle1)*sin(Init.Angle2)
                          - Init.ScaleY*y*cos(Init.Angle1)*sin(Init.Angle2)
                          - Init.ScaleZ*z*cos(Init.Angle2);
  return(Y);
}

inline int FrontPlane(  long double x1, long double y1,
                        long double x2, long double y2,
                        long double x3, long double y3,
                        long double x4, long double y4) {
  return(1);
}
