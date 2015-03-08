
/* $Id$ */

#include "MatlabStructures.h"
#include "matutils.h"

RadialCoordinate::RadialCoordinate(const mxArray *mx) :
  n(*(int *) mxGetData(mx, "n")),
  dr(*(double *) mxGetData(mx, "dr")),
  mass(*(double *) mxGetData(mx, "mass"))
{ 
  r = RVec(n, (double *) mxGetData(mx, "r"));
}

AngleCoordinate::AngleCoordinate(const mxArray *mx) :
  n(*(int *) mxGetData(mx, "n")),
  m(*(int *) mxGetData(mx, "m"))
{
  x = RVec(n, (double *) mxGetData(mx, "x"));
  w = RVec(n, (double *) mxGetData(mx, "w"));
}
  
