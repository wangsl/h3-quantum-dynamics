
/* $Id$ */

#include "MatlabStructures.h"
#include "matutils.h"
#include "fftwinterface.h"

RadialCoordinate::RadialCoordinate(const mxArray *mx) :
  n(*(int *) mxGetData(mx, "n")),
  dr(*(double *) mxGetData(mx, "dr")),
  mass(*(double *) mxGetData(mx, "mass"))
{ 
  r = RVec(n, (double *) mxGetData(mx, "r"));
  
  psq2m.resize(n);
  FFTWInterface::get_momentum_for_fftw(psq2m, n*dr);
  for(int i = 0; i < n; i++) {
    psq2m[i] =  psq2m[i]*psq2m[i]/(2*mass);
  }

  one2mr2.resize(n);
  const double m2 = mass+mass;
  for(int i = 0; i < n; i++) { 
    one2mr2[i] = 1.0/(m2*r[i]*r[i]);
  }
}

AngleCoordinate::AngleCoordinate(const mxArray *mx) :
  n(*(int *) mxGetData(mx, "n")),
  m(*(int *) mxGetData(mx, "m"))
{
  x = RVec(n, (double *) mxGetData(mx, "x"));
  w = RVec(n, (double *) mxGetData(mx, "w"));
  
  double *p = (double *) mxGetData(mx, "legendre");
  insist(p);

  legendre = RMat(n, m, p);
}
  
EvolutionTime::EvolutionTime(const mxArray *mx) :
  total_steps(*(int *) mxGetData(mx, "total_steps")),
  steps(*(int *) mxGetData(mx, "steps")),
  time_step(*(double *) mxGetData(mx, "time_step"))
{ }

