
/* $Id$ */

#include <iostream>
#include <cstring>
#include <cmath>
#include <mex.h>
#include "matutils.h"
#include "MatlabStructures.h"
#include "fort.h"
#include "timeEvol.h"

extern "C" int FORT(myisnan)(const double &x)
{
  return isnan(x);
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  const int np = std::cout.precision();
  std::cout.precision(15);
  
  std::cout << " 3D Time evolotion" << std::endl;

  insist(nrhs > 0);

  RadialCoordinate r1(prhs[0]);
  
  RadialCoordinate r2(prhs[1]);
  
  AngleCoordinate theta(prhs[2]);
  
  MatlabArray<double> pot(prhs[3]);

  MatlabArray<Complex> psi(prhs[4]);

  EvolutionTime time(prhs[5]);

  Options options(prhs[6]);

  DumpFunction dump1(prhs[7]);
  DumpFunction dump2(prhs[8]);

  TimeEvolution time_evol(pot, psi, r1, r2, theta, time, options, dump1, dump2);

  //cout << " module: " << time_evol.module_for_psi() << endl;
  
  const int n1 = r1.n;
  const int n2 = r2.n;
  const int n3 = theta.n;
  
  const double n1n2 = n1*n2;

  time_evol.time_evolution();

  std::cout.flush();
  std::cout.precision(np);
}
