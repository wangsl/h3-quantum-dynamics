
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

extern "C" void FORT(psitest)(const double *psi, 
			      const int &n1, const int &n2, const int &nTheta,
			      const double &dr1, const double &dr2, 
			      const double *x, const double *w);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  const int np = std::cout.precision();
  std::cout.precision(15);
  
  std::cout << " 3D Time evolotion" << std::endl;

  insist(nrhs>0);

  RadialCoordinate r1(prhs[0]);
  
  RadialCoordinate r2(prhs[1]);
  
  AngleCoordinate theta(prhs[2]);
  
  double *pot = mxGetPr(prhs[3]);
  insist(pot);
  
  double *psi = mxGetPr(prhs[4]);
  insist(psi);
  
  EvolutionTime time(prhs[5]);
  
  TimeEvolution time_evol(pot, (Complex *) psi, r1, r2, theta, time);
  
  cout << " module: " << time_evol.module_for_psi() << endl;
  
  const int n1 = r1.n;
  const int n2 = r2.n;
  const int n3 = theta.n;
  
  const double n1n2 = n1*n2;
  
  for(int i = 0; i < 10; i++) {
    
    cout << i << endl;

    time_evol.calculate_energy();

    cout << time_evol.module_for_legendre_psi() << endl;
  }

  std::cout.precision(np);
}
