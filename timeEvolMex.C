
/* $Id$ */

#include <iostream>
#include <cstring>
#include <cmath>
#include <mex.h>
#include "matutils.h"
#include "MatlabStructures.h"
#include "fort.h"

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
  // cout << r1 << endl;

  RadialCoordinate r2(prhs[1]);
  // cout << r2 << endl;

  AngleCoordinate theta(prhs[2]);
  // cout << theta << endl;
  
  double *pot = mxGetPr(prhs[3]);
  insist(pot);

  double *psi = mxGetPr(prhs[4]);
  insist(psi);

  FORT(psitest)(psi, r1.n, r2.n, theta.n, r1.dr, r2.dr, theta.x, theta.w);

  std::cout.precision(np);
}
