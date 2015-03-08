
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
  // cout << r1 << endl;

  RadialCoordinate r2(prhs[1]);
  // cout << r2 << endl;

  AngleCoordinate theta(prhs[2]);
  cout << theta << endl;
  
  double *pot = mxGetPr(prhs[3]);
  insist(pot);

  double *psi = mxGetPr(prhs[4]);
  insist(psi);

  //FORT(psitest)(psi, r1.n, r2.n, theta.n, r1.dr, r2.dr, theta.x, theta.w);

  TimeEvolution time_evol(pot, (Complex *) psi, r1, r2, theta);
  
  cout << " module: " << time_evol.module() << endl;

  const int n1 = r1.n;
  const int n2 = r2.n;
  const int n3 = theta.n;
  
  const double n1n2 = sqrt(1.0*n1*n2);
  
  for(int i = 0; i<10; i++) {
    cout << i << endl;
    
    time_evol.forward_transform();
    
#pragma omp parallel for if(n1*n2*n3 > 100)	\
  default(shared) schedule(static, 1)                 
    for(int i1 = 0; i1<2*n1*n2*n3; i1++)
      psi[i1] /= n1n2;
    
    cout << " module: " << time_evol.module() << endl;
    
    time_evol.backward_transform();
    
#pragma omp parallel for if(n1*n2*n3 > 100)	\
  default(shared) schedule(static, 1)                 
    for(int i1 = 0; i1<2*n1*n2*n3; i1++)
      psi[i1] /= n1n2;
    
    cout << " module: " << time_evol.module() << endl;
  }

  std::cout.precision(np);
}
