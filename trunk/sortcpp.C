
/* $Id$ */

#include <iostream>
#include <algorithm>
#include <cassert>
#include "fort.h"

typedef std::pair<double, double> PairDouble;

void sort2_double(double *x, double *y, const int n)
{
  PairDouble *xy = new PairDouble [n];
  assert(xy);
  
  for(int i = 0; i < n; i++) {
    xy[i].first = x[i];
    xy[i].second = y[i];
  }
  
  std::sort(xy, xy+n,
	    [](const PairDouble &xy1, const PairDouble &xy2)
	    { return (xy1.first < xy2.first); }
	    );
  
  for(int i = 0; i < n; i++) {
    x[i] = xy[i].first;
    y[i] = xy[i].second;
  }

  if(xy) { delete [] xy; xy = 0; }
}

// Fortran version: Sort2Double
extern "C" void FORT(sort2double)(double *x, double *y, const int &n)
{
  sort2_double(x, y, n);
}


#if 0
int main()
{
  double A[] = {9,6,3,5,7.9};
  double B[] = {1,2,3,4,5.2};

  std::cout << "Before:" << std::endl;
  
  const int n = sizeof(A)/sizeof(double);

  for(int i = 0; i < n; i++) 
    std::cout << " " << i << " " << A[i] << " " << B[i] << std::endl;
  
  sort2_double(A, B, n);

  std::cout << std::endl;
  for(int i = 0; i < n; i++) 
    std::cout << " " << i << " " << A[i] << " " << B[i] << std::endl;
}
#endif



