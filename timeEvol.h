
/* $Id$ */

#ifndef TIMEEVOL_H
#define TIMEEVOL_H

#include <iostream>
using namespace std;
#include "complex.h"
#include "MatlabStructures.h"

class TimeEvolution
{
public:
  double *pot;
  Complex *psi;

  //RadialCoordinate r1;
  //RadialCoordinate r2;
  //AngleCoordinate theta;

  TimeEvolution(double *pot, Complex *psi);

  ~TimeEvolution();
  
private:


};

#endif /* TIMEEVOL_H */
