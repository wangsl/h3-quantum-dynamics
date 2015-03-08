
/* $Id$ */

#ifndef TIMEEVOL_H
#define TIMEEVOL_H

#include <iostream>
using namespace std;
#include "complex.h"
#include "MatlabStructures.h"
#include "fftwinterface.h"

class TimeEvolution
{
public:
  TimeEvolution(double *pot, Complex *psi, 
		const RadialCoordinate &r1, const RadialCoordinate &r2,
		const AngleCoordinate &theta);
  
  ~TimeEvolution();

  double module() const;

  void forward_transform();
  void backward_transform();

private:
  
  double *pot;
  Complex *psi;
  const RadialCoordinate &r1;
  const RadialCoordinate &r2;
  const AngleCoordinate &theta;

  Vec<FFTWInterface *> fftw;

  void setup_fftw();
  void destroy_fftw();
};

#endif /* TIMEEVOL_H */
