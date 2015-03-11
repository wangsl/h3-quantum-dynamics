
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
		const AngleCoordinate &theta,
		EvolutionTime &time);
  
  ~TimeEvolution();

  double module_for_psi() const;
  double module_for_legendre_psi();
  
  void forward_fft_for_psi();
  void backward_fft_for_psi();
  
  void forward_fft_for_legendre_psi();
  void backward_fft_for_legendre_psi();

  void forward_legendre_transform();
  void backward_legendre_transform();

  void calculate_energy();
  
private:
  
  double *pot;
  Complex *psi;
  const RadialCoordinate &r1;
  const RadialCoordinate &r2;
  const AngleCoordinate &theta;
  EvolutionTime &time;

  Complex *_legendre_psi;
  Complex * &legendre_psi();

  Vec<FFTWInterface *> fftw_for_psi;
  Vec<FFTWInterface *> fftw_for_legendre_psi;

  void setup_fftw_interface_for_psi();
  void destroy_fftw_interface_for_psi();

  void setup_fftw_interface_for_legendre_psi();
  void destroy_fftw_interface_for_legendre_psi();

  double potential_energy();
  double rotational_energy();

  double kinetic_energy_for_psi();
  double kinetic_energy_for_legendre_psi();
};

#endif /* TIMEEVOL_H */
