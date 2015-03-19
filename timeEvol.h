
/* $Id$ */

#ifndef TIMEEVOL_H
#define TIMEEVOL_H

#include <iostream>
using namespace std;
#include "complex.h"
#include "MatlabStructures.h"
#include "fftwinterface.h"
#include "matlabArray.h"

class TimeEvolution
{
public:
  TimeEvolution(const MatlabArray<double> &pot,
		MatlabArray<Complex> &psi,
		const RadialCoordinate &r1,
		const RadialCoordinate &r2,
		const AngleCoordinate &theta,
		EvolutionTime &time,
		const Options &options);
  
  ~TimeEvolution();

  double module_for_psi() const;
  double module_for_legendre_psi();
  
  void calculate_energy();

  void test();

  void time_evolution();
  
private:
  
  double *pot;
  Complex *psi;

  const MatlabArray<double> &m_pot;
  MatlabArray<Complex> &m_psi;

  const RadialCoordinate &r1;
  const RadialCoordinate &r2;
  const AngleCoordinate &theta;
  EvolutionTime &time;
  const Options &options;

  Complex *_legendre_psi;
  Complex * &legendre_psi();

  Complex *exp_ipot_dt;
  Complex *exp_irot_dt_2;
  Complex *exp_ikin_dt;

  Vec<FFTWInterface *> fftw_for_psi;
  Vec<FFTWInterface *> fftw_for_legendre_psi;

  void setup_fftw_interface_for_psi();
  void destroy_fftw_interface_for_psi();

  void setup_fftw_interface_for_legendre_psi();
  void destroy_fftw_interface_for_legendre_psi();
  
  void forward_fft_for_psi();
  void backward_fft_for_psi();
  
  void forward_fft_for_legendre_psi();
  void backward_fft_for_legendre_psi();
  
  void forward_legendre_transform();
  void backward_legendre_transform();
  
  double potential_energy();
  double rotational_energy(const int do_legendre_transform = 1);
  
  double kinetic_energy_for_psi();
  double kinetic_energy_for_legendre_psi(const int do_fft = 1);

  void setup_exp_ipot_dt();
  void setup_exp_irot_dt_2();
  void setup_exp_ikin_dt();

  void pre_evolution_with_potential_dt_2();

  void evolution_with_potential_dt();
  void evolution_with_rotational_dt_2();
  void evolution_with_kinetic_dt();

  void evolution_dt();

};

#endif /* TIMEEVOL_H */
