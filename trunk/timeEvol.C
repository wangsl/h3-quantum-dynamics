
/* $Id$ */

#include "timeEvol.h"
#include "mat.h"

extern "C" {
  void FORT(forwardlegendretransform)(const Complex *CPsi, Complex *LPsi, 
				      const int &NR, const int &NTheta, const int &NLeg, 
				      const double *W, const double *LegP);
  
  void FORT(backwardlegendretransform)(Complex *CPsi, const Complex *LPsi, 
				       const int &NR, const int &NTheta, const int &NLeg, 
				       const double *W, const double *LegP);
}

TimeEvolution::TimeEvolution(double *pot_, Complex *psi_,
			     const RadialCoordinate &r1_, const RadialCoordinate &r2_,
			     const AngleCoordinate &theta_,
			     EvolutionTime &time_) :
  pot(pot_), psi(psi_), r1(r1_), r2(r2_), theta(theta_), time(time_),
  _legendre_psi(0)
{ } 

TimeEvolution::~TimeEvolution()
{
  pot = 0;
  psi = 0;
  
  destroy_fftw_interface_for_psi();
  destroy_fftw_interface_for_legendre_psi();
  
  if(_legendre_psi) { delete [] _legendre_psi; _legendre_psi = 0; }
}

Complex * &TimeEvolution::legendre_psi()
{
  if(!_legendre_psi) {
    const int &n1 = r1.n;
    const int &n2 = r2.n;
    const int m = theta.m + 1;
    _legendre_psi = new Complex [n1*n2*m];
    insist(_legendre_psi);
  }

  return _legendre_psi;
}

void TimeEvolution::destroy_fftw_interface_for_psi()
{
  Vec<FFTWInterface *> &fftw = fftw_for_psi;
  for(int i = 0; i < fftw.size(); i++) {
    if(fftw[i]) {
      delete fftw[i];
      fftw[i] = 0;
    }
  }
  fftw.resize(0);
}

void TimeEvolution::destroy_fftw_interface_for_legendre_psi()
{
  Vec<FFTWInterface *> &fftw = fftw_for_legendre_psi;
  for(int i = 0; i < fftw.size(); i++) {
    if(fftw[i]) {
      delete fftw[i];
      fftw[i] = 0;
    }
  }
  fftw.resize(0);
}

void TimeEvolution::setup_fftw_interface_for_psi()
{
  Vec<FFTWInterface *> &fftw = fftw_for_psi;

  const int &n_theta = theta.n;
  if(fftw.size() == n_theta) return;
  
  fftw_init_threads();
  
  fftw.resize(n_theta);
  fftw.zero();
  
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  double *p = (double *) psi;
  for(int k = 0; k < n_theta; k++) {
    fftw[k] = new FFTWInterface(p, n1, n2, FFTW_MEASURE, 1);
    insist(fftw[k]);
    p += 2*n1*n2;
  }
}

void TimeEvolution::setup_fftw_interface_for_legendre_psi()
{
  Vec<FFTWInterface *> &fftw = fftw_for_legendre_psi;

  const int m = theta.m + 1;
  if(fftw.size() == m) return;
  
  fftw_init_threads();
  
  fftw.resize(m);
  fftw.zero();
  
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  double *p = (double *) legendre_psi();
  for(int k = 0; k < m; k++) {
    fftw[k] = new FFTWInterface(p+k*2*n1*n2, n1, n2, FFTW_MEASURE, 1);
    insist(fftw[k]);
  }
}

void TimeEvolution::forward_fft_for_psi()
{
  setup_fftw_interface_for_psi();
  for(int i = 0; i < fftw_for_psi.size(); i++) 
    fftw_for_psi[i]->forward_transform();
}

void TimeEvolution::backward_fft_for_psi()
{
  setup_fftw_interface_for_psi();
  for(int i = 0; i < fftw_for_psi.size(); i++) 
    fftw_for_psi[i]->backward_transform();
}

void TimeEvolution::forward_fft_for_legendre_psi()
{
  setup_fftw_interface_for_legendre_psi();
  for(int i = 0; i < fftw_for_legendre_psi.size(); i++) 
    fftw_for_legendre_psi[i]->forward_transform();
}

void TimeEvolution::backward_fft_for_legendre_psi()
{
  setup_fftw_interface_for_legendre_psi();
  for(int i = 0; i < fftw_for_legendre_psi.size(); i++) 
    fftw_for_legendre_psi[i]->backward_transform();
}

void TimeEvolution::forward_legendre_transform()
{
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int &n_theta = theta.n;
  const int m = theta.m + 1;
  
  FORT(forwardlegendretransform)(psi, legendre_psi(), n1*n2, n_theta, m, 
				 theta.w, theta.legendre);
}

void TimeEvolution::backward_legendre_transform()
{
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int &n_theta = theta.n;
  const int m = theta.m + 1;
  
  FORT(backwardlegendretransform)(psi, legendre_psi(), n1*n2, n_theta, m, 
				  theta.w, theta.legendre);
}

double TimeEvolution::module_for_psi() const
{
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int &n_theta = theta.n;
  const RVec &w = theta.w;
  
  Mat<Complex> p(n1*n2, n_theta, psi);
  
  double s = 0.0;
#pragma omp parallel for if(p.columns() > 100)	      \
  default(shared) schedule(static, 1)                 \
  reduction(+:s)
  for(int k = 0; k < p.columns(); k++) {
    double sk = 0.0;
    for(int i = 0; i < p.rows(); i++) {
      sk += abs2(p(i,k));
    }
    s += w[k]*sk;
  }
  
  s *= r1.dr*r2.dr;
  return s;
}

void TimeEvolution::calculate_energy()
{
  double e_kin = kinetic_energy_for_psi();
  double e_rot = rotational_energy();
  double e_pot = potential_energy();

  double e_total = e_kin + e_rot + e_pot;

  cout << " Total energy: " << e_total << endl;

  kinetic_energy_for_legendre_psi();
}

double TimeEvolution::potential_energy()
{
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const double &dr1 = r1.dr;
  const double &dr2 = r2.dr;
  const int &n_theta = theta.n;
  const RVec &w = theta.w;
  
  const Mat<Complex> Psi(n1*n2, n_theta, psi);
  const RMat Pot(n1*n2, n_theta, pot);
  
  double s = 0.0;
#pragma omp parallel for if(Psi.columns() > 100)      \
  default(shared) schedule(static, 1)                 \
  reduction(+:s)
  for(int k = 0; k < Psi.columns(); k++) {
    double sk = 0.0;
    for(int i = 0; i < Psi.rows(); i++) {
      sk += abs2(Psi(i,k))*Pot(i,k);
    }
    s += w[k]*sk;
  }
  
  double e_pot = s*r1.dr*r2.dr;
  cout << " e_pot: " << e_pot << endl;

  return e_pot;
}

double TimeEvolution::kinetic_energy_for_psi()
{ 
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int n_theta = theta.n;
  const double &dr1 = r1.dr;
  const double &dr2 = r2.dr;
  const RVec &w = theta.w;
  
  const RVec &kin1 = r1.psq2m;
  const RVec &kin2 = r2.psq2m;
  const double n1n2 = n1*n2;

  forward_fft_for_psi();
  
  double e = 0.0;
#pragma omp parallel for			\
  default(shared) schedule(static, 1)		\
  reduction(+:e)
  for(int k = 0; k < n_theta; k++) {
    Mat<Complex> Psi(n1, n2, psi+k*n1*n2);
    double ek = 0.0;
    for(int j = 0; j < n2; j++) {
      for(int i = 0; i < n1; i++) {
	ek += abs2(Psi(i,j))*(kin1[i] + kin2[j]);
	Psi(i,j) /= n1n2;
      }
    }
    e += w[k]*ek;
  }
  
  backward_fft_for_psi();
  
  double e_kin = e*dr1*dr2/(n1*n2);

  cout << " e_kin: " << e_kin << endl;

  return e_kin;
}

double TimeEvolution::rotational_energy()
{ 
  forward_legendre_transform();

  const double &dr1 = r1.dr;
  const double &dr2 = r2.dr;
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int m = theta.m + 1;
  
  const RVec &I1 = r1.one2mr2;
  const RVec &I2 = r2.one2mr2;
  
  Complex *p = legendre_psi();
  
  double s = 0.0;
#pragma omp parallel for			\
  default(shared) schedule(static, 1)		\
  reduction(+:s)
  for(int l = 0; l < m; l++) {
    const Mat<Complex> LPsi(n1, n2, p+l*n1*n2);
    double sl = 0.0;
    for(int j = 0; j < n2; j++) {
      for(int i = 0; i < n1; i++) {
	sl += abs2(LPsi(i,j))*(I1[i]+I2[j]);
      }
    }
    s += sl*l*(l+1)/(l+0.5);
  }
  
  double e_rot = s*dr1*dr2;

  cout << " e_rot: " << e_rot << endl;

  return e_rot;
} 

double TimeEvolution::module_for_legendre_psi()
{ 
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int m = theta.m + 1;
  
  const double &dr1 = r1.dr;
  const double &dr2 = r2.dr;
  
  Mat<Complex> LPsi(n1*n2, m, legendre_psi());
  
  double s = 0.0;
#pragma omp parallel for			      \
  default(shared) schedule(static, 1)                 \
  reduction(+:s)
  for(int l = 0; l < LPsi.columns(); l++) {
    double sl = 0.0;
    for(int i = 0; i < LPsi.rows(); i++) {
      sl += abs2(LPsi(i,l));
    }
    s += sl/(l+0.5);
  }
  
  s *= dr1*dr2;

  cout << " module from Legendre Psi: " << s << endl;
  
  return s;
}

double TimeEvolution::kinetic_energy_for_legendre_psi()
{ 
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const double &dr1 = r1.dr;
  const double &dr2 = r2.dr;
  
  const RVec &kin1 = r1.psq2m;
  const RVec &kin2 = r2.psq2m;

  const int m = theta.m + 1;

  const double n1n2 = n1*n2;
  
  forward_fft_for_legendre_psi();

  Complex *p = legendre_psi();

  double s = 0.0;
#pragma omp parallel for			\
  default(shared) schedule(static, 1)		\
  reduction(+:s)
  for(int l = 0; l < m; l++) {
    Mat<Complex> Psi(n1, n2, p+l*n1*n2);
    double sl = 0.0;
    for(int j = 0; j < n2; j++) {
      for(int i = 0; i < n1; i++) {
	sl += abs2(Psi(i,j))*(kin1[i] + kin2[j]);
	Psi(i,j) /= n1n2;
      }
    }
    s += sl/(l+0.5);
  }
  
  backward_fft_for_legendre_psi();
  
  double e_kin = s*dr1*dr2/(n1*n2);

  cout << " e_kin from Legendre Psi: " << e_kin << endl;

  return e_kin;
}

