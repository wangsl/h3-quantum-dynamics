
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
			     const AngleCoordinate &theta_) :
  pot(pot_), psi(psi_), r1(r1_), r2(r2_), theta(theta_)
{ }

TimeEvolution::~TimeEvolution()
{
  pot = 0;
  psi = 0;
  destroy_fftw_interface();
}

void TimeEvolution::destroy_fftw_interface()
{
  for(int i = 0; i < fftw.size(); i++) {
    if(fftw[i]) {
      delete fftw[i];
      fftw[i] = 0;
    }
  }
  fftw.resize(0);
}

void TimeEvolution::setup_fftw_interface()
{
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

void TimeEvolution::forward_fft_transform()
{
  setup_fftw_interface();
  for(int i = 0; i < fftw.size(); i++) 
    fftw[i]->forward_transform();
}

void TimeEvolution::backward_fft_transform()
{
  setup_fftw_interface();
  for(int i = 0; i < fftw.size(); i++) 
    fftw[i]->backward_transform();
}

void TimeEvolution::forward_legendre_transform()
{
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int &n_theta = theta.n;
  const int m = theta.m + 1;
  
  legendre_psi.resize(n1*n2, m);
  
  FORT(forwardlegendretransform)(psi, legendre_psi, n1*n2, n_theta, m, 
				 theta.w, theta.legendre);
}

void TimeEvolution::backward_legendre_transform()
{
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int &n_theta = theta.n;
  const int m = theta.m + 1;
  
  FORT(backwardlegendretransform)(psi, legendre_psi, n1*n2, n_theta, m, 
				  theta.w, theta.legendre);
}

double TimeEvolution::module() const
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
  e_pot = 0.0;
  e_kin = 0.0;
  e_rot = 0.0;
  
  calculate_potential_energy();
  calculate_kinetic_energy();
  calculate_rotational_energy();

  double e = e_kin + e_rot + e_pot;

  cout << " Total energy: " << e << endl;
}

void TimeEvolution::calculate_potential_energy()
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
  
  e_pot = s*r1.dr*r2.dr;
  cout << " e_pot: " << e_pot << endl;
}

void TimeEvolution::calculate_kinetic_energy()
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

  forward_fft_transform();
  
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
  
  backward_fft_transform();
  
  e_kin = e*dr1*dr2/(n1*n2);
  cout << " e_kin: " << e_kin << endl;
}

void TimeEvolution::calculate_rotational_energy()
{ 
  forward_legendre_transform();

  const double &dr1 = r1.dr;
  const double &dr2 = r2.dr;
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int m = theta.m + 1;
  
  const RVec &I1 = r1.one2mr2;
  const RVec &I2 = r2.one2mr2;
  
  Complex *p = legendre_psi;
  
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
  
  e_rot = s*dr1*dr2;
  cout << " e_rot: " << e_rot << endl;
} 

void TimeEvolution::legendre_transform_test() const
{ 
  const double &dr1 = r1.dr;
  const double &dr2 = r2.dr;
  
  const Mat<Complex> &LPsi = legendre_psi;
  
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
  cout << " Norm test after forward Legendre transform: " << s << endl;
}
