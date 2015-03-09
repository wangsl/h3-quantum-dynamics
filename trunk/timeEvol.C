
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
  destroy_fftw();
}

void TimeEvolution::destroy_fftw()
{
  for(int i = 0; i < fftw.size(); i++) {
    if(fftw[i]) {
      delete fftw[i];
      fftw[i] = 0;
    }
  }
  fftw.resize(0);
}

void TimeEvolution::setup_fftw()
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
  setup_fftw();
  for(int i = 0; i < fftw.size(); i++) 
    fftw[i]->forward_transform();
}

void TimeEvolution::backward_fft_transform()
{
  setup_fftw();
  for(int i = 0; i < fftw.size(); i++) 
    fftw[i]->backward_transform();
}

void TimeEvolution::forward_legendre_transform()
{
  cout << " TimeEvolution::forward_legendre_transform" << endl;
  
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
  cout << " TimeEvolution::backward_legendre_transform" << endl;

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
