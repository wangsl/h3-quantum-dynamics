
/* $Id$ */

#include "timeEvol.h"
#include "mat.h"

//#define _PrintFunction_ { cout << " " << __func__ << endl; }

#define _PrintFunction_ { }

extern "C" {
#if 0
  void FORT(forwardlegendretransform)(const Complex *CPsi, Complex *LPsi, 
				      const int &NR, const int &NTheta, const int &NLeg, 
				      const double *W, const double *LegP);
#endif
  
  void FORT(forwardlegendretransform)(const Complex *CPsi, Complex *LPsi, 
				      const int &NR, const int &NTheta, const int &NLeg, 
				      const double *WLegP);
  
#if 0
  void FORT(backwardlegendretransform)(Complex *CPsi, const Complex *LPsi, 
				       const int &NR, const int &NTheta, const int &NLeg, 
				       const double *W, const double *LegP);
#endif
  
  void FORT(backwardlegendretransform)(Complex *CPsi, const Complex *LPsi, 
				       const int &NR, const int &NTheta, const int &NLeg, 
				       const double *LegP);
}

TimeEvolution::TimeEvolution(const MatlabArray<double> &m_pot_,
			     MatlabArray<Complex> &m_psi_,
			     const RadialCoordinate &r1_, 
			     const RadialCoordinate &r2_,
			     const AngleCoordinate &theta_,
			     EvolutionTime &time_,
			     const Options &options_,
			     const DumpFunction &dump1_,
			     const DumpFunction &dump2_) :
  m_pot(m_pot_), 
  m_psi(m_psi_), 
  r1(r1_),
  r2(r2_), 
  theta(theta_), 
  time(time_),
  options(options_),
  dump1(dump1_),
  dump2(dump2_),
  _legendre_psi(0), 
  exp_ipot_dt(0), 
  exp_irot_dt_2(0), 
  exp_ikin_dt(0),
  weight_legendre(0),
  dump(0)
{ 
  pot = m_pot.data;
  insist(pot);

  psi = m_psi.data;
  insist(psi);

  _PrintFunction_;
} 

TimeEvolution::~TimeEvolution()
{
  _PrintFunction_;
  
  pot = 0;
  psi = 0;
  
  destroy_fftw_interface_for_psi();
  destroy_fftw_interface_for_legendre_psi();
  
  if(_legendre_psi) { delete [] _legendre_psi; _legendre_psi = 0; }
  if(exp_ipot_dt) { delete [] exp_ipot_dt; exp_ipot_dt = 0; }
  if(exp_irot_dt_2) { delete [] exp_irot_dt_2; exp_irot_dt_2 = 0; }
  if(exp_ikin_dt) { delete [] exp_ikin_dt; exp_ikin_dt = 0; }
  if(dump) { delete [] dump; dump = 0; }
  if(weight_legendre) { delete [] weight_legendre; weight_legendre = 0; }
}

Complex * &TimeEvolution::legendre_psi()
{
  if(!_legendre_psi) {

    _PrintFunction_;

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
  _PrintFunction_;

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
  _PrintFunction_;

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

  _PrintFunction_;
  
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

  _PrintFunction_;
  
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
  _PrintFunction_;

  setup_fftw_interface_for_psi();
  for(int i = 0; i < fftw_for_psi.size(); i++) 
    fftw_for_psi[i]->forward_transform();
}

void TimeEvolution::backward_fft_for_psi()
{
  _PrintFunction_;

  setup_fftw_interface_for_psi();
  for(int i = 0; i < fftw_for_psi.size(); i++) 
    fftw_for_psi[i]->backward_transform();
}

void TimeEvolution::forward_fft_for_legendre_psi()
{
  _PrintFunction_;

  setup_fftw_interface_for_legendre_psi();
  for(int i = 0; i < fftw_for_legendre_psi.size(); i++) 
    fftw_for_legendre_psi[i]->forward_transform();
}

void TimeEvolution::backward_fft_for_legendre_psi()
{
  _PrintFunction_;
  
  setup_fftw_interface_for_legendre_psi();
  for(int i = 0; i < fftw_for_legendre_psi.size(); i++) 
    fftw_for_legendre_psi[i]->backward_transform();
}

void TimeEvolution::forward_legendre_transform()
{
  _PrintFunction_;
  
  setup_weight_legendre();

  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int &n_theta = theta.n;
  const int m = theta.m + 1;
  
  FORT(forwardlegendretransform)(psi, legendre_psi(), n1*n2, n_theta, m, 
				 weight_legendre);
}

void TimeEvolution::backward_legendre_transform()
{
  _PrintFunction_;
  
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int &n_theta = theta.n;
  const int m = theta.m + 1;
  
  FORT(backwardlegendretransform)(psi, legendre_psi(), n1*n2, n_theta, m, 
				  theta.legendre);
}

double TimeEvolution::module_for_psi() const
{
  _PrintFunction_;

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
  _PrintFunction_;

  double e_kin = kinetic_energy_for_psi();
  double e_rot = rotational_energy();
  double e_pot = potential_energy();

  double e_total = e_kin + e_rot + e_pot;

  cout << " Total energy: " << e_total << endl;
}

double TimeEvolution::potential_energy()
{
  _PrintFunction_;

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

  return e_pot;
}

double TimeEvolution::kinetic_energy_for_psi()
{ 
  _PrintFunction_;

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

  return e_kin;
}

double TimeEvolution::rotational_energy(const int do_legendre_transform)
{ 
  _PrintFunction_;

  if(do_legendre_transform)
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

  if(do_legendre_transform)
    backward_legendre_transform();
  
  double e_rot = s*dr1*dr2;

  return e_rot;
} 

double TimeEvolution::module_for_legendre_psi()
{ 
  _PrintFunction_;

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

  return s;
}

double TimeEvolution::kinetic_energy_for_legendre_psi(const int do_fft)
{ 
  _PrintFunction_;

  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const double &dr1 = r1.dr;
  const double &dr2 = r2.dr;
  
  const RVec &kin1 = r1.psq2m;
  const RVec &kin2 = r2.psq2m;

  const int m = theta.m + 1;

  const double n1n2 = n1*n2;

  if(do_fft)
    forward_fft_for_legendre_psi();
  
  Complex *p = legendre_psi();
  
  double s = 0.0;
#pragma omp parallel for			\
  default(shared) schedule(static, 1)		\
  reduction(+:s)
  for(int l = 0; l < m; l++) {
    Mat<Complex> LPsi(n1, n2, p+l*n1*n2);
    double sl = 0.0;
    for(int j = 0; j < n2; j++) {
      for(int i = 0; i < n1; i++) {
	sl += abs2(LPsi(i,j))*(kin1[i] + kin2[j]);
	if(do_fft)
	  LPsi(i,j) /= n1n2;
      }
    }
    s += sl/(l+0.5);
  }

  double e_kin = s*dr1*dr2;
  
  if(do_fft) {
    backward_fft_for_legendre_psi();
    e_kin /= n1*n2;
  } else {
    e_kin *= n1*n2;
  }
  
  return e_kin;
}

void TimeEvolution::setup_exp_ipot_dt()
{
  if(exp_ipot_dt) return;

  _PrintFunction_;
  
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int &n_theta = theta.n;
  
  exp_ipot_dt = new Complex [n1*n2*n_theta];
  insist(exp_ipot_dt);
  
  const double &dt = time.time_step;
  const Complex Idt(0.0, -dt);
  
#pragma omp parallel for			\
  default(shared) schedule(static, 1)		
  for(int k = 0; k < n_theta; k++) {
    RMat v(n1, n2, pot+k*n1*n2);
    Mat<Complex> p(n1, n2, exp_ipot_dt+k*n1*n2);
    for(int j = 0; j < n2; j++) {
      for(int i = 0; i < n1; i++) {
	p(i,j) = exp(Idt*v(i,j));
      }
    }
  }
}

void TimeEvolution::setup_exp_irot_dt_2()
{ 
  if(exp_irot_dt_2) return;

  _PrintFunction_;
  
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int m = theta.m + 1;
  
  exp_irot_dt_2 = new Complex [n1*n2*m];
  insist(exp_irot_dt_2);

  const double &dt = time.time_step;
  const Complex Idt2(0.0, -dt/2);

#pragma omp parallel for			\
  default(shared) schedule(static, 1)	
  for(int l = 0; l < m; l++) {
    Mat<Complex> p(n1, n2, exp_irot_dt_2+l*n1*n2);
    for(int j = 0; j < n2; j++) {
      const double &I2 = r2.one2mr2[j];
      for(int i = 0; i < n1; i++) {
	const double &I1 = r1.one2mr2[i];
	p(i,j) = exp(l*(l+1)*(I1+I2)*Idt2);
      }
    }
  }
}

void TimeEvolution::setup_exp_ikin_dt()
{ 
  if(exp_ikin_dt) return;

  const int &n1 = r1.n;
  const int &n2 = r2.n;
  
  exp_ikin_dt = new Complex [n1*n2];
  insist(exp_ikin_dt);
  
  const double &dt = time.time_step;
  const Complex Idt(0.0, -dt);
  
  Mat<Complex> p(n1, n2, exp_ikin_dt);
  
#pragma omp parallel for			\
  default(shared) schedule(static, 1)	
  for(int j = 0; j < n2; j++) {
    const double &T2 = r2.psq2m[j];
    for(int i = 0; i < n1; i++) {
      const double &T1 = r1.psq2m[i];
      p(i,j) = exp((T1+T2)*Idt)/n1/n2;
    }
  }
}

void TimeEvolution::pre_evolution_with_potential_dt_2()
{ 
  _PrintFunction_;
  
  setup_exp_ipot_dt();
  
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int &n_theta = theta.n;
  
#pragma omp parallel for			\
  default(shared) schedule(static, 1)	    
  for(int k = 0; k < n_theta; k++) {
    const Mat<Complex> v(n1, n2, exp_ipot_dt+k*n1*n2);
    Mat<Complex> Psi(n1, n2, psi+k*n1*n2);
    for(int j = 0; j < n2; j++) {
      for(int i = 0; i < n1; i++) {
	Psi(i,j) /= sqrt(v(i,j));
      }
    }
  }
}

void TimeEvolution::evolution_with_potential_dt()
{ 
  _PrintFunction_;

  setup_exp_ipot_dt();
  
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int &n_theta = theta.n;
  
#pragma omp parallel for			\
  default(shared) schedule(static, 1)	    
  for(int k = 0; k < n_theta; k++) {
    const Mat<Complex> v(n1, n2, exp_ipot_dt+k*n1*n2);
    Mat<Complex> Psi(n1, n2, psi+k*n1*n2);
    for(int j = 0; j < n2; j++) {
      for(int i = 0; i < n1; i++) {
	Psi(i,j) *= v(i,j);
      }
    }
  }
}

void TimeEvolution::evolution_with_rotational_dt_2()
{
  _PrintFunction_;

  setup_exp_irot_dt_2();
  
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int m = theta.m + 1;
  
  Complex *lpsi = legendre_psi();
  
#pragma omp parallel for			\
  default(shared) schedule(static, 1)	
  for(int l = 0; l < m; l++) {
    Mat<Complex> rot(n1, n2, exp_irot_dt_2+l*n1*n2);
    Mat<Complex> LPsi(n1, n2, lpsi+l*n1*n2);
    for(int j = 0; j < n2; j++) {
      for(int i = 0; i < n1; i++) {
	LPsi(i,j) *= rot(i,j);
      }
    }
  }
}

void TimeEvolution::evolution_with_kinetic_dt()
{ 
  setup_exp_ikin_dt();

  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int m = theta.m + 1;
  
  Mat<Complex> kin(n1, n2, exp_ikin_dt);

  Complex *lpsi = legendre_psi();

#pragma omp parallel for			\
  default(shared) schedule(static, 1)	
  for(int l = 0; l < m; l++) {	       
    Mat<Complex> LPsi(n1, n2, lpsi+l*n1*n2);
    for(int j = 0; j < n2; j++) {
      for(int i = 0; i < n1; i++) {
	LPsi(i,j) *= kin(i,j);
      }
    }
  }
}

void TimeEvolution::test()
{
  _PrintFunction_;

  kinetic_energy_for_psi();

  forward_legendre_transform();
  kinetic_energy_for_legendre_psi();

  return;

  setup_exp_ipot_dt();
  setup_exp_irot_dt_2();
  setup_exp_ikin_dt();
}

void TimeEvolution::evolution_dt()
{
  _PrintFunction_;

  evolution_with_potential_dt();

  forward_legendre_transform();

  evolution_with_rotational_dt_2();

  forward_fft_for_legendre_psi();

  evolution_with_kinetic_dt();

  const double e_kin = kinetic_energy_for_legendre_psi(0);

  backward_fft_for_legendre_psi();

  evolution_with_rotational_dt_2();

  const double e_rot = rotational_energy(0);

  backward_legendre_transform();

  const double e_pot = potential_energy();

  cout << " e_kin: " << e_kin << "\n"
       << " e_rot: " << e_rot << "\n"
       << " e_pot: " << e_pot << "\n"
       << " e_tot: " << e_kin + e_rot + e_pot << endl;
}

void TimeEvolution::time_evolution()
{
  const int &total_steps = time.total_steps;
  int &steps = time.steps;

  for(int i_step = 0; i_step < total_steps; i_step++) {
    
    cout << "\n Step: " << i_step << endl;

    cout << " Module: " << module_for_psi() << endl;
    
    if(i_step == 0 && steps == 0)
      pre_evolution_with_potential_dt_2();
    
    evolution_dt();

    dump_psi();

    if(options.wave_to_matlab) {
      mxArray *mx[] = { (mxArray *) r1.mx, (mxArray *) r2.mx, (mxArray *) theta.mx,
			(mxArray *) m_pot.mx, (mxArray *) m_psi.mx, (mxArray *) time.mx,
			(mxArray *) options.mx 
      };
      
      const int n = sizeof(mx)/sizeof(mxArray *);
      wavepacket_to_matlab(options.wave_to_matlab, n, mx);
    }
    
    steps++;
    cout.flush();
  }
}

void TimeEvolution::setup_dump()
{
  if(dump) return;
  if(!apply_dump()) return;
  
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  
  if(!dump) { 
    dump = new double [n1*n2];
    assert(dump);
  }
  
  RMat d(n1, n2, dump);
  const double *d1 = dump1.dump;
  const double *d2 = dump2.dump;
  
#pragma omp parallel for			\
  default(shared) schedule(static, 1)
  for(int j = 0; j < n2; j++) {
    for(int i = 0; i < n1; i++) {
      d(i,j) = d1[i]*d2[j];
    } 
  }
}

void TimeEvolution::dump_psi()
{
  if(!apply_dump()) return;
  
  setup_dump();
  
  const int &n1 = r1.n;
  const int &n2 = r2.n;
  const int &n_theta = theta.n;

  const RMat d(n1, n2, dump);

#pragma omp parallel for			\
  default(shared) schedule(static, 1)	    
  for(int k = 0; k < n_theta; k++) {
    Mat<Complex> Psi(n1, n2, psi+k*n1*n2);
    for(int j = 0; j < n2; j++) {
      for(int i = 0; i < n1; i++) {
	Psi(i,j) *= d(i,j);
      }
    }
  }
}

void  TimeEvolution::setup_weight_legendre()
{
  if(weight_legendre) return;

  const int &n_theta = theta.n;
  const int m = theta.m + 1;

  weight_legendre = new double[n_theta*m];
  insist(weight_legendre);
  
  const double *w = theta.w;
  const RMat &P = theta.legendre;
  
  RMat wp(n_theta, m, weight_legendre);

#pragma omp parallel for			\
  default(shared) schedule(static, 1)
  for(int l = 0; l < m; l++) {
    const double f = l+0.5;
    for(int k = 0; k < n_theta; k++) {
      wp(k,l) = f*w[k]*P(k,l);
    }
  }
}
