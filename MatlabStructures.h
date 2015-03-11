
/* $Id$ */

#ifndef MATLAB_STRUCTURES_H
#define MATLAB_STRUCTURES_H

#include <iostream>
using namespace std;
#include <mex.h>
#include "rmat.h"

class RadialCoordinate
{
public:
  const int &n; // out
  RVec r;
  const double &dr; // out
  const double &mass; // out

  RVec psq2m;
  RVec one2mr2;

  RadialCoordinate(const mxArray *mx);
  
private:
  // to prevent assigment and copy operation
  RadialCoordinate(const RadialCoordinate &);
  RadialCoordinate & operator =(const RadialCoordinate &);
  
  /* IO */
  friend ostream & operator <<(ostream &s, const RadialCoordinate &c);
  void write_fields(ostream &s) const;
};

class AngleCoordinate
{
public:
  const int &n; // out
  const int &m; // out
  RVec x;
  RVec w;
  RMat legendre; 

  AngleCoordinate(const mxArray *mx);

private:
  // to prevent assigment and copy operation
  AngleCoordinate(const AngleCoordinate &);
  AngleCoordinate & operator =(const AngleCoordinate &);

  /* IO */
  friend ostream & operator <<(ostream &s, const AngleCoordinate &c);
  void write_fields(ostream &s) const;
};

class EvolutionTime
{
public:
  const int &total_steps; // out
  const double &time_step; // out
  int &steps; // out

  EvolutionTime(const mxArray *mx);

private:
  EvolutionTime(const EvolutionTime &);
  EvolutionTime & operator =(const EvolutionTime &);
  
  /* IO */
  friend ostream & operator <<(ostream &s, const EvolutionTime &c);
  void write_fields(ostream &s) const;
};

#endif /* MATLAB_STRUCTURES_H */
