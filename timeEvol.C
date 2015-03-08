
/* $Id$ */

#include "timeEvol.h"

TimeEvolution::TimeEvolution(double *pot_, Complex *psi_) :
  pot(pot_), psi(psi_)
{ }

TimeEvolution::~TimeEvolution()
{
  pot = 0;
  psi = 0;
}

