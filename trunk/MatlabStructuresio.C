
/* created at: 2015-03-10 23:05:35 */

#include <iostream>
using namespace std;
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "indent.h"
#include "MatlabStructures.h"
#include "die.h"

ostream & operator <<(ostream &s, const RadialCoordinate &c)
{
  s << " {\n";
  IndentPush();
  c.write_fields(s);
  IndentPop();
  return s << Indent() << " }";
}

void RadialCoordinate::write_fields(ostream &s) const
{
  s << Indent() << "n " << n << "\n";
  s << Indent() << "dr " << dr << "\n";
  s << Indent() << "mass " << mass << "\n";
}

ostream & operator <<(ostream &s, const AngleCoordinate &c)
{
  s << " {\n";
  IndentPush();
  c.write_fields(s);
  IndentPop();
  return s << Indent() << " }";
}

void AngleCoordinate::write_fields(ostream &s) const
{
  s << Indent() << "n " << n << "\n";
  s << Indent() << "m " << m << "\n";
}

ostream & operator <<(ostream &s, const EvolutionTime &c)
{
  s << " {\n";
  IndentPush();
  c.write_fields(s);
  IndentPop();
  return s << Indent() << " }";
}

void EvolutionTime::write_fields(ostream &s) const
{
  s << Indent() << "total_steps " << total_steps << "\n";
  s << Indent() << "time_step " << time_step << "\n";
  s << Indent() << "steps " << steps << "\n";
}

