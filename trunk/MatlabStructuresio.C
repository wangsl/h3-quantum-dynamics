
/* created at: 2015-03-13 10:08:40 */

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

ostream & operator <<(ostream &s, const Options &c)
{
  s << " {\n";
  IndentPush();
  c.write_fields(s);
  IndentPop();
  return s << Indent() << " }";
}

void Options::write_fields(ostream &s) const
{
  if (wave_to_matlab)
    s << Indent() << "wave_to_matlab " << wave_to_matlab << "\n";
  if (test_name)
    s << Indent() << "test_name " << test_name << "\n";
}

