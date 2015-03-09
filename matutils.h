
/* $Id$ */

#ifndef MATUTILS_H
#define MATUTILS_H

#include <mex.h>

#define MCrash(x) MatlabCrashLoc(x, __FILE__, __LINE__)
#define MatCrash(x) MatlabCrashLoc(x, __FILE__, __LINE__)
#define MatlabCrash(x) MatlabCrashLoc(x, __FILE__, __LINE__)

#define insist(x) if (!(x)) MatlabCrashLoc("insist failed: " #x, __FILE__, __LINE__)

void MatlabCrashLoc(const char *message, const char *file_name, const int line);

inline void *mxGetData(const mxArray *mx, const char *field)
{
  insist(mx);
  mxArray *mxPtr = mxGetField(mx, 0, field);
  insist(mxPtr);
  void *ptr = mxGetData(mxPtr);
  insist(ptr);
  return ptr;
}

#endif /* MATUTILS_H */
