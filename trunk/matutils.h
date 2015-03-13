
/* $Id$ */

#ifndef MATUTILS_H
#define MATUTILS_H

#include <mex.h>
#include <cstring>

#define MCrash(x) MatlabCrashLoc(x, __FILE__, __LINE__)
#define MatCrash(x) MatlabCrashLoc(x, __FILE__, __LINE__)
#define MatlabCrash(x) MatlabCrashLoc(x, __FILE__, __LINE__)

#define insist(x) if (!(x)) MatlabCrashLoc("insist failed: " #x, __FILE__, __LINE__)

void MatlabCrashLoc(const char *message, const char *file_name, const int line);

void wavepacket_to_matlab(const char *script, const int nrhs = 0, mxArray *prhs[] = 0);

inline void *mxGetData(const mxArray *mx, const char *field)
{
  insist(mx);
  mxArray *mxPtr = mxGetField(mx, 0, field);
  insist(mxPtr);
  void *ptr = mxGetData(mxPtr);
  insist(ptr);
  return ptr;
}

inline char *mxGetString(const mxArray *mx, const char *field)
{
  insist(mx);
  char *tmp = mxArrayToString(mxGetField(mx, 0, field));
  if(!tmp) return 0;
  
  char *string = new char [strlen(tmp) + 1];
  insist(string);
  memcpy(string, tmp, strlen(tmp)*sizeof(char));
  string[strlen(tmp)] = '\0';
  if(tmp) { mxFree(tmp); tmp = 0; }
  return string;
}

#endif /* MATUTILS_H */
