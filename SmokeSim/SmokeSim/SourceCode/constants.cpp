// Written by Peter Kutz.

#include "constants.h"

// TODO: Adjust the various settings below, and add more if you want.

const int theMillisecondsPerFrame = 10;

#ifdef _DEBUG
const int theDim[3] = {4, 10, 4};
#else
const int theDim[3] = {2, 2, 1};//12, 12, 4
#endif

const double theCellSize = 0.5;

