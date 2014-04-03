// Written by Peter Kutz.

#include "constants.h"

// TODO: Adjust the various settings below, and add more if you want.

const int theMillisecondsPerFrame = 10;

#ifdef _DEBUG
const int theDim[3] = {9, 20, 1};
#else
//const int theDim[3] = {31, 31, 31};
const int theDim[3] = {9, 20, 1};
#endif

const double theCellSize = 0.5;

const double AMBIENT_TEMP = 270.0f;

// sources
const double SOURCE_DENSITY = 1.0f;
const double SOURCE_TEMPERATURE = 280.0f;
const double SOURCE_VELOCITY = 2.0f;

const double FLUID_DENSITY = 1.0f;
const double INITIAL_TEMPERATURE = 0.0f;
const double PARTICLE_MASS = 1.0f;
const double FLUID_VISCOSITY = 3.5f;

// external forces
const bool BUOYANCY_FORCE = true;
const bool VORTICITY_CONFINEMENT_FORCE = true;
const bool VISCOSITY_FORCE = false;

// buoyancy force
// target.mTemp( i, j, k ) = -1.0f * alpha * mD( i, j, k ) + beta * ( mT( i, j, k )  - ambient_temp );
const double BUOYANCY_ALPHA = 0.5f;
const double BUOYANCY_BETA = 0.01f;


// to bypass preconditioner in MACGrid::conjugateGradient() method
const bool USE_PRECONDITIONER = true;