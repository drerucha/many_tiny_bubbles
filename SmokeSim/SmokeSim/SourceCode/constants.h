// Written by Peter Kutz.

#ifndef CONSTANTS_H
#define CONSTANTS_H



#define LERP(a,b,t) (1-t)*a + t*b


// Don't try to modify the values of these here.
// Modify the values of these in constants.cpp instead.
extern const int theMillisecondsPerFrame;
extern const int theDim[3];
extern const double theCellSize;

extern const double AMBIENT_TEMP;

// sources
extern const double SOURCE_DENSITY;
extern const double SOURCE_TEMPERATURE;
extern const double SOURCE_VELOCITY;

extern const double FLUID_DENSITY;
extern const double INITIAL_TEMPERATURE;
extern const double PARTICLE_MASS;
extern const double FLUID_VISCOSITY;

// external forces
extern const bool BUOYANCY_FORCE;
extern const bool VORTICITY_CONFINEMENT_FORCE;
extern const bool VISCOSITY_FORCE;

// buoyancy force
extern const double BUOYANCY_ALPHA;
extern const double BUOYANCY_BETA;

// to bypass preconditioner in MACGrid::conjugateGradient() method
extern const bool USE_PRECONDITIONER;

#endif