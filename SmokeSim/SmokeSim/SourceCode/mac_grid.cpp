// Modified by Peter Kutz, 2011 and 2012.

#include "mac_grid.h"
#include "open_gl_headers.h"
#include "camera.h"
#include "custom_output.h"
#include "constants.h"
#include <math.h>
#include <map>
#include <stdio.h>
#undef max
#undef min
#include <fstream>


////////////////////////////////////////////////////
// globals
////////////////////////////////////////////////////

MACGrid target;

// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = false;

#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_CELL_REVERSE \
   for(int k = theDim[MACGrid::Z] - 1; k >= 0; k--)  \
      for(int j = theDim[MACGrid::Y] - 1; j >= 0; j--) \
         for(int i = theDim[MACGrid::X] - 1; i >= 0; i--) 

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) 


////////////////////////////////////////////////////
// default constructor
////////////////////////////////////////////////////
MACGrid::MACGrid()
{
   initialize();
}


////////////////////////////////////////////////////
// copy constructor
////////////////////////////////////////////////////
MACGrid::MACGrid( const MACGrid& orig )
{
	mU = orig.mU;
	mV = orig.mV;
	mW = orig.mW;
	mP = orig.mP;
	mD = orig.mD;
	mT = orig.mT;
	mTemp = orig.mTemp;
	mConfForceX = orig.mConfForceX;
	mConfForceY = orig.mConfForceY;
	mConfForceZ = orig.mConfForceZ;
}


////////////////////////////////////////////////////
// assignment operator
////////////////////////////////////////////////////
MACGrid& MACGrid::operator=( const MACGrid& orig )
{
	if (&orig == this)
	{
		return *this;
	}
	mU = orig.mU;
	mV = orig.mV;
	mW = orig.mW;
	mP = orig.mP;
	mD = orig.mD;
	mT = orig.mT;   
	mTemp = orig.mTemp;
	mConfForceX = orig.mConfForceX;
	mConfForceY = orig.mConfForceY;
	mConfForceZ = orig.mConfForceZ;
   return *this;
}


////////////////////////////////////////////////////
// destructor
////////////////////////////////////////////////////
MACGrid::~MACGrid()
{
}


////////////////////////////////////////////////////
// reset()
////////////////////////////////////////////////////
void MACGrid::reset()
{
	mU.initialize();
	mV.initialize();
	mW.initialize();
	mP.initialize();
	mD.initialize();
	mT.initialize( INITIAL_TEMPERATURE );
	mTemp.initialize();
	mConfForceX.initialize();
	mConfForceY.initialize();
	mConfForceZ.initialize();
	setUpAMatrix();
}


////////////////////////////////////////////////////
// initialize()
////////////////////////////////////////////////////
void MACGrid::initialize()
{
   reset();
}


////////////////////////////////////////////////////
// updateSources()
// set mT, mD, mV to "push" smoke into the container
////////////////////////////////////////////////////
void MACGrid::updateSources()
{
	// TODO: create intuitive constants to control size and velocity of sources

	// 1x1 source in middle of 9x1 x/z grid

	mT( 4, 0, 0 ) = SOURCE_TEMPERATURE;
	mD( 4, 0, 0 ) = SOURCE_DENSITY;
	mV( 4, 1, 0 ) = SOURCE_VELOCITY;


	// 3x3 source in middle of 30x30 x/z grid

	//mT( 14, 0, 14 ) = 280.0f;
	//mD( 14, 0, 14 ) = 1.0f;
	//mV( 14, 1, 14 ) = 2.0f;

	//mT( 14, 0, 15 ) = 280.0f;
	//mD( 14, 0, 15 ) = 1.0f;
	//mV( 14, 1, 15 ) = 2.0f;

	//mT( 14, 0, 16 ) = 280.0f;
	//mD( 14, 0, 16 ) = 1.0f;
	//mV( 14, 1, 16 ) = 2.0f;

	//mT( 15, 0, 14 ) = 280.0f;
	//mD( 15, 0, 14 ) = 1.0f;
	//mV( 15, 1, 14 ) = 2.0f;

	//mT( 15, 0, 15 ) = 280.0f;
	//mD( 15, 0, 15 ) = 1.0f;
	//mV( 15, 1, 15 ) = 2.0f;

	//mT( 15, 0, 16 ) = 280.0f;
	//mD( 15, 0, 16 ) = 1.0f;
	//mV( 15, 1, 16 ) = 2.0f;

	//mT( 16, 0, 14 ) = 280.0f;
	//mD( 16, 0, 14 ) = 1.0f;
	//mV( 16, 1, 14 ) = 2.0f;

	//mT( 16, 0, 15 ) = 280.0f;
	//mD( 16, 0, 15 ) = 1.0f;
	//mV( 16, 1, 15 ) = 2.0f;

	//mT( 16, 0, 16 ) = 280.0f;
	//mD( 16, 0, 16 ) = 1.0f;
	//mV( 16, 1, 16 ) = 2.0f;
}


////////////////////////////////////////////////////
// generateBubbles()
////////////////////////////////////////////////////
void MACGrid::generateBubbles()
{
	vec3 position( theDim[MACGrid::X] * theCellSize / 2.0f, theCellSize, theCellSize / 2.0f );
	bubblePos.push_back( position );
}


////////////////////////////////////////////////////
// advectVelocity()
// use traceback method to update velocities at each cell face
////////////////////////////////////////////////////
void MACGrid::advectVelocity( double dt )
{
	// loop through every MacGrid face
	// i, j, k are defined by FOR_EACH_FACE constant
	FOR_EACH_FACE {
		// mU = GridDataX, mV = GridDataY, mW = GridDataZ
		double velU = mU( i, j, k );
		double velV = mV( i, j, k );
		double velW = mW( i, j, k );

		// compute face center positions
		vec3 centerPosition = getCenter( i, j, k );
		vec3 grid_x_bottom_border_pos = centerPosition - vec3( theCellSize * 0.5f, 0.0f, 0.0f );
		vec3 grid_y_bottom_border_pos = centerPosition - vec3( 0.0f, theCellSize * 0.5f, 0.0f );
		vec3 grid_z_bottom_border_pos = centerPosition - vec3( 0.0f, 0.0f, theCellSize * 0.5f );


		////////////////////////////////////////////////////
		// compute velocities at face centers
		////////////////////////////////////////////////////

		vec3 grid_x_bottom_border_vel, grid_y_bottom_border_vel, grid_z_bottom_border_vel;

		// compute velocity at grid_x_bottom_border_pos
		if ( i == 0 ) {
			// low boundary cell
			grid_x_bottom_border_vel[VX] = velU;
			grid_x_bottom_border_vel[VY] = 0.5f * ( velV + mV( i, j+1, k ) );
			grid_x_bottom_border_vel[VZ] = 0.5f * ( velW + mW( i, j, k+1 ) );
		}
		else if ( i == theDim[MACGrid::X] ) {
			// high boundary cell
			grid_x_bottom_border_vel[VX] = velU;
			grid_x_bottom_border_vel[VY] = 0.5f * ( mV( i-1, j, k ) + mV( i-1, j+1, k ) );
			grid_x_bottom_border_vel[VZ] = 0.5f * ( mW( i-1, j, k ) + mW( i-1, j, k+1 ) );
		}
		else {
			// not boundary cell - cell is somehwere in middle of container
			grid_x_bottom_border_vel[VX] = velU;
			grid_x_bottom_border_vel[VY] = 0.25f * ( velV + mV( i-1, j, k ) + mV( i, j+1, k ) + mV( i-1, j+1, k ) );
			grid_x_bottom_border_vel[VZ] = 0.25f * ( velW + mW( i-1, j, k ) + mW( i, j, k+1 ) + mW( i-1, j, k+1 ) );
		}

		// compute velocity at grid_y_bottom_border_pos
		if ( j == 0 ) {
			// low boundary cell
			grid_y_bottom_border_vel[VX] = 0.5f * ( velU + mU( i+1, j, k ) );
			grid_y_bottom_border_vel[VY] = velV;
			grid_y_bottom_border_vel[VZ] = 0.5f * ( velW + mW( i, j, k+1 ) );
		}
		else if ( j == theDim[MACGrid::Y] ) {
			// high boundary cell
			grid_y_bottom_border_vel[VX] = 0.5f * (  mU( i, j-1, k ) + mU( i+1, j-1, k ) );
			grid_y_bottom_border_vel[VY] = velV;
			grid_y_bottom_border_vel[VZ] = 0.5f * ( mW( i, j-1, k ) + mW( i, j-1, k+1 ) );
		}
		else {
			// not boundary cell - cell is somehwere in middle of container
			grid_y_bottom_border_vel[VX] = 0.25f * ( velU + mU( i, j-1, k ) + mU( i+1, j, k ) + mU( i+1, j-1, k ) );
			grid_y_bottom_border_vel[VY] = velV;
			grid_y_bottom_border_vel[VZ] = 0.25f * ( velW + mW( i, j-1, k ) + mW( i, j, k+1 ) + mW( i, j-1, k+1 ) );
		}

		// compute velocity at grid_z_bottom_border_pos
		if ( k == 0 ) {
			// low boundary cell
			grid_z_bottom_border_vel[VX] = 0.5f * ( velU + mU( i+1, j, k ) );
			grid_z_bottom_border_vel[VY] = 0.5f * ( velV + mV( i, j+1, k ) );
			grid_z_bottom_border_vel[VZ] = velW;
		}
		else if ( k == theDim[MACGrid::Z] ) {
			// high boundary cell
			grid_z_bottom_border_vel[VX] = 0.5f * ( mU( i, j, k-1 ) + mU( i+1, j, k-1 ) );
			grid_z_bottom_border_vel[VY] = 0.5f * ( mV( i, j, k-1 ) + mV( i, j+1, k-1 ) );
			grid_z_bottom_border_vel[VZ] = velW;
		}
		else {
			// not boundary cell - cell is somehwere in middle of container
			grid_z_bottom_border_vel[VX] = 0.25f * ( velU + mU( i, j, k-1 ) + mU( i+1, j, k ) + mU( i+1, j, k-1 ) );
			grid_z_bottom_border_vel[VY] = 0.25f * ( velV + mV( i, j, k-1 ) + mV( i, j+1, k ) + mV( i, j+1, k-1 ) );
			grid_z_bottom_border_vel[VZ] = velW;
		}

		// solve for advection
		velU = getVelocityX( grid_x_bottom_border_pos - ( dt * grid_x_bottom_border_vel ) );
		velV = getVelocityY( grid_y_bottom_border_pos - ( dt * grid_y_bottom_border_vel ) );
		velW = getVelocityZ( grid_z_bottom_border_pos - ( dt * grid_z_bottom_border_vel ) );

		// store in target, checking boundaries
		if ( i > 0 && i < theDim[MACGrid::X] ) {
			target.mU( i, j, k ) = velU;
		}
		else {
			target.mU( i, j, k ) = 0.0f;
		}

		if ( j > 0 && j < theDim[MACGrid::Y] ) {
			target.mV( i, j, k ) = velV;
		}
		else {
			target.mV( i, j, k ) = 0.0f;
		}

		if ( k != 0 && k != theDim[MACGrid::Z] ) {
			target.mW( i, j, k ) = velW;
		}
		else {
			target.mW( i, j, k ) = 0.0f;
		}
	}

    // save result to object
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
}


////////////////////////////////////////////////////
// advectTemperature()
// use traceback method to update temperatures at each cell center
////////////////////////////////////////////////////
void MACGrid::advectTemperature( double dt )
{
	FOR_EACH_CELL {
		// velocity at cell center
		vec3 velocity( ( mU( i, j, k ) + mU( i+1, j, k ) ) / 2.0f,
					   ( mV( i, j, k ) + mV( i, j+1, k ) ) / 2.0f,
					   ( mW( i, j, k ) + mW( i, j, k+1 ) ) / 2.0f );

		// trace back particle position using known velocity
		vec3 pos = getCenter( i, j, k );
		pos -= dt * velocity;

		// interpolate temperature for passed-in position, and store in target
		target.mT( i, j, k ) = getTemperature( pos );
	}

	// save result to object
	mT = target.mT;
}


////////////////////////////////////////////////////
// advectDensity()
// use traceback method to update temperatures at each cell center
////////////////////////////////////////////////////
void MACGrid::advectDensity( double dt )
{
	// use an identical trace back method to the one used in MACGrid::advectTemperature()
	FOR_EACH_CELL {
		vec3 velocity( ( mU( i, j, k ) + mU( i+1, j, k ) ) / 2.0f,
					   ( mV( i, j, k ) + mV( i, j+1, k ) ) / 2.0f,
					   ( mW( i, j, k ) + mW( i, j, k+1 ) ) / 2.0f );

		vec3 pos = getCenter( i, j, k );
		pos -= dt * velocity;
		target.mD( i, j, k ) = getDensity( pos );
	}

	mD = target.mD;
}


////////////////////////////////////////////////////
// addExternalForces()
////////////////////////////////////////////////////
void MACGrid::addExternalForces( double dt )
{
	// to activate/deactivate these external forces, set the appropriate constants in constants.cpp

	if ( BUOYANCY_FORCE ) {
		computeBouyancy( dt );
	}
	if ( VORTICITY_CONFINEMENT_FORCE ) {
		computeVorticityConfinement( dt );
	}
	if ( VISCOSITY_FORCE ) {
		computeViscosityForce( dt );
	}
}


////////////////////////////////////////////////////
// computeBouyancy()
// add external force influenced by liquid's temperature
////////////////////////////////////////////////////
void MACGrid::computeBouyancy( double dt )
{
	// set BUOYANCY_ALPHA, BUOYANCY_BETA, AMBIENT_TEMP in constants.cpp

	FOR_EACH_CELL {
		target.mTemp( i, j, k ) = -1.0f * BUOYANCY_ALPHA * mD( i, j, k ) + BUOYANCY_BETA * ( mT( i, j, k )  - AMBIENT_TEMP );
	}

	FOR_EACH_CELL {
		if ( j != 0 ) {
			// get cell center position
			vec3 cell_center_pos = getCenter( i, j, k );

			// get y face center position
			vec3 y_face_center_pos = cell_center_pos - vec3( 0.0f, theCellSize * 0.5f, 0.0f );

			// interpolate buoyancy force from cell center to face center
			double buoyancy_force = target.mTemp.interpolate( y_face_center_pos );
		
			// acceleration = force on particle / mass of particle
			double acceleration = buoyancy_force / PARTICLE_MASS;

			// perform explicit Euler integration to update existing velocity with applied force
			// v' = v + at
			target.mV( i, j, k ) = mV( i, j, k ) + dt * acceleration;
		}
	}

	mV = target.mV;
}


////////////////////////////////////////////////////
// computeVorticityConfinement()
// add external force to preserve localized vortices
////////////////////////////////////////////////////
void MACGrid::computeVorticityConfinement( double dt )
{
	//double epsilon = 1.0f;
	double epsilon = 0.1f;
	double two_cell_widths = 2.0f * theCellSize;
	double very_small = pow( 10.0f, -20.0f );

	FOR_EACH_CELL {

		// TODO: ask about boundary conditions here
		// index out of bounds when i == 0, 1, theDim[MACGrid::X]
		// index out of bounds when j == 0, 1, theDim[MACGrid::Y]
		// index out of bounds when k == 0, 1, theDim[MACGrid::Z]
		// velocities return 0 by default for these cases

		vec3 omegaGradient( ( getOmegaVector( i+1, j, k ).Length() - getOmegaVector( i-1, j, k ).Length() ) / two_cell_widths,
							( getOmegaVector( i, j+1, k ).Length() - getOmegaVector( i, j-1, k ).Length() ) / two_cell_widths,
							( getOmegaVector( i, j, k+1 ).Length() - getOmegaVector( i, j, k-1 ).Length() ) / two_cell_widths );
	
		// add very_small to prevent divide by zero
		vec3 normal = omegaGradient / ( omegaGradient.Length() + very_small );

		//vec3 temp = getOmegaVector( i, j, k );
		vec3 confinement = epsilon * theCellSize * normal.Cross( getOmegaVector( i, j, k ) );
		target.mConfForceX(i, j, k) = confinement[0];
		target.mConfForceY(i, j, k) = confinement[1];
		target.mConfForceZ(i, j, k) = confinement[2];
	}

	// take vorticity confinement forces computed at cell centers and approximate at faces for velocity field 
	FOR_EACH_CELL {

		// TODO: ask about boundary conditions here
		// currently, faces with 0 indices are being ignored

		// first, approximate vorticity force at face center with linear interpolation
		// next, acceleration = force on particle / mass of particle
		// finally, perform explicit Euler integration to update existing velocity with applied force
		// v' = v + at
		if ( i != 0 ) {
			double vorticity_force_x = ( target.mConfForceX( i, j, k ) + target.mConfForceX( i-1, j, k ) ) / 2.0f;
			double acceleration_x = vorticity_force_x / PARTICLE_MASS;
			target.mU( i, j, k ) = mU( i, j, k ) + dt * acceleration_x;
		}
		if ( j != 0 ) {
			double vorticity_force_y = ( target.mConfForceY( i, j, k ) + target.mConfForceY( i, j-1, k ) ) / 2.0f;
			double acceleration_y = vorticity_force_y / PARTICLE_MASS;
			target.mV( i, j, k ) = mV( i, j, k ) + dt * acceleration_y;
		}
		if ( k != 0 ) {
			double vorticity_force_z = ( target.mConfForceZ( i, j, k ) + target.mConfForceZ( i, j, k-1 ) ) / 2.0f;
			double acceleration_z = vorticity_force_z / PARTICLE_MASS;
			target.mW( i, j, k ) = mW( i, j, k ) + dt * acceleration_z;	
		}
	}

	// save result to object
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}


////////////////////////////////////////////////////
// getOmegaVector()
// helper method for computing vorticity confinement
////////////////////////////////////////////////////
vec3 MACGrid::getOmegaVector( int i, int j, int k )
{
	double half_cell_width = 1.0f / ( 2.0f * theCellSize );

	// return omega vector
	return vec3( ( mW( i, j+1, k ) - mW( i, j-1, k ) ) * half_cell_width - ( mV( i, j, k+1 ) - mV( i, j, k-1 ) ) * half_cell_width,
				 ( mU( i, j, k+1 ) - mU( i, j, k-1 ) ) * half_cell_width - ( mW( i+1, j, k ) - mW( i-1, j, k ) ) * half_cell_width,
				 ( mV( i+1, j, k ) - mV( i-1, j, k ) ) * half_cell_width - ( mU( i, j+1, k ) - mU( i, j-1, k ) ) * half_cell_width );
}


////////////////////////////////////////////////////
// computeViscosityForce()
// add external force to control fluid's thickness (higher viscosity == more like honey than water)
////////////////////////////////////////////////////
void MACGrid::computeViscosityForce( double dt )
{
	FOR_EACH_CELL {
		// compute velocity differences between adjacent cells
		double velDiff_I_IMinus1, velDiff_IPlus1_I,
			   velDiff_J_JMinus1, velDiff_JPlus1_J,
			   velDiff_K_KMinus1, velDiff_KPlus1_K;

		// velDiff_I_IMinus1
		if ( i > 0 ) {
			velDiff_I_IMinus1 = mU( i, j, k ) - mU( i-1, j, k );
		}
		else {
			velDiff_I_IMinus1 = 0.0f;
		}

		// velDiff_IPlus1_I
		if ( i < theDim[MACGrid::X] ) {
			velDiff_IPlus1_I  = mU( i+1, j, k ) - mU( i, j, k );
		}
		else {
			velDiff_IPlus1_I = 0.0f;
		}

		// velDiff_J_JMinus1
		if ( j > 0 ) {
			velDiff_J_JMinus1  = mV( i, j, k ) - mV( i, j-1, k );
		}
		else {
			velDiff_J_JMinus1  = 0.0f;
		}

		// velDiff_JPlus1_J
		if ( i < theDim[MACGrid::Y] ) {
			velDiff_JPlus1_J  = mV( i, j+1, k ) - mV( i, j, k );
		}
		else {
			velDiff_JPlus1_J = 0.0f;
		}

		// velDiff_K_KMinus1
		if ( k > 0 ) {
			velDiff_K_KMinus1 = mW( i, j, k ) - mW( i, j, k-1 );
		}
		else {
			velDiff_K_KMinus1 = 0.0f;
		}

		// velDiff_KPlus1_K
		if ( k < theDim[MACGrid::Z] ) {
			velDiff_KPlus1_K = mW( i, j, k+1 ) - mW( i, j, k );
		}
		else {
			velDiff_KPlus1_K = 0.0f;
		}

		double one_over_cell_width_squared = 1.0f / ( theCellSize * theCellSize );

		// acceleration = force on particle / mass of particle
		// finally, perform explicit Euler integration to update existing velocity with applied force
		// v' = v + at
		if ( i != 0 ) {
			double viscosityForceX = FLUID_VISCOSITY * ( velDiff_IPlus1_I - velDiff_I_IMinus1 ) * one_over_cell_width_squared;
			double acceleration_x = viscosityForceX / PARTICLE_MASS;
			target.mU( i, j, k ) = mU( i, j, k ) + dt * acceleration_x;
		}
		if ( j != 0 ) {
			double viscosityForceY = FLUID_VISCOSITY * ( velDiff_JPlus1_J - velDiff_J_JMinus1 ) * one_over_cell_width_squared;
			double acceleration_y = viscosityForceY / PARTICLE_MASS;
			target.mV( i, j, k ) = mV( i, j, k ) + dt * acceleration_y;
		}
		if ( k != 0 ) {
			double viscosityForceZ = FLUID_VISCOSITY * ( velDiff_KPlus1_K - velDiff_K_KMinus1 ) * one_over_cell_width_squared;
			double acceleration_z = viscosityForceZ / PARTICLE_MASS;
			target.mW( i, j, k ) = mW( i, j, k ) + dt * acceleration_z;	
		}

	}

	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}


////////////////////////////////////////////////////
// project()
// solve for pressure in a manner that creates a divergence free system
////////////////////////////////////////////////////
void MACGrid::project( double dt )
{
	////////////////////////////////////////////////////
	// solve for pressure
	////////////////////////////////////////////////////

	double constant = -1.0f * theCellSize * theCellSize / dt;

	// AMatrix is a MACGrid member and the pressure coefficient matrix
	// setUpAMatrix() initializes AMatrix given the current grid structure
	setUpAMatrix();

	// d is divergence vector; p is pressure vector
	GridData d, p;
	d.initialize( 0.0f );
	p.initialize( 0.0f );

	FOR_EACH_CELL {
		// fill divergence vector - ( constant * density * velocity gradient )
		d( i, j, k ) = constant * FLUID_DENSITY * ( ( mU( i+1, j, k ) - mU( i, j, k ) ) / theCellSize + 
													( mV( i, j+1, k ) - mV( i, j, k ) ) / theCellSize +
													( mW( i, j, k+1 ) - mW( i, j, k ) ) / theCellSize );
	}

	// solve pressure vector by solving system of linear equations
	int max_num_iterations = 100;
	double tolerance = 0.00001f;

	// conjugateGradient method includes preconditioner
	conjugateGradient( AMatrix, p, d, max_num_iterations, tolerance );

	// store in target
	target.mP = p;
	

	////////////////////////////////////////////////////
	// apply computed pressures to velocity field
	////////////////////////////////////////////////////

	FOR_EACH_FACE {
		double velU = mU( i, j, k );
		double velV = mV( i, j, k );
		double velW = mW( i, j, k );

		double current_pressure = target.mP( i, j, k );

		// to check boundary conditions
		double pressure_i_minus_1, pressure_j_minus_1, pressure_k_minus_1;

		// set pressure_i_minus_1
		if ( i-1 < 0 ) {
			pressure_i_minus_1 = 0.0f;
		}
		else {
			pressure_i_minus_1 = target.mP( i-1, j, k );
		}

		// set pressure_j_minus_1
		if ( j-1 < 0 ) {
			pressure_j_minus_1 = 0.0f;
		}
		else {
			pressure_j_minus_1 = target.mP( i, j-1, k );
		}

		// set pressure_k_minus_1
		if ( k-1 < 0 ) {
			pressure_k_minus_1 = 0.0f;
		}
		else {
			pressure_k_minus_1 = target.mP( i, j, k-1 );
		}

		// apply computed pressures to velocity field
		// if face is an outside border of container, then velocity at that face must be 0
		if ( i == 0 || i == theDim[MACGrid::X] ) {
			velU = 0.0f;
		}
		else {
			velU -= dt * ( 1.0f / FLUID_DENSITY ) * ( ( target.mP( i, j, k ) - pressure_i_minus_1 ) / theCellSize );
		}
		
		if ( j == 0 || j == theDim[MACGrid::Y] ) {
			velV = 0.0f;
		}
		else {		
			velV -= dt * ( 1.0f / FLUID_DENSITY ) * ( ( target.mP( i, j, k ) - pressure_j_minus_1 ) / theCellSize );
		}
		
		if ( k == 0 || k == theDim[MACGrid::Z] ) {
			velW = 0.0f;
		}
		else {
			velW -= dt * ( 1.0f / FLUID_DENSITY ) * ( ( target.mP( i, j, k ) - pressure_k_minus_1 ) / theCellSize );
		}

		// store in target with additional container boundary checking
		// ex: for a 2x2 matrix, the j value for mU can be only 0 or 1, but j can have values of 0, 1, or 2 for a 2x2 matrix
		if ( j < theDim[MACGrid::Y] && k < theDim[MACGrid::Z] ) {
			target.mU( i, j, k ) = velU;
		}
		if ( i < theDim[MACGrid::X] && k < theDim[MACGrid::Z] ) {
			target.mV( i, j, k ) = velV;
		}
		if ( i < theDim[MACGrid::X] && j < theDim[MACGrid::Y] ) {
			target.mW( i, j, k ) = velW;
		}
	}
	
	// save result to object
	mP = target.mP;
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;

	// debug
	if ( !testForDivergenceFreeSystem() ) {
		std::cout << "it appears the system is not divergence free" << std::endl;
		std::cin.ignore();
	}
}


////////////////////////////////////////////////////
// testForDivergenceFreeSystem()
// method used for debugging purposes
////////////////////////////////////////////////////
bool MACGrid::testForDivergenceFreeSystem() const
{
	FOR_EACH_CELL {
		double sum = ( mU( i+1, j, k ) - mU( i,j,k ) ) + 
					 ( mV( i, j+1, k ) - mV( i,j,k ) ) +
					 ( mW( i, j, k+1 ) - mW( i,j,k ) );

		if ( abs( sum ) > 0.01 ) {
			return false;
		}
	}

	return true;
}


////////////////////////////////////////////////////
// conjugateGradient()
// solves Ap = d for p, includes Modified Incomplete Cholesky preconditioner
////////////////////////////////////////////////////
bool MACGrid::conjugateGradient( const GridDataMatrix& A, GridData& p, const GridData& d, int maxIterations, double tolerance ) {
	FOR_EACH_CELL {
		p( i, j, k ) = 0.0; // initial guess p = 0
	}

	GridData r = d; // residual vector

	GridData z; 
	z.initialize();

	GridData mPrecon;
	mPrecon.initialize();

	GridData mQ;
	mQ.initialize();

	GridData s; // search vector

	// apply preconditioner
	if ( USE_PRECONDITIONER ) {
		applyPreconditioner( A, r, z, mPrecon, mQ, true );
		s = z;
	}
	else {
		z = r;
		s = r;
	}

	double sigma = dotProduct( z, r );

	for ( int iteration = 0; iteration < maxIterations; ++iteration ) {

		double rho = sigma;

		apply(A, s, z);

		double alpha = rho/dotProduct(z, s);

		GridData alphaTimesS; alphaTimesS.initialize();
		multiply(alpha, s, alphaTimesS);
		add(p, alphaTimesS, p);

		GridData alphaTimesZ; alphaTimesZ.initialize();
		multiply(alpha, z, alphaTimesZ);
		subtract(r, alphaTimesZ, r);

		if (maxMagnitude(r) <= tolerance) {
			//PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations.");
			return true;
		}

		// apply preconditioner
		if ( USE_PRECONDITIONER ) {
			applyPreconditioner( A, r, z, mPrecon, mQ, false );
		}
		else {
			z = r;
		}

		double sigmaNew = dotProduct( z, r );

		double beta = sigmaNew / rho;

		GridData betaTimesS; betaTimesS.initialize();
		multiply(beta, s, betaTimesS);
		add(z, betaTimesS, s);
		//s = z + beta * s;

		sigma = sigmaNew;
	}

	PRINT_LINE( "PCG didn't converge!" );
	return false;
}


////////////////////////////////////////////////////
// applyPreconditioner()
// Modified Incomplete Cholesky preconditioner
////////////////////////////////////////////////////
void MACGrid::applyPreconditioner( const GridDataMatrix& A, const GridData& r, GridData& z,
								   GridData& mPrecon, GridData& mQ, const bool& compute_mPrecon )
{
	// TODO: preconditioner method needs to be modified if fluid does not fill entire grid, i.e. water instead of smoke
	// see pg 35 of Bridson's 2007 SIGGRAPH course notes

	if ( compute_mPrecon ) {
		double torque = 0.9;
		FOR_EACH_CELL {
			double E_i_j_k = A.diag( i, j, k ) - pow( A.plusI( i-1, j, k ) * mPrecon( i-1, j, k ), 2 )
											   - pow( A.plusJ( i, j-1, k ) * mPrecon( i, j-1, k ), 2 )
											   - pow( A.plusK( i, j, k-1 ) * mPrecon( i, j, k-1 ), 2 )
											   - torque * ( A.plusI( i-1, j, k ) * ( A.plusJ( i-1, j, k ) + A.plusK( i-1, j, k ) ) * pow( mPrecon( i-1, j, k ), 2 ) 
													+ A.plusJ( i, j-1, k ) * ( A.plusI( i, j-1, k ) + A.plusK( i, j-1, k ) ) * pow( mPrecon( i, j-1, k ), 2 ) 
													+ A.plusK( i, j, k-1 ) * ( A.plusI( i, j, k-1 ) + A.plusJ( i, j, k-1 ) ) * pow( mPrecon( i, j, k-1 ), 2 ) ); 
			mPrecon( i, j, k ) = 1 / sqrt( E_i_j_k + pow( 10.0f, -30 ) );
		}
	}

	FOR_EACH_CELL {
		double t = r( i, j, k ) - A.plusI( i-1, j, k ) * mPrecon( i-1, j, k ) * mQ( i-1, j, k )
								- A.plusJ( i, j-1, k ) * mPrecon( i, j-1, k ) * mQ( i, j-1, k )
								- A.plusK( i, j, k-1 ) * mPrecon( i, j, k-1 ) * mQ( i, j, k-1 );
		mQ( i, j, k ) = t * mPrecon( i, j, k );
	}

	FOR_EACH_CELL_REVERSE {
		double t = mQ( i, j, k ) - A.plusI( i, j, k ) * mPrecon( i, j, k ) * z( i+1, j, k )
								 - A.plusJ( i, j, k ) * mPrecon( i, j, k ) * z( i, j+1, k )
								 - A.plusK( i, j, k ) * mPrecon( i, j, k ) * z( i, j, k+1 );
		z( i, j, k ) = t * mPrecon( i, j, k );
	}
}



/*********** METHODS FOR BUBBLES ***********/


////////////////////////////////////////////////////
// advectBubbles()
////////////////////////////////////////////////////
void MACGrid::advectBubbles( double dt )
{
	int i = 0;
	for(std::vector<vec3>::iterator iter = bubblePos.begin(); iter != bubblePos.end(); ++iter)
	{
		vec3 position = bubblePos[i];
		vec3 velocity = getVelocity(position);
		position += dt * velocity;

		// particle goes outside the container
		if(position[0] < 0 || position[0] > theDim[MACGrid::X] * theCellSize
		|| position[1] < 0 || position[1] > theDim[MACGrid::Y] * theCellSize
		|| position[2] < 0 || position[2] > theDim[MACGrid::Z] * theCellSize)
		{
			bubblePos.erase(iter);
		}
		else
			bubblePos[i] = position;

		i++;
	}
}


////////////////////////////////////////////////////
// getBubblePosition()
////////////////////////////////////////////////////
float* MACGrid::getBubblePosition( int *size )
{
	float *bubblePositionArray = new float[bubblePos.size() * 3];
	*size = bubblePos.size();
	for( unsigned int i = 0; i < bubblePos.size(); i++ ) {
		vec3 position = bubblePos[i];
		bubblePositionArray[i * 3] = position[0];
		bubblePositionArray[i * 3 + 1] = position[1];
		bubblePositionArray[i * 3 + 2] = position[2];
	}

	return bubblePositionArray;
}



/*********** GETTERS ***********/


vec3 MACGrid::getVelocity(const vec3& pt)
{
   vec3 vel;
   vel[0] = getVelocityX(pt); 
   vel[1] = getVelocityY(pt); 
   vel[2] = getVelocityZ(pt); 
   return vel;
}

double MACGrid::getVelocityX(const vec3& pt)
{
   return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt)
{
   return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt)
{
   return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt)
{
   return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt)
{
   return mD.interpolate(pt);
}

vec3 MACGrid::getCenter(int i, int j, int k)
{
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return vec3(x, y, z);
}



/*********** HELPER METHODS ***********/
/*********** THE FOLLOWING METHODS WERE GIVEN TO US BY OUR INSTRUCTOR IN THE FRAMEWORK ***********/


bool MACGrid::isValidCell(int i, int j, int k)
{
	if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
		return false;
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}

void MACGrid::setUpAMatrix() {

	FOR_EACH_CELL {

		int numFluidNeighbors = 0;
		if (i-1 >= 0) {
			AMatrix.plusI(i-1,j,k) = -1;
			numFluidNeighbors++;
		}
		if (i+1 < theDim[MACGrid::X]) {
			AMatrix.plusI(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (j-1 >= 0) {
			AMatrix.plusJ(i,j-1,k) = -1;
			numFluidNeighbors++;
		}
		if (j+1 < theDim[MACGrid::Y]) {
			AMatrix.plusJ(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (k-1 >= 0) {
			AMatrix.plusK(i,j,k-1) = -1;
			numFluidNeighbors++;
		}
		if (k+1 < theDim[MACGrid::Z]) {
			AMatrix.plusK(i,j,k) = -1;
			numFluidNeighbors++;
		}
		// Set the diagonal:
		AMatrix.diag(i,j,k) = numFluidNeighbors;
	}
}

double MACGrid::dotProduct(const GridData & vector1, const GridData & vector2) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		result += vector1(i,j,k) * vector2(i,j,k);
	}

	return result;
}

void MACGrid::add(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
	}

}

void MACGrid::subtract(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
	}

}

void MACGrid::multiply(const double scalar, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = scalar * vector(i,j,k);
	}

}

double MACGrid::maxMagnitude(const GridData & vector) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		if (abs(vector(i,j,k)) > result) result = abs(vector(i,j,k));
	}

	return result;
}

void MACGrid::apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL { // For each row of the matrix.

		double diag = 0;
		double plusI = 0;
		double plusJ = 0;
		double plusK = 0;
		double minusI = 0;
		double minusJ = 0;
		double minusK = 0;

		diag = matrix.diag(i,j,k) * vector(i,j,k);
		if (isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
		if (isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
		if (isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
		if (isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
		if (isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
		if (isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);

		result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
	}

}




/////////////////////////////////////////////////////////////////////

void MACGrid::saveSmoke(const char* fileName) {
	std::ofstream fileOut(fileName);
	if (fileOut.is_open()) {
		FOR_EACH_CELL {
			fileOut << mD(i,j,k) << std::endl;
		}
		fileOut.close();
	}
}





/////////////////////////////////////////////////////////////////////

void MACGrid::draw(const Camera& c)
{   
   drawWireGrid();
   if (theDisplayVel) drawVelocities();   
   if (theRenderMode == CUBES) drawSmokeCubes(c);
   else drawSmoke(c);
}

void MACGrid::drawVelocities()
{
   // Draw line at each center
   //glColor4f(0.0, 1.0, 0.0, 1.0); // Use this if you want the lines to be a single color.
   glBegin(GL_LINES);
      FOR_EACH_CELL
      {
         vec3 pos = getCenter(i,j,k);
         vec3 vel = getVelocity(pos);
         if (vel.Length() > 0.0001)
         {
		   // Un-comment the line below if you want all of the velocity lines to be the same length.
           //vel.Normalize();
           vel *= theCellSize/2.0;
           vel += pos;
		   glColor4f(1.0, 1.0, 0.0, 1.0);
           glVertex3dv(pos.n);
		   glColor4f(0.0, 1.0, 0.0, 1.0);
           glVertex3dv(vel.n);
         }
      }
   glEnd();
}

vec4 MACGrid::getRenderColor(int i, int j, int k)
{

	// Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = mD(i, j, k); 
    return vec4(1.0, 1.0, 1.0, value);

}

vec4 MACGrid::getRenderColor(const vec3& pt)
{

	// TODO: Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = getDensity(pt); 
    return vec4(1.0, 1.0, 1.0, value);

}

void MACGrid::drawZSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double startk = back - stepsize;
   double endk = 0;
   double stepk = -theCellSize;

   if (!backToFront)
   {
      startk = 0;
      endk = back;   
      stepk = theCellSize;
   }

   for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double i = 0.0; i <= right; i += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double i = right; i >= 0.0; i -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawXSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double starti = right - stepsize;
   double endi = 0;
   double stepi = -theCellSize;

   if (!backToFront)
   {
      starti = 0;
      endi = right;   
      stepi = theCellSize;
   }

   for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double k = 0.0; k <= back; k += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double k = back; k >= 0.0; k -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawSmoke(const Camera& c)
{
   vec3 eyeDir = c.getBackward();
   double zresult = fabs(Dot(eyeDir, vec3(1,0,0)));
   double xresult = fabs(Dot(eyeDir, vec3(0,0,1)));
   //double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

   if (zresult < xresult)
   {      
      drawZSheets(c.getPosition()[2] < 0);
   }
   else 
   {
      drawXSheets(c.getPosition()[0] < 0);
   }
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
   std::multimap<double, MACGrid::Cube, std::greater<double> > cubes;
   FOR_EACH_CELL
   {
      MACGrid::Cube cube;
      cube.color = getRenderColor(i,j,k);
      cube.pos = getCenter(i,j,k);
      cube.dist = DistanceSqr(cube.pos, c.getPosition());
      cubes.insert(make_pair(cube.dist, cube));
   } 

   // Draw cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double> >::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawCube(it->second);
   }
}

void MACGrid::drawWireGrid()
{
   // Display grid in light grey, draw top & bottom

   double xstart = 0.0;
   double ystart = 0.0;
   double zstart = 0.0;
   double xend = theDim[0]*theCellSize;
   double yend = theDim[1]*theCellSize;
   double zend = theDim[2]*theCellSize;

   glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);
      glColor3f(0.25, 0.25, 0.25);

      glBegin(GL_LINES);
      for (int i = 0; i <= theDim[0]; i++)
      {
         double x = xstart + i*theCellSize;
         glVertex3d(x, ystart, zstart);
         glVertex3d(x, ystart, zend);

         glVertex3d(x, yend, zstart);
         glVertex3d(x, yend, zend);
      }

      for (int i = 0; i <= theDim[2]; i++)
      {
         double z = zstart + i*theCellSize;
         glVertex3d(xstart, ystart, z);
         glVertex3d(xend, ystart, z);

         glVertex3d(xstart, yend, z);
         glVertex3d(xend, yend, z);
      }

      glVertex3d(xstart, ystart, zstart);
      glVertex3d(xstart, yend, zstart);

      glVertex3d(xend, ystart, zstart);
      glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zend);
      glVertex3d(xstart, yend, zend);

      glVertex3d(xend, ystart, zend);
      glVertex3d(xend, yend, zend);
      glEnd();
   glPopAttrib();

   glEnd();
}

#define LEN 0.5
void MACGrid::drawFace(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);
      glEnd();
   glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0, -1.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN, -LEN);         

         glNormal3d( 0.0,  0.0, -0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN,  LEN, -LEN);
         glVertex3d( LEN,  LEN, -LEN);
         glVertex3d( LEN, -LEN, -LEN);

         glNormal3d(-1.0,  0.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d(-LEN,  LEN,  LEN);
         glVertex3d(-LEN,  LEN, -LEN);

         glNormal3d( 0.0, 1.0,  0.0);
         glVertex3d(-LEN, LEN, -LEN);
         glVertex3d(-LEN, LEN,  LEN);
         glVertex3d( LEN, LEN,  LEN);
         glVertex3d( LEN, LEN, -LEN);

         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);

         glNormal3d(1.0,  0.0,  0.0);
         glVertex3d(LEN, -LEN, -LEN);
         glVertex3d(LEN, -LEN,  LEN);
         glVertex3d(LEN,  LEN,  LEN);
         glVertex3d(LEN,  LEN, -LEN);
      glEnd();
   glPopMatrix();
}
