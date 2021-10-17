#include <cmath>
#include "pos_vel.h"


// Check periodic boundary conditions
// Note that this implies 3D PBC
void check_pbc( double* pos, const int natoms, const double* const boxlen ) 
{
  const double boxlen0 = boxlen[0];
  const double boxlen1 = boxlen[1];
  const double boxlen2 = boxlen[2];

  for ( int i = 0; i < natoms; i++ ) {
    double x = pos[ 3 * i + 0 ];
    if ( x >= boxlen0 ) { x -= boxlen0; pos[ 3 * i + 0 ] = x; }
    if ( x < 0. )      { x += boxlen0; pos[ 3 * i + 0 ] = x; }
  
    double y = pos[ 3 * i + 1 ];
    if ( y >= boxlen1 ) { y -= boxlen1; pos[ 3 * i + 1 ] = y; }
    if ( y < 0. )      { y += boxlen1; pos[ 3 * i + 1 ] = y; }
  
    double z = pos[ 3 * i + 2 ];
    if ( z >= boxlen2 ) { z -= boxlen2; pos[ 3 * i + 2 ] = z; }
    if ( z < 0. )      { z += boxlen2; pos[ 3 * i + 2 ] = z; }
  }

  return;
}


// Update positions and check PBC meanwhile
void update_pos_pbc( double* pos, double* posraw, 
    const double* const vel, const double* const forc, const int natoms, 
    const double dt, const double forc_hdtsq_scale, 
    const double imass, const double* const boxlen ) 
{
  const double boxlen0 = boxlen[0];
  const double boxlen1 = boxlen[1];
  const double boxlen2 = boxlen[2];

  for ( int i = 0; i < natoms; i++ ) {
    double dx = vel[ 3 * i + 0 ] * dt + forc[ 3 * i + 0 ] * forc_hdtsq_scale * imass;
    double x = pos[ 3 * i + 0 ] + dx;
//    if ( x >= boxlen0 ) { x -= boxlen0; }
//    if ( x < 0. )      { x += boxlen0; }
    x -= floor( x / boxlen0 ) * boxlen0;
    pos[ 3 * i + 0 ] = x;
    posraw[ 3 * i + 0 ] += dx;
  
    double dy = vel[ 3 * i + 1 ] * dt + forc[ 3 * i + 1 ] * forc_hdtsq_scale * imass;
    double y = pos[ 3 * i + 1 ] + dy;
//    if ( y >= boxlen1 ) { y -= boxlen1; }
//    if ( y < 0. )      { y += boxlen1; }
    y -= floor( y / boxlen1 ) * boxlen1;
    pos[ 3 * i + 1 ] = y;
    posraw[ 3 * i + 1 ] += dy;
  
    double dz = vel[ 3 * i + 2 ] * dt + forc[ 3 * i + 2 ] * forc_hdtsq_scale * imass;
    double z = pos[ 3 * i + 2 ] + dz;
//    if ( z >= boxlen2 ) { z -= boxlen2; }
//    if ( z < 0. )      { z += boxlen2; }
    z -= floor( z / boxlen2 ) * boxlen2;
    pos[ 3 * i + 2 ] = z;
    posraw[ 3 * i + 2 ] += dz;
  }

  return;
}


// Update velocities
// Note that this implies Velocity Verlet integrator
void update_vel( double* vel, 
    const double* const forcold, const double* const forc, 
    const int natoms, const double forc_hdt_scale, const double imass ) 
{
  for ( int i = 0; i < natoms; i++ ) {
    vel[ 3 * i + 0 ] += ( forcold[ 3 * i + 0 ] + forc[ 3 * i + 0 ] ) * forc_hdt_scale * imass;
    vel[ 3 * i + 1 ] += ( forcold[ 3 * i + 1 ] + forc[ 3 * i + 1 ] ) * forc_hdt_scale * imass;
    vel[ 3 * i + 2 ] += ( forcold[ 3 * i + 2 ] + forc[ 3 * i + 2 ] ) * forc_hdt_scale * imass;
  }

  return;
}
