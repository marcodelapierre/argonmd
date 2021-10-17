#include <cmath>
#include "force.h"


// Compute forces and potential energy
// Note that this implies 3D PBC
// Test against LAMMPS successful! (forces and accelerations)
void compute_forc_epot( const double* const pos, const int natoms, 
    const int maxneigh, const int* const numneigh, const int* const neigh, 
    const double* const boxlen, const double* const boxhalf, 
    const double cutsq, const double sigma6, const double eps, 
    double* forc, double& epot )
{
  const double boxlen0 = boxlen[0];
  const double boxlen1 = boxlen[1];
  const double boxlen2 = boxlen[2];
  const double boxhalf0 = boxhalf[0];
  const double boxhalf1 = boxhalf[1];
  const double boxhalf2 = boxhalf[2];

  epot = 0.;
  for ( int i = 0; i < natoms; i++ ) {
    const int* const neighs = &neigh[ i * maxneigh ];
    const int numneighs = numneigh[ i ];
    const double x = pos[ 3 * i + 0 ];
    const double y = pos[ 3 * i + 1 ];
    const double z = pos[ 3 * i + 2 ];
    double fx = 0.;
    double fy = 0.;
    double fz = 0.;

    for ( int k = 0; k < numneighs; k++ ) {
      const int j = neighs[k];

      double dx = x - pos[ 3 * j + 0 ];
//      if ( dx > boxhalf0 ) { dx -= boxlen0; }
//      if ( dx < - boxhalf0 ) { dx += boxlen0; }
      dx -= floor( ( dx + boxhalf0 ) / boxlen0 ) * boxlen0;

      double dy = y - pos[ 3 * j + 1 ];
//      if ( dy > boxhalf1 ) { dy -= boxlen1; }
//      if ( dy < - boxhalf1 ) { dy += boxlen1; }
      dy -= floor( ( dy + boxhalf1 ) / boxlen1 ) * boxlen1;

      double dz = z - pos[ 3 * j + 2 ];
//      if ( dz > boxhalf2 ) { dz -= boxlen2; }
//      if ( dz < - boxhalf2 ) { dz += boxlen2; }
      dz -= floor( ( dz + boxhalf2 ) / boxlen2 ) * boxlen2;

      const double rsq = dx * dx + dy * dy + dz * dz;
      if ( rsq <= cutsq ) {
        const double irsq = 1.0 / rsq;
        const double isr6 = irsq * irsq * irsq * sigma6;

        const double force_factor = 48.0 * isr6 * (isr6 - 0.5) * irsq * eps;
        fx += dx * force_factor;
        fy += dy * force_factor;
        fz += dz * force_factor;
        epot += isr6 * (isr6 - 1.0) * eps;
      }
    }
    forc[ 3 * i + 0 ] = fx;
    forc[ 3 * i + 1 ] = fy;
    forc[ 3 * i + 2 ] = fz;
  }
  epot *= 2.0; // 4.0 * 0.5 [4.0 from LJ formula, 0.5 to account for double counting of contributes]
  
  return;
}
