#include "force.h"


// Compute forces and potential energy
// Note that this implies 3D PBC
// Test against LAMMPS successful! (forces and accelerations)
void compute_forc_epot( const double* const pos, const int natoms, 
                    const int maxneigh, const int* const numneigh, const int* const neigh, 
                    const double boxlen, const double boxhalf, const double cutsq, 
                    const double sigma6, const double eps, 
                    double* forc, double& epot )
{
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
    double epot_omp = 0.;
  
    # pragma omp parallel for reduction(+:fx) reduction(+:fy) reduction(+:fz) reduction(+:epot_omp)
    for ( int k = 0; k < numneighs; k++ ) {
      const int j = neighs[k];
  
      double dx = x - pos[ 3 * j + 0 ];
      if ( dx > boxhalf ) { dx -= boxlen; }
      if ( dx < - boxhalf ) { dx += boxlen; }
  
      double dy = y - pos[ 3 * j + 1 ];
      if ( dy > boxhalf ) { dy -= boxlen; }
      if ( dy < - boxhalf ) { dy += boxlen; }
  
      double dz = z - pos[ 3 * j + 2 ];
      if ( dz > boxhalf ) { dz -= boxlen; }
      if ( dz < - boxhalf ) { dz += boxlen; }
  
      const double rsq = dx * dx + dy * dy + dz * dz;
      if ( rsq <= cutsq ) {
        const double irsq = 1.0 / rsq;
        const double isr6 = irsq * irsq * irsq * sigma6;
  
        const double force_factor = 48.0 * isr6 * (isr6 - 0.5) * irsq * eps;
        fx += dx * force_factor;
        fy += dy * force_factor;
        fz += dz * force_factor;
        epot_omp += isr6 * (isr6 - 1.0) * eps;
      }
    }
    forc[ 3 * i + 0 ] = fx;
    forc[ 3 * i + 1 ] = fy;
    forc[ 3 * i + 2 ] = fz;
    epot += epot_omp;
  }
  epot *= 2.0; // 4.0 * 0.5 [4.0 from LJ formula, 0.5 to account for double counting of contributes]
  
  return;
}