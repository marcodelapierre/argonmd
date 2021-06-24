#include "neigh.h"


// Build full neighbour list
// Note that this implies 3D PBC
void compute_neigh( const double* const pos, const int natoms, 
                const double boxlen, const double boxhalf, const double cutskinsq, 
                const int maxneigh, int* numneigh, int* neigh ) 
{
  for ( int i = 0; i < natoms; i++ ) {
    numneigh[ i ] = 0;
  }

  for ( int i = 0; i < natoms; i++ ) {
    int num_nn = 0;
    for ( int j = 0; j < natoms; j++ ) {
      if ( i == j ) continue;

      double dx = pos[ 3 * i + 0 ] - pos[ 3 * j + 0 ];
      if ( dx > boxhalf )   { dx -= boxlen; }
      if ( dx < - boxhalf ) { dx += boxlen; }

      double dy = pos[ 3 * i + 1 ] - pos[ 3 * j + 1 ];
      if ( dy > boxhalf )   { dy -= boxlen; }
      if ( dy < - boxhalf ) { dy += boxlen; }

      double dz = pos[ 3 * i + 2 ] - pos[ 3 * j + 2 ];
      if ( dz > boxhalf )   { dz -= boxlen; }
      if ( dz < - boxhalf ) { dz += boxlen; }

      double rsq = dx * dx + dy * dy + dz * dz;
      if ( rsq <= cutskinsq ) {
        neigh[ i * maxneigh + num_nn++ ] = j;
      }
    }
    numneigh[ i ] = num_nn;
  }

  return;
}
