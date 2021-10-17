#include <cmath>
#include "neigh.h"


// Build full neighbour list
// Note that this implies 3D PBC
void compute_neigh( const double* const pos, const int natoms, 
    const double* const boxlen, const double* const boxhalf, 
    const double cutskinsq, const int maxneigh, 
    int* numneigh, int* neigh ) 
{
  const double boxlen0 = boxlen[0];
  const double boxlen1 = boxlen[1];
  const double boxlen2 = boxlen[2];
  const double boxhalf0 = boxhalf[0];
  const double boxhalf1 = boxhalf[1];
  const double boxhalf2 = boxhalf[2];

  for ( int i = 0; i < natoms; i++ ) {
    numneigh[ i ] = 0;
  }

  for ( int i = 0; i < natoms; i++ ) {
    int num_nn = 0;
    for ( int j = 0; j < natoms; j++ ) {
      if ( i == j ) continue;

      double dx = pos[ 3 * i + 0 ] - pos[ 3 * j + 0 ];
//      if ( dx > boxhalf0 )   { dx -= boxlen0; }
//      if ( dx < - boxhalf0 ) { dx += boxlen0; }
      dx -= floor( ( dx + boxhalf0 ) / boxlen0 ) * boxlen0;

      double dy = pos[ 3 * i + 1 ] - pos[ 3 * j + 1 ];
//      if ( dy > boxhalf1 )   { dy -= boxlen1; }
//      if ( dy < - boxhalf1 ) { dy += boxlen1; }
      dy -= floor( ( dy + boxhalf1 ) / boxlen1 ) * boxlen1;

      double dz = pos[ 3 * i + 2 ] - pos[ 3 * j + 2 ];
//      if ( dz > boxhalf2 )   { dz -= boxlen2; }
//      if ( dz < - boxhalf2 ) { dz += boxlen2; }
      dz -= floor( ( dz + boxhalf2 ) / boxlen2 ) * boxlen2;

      double rsq = dx * dx + dy * dy + dz * dz;
      if ( rsq <= cutskinsq ) {
        neigh[ i * maxneigh + num_nn++ ] = j;
      }
    }
    numneigh[ i ] = num_nn;
  }

  return;
}
