#include "setup.h"
#include "random.h"


// Define structure and initialise velocities
// Note that this implies 3D PBC
void setup_struc_vel( const int funits, const int box_units, 
                      const double cellpar, const double* unitpos, 
                      const int natoms, double* pos, double* posraw, double* vel ) 
{
  const int fd = funits * 3;
  const int bfd = box_units * fd;
  const int bbfd = box_units * bfd;

  double vxtmp = 0.;
  double vytmp = 0.;
  double vztmp = 0.;
  for ( int i = 0; i < box_units; i++ ) {
    for ( int j = 0; j < box_units; j++ ) {
      for ( int k = 0; k < box_units; k++ ) {
        for ( int l = 0; l < funits; l++ ) {
          const int idx = i * bbfd + j * bfd + k * fd + l * 3;
          // positions
          pos[ idx + 0 ] = cellpar*i + unitpos[ 3 * l + 0 ];
          pos[ idx + 1 ] = cellpar*j + unitpos[ 3 * l + 1 ];
          pos[ idx + 2 ] = cellpar*k + unitpos[ 3 * l + 2 ];
          posraw[ idx + 0 ] = cellpar*i + unitpos[ 3 * l + 0 ];
          posraw[ idx + 1 ] = cellpar*j + unitpos[ 3 * l + 1 ];
          posraw[ idx + 2 ] = cellpar*k + unitpos[ 3 * l + 2 ];
  
          // velocities
          int seed = idx;
          for ( int m = 0; m < 5; m++ ) random( &seed );
          vel[ idx + 0 ] = random( &seed );
          for ( int m = 0; m < 5; m++ ) random( &seed );
          vel[ idx + 1 ] = random( &seed );
          for ( int m = 0; m < 5; m++ ) random( &seed );
          vel[ idx + 2 ] = random( &seed );
  
          vxtmp += vel[ idx + 0 ];
          vytmp += vel[ idx + 1 ];
          vztmp += vel[ idx + 2 ];
        }
      }
    }
  }
  vxtmp /= natoms;
  vytmp /= natoms;
  vztmp /= natoms;

  // Zero centre-of-mass motion
  for ( int i = 0; i < natoms; i++ ) {
    vel[ 3 * i + 0 ] -= vxtmp;
    vel[ 3 * i + 1 ] -= vytmp;
    vel[ 3 * i + 2 ] -= vztmp;
  }
  
  return;
}
