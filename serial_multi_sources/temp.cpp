#include <cmath>
#include "temp.h"

using namespace std;


// Compute temperature and kinetic energy
void compute_temp_ekin( const double* const vel, const int natoms, 
    const double mass, const double temp_scale, const double ekin_scale, 
    double& temp, double& ekin ) 
{
  double tmp = 0.;
  for ( int i = 0; i < natoms; i++ ) {
    const double vx = vel[ 3 * i + 0 ];
    const double vy = vel[ 3 * i + 1 ];
    const double vz = vel[ 3 * i + 2 ];
    tmp += (vx * vx + vy * vy + vz * vz) * mass; // mass: having it here is more general
  }

  temp = tmp * temp_scale;
  ekin = tmp * ekin_scale * 0.5;

  return;
}


// Rescale to desired temperature
void rescale_temp( double* vel, const int natoms, const double temp_ini, 
    double& temp, double& ekin ) 
{
  const double t_factor = temp_ini / temp;
  const double t_factor_sqrt = sqrt( t_factor );
  
  for ( int i = 0; i < natoms; i++ ) {
    vel[ 3 * i + 0 ] *= t_factor_sqrt;
    vel[ 3 * i + 1 ] *= t_factor_sqrt;
    vel[ 3 * i + 2 ] *= t_factor_sqrt;
  }
  temp *= t_factor;
  ekin *= t_factor;
  
  return;
}
