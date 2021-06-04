#include <iostream>
#include <cmath>
//#include <iomanip> // test only

using namespace std;

// function headers
void setup_struc_vel( const int, const int, const double, const double*, const int, double*, double* );
void get_temp_ekin(double*, const int, const double, const double, const double, double&, double& );
void rescale_temp( double*, const int, const double, double&, double& );
void get_neigh( double*, const int, const double, const double, const int, int*, int* );
double get_epot( double*, const int, const double, const double );
double random( int* ); // this one is taken from Mantevo/miniMD
void print_arr( double*, int );
void print_info ( const double, const double, const int, const double, const double, const double );




int main() {
//cout<<"Hello World!"<<endl;

// Define parameters - using LAMMPS "metal" units convention
//
// Input parameters - might become editable by input
const int box_side = 4; // no of unit cells per dimension
const int nsteps = 200000;
const int nthermo = 1000; // print thermo info every these steps
const double temp_ini = 10.; // K
// pressure unit is bar
//
// Other parameters from here on
//
const double step = 0.001; // ps
//
// Crystal structure for Argon (fcc)
// Note that fcc implies 3D PBC
const int ndims = 3; // no of periodic dimensions // beware: some of the code below implies 3D PBC
const double cellpar = 5.256; // angstrom
const double boxlen = cellpar * box_side;
const int funits = 4;
const int natoms = funits * (int)pow( (double)box_side, (double)ndims );
const double unitpos[ funits * 3 ] = {
  0., 0., 0.,
  0.5*cellpar, 0.5*cellpar, 0.,
  0.5*cellpar, 0., 0.5*cellpar,
  0., 0.5*cellpar, 0.5*cellpar
};
//
// Some physical constants here
const double k_B = 8.617343e-05; // eV/K
const double N_av = 6.02214129e23; // mol-1
const double J_eV = 1.602177e-19; // this is q_e
//
// Model parameters
const double mass = 39.95; // gram/mol (also amu)
const double eps_kB = 117.7; // K
const double eps = eps_kB * k_B; // eV
const double sigma = 3.504; // angstrom
const double cut_fac = 2.5; // adimensional, multiplies sigma
const double skin_fac = 0.3; // adimensional, multiplies sigma
const double cut = cut_fac * sigma;
const double cutskin = cut + skin_fac * sigma;
const double cutsq = cut * cut;
const double cutskinsq = cutskin * cutskin;
const int maxneigh = 150; // with an fcc of side 5.256, and cut+skin of 9.8112, the real maxneigh is 86
//
const double N_dof = ( natoms * 3 - 3 ); // note that this implies 3D PBC (different expressions for lower ndims)
const double mvv2e = 1.036427e-04; // this factor is needed for energy when using metal units
const double temp_scale = mvv2e / ( N_dof * k_B );


// Allocate arrays
int* numneigh = new int [ natoms ];
int* neigh = new int [ natoms * maxneigh ];
double* pos = new double [ natoms * 3 ];
double* vel = new double [ natoms * 3 ];
double* forc = new double [ natoms * 3 ];

// Define simulation variables
double temp, ekin, epot, etot;
int istep = 0;


// Define structure and initialise velocities
setup_struc_vel( funits, box_side, cellpar, unitpos, natoms, pos, vel );

// Rescale to desired temperature
get_temp_ekin( vel, natoms, mass, temp_scale, mvv2e, temp, ekin );
rescale_temp( vel, natoms, temp_ini, temp, ekin );

// Build (full) neighbour list
get_neigh( pos, natoms, boxlen, cutskinsq, maxneigh, numneigh, neigh );

// Compute initial potential energy
//epot 





// Get debug prints
if ( 0 ) { print_arr( pos, natoms ); print_arr( vel, natoms ); }
if ( 1 ) { print_info( cellpar, boxlen, natoms, temp, ekin, epot ); }







// big loop: time evolution

// later on: conditional neighbour update

// PBC check
// compute forces
// integrate
// update velocities

// compute and print output when required







// deallocate arrays
delete [] forc;
delete [] vel;
delete [] pos;
delete [] neigh;
delete [] numneigh;

return 0;
}




// Define structure and initialise velocities
// note that this implies 3D PBC
void setup_struc_vel( const int funits, const int box_side, const double cellpar, const double* unitpos, const int natoms, double* pos, double* vel ) 
{
  const int fd = funits * 3;
  const int bfd = box_side * fd;
  const int bbfd = box_side * bfd;

  int idx;
  int seed;
  double vxtmp = 0.;
  double vytmp = 0.;
  double vztmp = 0.;
  for ( int i = 0; i < box_side; i++ ) {
    for ( int j = 0; j < box_side; j++ ) {
      for ( int k = 0; k < box_side; k++ ) {
        for ( int l = 0; l < funits; l++ ) {
          idx = i * bbfd + j * bfd + k * fd + l * 3;
          // positions
          pos[ idx + 0 ] = cellpar*i + unitpos[ 3 * l + 0 ];
          pos[ idx + 1 ] = cellpar*j + unitpos[ 3 * l + 1 ];
          pos[ idx + 2 ] = cellpar*k + unitpos[ 3 * l + 2 ];
          //cout << left << setw(8) << pos[idx+0] << setw(8) << pos[idx+1] << setw(8) << pos[idx+2] << endl ; // test only
          //cout << idx << endl; // test only
  
          // velocities
          seed = idx;
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
        //cout << endl; // test only
      }
    }
  }
  vxtmp /= natoms;
  vytmp /= natoms;
  vztmp /= natoms;

  // Adjust velocities
  // Zero centre-of-mass motion
  for (int i = 0; i < natoms; i++) {
    vel[ 3 * i + 0 ] -= vxtmp;
    vel[ 3 * i + 1 ] -= vytmp;
    vel[ 3 * i + 2 ] -= vztmp;
  }
  
  return;
}


// Compute temperature and kinetic energy
void get_temp_ekin(double* vel, const int natoms, const double mass, 
                   const double temp_scale, const double ekin_scale, 
                   double& temp, double& ekin) 
{
  double tmp = 0.;
  for (int i = 0; i < natoms; i++) {
    double vx = vel[ 3 * i + 0 ];
    double vy = vel[ 3 * i + 1 ];
    double vz = vel[ 3 * i + 2 ];
    tmp += (vx * vx + vy * vy + vz * vz) * mass; // mass: having it here is more general; for heteroatomic systems, this will become an array
  }
  //cout << "Ave vv : " << tmp / natoms << endl;  // debug

  temp = tmp * temp_scale;
  ekin = tmp * ekin_scale * 0.5;
  return;
}


// Rescale to desired temperature
void rescale_temp( double* vel, const int natoms, const double temp_ini, double& temp, double& ekin ) 
{
  double t_factor, t_factor_sqrt;
  t_factor = temp_ini / temp;
  t_factor_sqrt = sqrt( t_factor );
  
  for (int i = 0; i < natoms; i++) {
    vel[ 3 * i + 0 ] *= t_factor_sqrt;
    vel[ 3 * i + 1 ] *= t_factor_sqrt;
    vel[ 3 * i + 2 ] *= t_factor_sqrt;
  }
  temp *= t_factor;
  ekin *= t_factor;
  
  return;
}


// Build full neighbour list
// note that this implies 3D PBC
// NOTE: in this first implementation, this is never updated; 
// should work at low temperatures, where atoms are likely to stick around their starting positions
void get_neigh( double* pos, const int natoms, const double boxlen, 
                const double cutskinsq, const int maxneigh, 
                int* numneigh, int* neigh ) 
{
  const double boxhalf = boxlen * 0.5;
  for (int i = 0; i < natoms; i++) {
    numneigh[ i ] = 0;
  }

  // int tot_nn = 0;
  // int max_nn = 0;
  // int min_nn = 1000000;
  for (int i = 0; i < natoms; i++) {
    int num_nn = 0;
    for (int j = 0; j < natoms; j++) {
      if ( i == j ) continue;
      double dx = pos[ 3 * i + 0 ] - pos[ 3 * j + 0 ];
      if ( dx > boxhalf ) { dx -= boxlen; }
      if ( dx < - boxhalf ) { dx += boxlen; }
      double dy = pos[ 3 * i + 1 ] - pos[ 3 * j + 1 ];
      if ( dy > boxhalf ) { dy -= boxlen; }
      if ( dy < - boxhalf ) { dy += boxlen; }
      double dz = pos[ 3 * i + 2 ] - pos[ 3 * j + 2 ];
      if ( dz > boxhalf ) { dz -= boxlen; }
      if ( dz < - boxhalf ) { dz += boxlen; }
      double rsq = dx * dx + dy * dy + dz * dz;
      if ( rsq <= cutskinsq ) {
        neigh[ i * maxneigh + num_nn++ ] = j;
      } //else {
        //cout << i << ' ' << j << ' ' << pos[ 3 * i + 0 ] << ' ' << pos[ 3 * i + 1 ] << ' ' << pos[ 3 * i + 2 ] << ' ' << pos[ 3 * j + 0 ] << ' ' << pos[ 3 * j + 1 ] << ' ' << pos[ 3 * j + 2 ] << endl;
      //} // debug
    }
    numneigh[ i ] = num_nn;
    // tot_nn += num_nn;
    // max_nn = max( max_nn, num_nn );
    // min_nn = min( min_nn, num_nn);
  }

  // cout << "tot_nn = " << tot_nn << endl;
  // cout << "max_nn = " << max_nn << endl;
  // cout << "min_nn = " << min_nn << endl;
  return;
}


// This is taken from Mantevo/miniMD
/* Park/Miller RNG w/out MASKING, so as to be like f90s version */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double random(int* idum) 
{
  int k;
  double ans;

  k = (*idum) / IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;

  if(*idum < 0) *idum += IM;

  ans = AM * (*idum);
  return ans;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK


// Generic function to print arrays
void print_arr(double* arr, const int natoms) 
{
  printf("%16c %16c %16c\n", 'X', 'Y', 'Z');
  for ( int i = 0; i < natoms; i++) {
    printf("%+16.6E %+16.6E %+16.6E\n", arr[ 3 * i + 0 ], arr[ 3 * i + 1 ], arr[ 3 * i + 2 ] );
  }

return;
}


// Print info on simulation model
void print_info ( const double cellpar, const double boxlen, const int natoms, const double temp, const double ekin, const double epot ) 
{
  cout << "Cell_par[Ang] : " << cellpar << endl;
  cout << "Box_len[Ang] : " << boxlen << endl;
  cout << "N_atoms : " << natoms << endl;
  cout << "Temp[K] : " << temp << endl;
  cout << "E_kin[eV] : " << ekin << endl;
  cout << "E_pot : " << epot << endl;
  cout << "E_tot : " << ekin + epot << endl;

  return;
}
