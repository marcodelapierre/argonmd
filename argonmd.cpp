#include <iostream>
#include <cmath>
//#include <iomanip> // test only

using namespace std;

// function headers
void get_temp_ekin(double*, const int, const double, const double, const double, double&, double& );
double get_epot( double*, const int, const double, const double );
double random( int* ); // This one is taken from Mantevo/miniMD
void print_arr( double*, int );


int main() {
//cout<<"Hello World!"<<endl;

// define parameters - using LAMMPS "metal" convention
//
// first two might become editable by cli arguments
const int box_side = 2; // no of unit cells per dimension
const int nsteps = 200000;
const int nthermo = 1000; // print thermo info every these steps
const double step = 0.001; // ps
const double temp_ini = 10.; // K
// pressure unit is bar
//
// argon crystal structure (fcc)
const int ndims = 3; // no of spatial dimensions // don't change this, some of the code below implies a value of 3
const double cellpar = 5.256; // angstrom
const int funits = 4;
// derived structural parameters
const int fd = funits * ndims;
const int bfd = box_side * fd;
const int bbfd = box_side * bfd;
const double unitpos[ fd ] = {
  0., 0., 0.,
  0.5*cellpar, 0.5*cellpar, 0.,
  0.5*cellpar, 0., 0.5*cellpar,
  0., 0.5*cellpar, 0.5*cellpar
};
const int natoms = funits * box_side * box_side * box_side;
//
// some physical constants here
const double k_B = 8.617343e-05; // eV/K
const double N_av = 6.02214129e23; // mol-1
const double J_eV = 1.602177e-19; // this is q_e
//
// model parameters
const double mass = 39.95; // gram/mol (also amu)
const double eps_kB = 117.7; // K
const double eps = eps_kB * k_B; // eV
const double sigma = 3.504; // angstrom
const double cut_fac = 2.5; // adimensional, multiplies sigma
const double skin_fac = 0.3; // adimensional, multiplies sigma
//
const double N_dof = ( natoms * 3 - 3 );
const double mvv2e = 1.036427e-04; // this factor is needed for energy when using metal units
const double temp_scale = mvv2e / ( N_dof * k_B );

// allocate arrays
//
double* pos = new double [ natoms * ndims ];
double* vel = new double [ natoms * ndims ];
double* force = new double [ natoms * ndims ];

// define structure AND 
// initialise velocities
int idx;
int seed;
double vxtmp = 0.;
double vytmp = 0.;
double vztmp = 0.;
for ( int i = 0; i < box_side; i++ ) {
  for ( int j = 0; j < box_side; j++ ) {
    for ( int k = 0; k < box_side; k++ ) {
      for ( int l = 0; l < funits; l++ ) {
        idx = i * bbfd + j * bfd + k * fd + l * ndims;
        // positions
        pos[ idx + 0 ] = cellpar*i + unitpos[ l * ndims + 0 ];
        pos[ idx + 1 ] = cellpar*j + unitpos[ l * ndims + 1 ];
        pos[ idx + 2 ] = cellpar*k + unitpos[ l * ndims + 2 ];
        //cout << left << setw(8) << pos[idx+0] << setw(8) << pos[idx+1] << setw(8) << pos[idx+2] << endl ; // test only
        //cout << idx << endl; // test only

        // velocities
        seed = idx;
        for(int m = 0; m < 5; m++) random(&seed);
        vel[ idx + 0 ] = random(&seed);
        for(int m = 0; m < 5; m++) random(&seed);
        vel[ idx + 1 ] = random(&seed);
        for(int m = 0; m < 5; m++) random(&seed);
        vel[ idx + 2 ] = random(&seed);

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

// adjust velocities
// zero centre-of-mass motion
for (int i =0; i < natoms; i++) {
  vel[ i * ndims + 0 ] -= vxtmp;
  vel[ i * ndims + 1 ] -= vytmp;
  vel[ i * ndims + 2 ] -= vztmp;
}
// rescale to desired temperature
double temp, t_factor;
double ekin, epot, etot;
// set debug prints
int debug_arr = 0;
int debug_info = 1;
get_temp_ekin(vel, natoms, mass, temp_scale, mvv2e, temp, ekin);
t_factor = sqrt( temp_ini / temp );
for (int i =0; i < natoms; i++) {
  vel[ i * ndims + 0 ] *= t_factor;
  vel[ i * ndims + 1 ] *= t_factor;
  vel[ i * ndims + 2 ] *= t_factor;
}
ekin *= temp_ini / temp; // order matters for these two: ekin, then temp
temp *= temp_ini / temp; // order matters for these two: ekin, then temp

// build neighbour list
// NOTE: in this first implementation, this is never updated; 
// should work at low temperatures, where atoms are likely to stick around their starting positions

//epot = 
//etot = ekin + epot;




// get debug prints
if ( debug_arr ) {
  print_arr( pos, natoms);
  print_arr( vel, natoms);
}
if ( debug_info ) {
  cout << "Cell_par[Ang] : " << cellpar << endl;
  cout << "Box_side[Ang] : " << box_side * cellpar << endl;
  cout << "N_atoms : " << natoms << endl;
  cout << "Temp[K] : " << temp << endl;
  cout << "E_kin[eV] : " << ekin << endl;
  //cout << "E_pot : " << epot << endl;
  //cout << "E_tot : " << etot << endl;
}


// big loop: time evolution


// compute forces

// integrate

// update velocities

// compute and print output when required





// deallocate arrays
delete [] force;
delete [] vel;
delete [] pos;

return 0;
}




void get_temp_ekin(double* vel, const int natoms, const double mass, 
                   const double temp_scale, const double ekin_scale, 
                   double& temp, double& ekin)
{
  double tmp = 0.;
  for (int i =0; i < natoms; i++) {
    double vx = vel[ 3*i + 0 ];
    double vy = vel[ 3*i + 1 ];
    double vz = vel[ 3*i + 2 ];
    tmp += (vx * vx + vy * vy + vz * vz) * mass; // mass: having it here is more general; for heteroatomic systems, this will become an array
  }
  //cout << "Ave vv : " << tmp / natoms << endl;  // debug

  temp = tmp * temp_scale;
  ekin = tmp * ekin_scale * 0.5;
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


void print_arr(double* arr, const int natoms)
{
  printf("%16c %16c %16c\n", 'X', 'Y', 'Z');
  for ( int i = 0; i < natoms; i++) {
    printf("%+16.6E %+16.6E %+16.6E\n", arr[ 3*i+0 ], arr[ 3*i+1 ], arr[ 3*i+2 ] );
  }

return;
}
