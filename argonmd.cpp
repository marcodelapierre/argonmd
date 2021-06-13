#include <iostream>
#include <cmath>
#include <ctime>
//#include <iomanip> // test only

using namespace std;

// function headers
void setup_struc_vel( const int, const int, const double, const double*, const int, double*, double*, double* );
void get_temp_ekin( const double* const, const int, const double, const double, const double, double&, double& );
void rescale_temp( double*, const int, const double, double&, double& );
void check_pbc( double*, const int, const double );
void get_neigh( const double* const, const int, const double, const double, const int, int*, int* );
void get_forc_epot( const double* const, const int, const int, const int* const, const int* const, 
                    const double, const double, const double, const double, double*, double& );
//
double random( int* ); // this one is taken from Mantevo/miniMD
void print_arr( const double* const, const int, const int );
void print_info ( const double, const double, const int, const double, const double, const double );




int main() {
//cout<<"Hello World!"<<endl;

// Define parameters - using LAMMPS "metal" physical units convention
// pressure unit is bar
// force unit is eV/Ang
//
// Input parameters - might become editable by input
const int box_units = 4; // no of unit cells per dimension in the simulation box
const int nsteps = 200000;
const double temp_ini = 10.; // K [117.7: datum from LAMMPS LJ example]
const int nneighupd = 20; // update neighbour list every these steps [from LAMMPS LJ example]
const int nthermo = 1000; // print thermo info every these steps
const int ndump = 1000; // dump structure every these steps
//
// Crystal structure for Argon (fcc)
// Note that fcc implies 3D PBC
const int funits = 4;
const int natoms = funits * box_units * box_units * box_units; // note that this implies 3D PBC
const double cellpar = 5.795; // angstrom [datum from LAMMPS LJ example] [5.256: from real data]
const double boxlen = cellpar * box_units;
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
const double step = 0.001; // ps
const double mass = 39.95; // gram/mol (also amu)
const double eps_kB = 117.7; // K
const double eps = eps_kB * k_B; // eV
const double sigma = 3.504; // angstrom
const double sigma2 = sigma * sigma;
const double sigma6 = sigma2 * sigma2 * sigma2;
const double cut_fac = 2.5; // adimensional, multiplies sigma
const double skin_fac = 0.3; // adimensional, multiplies sigma
const double cut = cut_fac * sigma;
const double cutskin = cut + skin_fac * sigma;
const double cutsq = cut * cut;
const double cutskinsq = cutskin * cutskin;
const int maxneigh = 150; // with an fcc of side 5.256, and cut+skin of 9.8112, the real maxneigh is 86
//
const double N_dof = ( natoms * 3 - 3 ); // note that this implies 3D PBC (different expressions for lower dimensionalities)
const double ekin_scale = 1.036427e-04; // this factor is needed when using metal units ("mvv2e" in Mantevo/miniMD)
const double temp_scale = ekin_scale / ( N_dof * k_B );


// Allocate arrays
int* numneigh = new int [ natoms ];
int* neigh = new int [ natoms * maxneigh ];
double* pos = new double [ natoms * 3 ];
double* posraw = new double [ natoms * 3 ];
double* vel = new double [ natoms * 3 ];
double* forc = new double [ natoms * 3 ];

// Define simulation variables
double temp, ekin, epot;
int istep = 0;


// Define structure and initialise velocities
setup_struc_vel( funits, box_units, cellpar, unitpos, natoms, pos, posraw, vel );
// Rescale to desired temperature
get_temp_ekin( vel, natoms, mass, temp_scale, ekin_scale, temp, ekin );
rescale_temp( vel, natoms, temp_ini, temp, ekin );

// PBC check // not needed at startup with current input structure, yet here for generality
check_pbc( pos, natoms, boxlen );
// Build (full) neighbour list
get_neigh( pos, natoms, boxlen, cutskinsq, maxneigh, numneigh, neigh );

// Compute initial forces
get_forc_epot( pos, natoms, maxneigh, numneigh, neigh, 
               boxlen, cutsq, sigma6, eps, forc, epot );



// Get debug prints
if ( 1 ) { print_arr( pos, 0, natoms ); print_arr( vel, 0, natoms ); }
if ( 1 ) { print_info( cellpar, boxlen, natoms, temp, ekin, epot ); }



// big loop: time evolution #3



// integrate N-body system #3

// PBC check
//check_pbc( pos, natoms, boxlen );
// Update (full) neighbour list
//if( ( istep > 0 ) && ( istep%nneighupd == 0 ) ) { get_neigh( pos, natoms, boxlen, cutskinsq, maxneigh, numneigh, neigh ); }






// compute and print output when required #4
// add timer for time loop #5
//clock_t start, watch;
//start = clock();
//watch = clock() - start;
// to print time: ((float)watch)/CLOCKS_PER_SEC

// dump xyz (optional) #6





// deallocate arrays
delete [] forc;
delete [] vel;
delete [] posraw;
delete [] pos;
delete [] neigh;
delete [] numneigh;

return 0;
}




// Define structure and initialise velocities
// note that this implies 3D PBC
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
          //cout << left << setw(8) << pos[idx+0] << setw(8) << pos[idx+1] << setw(8) << pos[idx+2] << endl ; // test only
          //cout << idx << endl; // test only
  
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
        //cout << endl; // test only
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


// Compute temperature and kinetic energy
void get_temp_ekin( const double* const vel, const int natoms, const double mass, 
                    const double temp_scale, const double ekin_scale, 
                    double& temp, double& ekin ) 
{
  double tmp = 0.;
  for ( int i = 0; i < natoms; i++ ) {
    const double vx = vel[ 3 * i + 0 ];
    const double vy = vel[ 3 * i + 1 ];
    const double vz = vel[ 3 * i + 2 ];
    tmp += (vx * vx + vy * vy + vz * vz) * mass; // mass: having it here is more general
  }
  //cout << "Ave vv : " << tmp / natoms << endl;  // debug

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


// Build full neighbour list
// note that this implies 3D PBC
void get_neigh( const double* const pos, const int natoms, const double boxlen, 
                const double cutskinsq, const int maxneigh, 
                int* numneigh, int* neigh ) 
{
  const double boxhalf = boxlen * 0.5;
  for ( int i = 0; i < natoms; i++ ) {
    numneigh[ i ] = 0;
  }

  // int tot_nn = 0;
  // int max_nn = 0;
  // int min_nn = 1000000;
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


// Check periodic boundary conditions
// note that this implies 3D PBC
void check_pbc( double* pos, const int natoms, const double boxlen ) 
{
  for ( int i = 0; i < natoms; i++ ) {
    double x = pos[ 3 * i + 0 ];
    if ( x >= boxlen ) { x -= boxlen; pos[ 3 * i + 0 ] = x; }
    if ( x < 0. )      { x += boxlen; pos[ 3 * i + 0 ] = x; }
  
    double y = pos[ 3 * i + 1 ];
    if ( y >= boxlen ) { y -= boxlen; pos[ 3 * i + 1 ] = y; }
    if ( y < 0. )      { y += boxlen; pos[ 3 * i + 1 ] = y; }
  
    double z = pos[ 3 * i + 2 ];
    if ( z >= boxlen ) { z -= boxlen; pos[ 3 * i + 2 ] = z; }
    if ( z < 0. )      { z += boxlen; pos[ 3 * i + 2 ] = z; }
  }

  return;
}


// Compute forces and potential energy
// note that this implies 3D PBC
void get_forc_epot( const double* const pos, const int natoms, 
                    const int maxneigh, const int* const numneigh, const int* const neigh, 
                    const double boxlen, const double cutsq, 
                    const double sigma6, const double eps, 
                    double* forc, double& epot )
{
  const double boxhalf = boxlen * 0.5;
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
  
        const double force_fac = 48.0 * isr6 * (isr6 - 0.5) * irsq * eps;
        fx += dx * force_fac;
        fy += dy * force_fac;
        fz += dz * force_fac;
        epot += isr6 * (isr6 - 1.0) * eps;
      }
    }
    forc[ 3 * i + 0 ] = fx;
    forc[ 3 * i + 1 ] = fy;
    forc[ 3 * i + 2 ] = fz;
  }
  epot *= 2.0; // 4.0 * 0.5 [4.0 from LJ formula, 0.5 to account for double counting of contributes]
  // normalise by natoms to compare with LAMMPS
  
  return;
}


// Random number generator
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
void print_arr( const double* const arr, const int istart, const int istop ) 
{
  printf("%16c %16c %16c\n", 'X', 'Y', 'Z');
  for ( int i = istart; i < istop; i++) {
    printf("%+16.6E %+16.6E %+16.6E\n", arr[ 3 * i + 0 ], arr[ 3 * i + 1 ], arr[ 3 * i + 2 ] );
  }

  return;
}


// Print information on simulation model
void print_info ( const double cellpar, const double boxlen, 
                  const int natoms, const double temp, 
                  const double ekin, const double epot ) 
{
  cout << "Cell_par[Ang] : " << cellpar << endl;
  cout << "Box_len[Ang] : " << boxlen << endl;
  cout << "N_atoms : " << natoms << endl;
  cout << "Temp[K] : " << temp << endl;
  cout << "E_kin[eV] : " << ekin << endl;
  cout << "E_pot[eV] : " << epot << endl;
  cout << "E_tot[eV] : " << ekin + epot << endl;

  return;
}
