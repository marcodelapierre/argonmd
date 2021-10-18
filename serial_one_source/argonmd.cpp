#include <iostream>
#include <cstring>
#include <cmath>
#include <ctime>

using namespace std;

// struct definitions
struct InputParams {
  int box_units[ 3 ];
  int nsteps;
  double temp_ini;
  int nneighupd;
  int nthermo;
  int ndump;
};

// function headers
InputParams get_input_params( const int, char** );
//
void setup_struc_vel( const int, const int* const, const double, const double* const, const int, double*, double*, double* );
//
void compute_temp_ekin( const double* const, const int, const double, const double, const double, double&, double& );
void rescale_temp( double*, const int, const double, double&, double& );
//
void compute_neigh_verlet( const double* const, const int, const double* const, const double* const, const double, const int, int*, int* );
//
void compute_forc_epot( const double* const, const int, const int, const int* const, const int* const, 
    const double* const, const double* const, const double, const double, const double, double*, double& );
//
void check_pbc( double*, const int, const double* const );
void update_pos_pbc( double*, double*, const double* const, const double* const, const int, 
    const double, const double, const double, const double* const );
void update_vel( double*, const double* const, const double* const, const int, const double, const double );
//
double random( int* ); // this one is taken from Mantevo/miniMD
//
void print_arr( const double* const, const int, const int );
void print_info( const int* const, const int, const double, const int, const int, const int, 
    const double, const double, const double, const double* const, const int );
void print_thermo( const int, const double, const double, const double, const double, const double, const double );
void dump_pdb( FILE*, const int, const double* const, const double, const char* const, const double* const, const int );



int main( int argc, char** argv ) {
//cout<<"Ciao Mondo!"<<endl;

// Define parameters - using LAMMPS "metal" physical units convention
// pressure unit is bar
// force unit is eV/Ang
//
// Input parameters - editable by input
InputParams input_params = get_input_params( argc, argv );
const int box_units[ 3 ] = {
  input_params.box_units[0],
  input_params.box_units[1],
  input_params.box_units[2]
}; // no of unit cells per each dimension in the simulation box
const int nsteps = input_params.nsteps; // no of time steps in the simulation
const double temp_ini = input_params.temp_ini; // K [117.7: datum from LAMMPS LJ example]
const int nneighupd = input_params.nneighupd; // update neighbour list every these steps [from LAMMPS LJ example]
const int nthermo = input_params.nthermo; // print thermo info every these steps
const int ndump = input_params.ndump; // dump structure every these steps
//
const char* coorfile = "coord.pdb"; // filename for initial atomic coordinates
const char* trajfile = "traj.pdb"; // filename for trajectory atomic coordinates
//
// Crystal structure for Argon (fcc)
// Note that fcc implies 3D PBC
const int funits = 4;
const int natoms = funits * box_units[0] * box_units[1] * box_units[2]; // note that this implies 3D PBC // affected by input parameters
const double cellpar = 5.795; // angstrom [datum from LAMMPS LJ example] [5.256: from real data]
const double cellang = 90.0; // degrees
const double boxlen[ 3 ] = {
  cellpar * box_units[0],
  cellpar * box_units[1],
  cellpar * box_units[2]
}; // affected by input parameters
const double boxhalf[ 3 ] = {
  0.5 * boxlen[0],
  0.5 * boxlen[1],
  0.5 * boxlen[2]
}; // affected by input parameters
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
const double dt = 0.001; // ps
const double hdt = 0.5 * dt;
const double hdtsq = 0.5 * dt * dt;
const double mass = 39.95; // gram/mol (also amu)
const double imass = 1.0 / mass;
const int atno = 18; // atomic number
const char* elsym = "Ar"; // element symbol
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
const double N_dof = ( natoms * 3 - 3 ); // note that this implies 3D PBC (different expressions for lower dimensionalities) // affected by input parameters
const double ekin_scale = 10.0 / N_av / J_eV; // 1.036427e-04; // this factor is needed when using metal units ("mvv2e" in Mantevo/miniMD) [from my notes]
const double temp_scale = ekin_scale / ( N_dof * k_B ); // [from my notes] // affected by input parameters
const double forc_scale = 1.0 / ekin_scale; // J_eV * N_av * 0.1; // 9648.536; // this factor is needed when using metal units [from my notes]
const double forc_hdt_scale = forc_scale * hdt;
const double forc_hdtsq_scale = forc_scale * hdtsq;


// Allocate arrays
int* numneigh = new int [ natoms ];
int* neigh = new int [ natoms * maxneigh ];
double* pos = new double [ natoms * 3 ];
double* posraw = new double [ natoms * 3 ];
double* vel = new double [ natoms * 3 ];
double* forc = new double [ natoms * 3 ];
double* forcold = new double [ natoms * 3 ];

// Define variables and pointers
double temp, ekin, epot, clocktime;
int istep = 0;
double* forctmp;
clock_t start, watch;
FILE* coor; // file for initial atomic coordinates
FILE* traj; // file for trajectory atomic coordinates

// Print program header
printf( "\n** ArgonMD **\n" );


// Define structure and initialise velocities
setup_struc_vel( funits, box_units, cellpar, unitpos, natoms, pos, posraw, vel );
// Rescale to desired temperature
compute_temp_ekin( vel, natoms, mass, temp_scale, ekin_scale, temp, ekin );
rescale_temp( vel, natoms, temp_ini, temp, ekin );

// PBC check // not needed at startup with current input structure, yet here for generality
check_pbc( pos, natoms, boxlen );
// Build (full) neighbour list
compute_neigh_verlet( pos, natoms, boxlen, boxhalf, cutskinsq, maxneigh, numneigh, neigh );

// Compute initial forces
compute_forc_epot( pos, natoms, maxneigh, numneigh, neigh, 
                   boxlen, boxhalf, cutsq, sigma6, eps, forc, epot );

// Print simulation info
print_info( box_units, nsteps, temp_ini, nneighupd, nthermo, ndump, dt, cut, cellpar, boxlen, natoms );
//if ( 0 ) { print_arr( pos, 0, natoms ); print_arr( vel, 0, natoms ); } // debug print
//
// Print initial thermo output
printf( "\n %9s  %10s  %8s  %12s  %12s  %12s  %10s\n", "Step", "Time[ps]", "Temp[K]", "Ekin[eV]", "Epot[eV]", "Etot[eV]", "Clock[s]" );
print_thermo( istep, dt*istep, temp, ekin/natoms, epot/natoms, (ekin+epot)/natoms, 0.0 );

// Dump initial atomic coordinates
if ( ndump > 0 ) {
  coor = fopen( coorfile, "w" );
  dump_pdb( coor, istep, boxlen, cellang, elsym, pos, natoms );
  fclose( coor );
  traj = fopen( trajfile, "w" );
}


// Time evolution loop
start = clock();
for (istep = 1; istep <= nsteps; istep++) {

// This code block could become a routine "integrate"; leaving it here for algorithm readability
// Note that this implies Velocity Verlet integrator
  {
// Update positions and check PBC meanwhile
  update_pos_pbc( pos, posraw, vel, forc, natoms, dt, forc_hdtsq_scale, imass, boxlen );

// Update (full) neighbour list
  if( nneighupd > 0 && istep%nneighupd == 0 ) { 
    compute_neigh_verlet( pos, natoms, boxlen, boxhalf, cutskinsq, maxneigh, numneigh, neigh );
  }

// Store old forces and compute new forces
  forctmp = forcold;
  forcold = forc;
  forc = forctmp;
  compute_forc_epot( pos, natoms, maxneigh, numneigh, neigh, 
                     boxlen, boxhalf, cutsq, sigma6, eps, forc, epot );

// Update velocities
  update_vel( vel, forcold, forc, natoms, forc_hdt_scale, imass );
  }


  if ( nthermo > 0 && istep%nthermo == 0 ) {
// Compute temperature when required
    compute_temp_ekin( vel, natoms, mass, temp_scale, ekin_scale, temp, ekin );

// Get clock time when required
    watch = clock() - start;
    clocktime = ((float)watch)/CLOCKS_PER_SEC;

// Print thermo output when required
    print_thermo( istep, dt*istep, temp, ekin/natoms, epot/natoms, (ekin+epot)/natoms, clocktime );
  }

// Dump atomic coordinates when required
// Note that PDB files are large, so enable dumping only for demonstrations; ideally these should go in a binary format (eg DCD)
// Note also that ideally "posraw" should be used for positions; "pos" is used here instead, for convenience when demonstrating with VMD
  if ( ndump > 0 && istep%ndump == 0 ) {
    dump_pdb( traj, istep, boxlen, cellang, elsym, pos, natoms );
  }

}

// Get and print final clocktime
watch = clock() - start;
clocktime = ((float)watch)/CLOCKS_PER_SEC;
printf( "\n Loop Clock Time [s] : %10.3F\n" , clocktime );


// Close trajectory file if needed
if ( ndump > 0 ) { fclose( traj ); }

// Deallocate arrays
delete [] forcold;
delete [] forc;
delete [] vel;
delete [] posraw;
delete [] pos;
delete [] neigh;
delete [] numneigh;

return 0;
}




// Get input parameters
InputParams get_input_params( const int argc, char** argv ) 
{
  InputParams input_params;

  if ( argc > 1 ) {
// allowing requests for orthorhombic boxes
    char input_string[ 63 ];
    char* token_string;
    int count = 0;
    int box_units[ 3 ] = { 0, 0, 0 };
    strcpy ( input_string, argv[1] );
    token_string = strtok( input_string, "," );
    while ( token_string != NULL ) {
      box_units[count++] = atoi( token_string );
      token_string = strtok( NULL, "," );
    }
    if ( count == 1 ) {
      input_params.box_units[0] = box_units[0];
      input_params.box_units[1] = box_units[0];
      input_params.box_units[2] = box_units[0];
    } else if ( count == 2 ) {
      input_params.box_units[0] = box_units[0];
      input_params.box_units[1] = box_units[0];
      input_params.box_units[2] = box_units[1];
    } else if ( count == 3 ) {
      input_params.box_units[0] = box_units[0];
      input_params.box_units[1] = box_units[1];
      input_params.box_units[2] = box_units[2];
    } else {
      input_params.box_units[0] = 5;
      input_params.box_units[1] = 5;
      input_params.box_units[2] = 5;
    }
  } else {
    input_params.box_units[0] = 5;
    input_params.box_units[1] = 5;
    input_params.box_units[2] = 5;
  }

  if ( argc > 2 ) {
    input_params.nsteps = atoi( argv[2] );
  } else {
    input_params.nsteps = 10000;
  }
  if ( argc > 3 ) {
    input_params.temp_ini = atof( argv[3] );
  } else {
    input_params.temp_ini = 10.;
  }
  if ( argc > 4 ) {
    input_params.nneighupd = atoi( argv[4] );
  } else {
    input_params.nneighupd = 20;
  }
  if ( argc > 5 ) {
    input_params.nthermo = atoi( argv[5] );
  } else {
    input_params.nthermo = 1000;
  }
  if ( argc > 6 ) {
    input_params.ndump = atoi( argv[6] );
  } else {
    input_params.ndump = 0;
  }

  return input_params;
}


// Define structure and initialise velocities
// Note that this implies 3D PBC
void setup_struc_vel( const int funits, const int* const box_units, 
    const double cellpar, const double* const unitpos, const int natoms, 
    double* pos, double* posraw, double* vel ) 
{
  const int fd = funits * 3;
  const int bfd = box_units[2] * fd;
  const int bbfd = box_units[1] * bfd;

  double vxtmp = 0.;
  double vytmp = 0.;
  double vztmp = 0.;
  for ( int i = 0; i < box_units[0]; i++ ) {
    for ( int j = 0; j < box_units[1]; j++ ) {
      for ( int k = 0; k < box_units[2]; k++ ) {
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


// Build full neighbour list - Verlet algorithm
// Note that this implies 3D PBC
void compute_neigh_verlet( const double* const pos, const int natoms, 
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
      dx -= floor( ( dx + boxhalf0 ) / boxlen0 ) * boxlen0;

      double dy = pos[ 3 * i + 1 ] - pos[ 3 * j + 1 ];
      dy -= floor( ( dy + boxhalf1 ) / boxlen1 ) * boxlen1;

      double dz = pos[ 3 * i + 2 ] - pos[ 3 * j + 2 ];
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
      dx -= floor( ( dx + boxhalf0 ) / boxlen0 ) * boxlen0;

      double dy = y - pos[ 3 * j + 1 ];
      dy -= floor( ( dy + boxhalf1 ) / boxlen1 ) * boxlen1;

      double dz = z - pos[ 3 * j + 2 ];
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
    x -= floor( x / boxlen0 ) * boxlen0;
    pos[ 3 * i + 0 ] = x;
    posraw[ 3 * i + 0 ] += dx;
  
    double dy = vel[ 3 * i + 1 ] * dt + forc[ 3 * i + 1 ] * forc_hdtsq_scale * imass;
    double y = pos[ 3 * i + 1 ] + dy;
    y -= floor( y / boxlen1 ) * boxlen1;
    pos[ 3 * i + 1 ] = y;
    posraw[ 3 * i + 1 ] += dy;
  
    double dz = vel[ 3 * i + 2 ] * dt + forc[ 3 * i + 2 ] * forc_hdtsq_scale * imass;
    double z = pos[ 3 * i + 2 ] + dz;
    z -= floor( z / boxlen2 ) * boxlen2;
    pos[ 3 * i + 2 ] = z;
    posraw[ 3 * i + 2 ] += dz;
  }
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
  printf( "\n%16c %16c %16c\n", 'X', 'Y', 'Z' );
  for ( int i = istart; i < istop; i++) {
    printf( "%+16.6E %+16.6E %+16.6E\n", arr[ 3 * i + 0 ], arr[ 3 * i + 1 ], arr[ 3 * i + 2 ] );
  }

  return;
}


// Print information on simulation
void print_info ( const int* const box_units, const int nsteps, 
    const double temp_ini, const int nneighupd, 
    const int nthermo, const int ndump, 
    const double dt, const double cut, const double cellpar, 
    const double* const boxlen, const int natoms ) 
{
  printf( "\n Box Units : %-3i  %-3i  %-3i\n", 
    box_units[0], box_units[1], box_units[2] );
  printf( " No. Time Steps : %i\n", nsteps );
  printf( " Initial Temp [K] : %-6.1F\n", temp_ini );
  printf( " Neigh Update Freq : %i\n", nneighupd );
  printf( " Thermo Print Freq : %i\n", nthermo );
  printf( " Coord Dump Freq : %i\n", ndump );

  printf( "\n Time Step [ps] : %-5.3F\n", dt );
  printf( " Cutoff Dist [Ang] : %-5.3F\n", cut );
  printf( " Cell Par [Ang] : %-5.3F\n", cellpar );
  printf( " Box Length [Ang] : %-7.3F  %-7.3F  %-7.3F\n", 
    boxlen[0], boxlen[1], boxlen[2] );
  printf( " No. Atoms : %i\n", natoms );

  return;
}


// Print thermodynamic information
void print_thermo( const int istep, const double time, const double temp, 
    const double ekin, const double epot, const double etot, 
    const double clock ) 
{
  printf( " %9i  %10.3F  %8.3F  %+12.9F  %+12.9F  %+12.9F  %10.3F\n", 
         istep, time, temp, ekin, epot, etot, clock );

  return;
}


// Dump atomic coordinates
void dump_pdb( FILE* file, const int istep, 
    const double* const boxlen, const double boxang, 
    const char* const elsym, const double* const pos, const int natoms ) 
{
  fprintf( file, "REMARK --- frame: %-5i\n", istep );
  fprintf( file, "CRYST1%9.3F%9.3F%9.3F%7.2F%7.2F%7.2F\n", boxlen[0], boxlen[1], boxlen[2], boxang, boxang, boxang );
  if ( natoms < 100000 ) {
    for (int i = 0; i < natoms; i++ ) {
      fprintf( file, "ATOM  %5i %4s UNK  %-5i   %8.3F%8.3F%8.3F  1.00  0.00          %2s  \n", 
               i+1, elsym , i+1, pos[ 3 * i + 0 ], pos[ 3 * i + 1 ], pos[ 3 * i + 2 ], elsym );
    }
  } else {
    for (int i = 0; i < natoms; i++ ) {
      fprintf( file, "ATOM  %5X %4s UNK  %-5X   %8.3F%8.3F%8.3F  1.00  0.00          %2s  \n", 
               i+1, elsym , i+1, pos[ 3 * i + 0 ], pos[ 3 * i + 1 ], pos[ 3 * i + 2 ], elsym );
    }
  }
  fprintf( file, "END   \n" );
  return;
}
