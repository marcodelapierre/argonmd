#include <iostream>
#include <ctime>

#include "input.h"
#include "setup.h"
#include "temp.h"
#include "neigh.h"
#include "force.h"
#include "pos_vel.h"
#include "print.h"

using namespace std;


int main( int argc, char** argv ) {
//cout<<"Ciao Mondo!"<<endl;

// Input parameters - editable by input
InputPars input_pars = get_input_pars( argc, argv );

// Define all parameters

// Define parameters - using LAMMPS "metal" physical units convention
// pressure unit is bar
// force unit is eV/Ang
//
const int box_units = input_pars.box_units; // no of unit cells per dimension in the simulation box
const int nsteps = input_pars.nsteps; // no of time steps in the simulation
const double temp_ini = input_pars.temp_ini; // K [117.7: datum from LAMMPS LJ example]
const int nneighupd = input_pars.nneighupd; // update neighbour list every these steps [from LAMMPS LJ example]
const int nthermo = input_pars.nthermo; // print thermo info every these steps
const int ndump = input_pars.ndump; // dump structure every these steps
//
const char* coorfile = "coord.pdb"; // filename for initial atomic coordinates
const char* trajfile = "traj.pdb"; // filename for trajectory atomic coordinates
//
// Crystal structure for Argon (fcc)
// Note that fcc implies 3D PBC
const int funits = 4;
const int natoms = funits * box_units * box_units * box_units; // note that this implies 3D PBC // affected by input parameters
const double cellpar = 5.795; // angstrom [datum from LAMMPS LJ example] [5.256: from real data]
const double cellang = 90.0; // degrees
const double boxlen = cellpar * box_units; // affected by input parameters
const double boxhalf = boxlen * 0.5; // affected by input parameters
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
compute_neigh( pos, natoms, boxlen, boxhalf, cutskinsq, maxneigh, numneigh, neigh );

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
    compute_neigh( pos, natoms, boxlen, boxhalf, cutskinsq, maxneigh, numneigh, neigh );
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
