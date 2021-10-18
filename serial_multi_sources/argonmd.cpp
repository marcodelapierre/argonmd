#include <iostream>
#include <ctime>

#include "input.h"
#include "params.h"
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
const InputParams input_params = get_input_params( argc, argv );

// Define all parameters
// This structure definition is not ideal (in OOP parameters would be spread across multiple classes); using it here for algorithm readability
const Params params = get_params( input_params );

// Allocate arrays
int* numneigh = new int [ params.natoms ];
int* neigh = new int [ params.natoms * params.maxneigh ];
double* pos = new double [ params.natoms * 3 ];
double* posraw = new double [ params.natoms * 3 ];
double* vel = new double [ params.natoms * 3 ];
double* forc = new double [ params.natoms * 3 ];
double* forcold = new double [ params.natoms * 3 ];

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
setup_struc_vel( params.funits, params.box_units, params.cellpar, params.unitpos, params.natoms, pos, posraw, vel );
// Rescale to desired temperature
compute_temp_ekin( vel, params.natoms, params.mass, params.temp_scale, params.ekin_scale, temp, ekin );
rescale_temp( vel, params.natoms, params.temp_ini, temp, ekin );

// PBC check // not needed at startup with current input structure, yet here for generality
check_pbc( pos, params.natoms, params.boxlen );
// Build (full) neighbour list
compute_neigh_verlet( pos, params.natoms, params.boxlen, params.boxhalf, params.cutskinsq, params.maxneigh, numneigh, neigh );

// Compute initial forces
compute_forc_epot( pos, params.natoms, params.maxneigh, numneigh, neigh, 
                   params.boxlen, params.boxhalf, params.cutsq, params.sigma6, params.eps, forc, epot );

// Print simulation info
print_info( params.box_units, params.nsteps, params.temp_ini, params.nneighupd, params.nthermo, params.ndump, 
            params.dt, params.cut, params.cellpar, params.boxlen, params.natoms );
//if ( 0 ) { print_arr( pos, 0, natoms ); print_arr( vel, 0, natoms ); } // debug print
//
// Print initial thermo output
printf( "\n %9s  %10s  %8s  %12s  %12s  %12s  %10s\n", "Step", "Time[ps]", "Temp[K]", "Ekin[eV]", "Epot[eV]", "Etot[eV]", "Clock[s]" );
print_thermo( istep, params.dt*istep, temp, ekin/params.natoms, epot/params.natoms, (ekin+epot)/params.natoms, 0.0 );

// Dump initial atomic coordinates
if ( params.ndump > 0 ) {
  coor = fopen( params.coorfile, "w" );
  dump_pdb( coor, istep, params.boxlen, params.cellang, params.elsym, pos, params.natoms );
  fclose( coor );
  traj = fopen( params.trajfile, "w" );
}


// Time evolution loop
start = clock();
for (istep = 1; istep <= params.nsteps; istep++) {

// This code block could become a routine "integrate"; leaving it here for algorithm readability
// Note that this implies Velocity Verlet integrator
  {
// Update positions and check PBC meanwhile
  update_pos_pbc( pos, posraw, vel, forc, params.natoms, params.dt, params.forc_hdtsq_scale, params.imass, params.boxlen );

// Update (full) neighbour list
  if( params.nneighupd > 0 && istep%params.nneighupd == 0 ) { 
    compute_neigh_verlet( pos, params.natoms, params.boxlen, params.boxhalf, params.cutskinsq, params.maxneigh, numneigh, neigh );
  }

// Store old forces and compute new forces
  forctmp = forcold;
  forcold = forc;
  forc = forctmp;
  compute_forc_epot( pos, params.natoms, params.maxneigh, numneigh, neigh, 
                     params.boxlen, params.boxhalf, params.cutsq, params.sigma6, params.eps, forc, epot );

// Update velocities
  update_vel( vel, forcold, forc, params.natoms, params.forc_hdt_scale, params.imass );
  }


  if ( params.nthermo > 0 && istep%params.nthermo == 0 ) {
// Compute temperature when required
    compute_temp_ekin( vel, params.natoms, params.mass, params.temp_scale, params.ekin_scale, temp, ekin );

// Get clock time when required
    watch = clock() - start;
    clocktime = ((float)watch)/CLOCKS_PER_SEC;

// Print thermo output when required
    print_thermo( istep, params.dt*istep, temp, ekin/params.natoms, epot/params.natoms, (ekin+epot)/params.natoms, clocktime );
  }

// Dump atomic coordinates when required
// Note that PDB files are large, so enable dumping only for demonstrations; ideally these should go in a binary format (eg DCD)
// Note also that ideally "posraw" should be used for positions; "pos" is used here instead, for convenience when demonstrating with VMD
  if ( params.ndump > 0 && istep%params.ndump == 0 ) {
    dump_pdb( traj, istep, params.boxlen, params.cellang, params.elsym, pos, params.natoms );
  }

}

// Get and print final clocktime
watch = clock() - start;
clocktime = ((float)watch)/CLOCKS_PER_SEC;
printf( "\n Loop Clock Time [s] : %10.3F\n" , clocktime );


// Close trajectory file if needed
if ( params.ndump > 0 ) { fclose( traj ); }

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
