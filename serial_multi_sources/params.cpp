#include "params.h"
#include "input.h"


// Define parameters - using LAMMPS "metal" physical units convention
Params get_params( InputParams input_params )
{
  Params params;

  params.box_units = input_params.box_units; // no of unit cells per dimension in the simulation box
  params.nsteps = input_params.nsteps; // no of time steps in the simulation
  params.temp_ini = input_params.temp_ini; // K [117.7: datum from LAMMPS LJ example]
  params.nneighupd = input_params.nneighupd; // update neighbour list every these steps [from LAMMPS LJ example]
  params.nthermo = input_params.nthermo; // print thermo info every these steps
  params.ndump = input_params.ndump; // dump structure every these steps

  params.natoms = params.funits * params.box_units * params.box_units * params.box_units; // note that this implies 3D PBC // affected by input parameters
  params.boxlen = params.cellpar * params.box_units; // affected by input parameters
  params.boxhalf = params.boxlen * 0.5; // affected by input parameters
  params.N_dof = ( params.natoms * 3 - 3 ); // note that this implies 3D PBC (different expressions for lower dimensionalities) // affected by input parameters
  params.temp_scale = params.ekin_scale / ( params.N_dof * params.k_B ); // [from my notes] // affected by input parameters

  return params;
}
