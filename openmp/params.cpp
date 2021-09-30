#include "params.h"
#include "input.h"


// Define parameters that are affected by input parameters
Params get_params( InputParams input_params )
{
  Params params;

  params.box_units = input_params.box_units;
  params.nsteps = input_params.nsteps;
  params.temp_ini = input_params.temp_ini;
  params.nneighupd = input_params.nneighupd;
  params.nthermo = input_params.nthermo;
  params.ndump = input_params.ndump;

  params.natoms = params.funits * params.box_units * params.box_units * params.box_units;
  params.boxlen = params.cellpar * params.box_units;
  params.boxhalf = params.boxlen * 0.5;
  params.N_dof = ( params.natoms * 3 - 3 );
  params.temp_scale = params.ekin_scale / ( params.N_dof * params.k_B );

  return params;
}
