#include "params.h"
#include "input.h"


// Define parameters that are affected by input parameters
Params get_params( InputParams input_params )
{
  Params params;

  params.box_units[0] = input_params.box_units[0];
  params.box_units[1] = input_params.box_units[1];
  params.box_units[2] = input_params.box_units[2];
  params.nsteps = input_params.nsteps;
  params.temp_ini = input_params.temp_ini;
  params.nneighupd = input_params.nneighupd;
  params.nthermo = input_params.nthermo;
  params.ndump = input_params.ndump;

  params.natoms = params.funits * params.box_units[0] * params.box_units[1] * params.box_units[2];
  params.boxlen[0] = params.cellpar * params.box_units[0];
  params.boxlen[1] = params.cellpar * params.box_units[1];
  params.boxlen[2] = params.cellpar * params.box_units[2];
  params.boxhalf[0] = 0.5 * params.boxlen[0];
  params.boxhalf[1] = 0.5 * params.boxlen[1];
  params.boxhalf[2] = 0.5 * params.boxlen[2];
  params.N_dof = ( params.natoms * 3 - 3 );
  params.temp_scale = params.ekin_scale / ( params.N_dof * params.k_B );

  return params;
}
