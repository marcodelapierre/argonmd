#include <iostream>
#include "input.h"

using namespace std;


// Get input parameters
InputParams get_input_params( const int argc, char** argv ) 
{
  InputParams input_params;

  if ( argc > 1 ) {
    input_params.box_units = atoi(argv[1]);
  } else {
    input_params.box_units = 5; 
  }
  if ( argc > 2 ) {
    input_params.nsteps = atoi(argv[2]);
  } else {
    input_params.nsteps = 10000;
  }
  if ( argc > 3 ) {
    input_params.temp_ini = atof(argv[3]);
  } else {
    input_params.temp_ini = 10.;
  }
  if ( argc > 4 ) {
    input_params.nneighupd = atoi(argv[4]);
  } else {
    input_params.nneighupd = 20;
  }
  if ( argc > 5 ) {
    input_params.nthermo = atoi(argv[5]);
  } else {
    input_params.nthermo = 1000;
  }
  if ( argc > 6 ) {
    input_params.ndump = atoi(argv[6]);
  } else {
    input_params.ndump = 0;
  }

  return input_params;
}
