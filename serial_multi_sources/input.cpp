#include <iostream>
#include "input.h"

using namespace std;


// Get input parameters
InputPars get_input_pars( const int argc, char** argv ) 
{
  InputPars input_pars;

  if ( argc > 1 ) {
    input_pars.box_units = atoi(argv[1]);
  } else {
    input_pars.box_units = 5; 
  }
  if ( argc > 2 ) {
    input_pars.nsteps = atoi(argv[2]);
  } else {
    input_pars.nsteps = 10000;
  }
  if ( argc > 3 ) {
    input_pars.temp_ini = atof(argv[3]);
  } else {
    input_pars.temp_ini = 10.;
  }
  if ( argc > 4 ) {
    input_pars.nneighupd = atoi(argv[4]);
  } else {
    input_pars.nneighupd = 20;
  }
  if ( argc > 5 ) {
    input_pars.nthermo = atoi(argv[5]);
  } else {
    input_pars.nthermo = 1000;
  }
  if ( argc > 6 ) {
    input_pars.ndump = atoi(argv[6]);
  } else {
    input_pars.ndump = 0;
  }

  return input_pars;
}
