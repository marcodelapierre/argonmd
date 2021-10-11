#include <iostream>
#include <cstring>
#include "input.h"

using namespace std;


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
