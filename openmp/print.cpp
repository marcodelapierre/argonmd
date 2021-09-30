#include <iostream>
#include "print.h"

using namespace std;


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
void print_info ( const int box_units, const int nsteps, const double temp_ini, 
                  const int nneighupd, const int nthermo, const int ndump, 
                  const double dt, const double cut, 
                  const double cellpar, const double boxlen, const int natoms ) 
{
  printf( "\n Box Units : %i\n", box_units );
  printf( " No. Time Steps : %i\n", nsteps );
  printf( " Initial Temp [K] : %-6.1F\n", temp_ini );
  printf( " Neigh Update Freq : %i\n", nneighupd );
  printf( " Thermo Print Freq : %i\n", nthermo );
  printf( " Coord Dump Freq : %i\n", ndump );

  printf( "\n Time Step [ps] : %-5.3F\n", dt );
  printf( " Cutoff Dist [Ang] : %-5.3F\n", cut );
  printf( " Cell Par [Ang] : %-5.3F\n", cellpar );
  printf( " Box Length [Ang] : %-7.3F\n", boxlen );
  printf( " No. Atoms : %i\n", natoms );

  return;
}


// Print thermodynamic information
void print_thermo( const int istep, const double time, 
                   const double temp, const double ekin, const double epot, 
                   const double etot, const double clock ) 
{
  printf( " %9i  %10.3F  %8.3F  %+12.9F  %+12.9F  %+12.9F  %10.3F\n", 
         istep, time, temp, ekin, epot, etot, clock );

  return;
}


// Dump atomic coordinates
void dump_pdb( FILE* file, const int istep, 
               const double boxlen, const double boxang, 
               const char* elsym, const double* const pos, const int natoms ) 
{
  fprintf( file, "REMARK --- frame: %-5i\n", istep );
  fprintf( file, "CRYST1%9.3F%9.3F%9.3F%7.2F%7.2F%7.2F\n", boxlen, boxlen, boxlen, boxang, boxang, boxang );
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
