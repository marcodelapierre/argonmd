#include <iostream>
//#include <iomanip> // test only

using namespace std;


// service routines
double random(int*); // This one is taken from Mantevo/miniMD


int main() {
//cout<<"Hello World!"<<endl;

// define parameters
// first three might become editable by cli arguments
const int box_side = 2; // no of unit cells
const int nsteps = 1000000;
const double temp = 50.; // K
//
const double eps = 0.34; // nm
const double sigma = 120.; // K*k_B
const double mass = 39.9; // amu
const double step = 0.001; // ps
// argon crystal structure (fcc)
const int ndims = 3; // no of spatial dimensions
const double cellpar = 0.5256; // nm
const int funits = 4;
const int fd = funits * ndims;
const double unitpos[ fd ] = {
  0., 0., 0.,
  0.5*cellpar, 0.5*cellpar, 0.,
  0.5*cellpar, 0., 0.5*cellpar,
  0., 0.5*cellpar, 0.5*cellpar
};

// physical constants here


// allocate arrays
const int natoms = funits * box_side * box_side * box_side;
//
double* pos = new double [ natoms * ndims ];
double* vel = new double [ natoms * ndims ];
double* acc = new double [ natoms * ndims ];


// define structure AND 
// initialise velocities
int idx;
int seed;
int bfd = box_side * fd;
int bbfd = box_side * bfd;
double vxtmp = 0.;
double vytmp = 0.;
double vztmp = 0.;
for ( int i = 0; i < box_side; i++ ) {
  for ( int j = 0; j < box_side; j++ ) {
    for ( int k = 0; k < box_side; k++ ) {
      for ( int l = 0; l < funits; l++ ) {
        idx = i * bbfd + j * bfd + k * fd + l * ndims;
        // positions
        pos[ idx + 0 ] = cellpar*i + unitpos[ l * ndims + 0 ];
        pos[ idx + 1 ] = cellpar*j + unitpos[ l * ndims + 1 ];
        pos[ idx + 2 ] = cellpar*k + unitpos[ l * ndims + 2 ];
        //cout << left << setw(8) << pos[idx+0] << setw(8) << pos[idx+1] << setw(8) << pos[idx+2] << endl ; // test only
        //cout << idx << endl; // test only

        // velocities
        seed = idx;
        for(int m = 0; m < 5; m++) random(&seed);
        vel[ idx + 0 ] = 0; //random(&seed);
        for(int m = 0; m < 5; m++) random(&seed);
        vel[ idx + 1 ] = random(&seed);
        for(int m = 0; m < 5; m++) random(&seed);
        vel[ idx + 2 ] = random(&seed);

        vxtmp += vel[ idx + 0 ];
        vytmp += vel[ idx + 1 ];
        vztmp += vel[ idx + 2 ];
      }
      //cout << endl; // test only
    }
  }
}


// adjust velocities
// zero centre-of-mass motion
for (int i =0; i < natoms; i++) {
  vel[ i * ndims + 0 ] -= vxtmp;
  vel[ i * ndims + 1 ] -= vytmp;
  vel[ i * ndims + 2 ] -= vztmp;
}
// rescale to desired temperature




// big loop: time evolution

// compute forces

// integrate

// update velocities




// compute and print output when required



// deallocate arrays
delete [] acc;
delete [] vel;
delete [] pos;


return 0;
}


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
