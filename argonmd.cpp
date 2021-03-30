#include <iostream>
//#include <iomanip> // test only

using namespace std;


double ** arrayAlloc (int rows, int cols){
    double** arr = new double* [rows];
    double* pool = new double [cols];
    for(int i = 0;i < rows; i++, pool += cols){
           arr[i] = pool;
    }
    return arr;
}
void arrayDealloc(double **arr){
        delete [] arr[0];
        delete [] arr;
    }


int main() {
cout<<"Hello World!"<<endl;

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
const double cellpar = 0.5256; // nm
const int funits = 4;
const double unitpos[4][3] = {
  0., 0., 0.,
  0.5*cellpar, 0.5*cellpar, 0.,
  0.5*cellpar, 0., 0.5*cellpar,
  0., 0.5*cellpar, 0.5*cellpar
};

// physical constants here


// allocate arrays
int natoms = funits * box_side * box_side * box_side;
const int dims = 3;
//
double** pos = arrayAlloc( natoms, 3);
double** vel = arrayAlloc( natoms, 3);
double** acc = arrayAlloc( natoms, 3);


// define structure
int idx;
int bf = box_side * funits;
int bbf = box_side * bf;
for ( int i = 0; i < box_side; i++ ) {
  for ( int j = 0; j < box_side; j++ ) {
    for ( int k = 0; k < box_side; k++ ) {
      for ( int l = 0; l < funits; l++ ) {
        idx = i * bbf + j * bf + k * funits;
        pos[idx][0] = cellpar*i + unitpos[l][0];
        pos[idx][1] = cellpar*j + unitpos[l][1];
        pos[idx][2] = cellpar*k + unitpos[l][2];
        //cout << left << setw(8) << pos[idx][0] << setw(8) << pos[idx][1] << setw(8) << pos[idx][2] << endl; // test only
      }
      //cout << endl; // test only
    }
  }
}


// initialise velocities



// big loop: time evolution

// compute forces

// integrate

// update velocities




// compute and print output when required



// deallocate arrays
arrayDealloc( acc );
arrayDealloc( vel );
arrayDealloc( pos );


return 0;
}
