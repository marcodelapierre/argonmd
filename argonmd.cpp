#include <iostream>
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
const int box_side = 10; // no of unit cells
const double temp = 50.; // K
const double eps = 0.34; // nm
const double sigma = 120.; // K*k_B
const double mass = 39.9; // amu
const double step = 0.001; // ps
const int nsteps = 1000000;
// argon crystal structure (fcc)
const double cellpar = 0.5256; // nm
const int funits = 4;
// physical constants here


// allocate arrays
int natoms = funits * box_side * box_side * box_side;
const int dims = 3;
//
double** pos = arrayAlloc( natoms, 3);
double** vel = arrayAlloc( natoms, 3);
double** acc = arrayAlloc( natoms, 3);


// define structure



// initialise velocities



// big loop: time evolution

// compute forces

// integrate

// update velocities

// print output when required



// deallocate arrays
arrayDealloc( acc );
arrayDealloc( vel );
arrayDealloc( pos );


return 0;
}
