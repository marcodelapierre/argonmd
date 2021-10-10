## ArgonMD

Training code: Molecular Dynamics  
Lennard Jones solid with velocity Verlet integrator  

Written in C++  


### Assumptions

* Physical model
  * Homo-atomic Argon system
    * Implies single mass value, single set of force-field parameters
  * Condensed phase
  * Cubic simulation box as a supercell of the crystal faced-centered cubic (*fcc*) unit cell
  * Three-dimensional (3D) periodic boundary conditions (PBC)
  * NVE ensemble (constant Number of particles, Volume and total Energy)
  * Lennard-Jones (LJ) pairwise interactions, with distance cut-off

* Algorithmic aspects
  * Velocity Verlet integrator
  * Full Verlet neighbour list with regular updates
    * The cell list implementation also exists, with better performance (linear rather than quadratic scaling)
  * *Metal* physical units as in LAMMPS
  * Starting structure is the equilibrium *fcc* lattice for solid Argon
  * Velocities initialised with random number generator, deterministic seed
  * Saving both raw and PBC-wrapped atomic coordinates
  * Optional dumping of atomic coordinates
    * PDB format - large files, for demonstration only
    * PDB outputs can be opened by VMD


### Available Procedural Implementations
1. `serial_one_source`: serial, single source file
2. `serial_multi_sources`: serial, multiple source files


### Usage

* Compile with `make` (parent directory, or any implementation sub-directory)

* Run with `./argonmd.x` (any implementation sub-directory)

* Output of default run: 

  ```
  
  ** ArgonMD **
  
   Box Units : 5
   No. Time Steps : 10000
   Initial Temp [K] : 10.0
   Neigh Update Freq : 20
   Thermo Print Freq : 1000
   Coord Dump Freq : 0
  
   Time Step [ps] : 0.001
   Cutoff Dist [Ang] : 8.760
   Cell Par [Ang] : 5.795
   Box Length [Ang] : 28.975
   No. Atoms : 500
  
        Step    Time[ps]   Temp[K]      Ekin[eV]      Epot[eV]      Etot[eV]    Clock[s]
           0       0.000    10.000  +0.001290016  -0.072045884  -0.070755867       0.000
        1000       1.000     5.973  +0.000770507  -0.071526372  -0.070755865       0.458
        2000       2.000     4.778  +0.000616321  -0.071376488  -0.070760167       0.938
        3000       3.000     4.492  +0.000579534  -0.071339370  -0.070759837       1.393
        4000       4.000     4.805  +0.000619904  -0.071377424  -0.070757520       1.845
        5000       5.000     5.007  +0.000645927  -0.071402785  -0.070756858       2.306
        6000       6.000     4.770  +0.000615386  -0.071373898  -0.070758513       2.776
        7000       7.000     5.101  +0.000658049  -0.071414245  -0.070756196       3.289
        8000       8.000     4.975  +0.000641736  -0.071398594  -0.070756858       3.743
        9000       9.000     4.839  +0.000624277  -0.071381135  -0.070756857       4.209
       10000      10.000     4.879  +0.000629441  -0.071386299  -0.070756857       4.670
  
   Loop Clock Time [s] :      4.670
  ```

* Editable input parameters

  The first six parameters, as per program output, are customisable through optional command line arguments.  
  Add a first argument to edit the number of box units, a second argument to also edit the number of time steps, and so on.  

  Example of run with custom box units (`4`) and time steps (`50000`):
  ```
  ./argonmd.x 4 50000
  ```

  Example of fully customised run:
  ```
  ./argonmd.x 10 1000 20. 20 100 100
  ```

* Number of box units
  
  The first optional parameter specifies how many times the Argon unit cell is replicated along *x*, *y* and *z*, to create the simulation box.  As each unit cell contains 4 atoms, this parameter also directly determines the total number of atoms in the simulation.  
  For instance, the default value of 5 means that the simulation box contains 5 X 5 X 5 = 125 unit cells, and 125 X 4 = 500 atoms.


### Extras
* The directory `Sample_runs/` contains sample output files produced by *ArgonMD* with various sets of input parameters
* The directory `Examples_lammps/` contains LAMMPS input/output files that come from the [LAMMPS source](https://github.com/lammps/lammps), `examples/UNITS/`;  they allow to compare outputs from *ArgonMD* with those from LAMMPS


### Credits
* Some parts of the code were inspired by code in [Mantevo miniMD](https://github.com/Mantevo/miniMD)
* Parameter values were taken by `examples/UNITS/` in [LAMMPS](https://github.com/lammps/lammps)
