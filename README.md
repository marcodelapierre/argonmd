## Argonmd

Training code: Molecular Dynamics  
Lennard Jones solid with velocity Verlet integrator


### Assumptions

* Physical model
  * *Metal* units from LAMMPS
  * Homo-atomic Argon system
    * implies single mass value, single set of force-field parameters
  * Three-dimensional (3D) periodic boundary conditions (PBC)
  * Condensed phase
  * Starting structure is the equilibrium faced-centered cubic (*fcc*) lattice
  * NVE ensemble (constant Number of particles, Volume and total Energy)
  * Lennard-Jones (*lj*) pairwise interactions, with distance cut-off

* Algorithms
  * Velocity Verlet integrator
  * Full neighbour list


### To-do List

* Multi-file code base
* Parallel implementations
  * OpenMP multi-threaded
  * OpenMP GPU offloading
  * CUDA
  * HIP
* Object Oriented (far future)


### Credits
* Some parts of the code were inspired by code in [Mantevo miniMD](https://github.com/Mantevo/miniMD)
* Parameter values were taken by the `examples/UNITS` in [LAMMPS](https://github.com/lammps/lammps)
