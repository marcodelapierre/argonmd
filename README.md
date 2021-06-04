## Argonmd

Training code: Molecular Dynamics  
Lennard Jones solid with velocity Verlet integrator


### Assumptions

* Physical model
  * Homo-atomic Argon system
    * implies single mass value, single set of force-field parameters
  * Three-dimensional (3D) with periodic boundary conditions
  * Solid state (**for now** - Initial temperature (T) below 40 K to be safe)
  * Starting structure is the equilibrium faced-centered cubic (*fcc*) lattice
  * NVE ensemble (constant Number of particles, Volume and total Energy)
  * Lennard-Jones (*lj*) pairwise interactions, with distance cut-off

* Algorithms
  * Velocity Verlet integrator
  * Full neighbour list
    * computed once **right now** (requires solid state at low T not to crash!)


### Credits
* Some parts of the code were inspired by code in [Mantevo miniMD](https://github.com/Mantevo/miniMD)
* Parameter values were taken by the `examples/UNITS` in [LAMMPS](https://github.com/lammps/lammps)
