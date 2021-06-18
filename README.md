## ArgonMD

Training code: Molecular Dynamics  
Lennard Jones solid with velocity Verlet integrator


### Assumptions

* Physical model
  * Homo-atomic Argon system
    * Implies single mass value, single set of force-field parameters
  * Condensed phase
  * Starting structure is the equilibrium faced-centered cubic (*fcc*) lattice
  * Cubic box as supercell of crystal *fcc* unit cell
  * Three-dimensional (3D) periodic boundary conditions (PBC)
  * NVE ensemble (constant Number of particles, Volume and total Energy)
  * Lennard-Jones (LJ) pairwise interactions, with distance cut-off

* Algorithmic aspects
  * *Metal* physical units as in LAMMPS
  * Velocity Verlet integrator
  * Full neighbour list with regular updates
  * Saving both raw and PBC-wrapped atomic coordinates


### Available Implementations
* Serial, single source file


### To-Do
* Optional dumping of atomic coordinates


### Extras
* The `tests` directory contains sample output files produced by *ArgonMD* with various sets of input parameters
* The `examples_lammps` directory contains LAMMPS input/output files that come from the [LAMMPS source](https://github.com/lammps/lammps), `examples/UNITS/`;  they allow to compare outputs from *ArgonMD* with those from LAMMPS


### Credits
* Some parts of the code were inspired by code in [Mantevo miniMD](https://github.com/Mantevo/miniMD)
* Parameter values were taken by `examples/UNITS/` in [LAMMPS](https://github.com/lammps/lammps)
