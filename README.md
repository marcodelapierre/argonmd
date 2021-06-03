## Argonmd

Training code: Molecular Dynamics  
Lennard Jones solid with velocity Verlet integrator

### Assumptions
* Argon chosen as model system
  - *fcc* atomic lattice
  - lattice constant
  - *LJ* parameters, with cut-off
* 3D cubic box
  - defined as supercell of the unit cell of Argon
* Solid state (at least to begin with)
  - obtained using a low temperature value
  - simplifies handling of neighbours
* NVE ensemble
* Full neighbour list

### Credits
* Some parts of the code were inspired by code in [Mantevo miniMD](https://github.com/Mantevo/miniMD)
* Parameter values were taken by the `examples/UNITS` in [LAMMPS](https://github.com/lammps/lammps)
