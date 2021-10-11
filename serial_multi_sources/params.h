#ifndef PARAMS_H
  #define PARAMS_H
  #include "input.h"

// This structure definition is not ideal (in OOP parameters would be spread across multiple classes); using it here for algorithm readability
  struct Params {

// Define parameters - using LAMMPS "metal" physical units convention
// pressure unit is bar
// force unit is eV/Ang
//
    int box_units[ 3 ]; // no of unit cells per each dimension in the simulation box
    int nsteps; // no of time steps in the simulation
    double temp_ini; // K [117.7: datum from LAMMPS LJ example]
    int nneighupd; // update neighbour list every these steps [from LAMMPS LJ example]
    int nthermo; // print thermo info every these steps
    int ndump; // dump structure every these steps
//
    const char* coorfile = "coord.pdb"; // filename for initial atomic coordinates
    const char* trajfile = "traj.pdb"; // filename for trajectory atomic coordinates
//
// Crystal structure for Argon (fcc)
// Note that fcc implies 3D PBC
    static const int funits = 4;
    int natoms; // note that this implies 3D PBC // affected by input parameters
    const double cellpar = 5.795; // angstrom [datum from LAMMPS LJ example] [5.256: from real data]
    const double cellang = 90.0; // degrees
    double boxlen[ 3 ]; // affected by input parameters
    double boxhalf[ 3 ]; // affected by input parameters
    const double unitpos[ funits * 3 ] = {
      0., 0., 0.,
      0.5*cellpar, 0.5*cellpar, 0.,
      0.5*cellpar, 0., 0.5*cellpar,
      0., 0.5*cellpar, 0.5*cellpar
    };
//
// Some physical constants here
    const double k_B = 8.617343e-05; // eV/K
    const double N_av = 6.02214129e23; // mol-1
    const double J_eV = 1.602177e-19; // this is q_e
//
// Model parameters
    const double dt = 0.001; // ps
    const double hdt = 0.5 * dt;
    const double hdtsq = 0.5 * dt * dt;
    const double mass = 39.95; // gram/mol (also amu)
    const double imass = 1.0 / mass;
    const int atno = 18; // atomic number
    const char* elsym = "Ar"; // element symbol
    const double eps_kB = 117.7; // K
    const double eps = eps_kB * k_B; // eV
    const double sigma = 3.504; // angstrom
    const double sigma2 = sigma * sigma;
    const double sigma6 = sigma2 * sigma2 * sigma2;
    const double cut_fac = 2.5; // adimensional, multiplies sigma
    const double skin_fac = 0.3; // adimensional, multiplies sigma
    const double cut = cut_fac * sigma;
    const double cutskin = cut + skin_fac * sigma;
    const double cutsq = cut * cut;
    const double cutskinsq = cutskin * cutskin;
    const int maxneigh = 150; // with an fcc of side 5.256, and cut+skin of 9.8112, the real maxneigh is 86
//
    double N_dof; // note that this implies 3D PBC (different expressions for lower dimensionalities) // affected by input parameters
    const double ekin_scale = 10.0 / N_av / J_eV; // 1.036427e-04; // this factor is needed when using metal units ("mvv2e" in Mantevo/miniMD) [from my notes]
    double temp_scale; // [from my notes] // affected by input parameters
    const double forc_scale = 1.0 / ekin_scale; // J_eV * N_av * 0.1; // 9648.536; // this factor is needed when using metal units [from my notes]
    const double forc_hdt_scale = forc_scale * hdt;
    const double forc_hdtsq_scale = forc_scale * hdtsq;

  };

  Params get_params( InputParams );

#endif
