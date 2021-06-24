#ifndef INPUT_H
  #define INPUT_H

// struct definitions
  struct InputPars {
    int box_units;
    int nsteps;
    double temp_ini;
    int nneighupd;
    int nthermo;
    int ndump;
  };

  InputPars get_input_pars( const int, char** );

#endif
