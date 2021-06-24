#ifndef INPUT_H
  #define INPUT_H

  struct InputParams {
    int box_units;
    int nsteps;
    double temp_ini;
    int nneighupd;
    int nthermo;
    int ndump;
  };

  InputParams get_input_params( const int, char** );

#endif
