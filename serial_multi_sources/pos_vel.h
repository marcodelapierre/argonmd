#ifndef POS_VEL_H
  #define POS_VEL_H

  void check_pbc( double*, const int, const double* const );
  void update_pos_pbc( double*, double*, const double* const, const double* const, const int, 
      const double, const double, const double, const double* const );
  void update_vel( double*, const double* const, const double* const, const int, const double, const double );

#endif
