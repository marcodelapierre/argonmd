#ifndef PRINT_H
  #define PRINT_H

  void print_arr( const double* const, const int, const int );
  void print_info( const int* const, const int, const double, const int, const int, const int, 
                   const double, const double, const double, const double* const, const int );
  void print_thermo( const int, const double, const double, const double, const double, const double, const double );
  void dump_pdb( FILE*, const int, const double* const, const double, const char* const, const double* const, const int );

#endif
