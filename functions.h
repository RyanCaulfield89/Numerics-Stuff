///////////////////////////////////////////////////////////////////////////////
//
//    File: functions.h
//    Description: Header file for functions.h
//                 This is a collection of functions that we maybe care about.
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To do list:
//      1 - Nothing yet.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef FUNCTIONS_H
#define FUNCTIONS_H //To avoid including this file twice
//Headers

//structures for the function parameters

//For quadratic_polynomial
typedef struct
{
  std::complex<double> a;
  std::complex<double> b;
  std::complex<double> c;
}
quadratic_parameters;

///////////////////////////////////////////////////////////////////////////////
std::complex<double> coulomb_potential(double x, void *params);
std::complex<double> harmonic_potential(double x, void *params);
std::complex<double> nonlocal_coulomb_potential(double x, double y, void *params);
std::complex<double> nonlocal_harmonic_potential(double x, double y, void *params);
std::complex<double> nonlocal_potential(double x, double y, void *params);
std::complex<double> seperable_nonlocal_potential(double x, double y, void *params);
std::complex<double> quadratic_polynomial(double x, void *params);

#endif
