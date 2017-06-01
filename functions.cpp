///////////////////////////////////////////////////////////////////////////////
//
//    File: functions.cpp
//    Description: This is a collection of functions that we care about.
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To do list:
//      1 - Nothing yet.
//
///////////////////////////////////////////////////////////////////////////////

//Includes
#include <cmath>
#include <complex>
#include "functions.h"
using namespace std;

///////////////////////////////////////////////////////////////////////////////

complex<double> coulomb_potential(double x, void *params)
{
  complex<double> charge;
  charge = *(complex<double> *) params;
  complex<double> return_value;
  return_value = charge / x;
  return (return_value);
}

complex<double> harmonic_potential(double x, void *params)
{
  complex<double> k;
  k = *(complex<double> *) params;
  complex<double> return_value = 0.5*k*x*x;
  return (return_value);
}

complex<double> nonlocal_coulomb_potential(double x, double y, void *params)
{
  complex<double> beta;
  beta = *(complex<double> *) params;
  complex<double> return_value;
  return_value =  1 / fabs(x - y) * exp(-(x - y)*(x - y) / beta);
  return (return_value);
}

complex<double> nonlocal_harmonic_potential(double x, double y, void *params)
{
  complex<double> beta;
  beta = *(complex<double> *) params;
  complex<double> return_value;
  return_value = (x - y)*(x - y) * exp(-(x - y)*(x - y) / beta);
  return (return_value);
}

complex<double> nonlocal_potential(double x, double y, void *params)
{
  complex<double> beta;
  beta = *(complex<double> *) params;
  complex<double> return_value;
  return_value = - (x + y) * exp(-(abs(x - y) - 1)*(abs(x - y) - 1) / beta);
  return (return_value);
}

complex<double> seperable_nonlocal_potential(double x, double y, void *params){
  complex<double> beta;
  beta = *(complex<double> *) params;
  complex<double> return_value;
  return_value = x*x*exp(-beta*y*y);
  return return_value;
}

complex<double> quadratic_polynomial(double x, void *params){
  complex<double> a = ((quadratic_parameters *) params)->a;
  complex<double> b = ((quadratic_parameters *) params)->b;
  complex<double> c = ((quadratic_parameters *) params)->c;
  complex<double> return_value;
  return_value = a*x*x + b*x + c;
  return return_value;
}

std::complex<double> gaussian(double x, void *params){
  complex<double> beta;
  beta = *(complex<double> *) params;
  complex<double> return_value;
  return_value = exp(-beta*x*x);
  return return_value;
}

std::complex<double> rational_function(double x, void *params){
  complex<double> beta;
  beta = *(complex<double> *) params;
  complex<double> return_value;
  return_value = 1.0 / (1.0 + beta * x * x);
  return return_value;
}

std::complex<double> woods_saxon_potential(double x, void *params){
  double R = ((woods_saxon_parameters *) params)->R;
  double a = ((woods_saxon_parameters *) params)->a;
  complex<double> return_value;
  return_value = 1.0 / (1 + exp((x - R) / a));
  return return_value;
}

std::complex<double> diff_woods_saxon_potential(double x, void *params){
  double R = ((woods_saxon_parameters *) params)->R;
  double a = ((woods_saxon_parameters *) params)->a;
  complex<double> return_value;
  return_value = - 1.0 / (a * pow(1 + exp((x - R) / a),2));
  return return_value;
}

std::complex<double> optical_potential(double x, void *params){
  //Parameters for the real part of the potential
  woods_saxon_parameters real_parameters =
    ((optical_potential_parameters *) params)->real_parameters;
  void *real_parameters_ptr;
  real_parameters_ptr = &real_parameters;
  //Parameters for the imginary inelastic part of the potential.
  //Note that the same parameters will be used for the imaginary part U_I
  //and the boundary part U_D.
  woods_saxon_parameters imaginary_parameters =
    ((optical_potential_parameters *) params)->imaginary_parameters;
  void *imaginary_parameters_ptr;
  imaginary_parameters_ptr = &imaginary_parameters;
  //These are the strengths of the various components of the potential
  double V = ((optical_potential_parameters *) params)->V;
  double W = ((optical_potential_parameters *) params)->W;
  double W_D = ((optical_potential_parameters *) params)->W_D;
  complex<double> return_value;
  //The value is given by U(r) = U_R(r) + U_I(r) + U_D(r)
  // with U_R(r) = -V * woods_saxon_potential(r)
  //      U_I(r) = -W * woods_saxon_potential(r)
  //      U_D(r) = 4ia * W_D * diff_woods_saxon_potential(r)
  //I'll add spin-orbit and coulomb terms later
  const   complex<double> i(0.0,1.0);
  return_value = - 1.0 * V * woods_saxon_potential(x,real_parameters_ptr)
                 - W * i * woods_saxon_potential(x,imaginary_parameters_ptr)
                 + 4.0 * i *  imaginary_parameters.a * W_D
                 * diff_woods_saxon_potential(x,imaginary_parameters_ptr);
  return return_value;
}
