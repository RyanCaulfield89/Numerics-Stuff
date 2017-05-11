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
}
