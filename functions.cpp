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

std::complex<double> coulomb_potential(double x, void *params)
{
  std::complex<double> charge;
  charge = *(std::complex<double> *) params;
  std::complex<double> return_value;
  return_value = charge / x;
  return (return_value);
}

std::complex<double> harmonic_potential(double x, void *params)
{
  std::complex<double> k;
  k = *(std::complex<double> *) params;
  std::complex<double> return_value = 0.5*k*x*x;
  return (return_value);
}

std::complex<double> nonlocal_coulomb_potential(double x, double y, void *params)
{
  std::complex<double> beta;
  beta = *(std::complex<double> *) params;
  std::complex<double> return_value;
  return_value =  1 / fabs(x - y) * exp(-(x - y)*(x - y) / beta);
  return (return_value);
}

std::complex<double> nonlocal_harmonic_potential(double x, double y, void *params)
{
  std::complex<double> beta;
  beta = *(std::complex<double> *) params;
  std::complex<double> return_value;
  return_value = (x - y)*(x - y) * exp(-(x - y)*(x - y) / beta);
  return (return_value);
}

std::complex<double> nonlocal_potential(double x, double y, void *params)
{
  std::complex<double> beta;
  beta = *(std::complex<double> *) params;
  std::complex<double> return_value;
  return_value = - (x + y) * exp(-(abs(x - y) - 1)*(abs(x - y) - 1) / beta);
  return (return_value);
}
