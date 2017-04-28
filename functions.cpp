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
using namespace std;

///////////////////////////////////////////////////////////////////////////////

double coulomb_potential(double x, void *params)
{
  double charge;
  charge = *(double *) params.q;
  double return_value;
  return_value = charge / x;
  return (return_value);
}

double harmonic_potential(double x, void *params)
{
  double k;
  k = *(double*) params;
  double return_value = 0.5*k*x*x;
}
