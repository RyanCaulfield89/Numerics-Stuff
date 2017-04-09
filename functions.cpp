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

//Headers

//Coulomb potential. The parameter pointer needs to have an attribute q, which
//is the charge of the particle in units of 4*Pi*epsilon_0.
double coulomb_potential(double x, void *params);

///////////////////////////////////////////////////////////////////////////////

double coulomb_potential(double x, void *params)
{
  double charge;
  charge = *(double *) params.q;
  double return_value;
  return_value = charge / x;
  return (return_value);
}
