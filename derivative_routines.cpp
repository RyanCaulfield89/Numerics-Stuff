//////////////////////////////////////////////////////////////////////////////
//
//    File: derivative_routines.cpp
//    Description: This is a set of tools for numerically differentiating
//    a function.
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To do list:
//      1 - This needs to be tested.
//      2 - Find other derivative methods that will be useful and add them.
//      3 - make a seperate header file if this starts getting long.
//
///////////////////////////////////////////////////////////////////////////////

//Includes
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_diff.h>
using namespace std;

//Headers

//Wraps a gsl adaptive central difference method for taking a derivative. This
//adaptive so need for a step size h.
double gsl_adaptive_derivative(double x, void *params_ptr,
  double (*function) (double x,  void *params));

///////////////////////////////////////////////////////////////////////////////

double gsl_adaptive_derivative(double x, void *params_ptr,
  double (*function) (double x,  void *params))
{
  double derivative_value;  //This will be the return value
  double abserr;            //This is the absolute error. It's needed for the
                            //gsl routine. We won't use it, but could implement
                            //something later that uses it.
  //The function needs to be defined this way for the gsl routine
  gsl_function My_function;
  My_function.function = function;
  My_function.params = params_ptr;

  //Here we actually do the calculation
  gsl_diff_central (&My_function, x, &derivative_value, &abserr);

  return derivative_value;
}
