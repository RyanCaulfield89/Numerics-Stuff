///////////////////////////////////////////////////////////////////////
//
//    File: integration.cpp
//    Description: This is a set of tools for numerically integrating
//    a function.
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To do list:
//      1 - Find other integration methods will be useful and add them.
//      2 - Is 1000 the right size for workspace?
///////////////////////////////////////////////////////////////////////

//includes
#include <cmath>
#include <gsl/gsl_integration.h>
using namespace std;

//Headers

//Wraps a gsl integration method so that our code looks cleaner.
//This will use an adaptive gsl method so no need for a specified
//number of points.
double my_gsl_integration(double lower, double upper,
    double (*integrand) (double x,  void *params));

///////////////////////////////////////////////////////////////////////

double my_gsl_integration(double lower, double upper, void *params_ptr
    double (*integrand) (double x,  void *params))
{
  double abs_error = 1.0e-12;
  double rel_error = 1.0e-12;
  double result;		// This will be the return value
  double error;

  //We need to allocate a gsl workspace for the gsl method
  gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc(1000);

  //gsl take functions in a very specific format. It needs a gsl_function
  //object with a function and a void pointer of parameters as attributes.
  gsl_function My_function;
  My_function.function = integrand;
  My_function.params = params_ptr;

  //This is where the integration is happening. Google gsl qags if
  //you want to understand more about the method and what it does.
  gsl_integration_qags (&My_function, lower, upper, abs_error,
    rel_error, 1000, work_ptr, &result, &error);

  //This frees up the memory used by the gsl workspace.
  gsl_integration_workspace_free(work_ptr);

  //Finally, the result of the integration was stored in result.
  return (result);
}
