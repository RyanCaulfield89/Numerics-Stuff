//////////////////////////////////////////////////////////////////////////////
//
//    File: ode_solver.cpp
//    Description: This is a collection of ode solvers
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To do list:
//      1 - This needs to be tested.
//      2 - Add a switch statement which allows for different algorithyms.
//      3 - make a seperate header file if this starts getting long.
//
///////////////////////////////////////////////////////////////////////////////

//Includes
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <cmath>
using namespace std;

//Headers

// ode solver that uses the rk45 method
double[] ode_solve(double x_initial, y_initial[], double x, int dimension,
  double step_size, void *params_ptr,
  int (*rhs) (double , const double y[], double *dfdy,
  	  double dfdt[], void *params_ptr),
  int (*jacobian) (double , const double y[], double *dfdy,
  	  double dfdt[], void *params_ptr));

///////////////////////////////////////////////////////////////////////////////
double[] ode_solve(double x_initial, y_initial[], double x, int dimension,
  double step_size, void *params_ptr,
  int (*rhs) (double , const double y[], double *dfdy,
  	  double dfdt[], void *params_ptr),
  int (*jacobian) (double , const double y[], double *dfdy,
  	  double dfdt[], void *params_ptr))
{
  //This will be the return value.
  double y[dimension];
  //Set the initial values
  for(int i = 0; i < dimension; i++){
    y[i] = y_initial[i];
  }

  //These are needed for the gsl routine
  const double eps_abs = 1.e-8;	        // absolute error requested
  const double eps_rel = 1.e-10;	// relative error requested

  //This is the gsl type pointer that specifies the ode algorithym to use.
  //Later we can add a switch statment to allow for different algorithyms.
  const gsl_odeiv_step_type *type_ptr = gsl_odeiv_step_rkf45;

  //gsl objects
  gsl_odeiv_step *step_ptr = gsl_odeiv_step_alloc (type_ptr, dimension);
  gsl_odeiv_control *control_ptr = gsl_odeiv_control_y_new (eps_abs, eps_rel);
  gsl_odeiv_evolve *evolve_ptr = gsl_odeiv_evolve_alloc (dimension);

  //This is the ode system that will be solved. It needs the parameters,
  //the rhs, the number of dimensions and the jacobian.
  gsl_odeiv_system my_system;
  my_system.function = rhs;	// the right-hand-side functions dy[i]/dt
  my_system.jacobian = jacobian;	// the Jacobian df[i]/dy[j]
  my_system.dimension = dimension;	// number of diffeq's
  my_system.params = params_ptr;	// parameters to pass to rhs and jacobian

  //This is the step size in x
  double delta_x = (x - x_initial) / 10000;

  //This is the x we will itterate ode_solver, initially x_initial
  double x_itt = x_initial;

  //Here's where we actually do the work
  for (double x_next = x_initial + delta_x; x_next <= x; x_next += delta_x)
  {
    while (x_itt < x_next)	// evolve from t to t_next
	  {
	    gsl_odeiv_evolve_apply (evolve_ptr, control_ptr, step_ptr,
			                    	  &my_system, &x_itt, x_next, &step_size, y);
	  }
  }

  //Free the gsl objects
  gsl_odeiv_evolve_free (evolve_ptr);
  gsl_odeiv_control_free (control_ptr);
  gsl_odeiv_step_free (step_ptr);

  //Finally, return y[]
  return (y[]);
}
