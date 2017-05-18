///////////////////////////////////////////////////////////////////////////////
//
//    File: CCquadrature_examples.cpp
//    Description: This is a series of examples showing what can be done with
//    the CCquadrature class. Use the make file make_CCquadrature_examples, but
//    be sure to change the path to the armadillo library to match yours.
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//
///////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include "functions.h"
#include "CCquadrature.h"
using namespace arma;
using namespace std;

///////////////////////////////////////////////////////////////////////////////

//Function headers
cx_vec polynomial_fit_coefficients(int numpoints, double lower, double upper,
  complex<double>(*my_function)(double x, void *params), void *params);

  complex<double> interpolating_polynomial(int numpoints, cx_vec coeffs,
    double point);

///////////////////////////////////////////////////////////////////////////////
int main(){
  //First, lets test that we are calculating the chebyshev polynomials correctly.
  cout << "First, let's calculate the chebyshev polynomials up to n = 4." << endl;
  cout << "Calculating. Saved to cheb_pol_test.dat. Run cheb_pol_test.plt" << endl;

  //These parameters aren't important for this first part but we need them to
  //create the object

  //Pick lower and upper bounds
  double lower = -1.0;
  double upper = 1.0;
  //parameter to pass to the gaussian
  complex<double> gamma = 1.0;
  void *test_params = &gamma;
  //Number of points
  int test_numpoints = 5;

  CCquadrature test_CCquadrature(test_numpoints, lower, upper, &gaussian, test_params);
  ofstream cheb_test_out("cheb_pol_test.dat");
  cheb_test_out << "#  x    T_0    T_1    T_2    T_3    T_4  " << endl;
  double z;
  for(int i = 0; i <= 2000; i++){
    z = lower + 0.001*i;
    cheb_test_out << z << "  ";
    for(int n = 0; n < 5; n++){
      cheb_test_out << test_CCquadrature.nth_Tchebyshev_polynomial(n,z) << "  ";
    }
    cheb_test_out << endl;
  }
  cheb_test_out.close();

  //Next, let's just take a gaussian at chebyshev nodes and try to interpolate
  //to get back the gaussian

  //parameter to pass to the gaussian
  complex<double> beta = 1.0;
  void *params = &beta;

  int numpoints;

  cout << "Let's take a gaussian sampled at Chebyshev nodes and try to interpolate it."
       << endl;
  cout << "We can compare it one sampled at evenly spaced points" << endl;
  cout << "Number of points to use?" << endl;
  cin >> numpoints;

  CCquadrature gaussian_CCquadrature(numpoints, lower, upper, &gaussian, params);
  //Once we make the object, the expansian coefficients are already calculated
  //Here we can calculate the coefficients for the evenly sampled polynomial fit.
  cx_vec coefficients = polynomial_fit_coefficients(numpoints, lower, upper, &gaussian, params);
  //Let's make an output file comparing the real gaussian to the interpolated one
  ofstream my_out("interpolated_gaussian.dat");
  my_out << "#  x    actual    cheb_interpolated    poly_interpolated" << endl;
  double x = 0.0;
  for(int i = 0; i < 1000; i++){
    x = lower + double(i)*(upper - lower)/1000.0;
    my_out << setprecision(5) << x << "  " << setprecision(16)
    << real(gaussian(x, params)) << "  "
    << real(gaussian_CCquadrature.interpolate(x)) << "  "
    << real(interpolating_polynomial(numpoints, coefficients, x)) << endl;
  }
  my_out.close();
  //Now, run interpolated_gaussian.plt to see a comparison.
  cout << "Saved data to interpolated_gaussian.dat" << endl;
  cout << "Now, run interpolated_gaussian.plt to see a comparison." << endl;

  //Now, lets do an integration with a quadratic polynomial
  cout << "Let's integrate a quadratic polynomial using CCquadrature" << endl;
  cout << "What should the coefficients a,b,c be?" << endl;

  //parameters for the quadratic polynomial
  quadratic_parameters Q_parameters;
  cout << "a = ";
  cin >> Q_parameters.a;
  cout << endl;
  cout << "b = ";
  cin >> Q_parameters.b;
  cout << endl;
  cout << "c = ";
  cin >> Q_parameters.c;
  cout << endl;
  void *quadratic_params;
  quadratic_params = &Q_parameters;

  //We only need 4 points to do this exactly
  numpoints = 4;
  CCquadrature quadratic_CCquadrature(numpoints, lower, upper,
                                      &quadratic_polynomial, quadratic_params);
  //print out points and weights
  cout << "Lets see the points and weights. We only need 4." << endl;
  for(int i = 0; i < numpoints; i++){
    cout << "x_" << i << " = " << quadratic_CCquadrature.get_points()(i)
         << "  w_" << i << " = " << quadratic_CCquadrature.get_weights()(i) <<endl;
  }

  complex<double> value = quadratic_CCquadrature.evaluate_integral();
  complex<double> exact = 2.0/3.0*Q_parameters.a + 2.0*Q_parameters.c;
  //report results to the user
  cout << "Integrating a quadratic from x=-1 to x=1 with" << endl;
  cout << "a = " << real(Q_parameters.a) << endl;
  cout << "b = " << real(Q_parameters.b) << endl;
  cout << "c = " << real(Q_parameters.c) << endl;
  cout << "This should yeild an exact value of" << endl;
  cout << "exact answer = " << setprecision(16) << real(exact) << endl;
  cout << "CCquadrature yields" << endl;
  cout << "CCquadrature value = " << setprecision(16) << real(value) << endl;

  cout << "On the next version, we will solve the schrodinger equation." << endl;
}

cx_vec polynomial_fit_coefficients(int numpoints, double lower, double upper,
  complex<double>(*my_function)(double x, void *params), void *params){
  //This just does a standard polynomial fit to a function at evenly spaced points.
  //It returns the coefficient of the polynomial which goes through the N points.

  //Set up evenly spaced points;
  vec points = vec(numpoints + 1);
  for(int i = 0; i <= numpoints; i++){
    points(i) = i * (upper - lower) / numpoints + lower;
  }
  //Set up a vector of function values at the points
  cx_vec function_values = cx_vec(numpoints + 1);
  for(int i = 0; i <= numpoints; i++){
    function_values(i) = my_function(points(i), params);
  }
  //This is the matrix A_ij = x_i^j
  mat interp_matrix = mat(numpoints + 1, numpoints + 1);
  for(int i = 0; i <= numpoints; i++){
    for(int j = 0; j <= numpoints; j++){
      interp_matrix(i,j) = pow(points(i),j);
    }
  }
  //Now, the coefficients are given by A^-1*function_values
  cx_vec coefficients = cx_vec(numpoints + 1);
  coefficients = inv(interp_matrix) * function_values;
  return coefficients;
}

complex<double> interpolating_polynomial(int numpoints, cx_vec coeffs,
  double point){
  complex<double> return_value = 0;
  for(int i = 0; i <= numpoints; i++){
    return_value += coeffs(i)*pow(point,i);
  }
  return return_value;
}
