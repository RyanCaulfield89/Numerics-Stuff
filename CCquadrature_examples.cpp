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

//Define a structure for the quadratic polynomial
typedef struct
{
  complex<double> a;
  complex<double> b;
  complex<double> c;
}
quadratic_parameters;

///////////////////////////////////////////////////////////////////////////////
int main(){
  //First, let's just take a gaussian at chebyshev nodes and try to interpolate
  //to get back the gaussian

  //Pick lower and upper bounds
  double lower = -1.0;
  double upper = 1.0;

  //parameter to pass to the gaussian
  complex<double> beta = 1.0;
  void *params = &beta;

  cout << "Let's take a gaussian sampled at Chebyshev nodes and try to interpolate it."
       << endl;
  cout << "Number of points to use?" << endl;
  cin >> numpoints;

  CCquadrature gaussian_CCquadrature(numpoints, lower, upper, &gaussian, params);
  //Once we make the object, the expansian coefficients are already calculated

  //Let's make an output file comparing the real gaussian to the interpolated one
  ofstream my_out("interpolated_gaussian.dat");
  my_out << "#  x    actual    interpolated" << endl;
  double x = 0.0;
  for(int i = 0; i < 1000; i++){
    x = lower + double(i)*(upper - lower)/1000.0;
    my_out << setprecision(5) << x << "  " << setprecision(16)
    << real(gaussian(x, params)) << "  "
    << real(gaussian_CCquadrature.interpolate(x)) << endl;
  }
  my_out.close();
  //Now, run interpolated_gaussian.plt to see a comparison.
  cout << "Saved data to interpolated_gaussian.dat" << endl;
  cout << "Now, run interpolated_gaussian.plt to see a comparison." << endl;

  //Now, lets do an integration with a quadratic polynomial
  cout << "Let's integrate a quadratic polynomial using CCquadrature" << endl;
  //parameters for the quadratic polynomial
  void *quadratic_params;
  quadratic_parameters Q_parameters;
  Q_parameters.a = 1.0;
  Q_parameters.b = 2.0;
  Q_parameters.c = 1.0;

  //We only need 3 points to do this exactly
  numpoints = 3;
  CCquadrature quadratic_CCquadrature(numpoints, lower, upper,
                                      &quadratic_polynomial, quadratic_params);
  complex<double> value = quadratic_CCquadrature.evaluate_integral();
  complex<double> exact = 2.0/3.0*Q_parameters.a + 2.0*Q_parameters.c;
  //report results to the user
  cout << "Integrating a quadratic from x=-1 to x=1 with" << endl;
  cout << "a = " << Q_parameters.a << endl;
  cout << "b = " << Q_parameters.b << endl;
  cout << "c = " << Q_parameters.c << endl;
  cout << "This should yeild an exact value of" << endl;
  cout << "exact answer = " << exact << endl;
  cout << "CCquadrature yields" << endl;
  cout << "CCquadrature value = " << value << endl;

  cout << "On the next version, we will solve the schrodinger equation." << endl;
}
