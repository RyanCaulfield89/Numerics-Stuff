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
int main(){
  //First, let's just take a gaussian at chebyshev nodes and try to interpolate
  //to get back the gaussian

  //Pick lower and upper bounds
  double lower = -1.0;
  double upper = 1.0;

  //parameter to pass to the gaussian
  complex<double> beta = 1.0;
  void *params = &beta;

  int numpoints;

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
