///////////////////////////////////////////////////////////////////////////////
//
//    File: GLquadrature.cpp
//    Description: See GLquadrature.h
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To Do List:
//      1 - This needs to be tested.
//
///////////////////////////////////////////////////////////////////////////////

//Includes
#include "GLquadrature.h"
using namespace arma;
using namespace std;

///////////////////////////////////////////////////////////////////////////////
GLquadrature::GLquadrature(int n, double a, double b,
  complex<double>(*my_function)(double x, void *params), void *params){
  numpoints = n-1;
  left_boundary = a;
  right_boundary = b;
  integrand = my_function;
  parameters = params;
  find_points();
  find_weights();
  find_coefficients();
}

GLquadrature::~GLquadrature(){

}

void GLquadrature::find_points(){

}

void GLquadrature::find_weights(){

}

void GLquadrature::find_coefficients(){
  //First, take the function values at the legendre nodes.
  cx_vec function_values = cx_vec(numpoints + 1);
  for(int i = 0; i <= numpoints; i++){
    function_values(i) = integrand(points(i), parameters);
  }
  //Now, this is the matrix A_ij = x_i^j
  mat interp_matrix = mat(numpoints + 1, numpoints + 1);
  for(int i = 0; i <= numpoints; i++){
    for(int j = 0; j <= numpoints; j++){
      interp_matrix(i,j) = pow(points(i),j);
    }
  }
  //Now, the coefficients are given by solving Ac = f
  //c = A^-1*f
  coefficients = cx_vec(numpoints + 1);
  coefficients = inv(interp_matrix)*function_values;
}

double GLquadrature::nth_legendre_polynomial(int n, double x){

}

double GLquadrature::diff_nth_legendre_polynomial(int n, double x){

}

complex<double> GLquadrature::evaluate_integral(){
  complex<double> sum = 0;
  double x = 0;
  for(int i = 0; i <= numpoints; i++){
    //The quadrature rule works for integrals from -1 to 1 so
    //we need to change variables to go from a to b.
    x = (right_boundary - left_boundary) * points(i) / 2.0
        + (right_boundary + left_boundary) / 2.0;
    //The sqrt(1-x^2) factor is the inverse of the weighting
    //function for chebyshev polynomials.
    sum += integrand(x, parameters) * weights(i);
  }
  sum *= (right_boundary - left_boundary) / 2.0;
  return sum;
}

complex<double> GLquadrature::interpolate(double x){
  complex<double> return_value = 0.0;
  for(int i = 0; i < numpoints; i++){
    return_value += coefficients(i) * pow(x,i);
  }
  return return_value;
}

int GLquadrature::get_numpoints(){
  return numpoints;
}

double GLquadrature::get_left_boundary(){
  return left_boundary;
}

double GLquadrature::get_right_boundary(){
  return right_boundary;
}

vec GLquadrature::get_points(){
  return points;
}

vec GLquadrature::get_weights(){
  return weights;
}

cx_vec GLquadrature::get_coefficients(){
  return coefficients;
}
