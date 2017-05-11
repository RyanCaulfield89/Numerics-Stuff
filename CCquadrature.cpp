///////////////////////////////////////////////////////////////////////////////
//
//    File: CCquadrature.h.cpp
//    Description: See CCquadrature.h
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To Do List:
//      1 - This needs to be tested. I should make an example file tht does a
//          bunch of stuff. It should include solving an IDE, doing an integral
//          of a polynomial exactly, doing an integral of another function
//          accurately and representing a function as a series of polynomials.
///////////////////////////////////////////////////////////////////////////////

//Includes
#include "CCquadrature.h"
using namespace arma;
using namespace std;

///////////////////////////////////////////////////////////////////////////////
CCquadrature::CCquadrature(int n, double a, double b,
  complex<double>(*my_function)(double x, void *params), void *params){
  numpoints = n;
  left_boundary = a;
  right_boundary = b;
  integrand = my_function;
  parameters = params;
  find_points();
  find_weights();
  find_coefficients();
  construct_diff_matrix();
}

void CCquadrature::construct_diff_matrix(){
  diff_matrix = mat(numpoints + 1, numpoints + 1);
  //This formula is complicated. See Trefethen spectral methods in matlab,
  //page 53 Theorem 7.
  //First, do the corners.
  diff_matrix(0,0) = (2.0*numpoints + 1.0) / 6.0;
  diff_matrix(numpoints,numpoints) = -(2.0*numpoints + 1.0) / 6.0;
  diff_matrix(0,numpoints) = pow(-1.0,numpoints) / 2.0;
  diff_matrix(numpoints,0) = -pow(-1.0,numpoints) / 2.0;
  //Now, do the first and last column and row not including the corners.
  for(int j = 1; j < numpoints; j++){
    diff_matrix(0,j) = 2.0 * pow(-1.0,j) / (1.0 - points(j));
    diff_matrix(numpoints,j) = -2.0 * pow(-1.0,numpoints + j)/(1.0 + points(j));
  }
  for(int i = 1; i < numpoints; i++){
    diff_matrix(i,0) = -pow(-1.0,i) / (2.0 * (1.0 - points(i)));
    diff_matrix(i,numpoints) = pow(-1.0,numpoints + i)/(2.0 * (1.0 + points(i)));
  }
  //Now, do the interior.
  for(int i = 1; i < numpoints; i++){
    for(int j = 1; j < numpoints; j++){
      //Diagonal elements are different.
      if(i==j){
        diff_matrix(i,j) = -points(i) / (2.0*(1.0 - points(i)*points(i)));
      }
      else{
        diff_matrix(i,j) = pow(-1.0, i + j) / (points(i) - points(j));
      }
    }
  }
}

void CCquadrature::find_points(){
  //use extreme points of the Nth chebyshev polynomial.
  //Given by x_i = cos(Pi*i/N)
  points = vec(numpoints + 1);
  for(int i = 0; i <= numpoints; i++){
    points(i) = cos(double(i) * M_PI / numpoints);
  }
}

void CCquadrature::find_weights(){
  weights = vec(numpoints + 1);
  for(int i = 0; i <= numpoints; i++){
    weights(i) = M_PI / (2 * (nth_Tchebyshev_polynomial(numpoints, points(i))
                * diff_nth_Tchebyshev_polynomial(numpoints+1,points(i))));
  }
}

void CCquadrature::find_coefficients(){
  //Uses the formula c_i = (2/N)*sum(j,1,N,f(x_j)*T_i(x_j))
  //By convention, c_0 gets an extra factor of 1/2.
  coefficients = cx_vec(numpoints + 1);
  for(int i = 0; i <= numpoints; i++){
    for(int j = 1; j <= numpoints; j++){
      coefficients(i) += integrand(points(j),params) *
      nth_Tchebyshev_polynomial(i,points(j));
    }
    coefficients(i) *= 2.0/numpoints;
  }
  coefficients(0) *= 0.5;
}

double CCquadrature::nth_Tchebyshev_polynomial(int n, double x){
  //Exact formula defined piecewise. Look it up.
  if(abs(x)<=1){
    return cos(n * acos(x));
  }
  else if(x > 1){
    return cosh(n * acosh(x));
  }
  else if(x < -1){
    return pow(-1,n) * cosh(n * acosh(-x));
  }
  else{
    //should never get here.
    return 0;
  }
}

double CCquadrature::nth_Uchebyshev_polynomial(int n, double x){
  //Exact formula. Look it up.
  return (pow(x+sqrt(x*x-1),n+1) - pow(x-sqrt(x*x-1),n+1))/(2*sqrt(x*x-1));
}

double CCquadrature::diff_nth_Tchebyshev_polynomial(int n, double x){
  //uses dT_n/dx(x) = nU_n-1(x)
  //identity, look it up
  return n*nth_Uchebyshev_polynomial(n-1, x);
}

double CCquadrature::diff_nth_Uchebyshev_polynomial(int n, double x){
  //uses dU_n/dx(x) = ((n+1)*T_n+1(x)-xU_n(x))/(x^2-1)
  //identity, look it up
  return ((n+1)*nth_Tchebyshev_polynomial(n+1,x)
  - x*nth_Uchebyshev_polynomial(n,x))
  / (x*x - 1);
}

complex<double> CCquadrature::evaluate_integral(){
  double sum = 0;
  double x = 0;
  for(int i = 0; i <= numpoints; i++){
    //The quadrature rule works for integrals from -1 to 1 so
    //we need to change variables to go from a to b.
    x = (right_boundary - left_boundary) * points(i) / 2.0
        + (right_boundary + left_boundary) / 2.0;
    //The sqrt(1-x^2) factor is the inverse of the weighting
    //function for chebyshev polynomials.
    sum += sqrt(1-x*x) * integrand(x, params) * weights(i);
  }
  return sum;
}

complex<double> interpolate(double x){
  complex<double> return_value = 0;
  for(int i = 0; i < numpoints; i++){
    return_value += coefficients(i) * nth_Tchebyshev_polynomial(i,x);
  }
  return return_value;
}

int CCquadrature::get_numpoints(){
  return numpoints;
}

double CCquadrature::get_left_boundary(){
  return left_boundary;
}

double CCquadrature::get_right_boundary(){
  return right_boundary;
}

vec CCquadrature::get_points(){
  return points;
}

vec CCquadrature::get_weights(){
  return weights;
}

cx_vec CCquadrature::get_coefficients(){
  return coefficients;
}

mat CCquadrature::get_diff_matrix(){
  return diff_matrix;
}
