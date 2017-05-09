///////////////////////////////////////////////////////////////////////////////
//
//    File: CCquadrature.h.cpp
//    Description: See CCquadrature.h
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To Do List:
//      1 - This needs to be tested
///////////////////////////////////////////////////////////////////////////////

//Includes
#include "CCquadrature.h"
using namespace arma;
using namespace std;

///////////////////////////////////////////////////////////////////////////////
CCquadrature::CCquadrature(int n, double a, double b){
  numpoints = n;
  left_boundary = a;
  right_boundary = b;
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
  //Now, do the first and last column and row.
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

double nth_Tchebyshev_polynomial(int n, double x){
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

double nth_Uchebyshev_polynomial(int n, double x){
  return (pow(x+sqrt(x*x-1),n+1) - pow(x-sqrt(x*x-1),n+1))/(2*sqrt(x*x-1));
}

double diff_nth_Tchebyshev_polynomial(int n, double x){
  //uses dT_n/dx(x) = nU_n-1(x)
  //identity, look it up
  return n*nth_Uchebyshev_polynomial(n-1, x);
}

double diff_nth_Uchebyshev_polynomial(int n, double x){
  //uses dU_n/dx(x) = ((n+1)*T_n+1(x)-xU_n(x))/(x^2-1)
  //identity, look it up
  return ((n+1)*nth_Tchebyshev_polynomial(n+1,x)
  - x*nth_Uchebyshev_polynomial(n,x))
  / (x*x - 1);
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

arma::vec CCquadrature::get_points(){
  return points;
}

arma::vec CCquadrature::get_weights(){
  return weights;
}

arma::mat CCquadrature::get_diff_matrix(){
  return diff_matrix;
}
