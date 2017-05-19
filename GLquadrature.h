///////////////////////////////////////////////////////////////////////////////
//
//    File: GLquadrature.h
//    Description: Header file for GLquadrature.cpp.
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
///////////////////////////////////////////////////////////////////////////////
#ifndef GLQUADRATURE_H
#define GLQUADRATURE_H //To avoid including this file twice

//Include files
#include <armadillo>  //This package does all the linear algebra stuff.
#include <complex>    //Complex number stuff.
#include <cmath>

class GLquadrature
{
  public:
    //Constructs the object and initializes the # of points and the boundaries.
    GLquadrature(int n, double a, double b,
      std::complex<double>(*my_function)(double x, void *params), void *params);
    //Destructor
    ~GLquadrature();


    //This finds the legendre nodes.
    void find_points();
    //This finds the weight for each node.
    void find_weights();
    //This find the expansion coefficients for the polynomial interpolation
    //You can then recreate the function using f(x) = sum(i,0,N, c_i*x^i)
    void find_coefficients();
    //Calculates the nth chebyshev polynomial of the 1st kind at x.
    double nth_legendre_polynomial(int n, double x);
    //Calculates the derivative of the nth chebyshev polynomial of
    //the 1st kind at x.
    double diff_nth_legendre_polynomial(int n, double x);
    //This evaluated the integral using the quadrature rule.
    std::complex<double> evaluate_integral();
    //This approximates the function at x by a polynomial.
    std::complex<double> interpolate(double x);

    //Getters for the important stuff.
    int get_numpoints();
    double get_left_boundary();
    double get_right_boundary();
    arma::vec get_points();
    arma::vec get_weights();
    arma::cx_vec get_coefficients();

  private:
    int numpoints;            //The number of points used in the quadrature
    double left_boundary;     //Left boundary or lower bound
    double right_boundary;    //Right boundary or upper bound
    arma::vec points;         //Quadrature points
    arma::vec weights;        //Quadrature weights
    arma::cx_vec coefficients;//Coefficients in the interpolation expansion
    std::complex<double>(*integrand)(double x, void *params);
                              //The function to be integrated or interpolated
    void *parameters;         //parameters for the potentials
};

#endif
