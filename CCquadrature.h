///////////////////////////////////////////////////////////////////////////////
//
//    File: CCquadrature.h
//    Description: Header file for CCquadrature.cpp.
//    This uses Clenshaw-Curtis dicretization to approximate integrals and
//    derivatives. Idealy, it should allow you to approximate integrals using
//    very few points and should return an exact value for polynomials of
//    degree N-1 or less where N is the number of points. It should also
//    be able to solve integro-differential equations, since you can discretize
//    an integral into a sum while also discretizing a derivative into a matrix.
//    The form of the IDE will be left to the user to implement in a seperate
//    file. In addition, this will also get the expansion coefficients for a
//    Chebyshev polynomial expansion of a function. This might be less useful if
//    you already have the function explicitly. This uses the Armadillo package.
//    You can find documentation at http://arma.sourceforge.net/docs.html
//    See CCquadrature_examples for examples of how this works.
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
///////////////////////////////////////////////////////////////////////////////
#ifndef CCQUADRATURE_H
#define CCQUADRATURE_H //To avoid including this file twice

//Include files
#include <armadillo>  //This package does all the linear algebra stuff.
#include <complex>    //Complex number stuff.
#include <cmath>

class CCquadrature
{
  public:
    //Constructs the object and initializes the # of points and the boundaries.
    CCquadrature(int n, double a, double b,
      std::complex<double>(*my_function)(double x, void *params), void *params);
    //Destructor
    ~CCquadrature();


    //This makes the differentiation matrix.
    void construct_diff_matrix();
    //This finds the chebyshev nodes.
    void find_points();
    //This finds the weight for each node.
    void find_weights();
    //This find the expansion coefficients for the polynomial interpolation
    //You can then recreate the function using f(x) = sum(i,0,N, c_i*x^i)
    void find_coefficients();
    //Calculates the nth chebyshev polynomial of the 1st kind at x.
    double nth_Tchebyshev_polynomial(int n, double x);
    //Calculates the nth chebyshev polynomial of the 2nd kind at x.
    double nth_Uchebyshev_polynomial(int n, double x);
    //Calculates the derivative of the nth chebyshev polynomial of
    //the 1st kind at x.
    double diff_nth_Tchebyshev_polynomial(int n, double x);
    //Calculates the derivative of the nth chebyshev polynomial of
    //the 2nd kind at x.
    double diff_nth_Uchebyshev_polynomial(int n, double x);
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
    arma::mat get_diff_matrix();

  private:
    int numpoints;            //The number of points used in the quadrature
    double left_boundary;     //Left boundary or lower bound
    double right_boundary;    //Right boundary or upper bound
    arma::vec points;         //Quadrature points
    arma::vec weights;        //Quadrature weights
    arma::cx_vec coefficients;//Coefficients in the interpolation expansion
    arma::mat diff_matrix;    //Differentiation matrix
    std::complex<double>(*integrand)(double x, void *params);
                              //The function to be integrated or interpolated
    void *parameters;         //parameters for the potentials
};

#endif
