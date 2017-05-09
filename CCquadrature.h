///////////////////////////////////////////////////////////////////////////////
//
//    File: CCquadrature.h
//    Description: Header file for CCquadrature.cpp.
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
///////////////////////////////////////////////////////////////////////////////
#ifndef ARMADILLOCCQUADRATURE_H
#define ARMADILLOCCQUADRATURE_H //To avoid including this file twice

//Include files
#include <armadillo>  //This package does all the linear algebra stuff.
#include <complex>    //Complex number stuff.
#include <cmath>

class CCquadrature
{
  public:
    //Constructs the object and initializes the # of points and the boundaries.
    CCquadrature(int n, double a, double b);
    //Destructor.
    ~CCquadrature();

    //This makes the differentiation matrix.
    void construct_diff_matrix();
    //This finds the chebyshev nodes.
    void find_points();
    //This finds the weight for each node.
    void find_weights();
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


    //Getters for the important stuff.
    int get_numpoints();
    double get_left_boundary();
    double get_right_boundary();
    arma::vec get_points();
    arma::vec get_weights();
    arma::mat get_diff_matrix();

  private:
    int numpoints             //The number of points used in the quadrature
    double left_boundary;     //Left boundary or lower bound
    double right_boundary;    //Right boundary or upper bound
    arma::vec points;         //Quadrature points
    arma::vec weights;        //Quadrature weights
    arma::mat diff_matrix;    //Differentiation matrix
}

#endif
