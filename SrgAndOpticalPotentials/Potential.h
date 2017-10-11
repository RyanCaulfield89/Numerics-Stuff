///////////////////////////////////////////////////////////////////////////////
//
//    File: Potential.h
//    Description:
//
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To Do List:
//      1 -
//
///////////////////////////////////////////////////////////////////////////////

#ifndef POTENTIAL_H
#define POTENTIAL_H //To avoid including this file twice

//Includes
#include <armadillo>  //This package does all the linear algebra stuff.
#include <complex>    //Complex number stuff.
#include <cmath>

Class Potential
{
  public:
    //Constructor for local potentials
    Potential(double q_min, double q_max, int num_points,
      std::complex<double>(*my_potential)(double q, void *params), void *params);

    //Overloaded constructor for nonlocal potentials
    Potential(double q_min, double q_max, int num_points,
      std::complex<double>(*my_potential)(double q1, double q2, void *params),
      void *params);

    //Destructor
    ~Potential();

    //Creates a matrix by discretizing the local potential. This will be a
    //diagonal matrix since the potential is local.
    void construct_local_matrix_rep();

    //Creates a matrix by discretizing the nonlocal potential. This will have
    //nonzero off-diagonal elements.
    void construct_nonlocal_matrix_rep();

  private:
    double q_min;             //The lower value of q used for discretization

    double q_max;             //The upper value of q used for discretization

    int num_points;           //The number of points used for discretization

    arma::mat potential_matrix;
                              //Matrix representation of the potential

    std::complex<double>(*local_potential)(double q, void *params);
                              //A local Potential

    std::complex<double>(*nonlocal_potential)(double q1, double q2, void *params);
                              //A nonlocal Potential

    void *parameters;         //parameters for the potentials

}
#endif
