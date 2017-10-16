///////////////////////////////////////////////////////////////////////////////
//
//    File: Hamiltonian.h
//    Description:
//
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To Do List:
//      1 -
//
///////////////////////////////////////////////////////////////////////////////

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H //To avoid including this file twice

//Includes
#include <armadillo>  //This package does all the linear algebra stuff.
#include <complex>    //Complex number stuff.
#include <cmath>      //Math stuff
#include <CCquadrature.h> //This is use for discretizing the Hamiltonian

Class Hamiltonian
{
  public:
    //Constructor for local potentials. The potential_type should be "x" or
    //"k", q_min and q_max the left and right bounaries, num_points is the number
    //of points used in discretization and my_potential is the potential energy.
    //Here, the potential is local.
    Hamiltonian(double q_min, double q_max, int num_points, const char potential_type,
      std::complex<double>(*my_potential)(double q, void *params), void *params);

    //Constructor for nonlocal potentials. The potential_type should be "x" or
    //"k", q_min and q_max the left and right bounaries, num_points is the number
    //of points used in discretization and my_potential is the potential energy.
    //Here, the potential is non-local.
    Hamiltonian(double q_min, double q_max, int num_points, const char potential_type,
      std::complex<double>(*my_potential)(double q1, double q2, void *params),
      void *params);

    //Destructor
    ~Hamiltonian();

    //Creates a matrix by discretizing the local potential. This will be a
    //diagonal matrix since the potential is local.
    void construct_local_xmatrix_rep();
    void construct_local_kmatrix_rep();

    //Creates a matrix by discretizing the nonlocal potential. This will have
    //nonzero off-diagonal elements.
    void construct_nonlocal_xmatrix_rep();
    void construct_nonlocal_kmatrix_rep();

    //Getters
    arma::mat get_hamiltonian_matrix();
    arma::mat get_kinetic_matrix();
    arma::mat get_potential_matrix();
    arma::vec get_points();
    arma::vec get_weights();


  private:
    double q_min;             //The lower value of q used for discretization

    double q_max;             //The upper value of q used for discretization

    int num_points;           //The number of points used for discretization

    arma::vec points;
                              //Quadrature points used in discretization

    arma::vec weights;        //Quadrature weights used in discretization

    arma::mat kinetic_matrix
                              //Matrix of the kinetic term in the Hamiltonian

    arma::mat potential_matrix;
                              //Matrix representation of the potential

    arma::mat hamiltonian_matrix;
                              //Matrix representation of the Hamiltonian

    std::complex<double>(*local_potential)(double q, void *params);
                              //A local Potential

    std::complex<double>(*nonlocal_potential)(double q1, double q2, void *params);
                              //A nonlocal Potential

    void *parameters;         //parameters for the potentials

}
#endif
