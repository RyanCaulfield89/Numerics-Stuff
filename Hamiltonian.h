///////////////////////////////////////////////////////////////////////////////
//
//    File: Hamiltonian.h
//    Description: Header file for Hamiltonian.cpp.
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
///////////////////////////////////////////////////////////////////////////////

#ifndef ARMADILLOHAMILTONIAN_H
#define ARMADILLOHAMILTONIAN_H //To avoid including this file twice

//Include files
#include <armadillo>  //This package does all the linear algebra stuff.
#include <complex>    //Complex number stuff.
#include <fstream>    //stuff for saving to files
#include <iostream>   //stuff for saving to files
#include <iomanip>    //stuff for saving to files

class Hamiltonian
{
  public:
    //A constructor that just initializes the dimension of H and sets all it'
    //elements to zero. You need to set the elements using set_element(i;j)
    Hamiltonian (const int dim);
    //A constructor which sets the dimensionality and the potential. It then
    //builds the matrix from the potential in xspace or kspace.
    //The potential_type should be "x" or "k".
    Hamiltonian (double Rmin, double Rmax, double h, const char potential_type,
      std::complex<double>(*potential)(double x, void *params), void *params);
    //This does the same as the previous constructor but for a nonlocal
    //potential in x-space or k-space. The potential_type should be "x" or "k".
    Hamiltonian (double Rmin, double Rmax, double h, const char potential_type,
      std::complex<double>(*potential)(double x1, double x2, void *params), void *params);
    //This constructor will be for 1-D scattering using the born approximation.
    //It constructs (E - H_0 + iepsilon)^(-1)V and has it act on a initial
    //plane wave state. Just as in the previous, the potential type should be
    //"x" or "k"
    Hamiltonian (double energy, double Rmin, double Rmax, double h,
      const char potential_type,
      std::complex<double>(*potential)(double x, void *params), void *params);

    ~Hamiltonian ();  // destructor

    //These construct the hamiltonian matrix for x or k both local and nonlocal.
    //These utilize dirichlet boundary conditions.
    void construct_localXmatrix();
    void construct_localKmatrix();
    void construct_nonlocalXmatrix();
    void construct_nonlocalKmatrix();

    //Constructs the scattering matrix for 1D scattering.
    void construct_scattering_matrix();
    //constructs a plane wave initial wave function.
    //used in 1D scattering.
    void construct_initial_wf();

    //This solves for the eigenvectors and eigenvalues. After, the can be
    //accesed with the appropriate getter function.
    void solve_eigensystem();

    //getters and setters
    void save_eigenvector(const int i, const char* filename);
    void set_element(const int i, const int j, const std::complex<double> value);
    std::complex<double> get_element(const int i, const int j);
    std::complex<double> get_eigenvalue(const int i);
    std::complex<double> get_eigenvector(const int i, const int j);
    double calculate_T();
    double calculate_R();

  private:
    double left_boundary;     // The minimum value of x
    double right_boundary;    // maximum value of x
    int dimension;            // Dimensionality
    double step_size;         // The step size used for FD approximation of
                              // the derivative.
    double initial_state_energy;//The energy of the initial wave function used
                              //in 1D scattering
    arma::cx_mat Hmatrix;     // The matrix form of the Hamiltonian
    arma::cx_mat scattering_matrix; // (E - H_0 + iepsilon)^(-1)V
                                    // used for 1D scattering
    arma::vec eigenvalues;    // vector of eigenvalues
    arma::cx_mat eigenvectors;// matrix of eigenvectors
    arma::cx_vec initial_wf;  // The initial wave function used in 1D scattering
    arma::cx_vec scattered_wf;// The scattered wave function used in 1D scattering
    std::complex<double>(*xpotential)(double x, void *params);
                              //potential in xspace
    std::complex<double>(*kpotential)(double k, void *params);
                              //potential in kspace
    std::complex<double>(*xnonLocalPotential)(double x1, double x2, void *params);
                              //Non-Local potnetial in x-space
    std::complex<double>(*knonLocalPotential)(double k1, double k2, void *params);
                              //Non-Local potnetial in k-space
    void *parameters;         //parameters for the potentials
};

#endif
