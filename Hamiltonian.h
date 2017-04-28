///////////////////////////////////////////////////////////////////////////////
//
//    File: Hamiltonian.h
//    Description: Header file for Hamiltonian.cpp. This contains information
//                 about potentials and finite difference matrix forms. The
//                 discretization is done with uniform spacing for a step
//                 size which is also stored as a parameter. Using these methods
//                 we can solve for eigenstates and eigenvalues.
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
///////////////////////////////////////////////////////////////////////////////

#ifndef ARMADILLOHAMILTONIAN_H
#define ARMADILLOHAMILTONIAN_H //To avoid including this file twice

//Include files
#include <armadillo>  //This package does all the linear algebra stuff.

class Hamiltonian
{
public:
  //A constructor that just initializes the dimension of H and sets all it'
  //elements to zero. You need to set the elements using set_element(i;j)
  Hamiltonian (const int dim);
  //A constructor which sets the dimensionality and the potential. It then
  //builds the matrix from the potential in xspace or kspace.
  //The potential_type should be "x" or "k".
  Hamiltonian (const int dim, , double h, const char potential_type,
    double(*potential)(double x, void *params), void *params);
  //This does the same as the previous constructor but for a nonlocal
  //potential in x-space or k-space. The potential_type should be "x" or "k".
  Hamiltonian::Hamiltonian (const int dim,, double h, const char potential_type,
    double(*potential)(double x1, double x2, void *params), void *params);

  ~Hamiltonian ();  // destructor

  //These construct the hamiltonian matrix for x or k both local and nonlocal.
  void construct_localXmatrix();
  void construct_localKmatrix();
  void construct_nonlocalXmatrix();
  void construct_nonlocalKmatrix();

  //This solves for the eigenvectors and eigenvalues. After, the can be
  //accesed with the appropriate getter function.
  void solve_eigensystem();

  //getters and setters
  void set_element(const int i, const int j, const double value);
  double get_element(const int i, const int j);
  double get_eigenvalue(const int i);
  double get_eigenvector(const int i, const int j);

private:
  int dimension;            // Dimensionality
  double step_size          // The step size used for FD approximation of
                            // the derivative.
  arma::mat Hmatrix;        // The matrix form of the Hamiltonian
  arma::vec eigenvalues;    // vector of eigenvalues
  arma::mat eigenvectors;   // matrix of eigenvectors
  double (*xpotential) (double x, void *params); //potential in xspace
  double (*kpotential) (double k, void *params) ; //potential in kspace
  double (*xnonLocalPotential) (double x1, double x2, void *params);
                            //Non-Local potnetial in x-space
  double (*knonLocalPotential) (double k1, double k2, void *params);
                            //Non-Local potnetial in k-space
  void *parameters              //parameters for the potentials
}

#endif
