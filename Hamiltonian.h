///////////////////////////////////////////////////////////////////////////////
//
//    File: Hamiltonian.h
//    Description: Header file for Hamiltonian.cpp
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
  //builds the matrix from the potential in position space and solves the
  //eigensystem.
  Hamiltonian (const int dim, double(*potential)(double x, void *params));
  //This does the same as the previous constructor but for a nonlocal
  //potential in x-space.
  Hamiltonian::Hamiltonian (const int dim,
    double(*potential)(double x1, double x2, void *params));

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
  arma::mat Hmatrix;           // The matrix form of the Hamiltonian
  arma::vec eigenvalues;    // vector of eigenvalues
  arma::mat eigenvectors;   // matrix of eigenvectors
  double (*potential) (double x, void *params) xpotential; //potential in xspace
  double (*potential) (double k, void *params) kpotential; //potential in kspace
  double (*potential) (double x1, double x2, void *params) xnonLocalPotential;
  //Non-Local potnetial in x-space
  double (*potential) (double k1, double k2, void *params) knonLocalPotential;
  //Non-Local potnetial in k-space
}

#endif
