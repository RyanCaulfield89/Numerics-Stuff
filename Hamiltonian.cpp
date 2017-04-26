///////////////////////////////////////////////////////////////////////////////
//
//    File: Hamiltonian.cpp
//    Description:
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To Do List:
//      1 -
///////////////////////////////////////////////////////////////////////////////

//Includes
#include <Hamiltonian.h>
using namespace arma;

///////////////////////////////////////////////////////////////////////////////
Hamiltonian::Hamiltonian (const int dim)
{
  dimension = dim;
  Hmatrix = mat(dimension,dimension);
  Hmatrix.zeroes();
}

Hamiltonian::~Hamiltonian () // Destructor for Hamiltonian
{
  //Doesn't do anything
}

//Constructor for an x-space local potential
Hamiltonian::Hamiltonian (const int dim,
  double(*potential)(double x, void *params))
{
  dimension = dim;
  xpotential = potential;
  construct_Hmatrix();
  solve_eigensystem();
}

void construct_localXmatrix()
{

}

void solve_eigensystem()
{
  eig_sym(eigenvalues, eigenvectors, Hmatrix);
}

void set_element(const int i, const int j, const double value)
{
  Hmatrix(i-1,j-1) = value;
}

double get_element(const int i, const int j)
{
  return Hmatrix(i-1,j-1);
}

double get_eigenvalue(const int i)
{
  return eigenvalues (i-1);
}

double get_eigenvector(const int i, const int j)
{
  return eigenvectors (j-1, i-1);
}
