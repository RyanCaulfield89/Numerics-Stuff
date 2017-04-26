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
Hamiltonian::Hamiltonian (const int dim, const char potential_type,
  double(*potential)(double x, void *params))
{
  dimension = dim;
  if(potential_type == "x")
  {
    xpotential = potential;
    construct_localXmatrix();
  }
  if(potential_type == "k")
  {
    kpotential = potential;
    construct_localKmatrix();
  }
  else
  {
    return();
  }
  solve_eigensystem();
}

//Constructor for an x-space nonlocal potential
Hamiltonian::Hamiltonian (const int dim, const char potential_type,
  double(*potential)(double x1, double x2, void *params))
{
  dimension = dim;
  if(potential_type == "x")
  {
    xnonLocalpotential = potential;
    construct_nonlocalXmatrix();
  }
  if(potential_type == "k")
  {
    knonLocalpotential = potential;
    construct_nonlocalKmatrix();
  }
  else
  {
    return();
  }
  solve_eigensystem();
}

void construct_localXmatrix()
{

}

void construct_localKmatrix()
{

}

void construct_nonlocalXmatrix()
{

}

void construct_nonlocalKmatrix()
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
