///////////////////////////////////////////////////////////////////////////////
//
//    File: Hamiltonian.cpp
//    Description:
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To Do List:
//      1 - Implement the other construct matrix functions
///////////////////////////////////////////////////////////////////////////////

//Includes
#include <Hamiltonian.h>
using namespace arma;

///////////////////////////////////////////////////////////////////////////////
Hamiltonian::Hamiltonian (const int dim, double h)
{
  dimension = dim;
  step_size = h;
  Hmatrix = mat(dimension,dimension);
  Hmatrix.zeroes();
}

Hamiltonian::~Hamiltonian () // Destructor for Hamiltonian
{
  //Doesn't do anything
}

//Constructor for an x-space local potential
Hamiltonian::Hamiltonian (const int dim, double h, const char potential_type,
  double(*potential)(double x, void *params), void *params)
{
  dimension = dim;
  step_size = h;
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
  //solve_eigensystem();
}

//Constructor for an x-space nonlocal potential
Hamiltonian::Hamiltonian (const int dim, double h, const char potential_type,
  double(*potential)(double x1, double x2, void *params), void *params)
{
  dimension = dim;
  step_size = h;
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
  //solve_eigensystem();
}

void construct_localXmatrix()
{
  Hmatrix = mat(dimension,dimension);
  for(int i = 1; i <= dimension; i++){
    for(int j = 1; j <= dimension; j++){
      if(i==j){
        set_element(i, j, 2.0/step_size + xpotential(double(i)*step_size));
      }
      else if(i==j-1){
        set_element(i, j, -1.0/step_size);
      }
      else if(i==j+1){
        set_element(i, j, -1.0/step_size);
      }
      else{
        set_element(i, j, 0);
      }
    }
  }
}

void construct_localKmatrix()
{
  Hmatrix = mat(dimension,dimension);
  Hmatrix.zeroes();
}

void construct_nonlocalXmatrix()
{
  Hmatrix = mat(dimension,dimension);
  Hmatrix.zeroes();
}

void construct_nonlocalKmatrix()
{
  Hmatrix = mat(dimension,dimension);
  Hmatrix.zeroes();
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
