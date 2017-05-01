///////////////////////////////////////////////////////////////////////////////
//
//    File: Hamiltonian.cpp
//    Description: See Hamiltonian.h
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To Do List:
//      1 - Implement the other construct matrix functions
///////////////////////////////////////////////////////////////////////////////

//Includes
#include "Hamiltonian.h"
using namespace arma;

///////////////////////////////////////////////////////////////////////////////
Hamiltonian::Hamiltonian (const int dim)
{
  dimension = dim;
  Hmatrix = zeros<cx_mat>(dimension,dimension);
}

Hamiltonian::~Hamiltonian () // Destructor for Hamiltonian
{
  //Doesn't do anything
}

//Constructor for an x-space local potential
Hamiltonian::Hamiltonian (double Rmin, double Rmax, double h, const char potential_type,
  std::complex<double>(*potential)(double x, void *params), void *params)
{
  left_boundary = Rmin;
  right_boundary = Rmax;
  step_size = h;
  dimension = (Rmax - Rmin) / h;
  parameters = params;
  if(potential_type == 'x')
  {
    xpotential = potential;
    construct_localXmatrix();
  }
  if(potential_type == 'k')
  {
    kpotential = potential;
    construct_localKmatrix();
  }
  //solve_eigensystem();
}

//Constructor for an x-space nonlocal potential
Hamiltonian::Hamiltonian (double Rmin, double Rmax, double h, const char potential_type,
  std::complex<double>(*potential)(double x1, double x2, void *params), void *params)
{
  left_boundary = Rmin;
  right_boundary = Rmax;
  step_size = h;
  dimension = (Rmax - Rmin) / h;
  parameters = params;
  if(potential_type == 'x')
  {
    xnonLocalPotential = potential;
    construct_nonlocalXmatrix();
  }
  if(potential_type == 'k')
  {
    knonLocalPotential = potential;
    construct_nonlocalKmatrix();
  }
  //solve_eigensystem();
}

void Hamiltonian::construct_localXmatrix()
{
  Hmatrix = cx_mat(dimension,dimension);
  for(int i = 1; i <= dimension; i++){
    for(int j = 1; j <= dimension; j++){
      double x = left_boundary + double(i) * step_size;
      if(i==j){
        set_element(i, j, 2.0 / step_size + xpotential(x, parameters));
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

void Hamiltonian::construct_localKmatrix()
{
  Hmatrix = cx_mat(dimension,dimension);
  for(int i = 1; i <= dimension; i++){
    for(int j = 1; j <= dimension; j++){
      double k = left_boundary + double(i) * step_size;
      if(i==j){
        set_element(i, j, k*k + kpotential(k, parameters));
      }
      else{
        set_element(i, j, 0);
      }
    }
  }
}

void Hamiltonian::construct_nonlocalXmatrix()
{
  Hmatrix = cx_mat(dimension,dimension);
  for(int i = 1; i <= dimension; i++){
    for(int j = 1; j <= dimension; j++){
      double x = left_boundary + double(i) * step_size;
      double y = right_boundary + double(j) * step_size;
      if(i==j){
        set_element(i, j, 2.0 / step_size + xnonLocalPotential(x, y, parameters));
      }
      else if(i==j-1){
        set_element(i, j, -1.0 / step_size + xnonLocalPotential(x, y, parameters));
      }
      else if(i==j+1){
        set_element(i, j, -1.0 / step_size + xnonLocalPotential(x, y, parameters));
      }
      else{
        set_element(i, j, xnonLocalPotential(x, y, parameters));
      }
    }
  }
}

void Hamiltonian::construct_nonlocalKmatrix()
{
  Hmatrix = cx_mat(dimension,dimension);
  for(int i = 1; i <= dimension; i++){
    for(int j = 1; j <= dimension; j++){
      double k = left_boundary + double(i) * step_size;
      double p = right_boundary + double(j) * step_size;
      if(i==j){
        set_element(i, j, k*k + knonLocalPotential(k, p, parameters));
      }
      else{
        set_element(i, j, knonLocalPotential(k, p, parameters));
      }
    }
  }
}

void Hamiltonian::solve_eigensystem()
{
  eig_sym(eigenvalues, eigenvectors, Hmatrix);
}

void Hamiltonian::set_element(const int i, const int j, const std::complex<double> value)
{
  Hmatrix(i-1,j-1) = value;
}

std::complex<double> Hamiltonian::get_element(const int i, const int j)
{
  return Hmatrix(i-1,j-1);
}

std::complex<double> Hamiltonian::get_eigenvalue(const int i)
{
  return eigenvalues (i-1);
}

std::complex<double> Hamiltonian::get_eigenvector(const int i, const int j)
{
  return eigenvectors (i-1, j-1);
}
