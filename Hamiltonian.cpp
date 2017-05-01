///////////////////////////////////////////////////////////////////////////////
//
//    File: Hamiltonian.cpp
//    Description: See Hamiltonian.h
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To Do List:
//      1 - make set_element more elegant so it can be used in
//          construct_scattering_matrix
///////////////////////////////////////////////////////////////////////////////

//Includes
#include "Hamiltonian.h"
using namespace arma;
using namespace std;

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
}

//Constructor for 1-D scattering
Hamiltonian::Hamiltonian (double energy, double Rmin, double Rmax, double h,
  const char potential_type,
  std::complex<double>(*potential)(double x, void *params), void *params)
{
  //First set all the important attributes
  left_boundary = Rmin;
  right_boundary = Rmax;
  step_size = h;
  dimension = (Rmax - Rmin) / h;
  parameters = params;
  initial_state_energy = energy;
  if(potential_type == 'x')
  {
    xpotential = potential;
  }
  if(potential_type == 'k')
  {
    kpotential = potential;
  }
  //Now construct (E - H_0 + iepsilon)^(-1)V
  construct_scattering_matrix();
  construct_initial_wf();
  scattered_wf = scattering_matrix * initial_wf;
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

void Hamiltonian::construct_scattering_matrix()
{
  scattering_matrix = cx_mat(dimension, dimension);
  //First construct (E - H_0 + iepsilon)
  cx_mat temp_matrix1;
  temp_matrix1 = cx_mat(dimension, dimension);
  for(int i = 1; i <= dimension; i++){
    for(int j = 1; j <= dimension; j++){
      if(i==j){
        temp_matrix1(i-1,j-1) = -2.0/step_size + initial_state_energy + 1i*10e-8;
      }
      else if(i==j-1){
        temp_matrix1(i-1,j-1) =  1.0 / step_size;
      }
      else if(i==j+1){
        temp_matrix1(i-1,j-1) =  1.0 / step_size;
      }
      else{
        temp_matrix1(i-1,j-1) =  0;
      }
    }
  }
  //Now construct V as a matrix.
  cx_mat temp_matrix2;
  temp_matrix2 = cx_mat(dimension, dimension);
  for(int i = 1; i <= dimension; i++){
    for(int j = 1; j <= dimension; j++){
      double x = left_boundary + double(i) * step_size;
      if(i==j){
        temp_matrix2(i-1,j-1) = xpotential(x, parameters);
      }
      else{
        temp_matrix2(i-1,j-1) =  0;
      }
    }
  }
  scattering_matrix = inv(temp_matrix1) * temp_matrix2;//(E-H_0+iepsilon)^(-1)*V
}

void Hamiltonian::construct_initial_wf()
{
  initial_wf = cx_vec(dimension);
  double k = sqrt(initial_state_energy);
  for(int i = 1; i < dimension; i++){
    double x = left_boundary + double(i) * step_size;
    complex<double> ikx = 1i * k * x;
    initial_wf(i) = exp(ikx);
  }
}

void Hamiltonian::solve_eigensystem()
{
  eig_sym(eigenvalues, eigenvectors, Hmatrix);
}

void Hamiltonian::save_eigenvector(const int i, const char* filename)
{
  //Set up an output file
  ofstream my_out(filename);
  my_out << "#  x       psi(x)     " << endl;
  my_out << setw(5) << right_boundary << "  " << 0.0 << endl;
  for(int j = 1; j <= dimension; j++){
    my_out << setw(5) << left_boundary + double(j) * step_size << "  "
                      << setprecision(16) << norm(get_eigenvector(i,j)) << endl;
  }
  my_out << setw(5) << left_boundary << "  " << 0.0 << endl;
  my_out.close();
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

double Hamiltonian::calculate_T()
{
  complex<double> t = 0;
  double k = sqrt(initial_state_energy);

  for(int i = 1; i < dimension; i++)
  {
    double x = left_boundary + double(i) * step_size;
    complex<double> ikx = 1i * k * x;
    t += exp(ikx)*scattered_wf(i);
  }
  t /= (right_boundary - left_boundary);
  return(norm(t));
}

double Hamiltonian::calculate_R()
{
  complex<double> r = 0;
  double k = sqrt(initial_state_energy);

  for(int i = 1; i < dimension; i++)
  {
    double x = left_boundary + double(i) * step_size;
    complex<double> ikx = 1i * k * x;
    r += exp(ikx)*scattered_wf(i);
  }
  r /= (right_boundary - left_boundary);
  return(norm(r));
}
