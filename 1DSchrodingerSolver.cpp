///////////////////////////////////////////////////////////////////////////////
//
//    File: 1DSchrodingerSolver.cpp
//    Description: Solves the 1-D eigenvalue problem H|psi> = E|psi>.
//                 This is mostly just to test the hamiltonian class for
//                 a 1D local potential in x-space.
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To Do List:
//      1 -
///////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include <iomanip>

//My classes
#include "Hamiltonian.h"
#include "functions.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
int main(){
  //Set up variables for discretizing the system
  double h = 0.01;
  double Rmax = 10.0;
  double dimension = Rmax / h;

  //These are the parameters to pass to the function.
  double alpha = 1.0;
  void *params = &alpha;

  //Create the Hamiltonian object and then solve for its eigenvalues and
  //eigenvectors.
  Hamiltonian my_hamiltonian(dimension, h, 'x', &harmonic_potential, params);
  my_hamiltonian.solve_eigensystem();

  //Set up an output file
  ofstream my_out("ground_state.dat");
  my_out << "#  x       psi(x)     " << endl;
  my_out << setw(5) << 0.0 << "  " << 0.0 << endl;
  for(int j = 1; j <= dimension; j++){
    my_out << setw(5) << double(j)*h << "  " << setprecision(16)
           << my_hamiltonian.get_eigenvector(1,j) << endl;
  }
  my_out.close();
  return(0);
}
