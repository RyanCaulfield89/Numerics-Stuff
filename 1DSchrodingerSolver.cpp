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

//My classes
#include "Hamiltonian.h"
#include "functions.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
int main(){
  //Set up variables for discretizing the system
  double h = 0.01;
  double Rmin = -2.0;
  double Rmax = 2.0;

  //These are the parameters to pass to the function.
  std::complex<double> alpha = 1.0;
  void *params = &alpha;

  //Create the Hamiltonian object and then solve for its eigenvalues and
  //eigenvectors. Save the ground state for plotting.
  Hamiltonian my_hamiltonian(Rmin, Rmax, h, 'x', &harmonic_potential, params);
  my_hamiltonian.solve_eigensystem();
  my_hamiltonian.save_eigenvector(1, "ground_state.dat");
  return(0);
}
