///////////////////////////////////////////////////////////////////////////////
//
//    File: main.cpp
//    Description:  This is the main executable for an SRG and Optical potential
//                  project. This takes potentials from some toy models, and
//                  numerically determines optical potentials for those models.
//                  It also can perform SRG transformations on potentials. It
//                  compares the results from SRG transforming potentials and
//                  then finding optical potentials to the results from SRG
//                  transforming the optical potentials directly. In principle,
//                  the physical observables, phase shifts, should be insensitive
//                  to this.
//
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To Do List:
//      1 - Implement classes
//      2 - Write the body
//
///////////////////////////////////////////////////////////////////////////////

//Classes that are included
#include "Functions.h"
#include "Hamiltonian.h"
#include "OpticalPotential.h"
#include "Scattering.h"

using namespace std;

int main(){

  //Decide which potential to use. To do this, we need to import a class which
  //has defined the potentials that we want to use. This is the function class,
  //which is included under Function.h.

  //Include a message here which tells the user which potentials we will use.

  //Now, we need to construct a Hamiltonian from the potential. This is done by
  //the Hamiltonian class included under Hamiltonian.h.

  //Before we SRG evolve, we should find the phase shifts. Later we will compare
  //the results to those obtained using the optical potential method. This can
  //be done using the scattering class.

  //Here, we need to SRG evolve the Hamiltonian. This is a unitary transformation
  //that is determined by solving a flow equation. This is also done by the
  //Hamiltonian class.

  //Inlude a message here that the Hamiltonian have been SRG evolved.

  //Now that we have our Hamiltonians, both evolved and unevolved, we need to
  //find the optical potentials for them. This involves finding the matrix
  //elements in the eigenbasis that you want to remove. In our toy models, this
  //will be the nuclear eigenbasis. We also need to find an inverse differential
  //operator, which means we need to find a Green's function kernal. All of this
  //is done by the OpticalPotential class, which is included under
  //OpticalPotential.h.

  //Include a message here that the optical potentials have been calculated.

  //Once we have the optical potentials, we can use them to determine phase
  //shifts. These are the physical observables that we care about. Since our
  //transformations are unitary, we should always get the same phase shifts.
  //This is done by the scattering class, which is included under Scattering.h.

  //Include a message here that the phase shifts have been calculated.

  //When all is said and done, we need a way to visualize the results. All of
  //potentials, evolved, unevolved and optical are saved in .dat files.
  //From there, plotting software can display the potentials in a way that
  //shows the changes in locality. The energy dependent phase shifts are also
  //save in .dat files, but these are much easier to plot. A simple gnuplot
  //script will plot the various phase shifts as a function of energy and the
  //relative error between the supposedly equvialent phase shifts.

  //Include instructions on how to plot the important data.
}
