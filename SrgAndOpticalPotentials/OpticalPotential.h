///////////////////////////////////////////////////////////////////////////////
//
//    File: OpticalPotential.h
//    Description:
//
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To Do List:
//      1 -
//
///////////////////////////////////////////////////////////////////////////////

#ifndef OPTICALPOT_H
#define OPTICALPOT_H //To avoid including this file twice

//Includes
#include <armadillo>  //This package does all the linear algebra stuff.
#include <complex>    //Complex number stuff.
#include <cmath>      //Math stuff
#include "Hamiltonian.h"

///////////////////////////////////////////////////////////////////////////////

Class OpticalPotential{
  Public:
    //Constructor that takes a Hamiltonian.
    OpticalPotential(Hamiltonian::Hamiltonian initial_Hamiltonian);

  Private:
    Hamiltonian::Hamiltonian initial_Hamiltonian;
    Hamiltonian::Hamiltonian final_Hamiltonian;
}

#endif
