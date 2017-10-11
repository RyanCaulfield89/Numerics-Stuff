///////////////////////////////////////////////////////////////////////////////
//
//    File: Optical_potential.cpp
//
//    Description: See Optical_potential.h
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To do list:
//      1 - Implement everything.
//
///////////////////////////////////////////////////////////////////////////////

//Includes
#include "Optical_potential.h"
using namespace arma;
using namespace std;

///////////////////////////////////////////////////////////////////////////////

Optical_potential::Optical_potential(int A, int Z){

}

double Optical_potential::real_term(double r){
  double return_value = -V * woods_saxon_potential(r, R_real, a_real);
  return return_value;
}

double Optical_potential::imaginary_term(double r){
  double return_value = -W * woods_saxon_potential(r, R_imag, a_imag);
  return return_value;
}

double Optical_potential::boundary_layer_term(double r){
  double return_value = 4.0*a_imag*W_D*diff_woods_saxon_potential(r,R_imag,a_imag);
  return return_value;
}

double Optical_potential::spin_orbit_term(double r){

}

double Optical_potential::coulomb_term(double r){

}

double Optical_potential::woods_saxon_potential(double r, double R, double a){
  double return_value;
  return_value = 1.0 / (1.0 + exp((r - R) / a));
  return return_value;
}

double Optical_potential::diff_woods_saxon_potential(double r, double R, double a){
  double return_value;
  return_value = - 1.0 / (a * pow(1.0 + exp((x - R) / a),2));
  return return_value;
}

std::complex<double> Optical_potential::calculate_optical_potential(double r){

}
