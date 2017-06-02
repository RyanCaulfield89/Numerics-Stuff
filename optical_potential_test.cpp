///////////////////////////////////////////////////////////////////////////////
//
//    File: optical_potential_test.cpp
//
//    Description:
//      File that plots the real and imaginary parts of the optical potential
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To do list:
//      1 - Nothing yet.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <complex>
#include "functions.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
int main() {
  //Set up parameters
  optical_potential_parameters op_parameters;
  op_parameters.E = 50.0;
  op_parameters.A = 56;
  op_parameters.Z = 26;
  double R = 1.25 * pow(op_parameters.A, 1.0/3.0);
  woods_saxon_parameters real_parameters;
  real_parameters.a = 0.65;
  real_parameters.R = R;
  woods_saxon_parameters imaginary_parameters;
  imaginary_parameters.a = 0.47;
  imaginary_parameters.R = R;
  op_parameters.real_parameters = real_parameters;
  op_parameters.imaginary_parameters = imaginary_parameters;
  void *op_parameters_ptr;
  op_parameters_ptr = &op_parameters;

  //Now, output the real and imaginary parts to a .dat file.
  ofstream my_out("optical_potential_test.dat");
  my_out << "#  x    Real part    Imaginary part" << endl;
  double x = 0.0;
  complex<double> value = 0.0;
  for(int i = 0; i < 2000; i++){
    x = 0.01 * i;
    value = optical_potential(x,op_parameters_ptr);
    my_out << setprecision(5) << x << "  " << setprecision(16)
           << real(value) << "  " << imag(value) << endl;
  }
  my_out.close();
}
