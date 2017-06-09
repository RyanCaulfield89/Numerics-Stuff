///////////////////////////////////////////////////////////////////////////////
//
//    File: Optical_potential.h
//
//    Description: Header file for the class Optical_potential. This should be
//    able to calculate optical potentials given a set of parameters, and solve
//    the relevent scattering problem.
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To do list:
//      1 - Nothing yet.
//
///////////////////////////////////////////////////////////////////////////////
#ifndef OPTICAL_POTENTIAL_H
#define OPTICAL_POTENTIAL_H //To avoid including this file twice

//Include files
#include <armadillo>  //This package does all the linear algebra stuff.
#include <complex>    //Complex number stuff.
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
class Optical_potential
{
  public:
    Optical_potential(int A, int Z);

    //The potential has 5 terms that need to be included. The forms of these
    //terms are all real valued functions, but some will have imaginary
    //coefficients when the potential is calculated. The following calculate
    //these terms at a distance r.
    double real_term(double r);
    double imaginary_term(double r);
    double boundary_layer_term(double r);
    double spin_orbit_term(double r);
    double coulomb_term(double r);

    //Woods-Saxon potential. This is used to calculate the real and imaginary
    //terms, and its derivative is used in the boundary term and the SO term.
    double woods_saxon_potential(double r, double R, double a);
    double diff_woods_saxon_potential(double r, double R, double a);

    //This calculates the value of the optical potential at a distance r.
    std::complex<double> calculate_optical_potential(double r);

  private:
    //Nucleus properties
    int A; //Total number of nucleons
    int Z; //Total number of protons
    double R_nucleus; //Radius of the nucleus

    //Projectile properties
    double E_projectile; //Energy of the projectile
    int Z_projectile; //Total charge of the projectile in fundamental charge units.
    int T_Z; //z-component of the isospin of projectle
    double S; //Spin angular momentum quantum number of the projectile
    double L; //Orbital angular momentum quantum number of the projectile
    double J; //Total angular momentum

    //Parameters for the potential
    //Real part of the potential
    double V; //Strength of the well
    double R_real; //Width of the well
    double a_real; //Width of the boundary layer

    //Imaginary/Absorptive part of the potential
    double W; //Strength of the well
    double R_imag; //Width of the well
    double a_imag; //Width of the boundary layer

    //Imaginary/Absorptive boundary term
    double W_D; //Strength of the boundary layer term

    //Spin Orbit coupling term
    double V_S; //Strength of the SO coupling

};

#endif
