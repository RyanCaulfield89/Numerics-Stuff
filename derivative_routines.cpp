////////////////////////////////////////////////////////////////////////////
//
//    File: derivative_routines.cpp
//    Description: This is a set of tools for numerically differentiating
//    a function.
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To do list:
//      1 - This needs to be tested.
//      2 - Find other derivative methods that will be useful and add them.
//      3 - Is 1000 the right size for workspace?
//      4 - make a seperate header file if this starts getting long.
//
//////////////////////////////////////////////////////////////////////////////

//includes
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_diff.h>
using namespace std;
