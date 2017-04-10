///////////////////////////////////////////////////////////////////////////////
//
//    File: Linalg.cpp
//    Description: This is a collection of linear algebra methods
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To Do List:
//      1 - This is not fully implemented
///////////////////////////////////////////////////////////////////////////////

//Includes
#include <gsl/gsl_eigen.h>
#include <cmath>
//using namespace std;

//Headers
class Linalg
{
  public:
    Linalg(int dim);
    void set_element(int i, int j, double value);
    double get_element(int i, int j);

  private:
    int dimension; //dimensionality of the matrix

    //This stuff is for the gsl routines
    gsl_matrix *Amat_ptr;   //This is the matrix that will be solved
    gsl_vector *Eigval_ptr;	//This will hold the eigenvalues
    gsl_matrix *Eigvec_ptr;	//Each column will be an eigenvector
    gsl_eigen_symmv_workspace *worksp;

}
///////////////////////////////////////////////////////////////////////////////
void eigen_system(int dimension, double **matrix)
{
  //These are needed for the gsl routine
  gsl_matrix *Amat_ptr;   //This is the matrix that will be solved
  gsl_vector *Eigval_ptr;	//This will hold the eigenvalues
  gsl_matrix *Eigvec_ptr;	//Each column will be an eigenvector
  gsl_eigen_symmv_workspace *worksp;

  //We need to allocate space for these gsl objects
  Amat_ptr = gsl_matrix_alloc (dimension, dimension);
  Eigval_ptr = gsl_vector_alloc (dimension);
  Eigvec_ptr = gsl_matrix_alloc (dimension, dimension);
  worksp = gsl_eigen_symmv_alloc (dimension);

  //This is the gsl method
  gsl_eigen_symmv (Amat_ptr, Eigval_ptr, Eigvec_ptr, worksp);
}
