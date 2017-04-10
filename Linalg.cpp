///////////////////////////////////////////////////////////////////////////////
//
//    File: Linalg.cpp
//    Description: This is a collection of linear algebra methods
//
//    Programmer: Ryan Caulfield Caulfield.16@osu.edu
//
//    To Do List:
//      1 - This needs to be tested.
//      2 - Find other derivative methods that will be useful and add them.
//      3 - Make a seperate header file.
///////////////////////////////////////////////////////////////////////////////

//Includes
#include <gsl/gsl_eigen.h>
#include <cmath>
//using namespace std;

//Class Header
class Linalg
{
  public:
    //Constructors and Destructors
    Linalg(int dim);
    ~Linalg();

    //Getters and Setters
    void set_element(int i, int j, double value);
    double get_element(int i, int j);
    double get_eigenvalue(int n);
    double get_eigenvector(int n, int vector_entry);
    void eigen_system();

  private:
    int dimension; //dimensionality of the matrix

    //This stuff is for the gsl routines
    gsl_matrix *Amat_ptr;   //This is the matrix that will be solved
    gsl_vector *Eigval_ptr;	//This will hold the eigenvalues
    gsl_matrix *Eigvec_ptr;	//Each column will be an eigenvector
    gsl_eigen_symmv_workspace *worksp; //workspace used by gsl
    gsl_vector *eigenvector_ptr;  // This is a single vector
}
///////////////////////////////////////////////////////////////////////////////
Linalg::Linalg(int dim)
{
  //Initialize everything and allocate space to stuff.
  dimension = dim;
  Amat_ptr = gsl_matrix_alloc (dimension, dimension);
  Eigval_ptr = gsl_vector_alloc (dimension);
  Eigvec_ptr = gsl_matrix_alloc (dimension, dimension);
  worksp = gsl_eigen_symmv_alloc (dimension);
  eigenvector_ptr = gsl_vector_alloc (dimension);
}

Linalg::~Linalg()
{
  //Free up space when you're done.
  gsl_matrix_free (Eigvec_ptr);
  gsl_vector_free (Eigval_ptr);
  gsl_matrix_free (Amat_ptr);
  gsl_vector_free (eigenvector_ptr);
  gsl_eigen_symmv_free (worksp);
}

void Linalg::set_element(int i, int j, double value)
{
    gsl_matrix_set (Amat_ptr, i-1, j-1, value);
}

double Linalg::get_element(int i, int j)
{
  return gsl_matrix_get(Amat_ptr, i - 1, j - 1);
}

double Linalg::get_eigenvalue(int n)
{
  return gsl_vector_get (Eigval_ptr, n - 1);
}

double Linalg::get_eigenvector(int n, int vector_entry)
{
  gsl_matrix_get_col (eigenvector_ptr, Eigvec_ptr, i-1);
  return gsl_vector_get (eigenvector_ptr, j-1);
}

void Linalg::eigen_system()
{
  //This is the gsl method.
  gsl_eigen_symmv (Amat_ptr, Eigval_ptr, Eigvec_ptr, worksp);
  //Sort the eigenvalues and eigenvectors.
  gsl_eigen_symmv_sort (Eigval_ptr, Eigvec_ptr, GSL_EIGEN_SORT_VAL_ASC);
}
