/**************************************************************************
 *
 * File name: myblas.h
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Version: 1.1
 * Last updated: 17.8.2009
 *
 * A collection of basic linear algebra functions, in the spirit of the
 * BLAS/LAPACK libraries.
 *
 *************************************************************************/



#ifndef __MY_BLAS_H__
#define __MY_BLAS_H__


#include "mex.h"
#include <math.h>



/**************************************************************************
 * Squared value.
 **************************************************************************/
#define SQR(X) ((X)*(X))



/**************************************************************************
 * Matrix-vector multiplication. 
 *
 * Computes an operation of the form:
 *
 *   y := alpha*A*x
 *
 * Parameters:
 *   A - matrix of size n X m
 *   x - vector of length m
 *   y - output vector of length n
 *   alpha - real constant
 *   n, m - dimensions of A
 *
 * Note: This function re-writes the contents of y.
 *
 **************************************************************************/
void mat_vec(double alpha, double A[], double x[], double y[], mwSize n, mwSize m);



/**************************************************************************
 * Matrix-transpose-vector multiplication. 
 *
 * Computes an operation of the form:
 *
 *   y := alpha*A'*x
 *
 * Parameters:
 *   A - matrix of size n X m
 *   x - vector of length n
 *   y - output vector of length m
 *   alpha - real constant
 *   n, m - dimensions of A
 *
 * Note: This function re-writes the contents of y.
 *
 **************************************************************************/
void matT_vec(double alpha, double A[], double x[], double y[], mwSize n, mwSize m);



/**************************************************************************
 * Sparse-matrix-vector multiplication. 
 *
 * Computes an operation of the form:
 *
 *   y := alpha*A*x
 *
 * where A is a sparse matrix.
 *
 * Parameters:
 *   pr,ir,jc - sparse representation of the matrix A, of size n x m
 *   x - vector of length m
 *   y - output vector of length n
 *   alpha - real constant
 *   n, m - dimensions of A
 *
 * Note: This function re-writes the contents of y.
 *
 **************************************************************************/
void mat_sp_vec(double alpha, double pr[], mwIndex ir[], mwIndex jc[], double x[], double y[], mwSize n, mwSize m);



/**************************************************************************
 * Sparse-matrix-transpose-vector multiplication. 
 *
 * Computes an operation of the form:
 *
 *   y := alpha*A'*x
 *
 * where A is a sparse matrix.
 *
 * Parameters:
 *   pr,ir,jc - sparse representation of the matrix A, of size n x m
 *   x - vector of length m
 *   y - output vector of length n
 *   alpha - real constant
 *   n, m - dimensions of A
 *
 * Note: This function re-writes the contents of y.
 *
 **************************************************************************/
void matT_sp_vec(double alpha, double pr[], mwIndex ir[], mwIndex jc[], double x[], double y[], mwSize n, mwSize m);



/**************************************************************************
 * Matrix-sparse-vector multiplication. 
 *
 * Computes an operation of the form:
 *
 *   y := alpha*A*x
 *
 * where A is a matrix and x is a sparse vector.
 *
 * Parameters:
 *   A - matrix of size n X m
 *   pr,ir,jc - sparse representation of the vector x, of length m
 *   y - output vector of length n
 *   alpha - real constant
 *   n, m - dimensions of A
 *
 * Note: This function re-writes the contents of y.
 *
 **************************************************************************/
void mat_vec_sp(double alpha, double A[], double pr[], mwIndex ir[], mwIndex jc[], double y[], mwSize n, mwSize m);



/**************************************************************************
 * Matrix-transpose-sparse-vector multiplication. 
 *
 * Computes an operation of the form:
 *
 *   y := alpha*A'*x
 *
 * where A is a matrix and x is a sparse vector.
 *
 * Parameters:
 *   A - matrix of size n X m
 *   pr,ir,jc - sparse representation of the vector x, of length n
 *   y - output vector of length m
 *   alpha - real constant
 *   n, m - dimensions of A
 *
 * Note: This function re-writes the contents of y.
 *
 **************************************************************************/
void matT_vec_sp(double alpha, double A[], double pr[], mwIndex ir[], mwIndex jc[], double y[], mwSize n, mwSize m);



/**************************************************************************
 * Sparse-matrix-sparse-vector multiplication. 
 *
 * Computes an operation of the form:
 *
 *   y := alpha*A*x
 *
 * where A is a sparse matrix and x is a sparse vector.
 *
 * Parameters:
 *   pr,ir,jc - sparse representation of the matrix A, of size n x m
 *   prx,irx,jcx - sparse representation of the vector x (of length m)
 *   y - output vector of length n
 *   alpha - real constant
 *   n, m - dimensions of A
 *
 * Note: This function re-writes the contents of y.
 *
 **************************************************************************/
void mat_sp_vec_sp(double alpha, double pr[], mwIndex ir[], mwIndex jc[], double prx[], mwIndex irx[], mwIndex jcx[], double y[], mwSize n, mwSize m);



/**************************************************************************
 * Sparse-matrix-transpose-sparse-vector multiplication. 
 *
 * Computes an operation of the form:
 *
 *   y := alpha*A'*x
 *
 * where A is a sparse matrix and x is a sparse vector.
 *
 * Importnant note: this function is provided for completeness, but is NOT efficient.
 * If possible, convert x to non-sparse representation and use matT_vec_sp instead.
 *
 * Parameters:
 *   pr,ir,jc - sparse representation of the matrix A, of size n x m
 *   prx,irx,jcx - sparse representation of the vector x (of length n)
 *   y - output vector of length n
 *   alpha - real constant
 *   n, m - dimensions of A
 *
 * Note: This function re-writes the contents of y.
 *
 **************************************************************************/
void matT_sp_vec_sp(double alpha, double pr[], mwIndex ir[], mwIndex jc[], double prx[], mwIndex irx[], mwIndex jcx[], double y[], mwSize n, mwSize m);



/**************************************************************************
 * Matrix-matrix multiplication. 
 *
 * Computes an operation of the form:
 *
 *   X := alpha*A*B
 *
 * Parameters:
 *   A - matrix of size n X m
 *   B - matrix of size m X k
 *   X - output matrix of size n X k
 *   alpha - real constant
 *   n, m, k - dimensions of A, B
 *
 * Note: This function re-writes the contents of X.
 *
 **************************************************************************/
void mat_mat(double alpha, double A[], double B[], double X[], mwSize n, mwSize m, mwSize k);



/**************************************************************************
 * Matrix-transpose-matrix multiplication. 
 *
 * Computes an operation of the form:
 *
 *   X := alpha*A*B
 *
 * Parameters:
 *   A - matrix of size n X m
 *   B - matrix of size m X k
 *   X - output matrix of size n X k
 *   alpha - real constant
 *   n, m, k - dimensions of A, B
 *
 * Note: This function re-writes the contents of X.
 *
 **************************************************************************/
void matT_mat(double alpha, double A[], double B[], double X[], mwSize n, mwSize m, mwSize k);



/**************************************************************************
 * Tensor-matrix multiplication. 
 *
 * This function accepts a 3-D tensor A of size n X m X k
 * and a 2-D matrix B of size l X k.
 * The function computes the 3-D tensor X of size n X m X l, where
 *
 *   X(i,j,:) = B*A(i,j,:)
 *
 * for all i,j.
 *
 * Parameters:
 *   A - tensor of size n X m X k
 *   B - matrix of size l X k
 *   X - output tensor of size n X m X l
 *   alpha - real constant
 *   n, m, k, l - dimensions of A, B
 *
 * Note: This function re-writes the contents of X.
 *
 **************************************************************************/
void tens_mat(double alpha, double A[], double B[], double X[], mwSize n, mwSize m, mwSize k, mwSize l);



/**************************************************************************
 * Tensor-matrix-transpose multiplication. 
 *
 * This function accepts a 3-D tensor A of size n X m X k
 * and a 2-D matrix B of size k X l.
 * The function computes the 3-D tensor X of size n X m X l, where
 *
 *   X(i,j,:) = B'*A(i,j,:)
 *
 * for all i,j.
 *
 * Parameters:
 *   A - tensor of size n X m X k
 *   B - matrix of size k X l
 *   X - output tensor of size n X m X l
 *   alpha - real constant
 *   n, m, k, l - dimensions of A, B
 *
 * Note: This function re-writes the contents of X.
 *
 **************************************************************************/
void tens_matT(double alpha, double A[], double B[], double X[], mwSize n, mwSize m, mwSize k, mwSize l);



/**************************************************************************
 * Vector-vector sum.
 *
 * Computes an operation of the form:
 *
 *   y := alpha*x + y
 *
 * Parameters:
 *   x - vector of length n
 *   y - output vector of length n
 *   alpha - real constant
 *   n - length of x,y
 *
 * Note: This function re-writes the contents of y.
 *
 **************************************************************************/
void vec_sum(double alpha, double x[], double y[], mwSize n);



/**************************************************************************
 * Triangular back substitution.
 *
 * Solve the set of linear equations
 *
 *   T*x = b
 *
 * where T is lower or upper triangular.
 *
 * Parameters:
 *   ul - 'U' for upper triangular, 'L' for lower triangular
 *   A  - matrix of size n x m containing T
 *   b  - vector of length k
 *   x  - output vector of length k
 *   n  - size of first dimension of A
 *   k  - the size of the equation set, k<=n,m
 *
 * Note:
 *   The matrix A can be of any size n X m, as long as n,m >= k. 
 *   Only the lower/upper triangle of the submatrix A(1:k,1:k) defines the
 *   matrix T (depending on the parameter ul).
 *
 **************************************************************************/
void backsubst(char ul, double A[], double b[], double x[], mwSize n, mwSize k);



/**************************************************************************
 * Solve a set of equations using a Cholesky decomposition.
 *
 * Solve the set of linear equations
 *
 *   M*x = b
 *
 * where M is positive definite with a known Cholesky decomposition:
 * either M=L*L' (L lower triangular) or M=U'*U (U upper triangular).
 *
 * Parameters:
 *   ul - 'U' for upper triangular, 'L' for lower triangular decomposition
 *   A  - matrix of size n x m with the Cholesky decomposition of M
 *   b  - vector of length k
 *   x  - output vector of length k
 *   n  - size of first dimension of A
 *   k  - the size of the equation set, k<=n,m
 *
 * Note:
 *   The matrix A can be of any size n X m, as long as n,m >= k. 
 *   Only the lower/upper triangle of the submatrix A(1:k,1:k) is used as
 *   the Cholesky decomposition of M (depending on the parameter ul).
 *
 **************************************************************************/
void cholsolve(char ul, double A[], double b[], double x[], mwSize n, mwSize k);



/**************************************************************************
 * Maximum absolute value.
 *
 * Returns the index of the coefficient with maximal absolute value in a vector.
 *
 * Parameters:
 *   x - vector of length n
 *   n - length of x
 *
 **************************************************************************/
mwIndex maxabs(double x[], mwSize n);



/**************************************************************************
 * Maximum vector element.
 *
 * Returns the index of the maximal coefficient in a vector.
 *
 * Parameters:
 *   x - vector of length n
 *   n - length of x
 *
 **************************************************************************/
mwIndex maxpos(double x[], mwSize n);



/**************************************************************************
 * Vector-vector dot product.
 *
 * Computes an operation of the form:
 *
 *   c = a'*b
 *
 * Parameters:
 *   a, b - vectors of length n
 *   n - length of a,b
 *
 * Returns: The dot product c.
 *
 **************************************************************************/
double dotprod(double a[], double b[], mwSize n);



/**************************************************************************
 * Indexed vector assignment.
 *
 * Perform a permutation assignment of the form
 *
 *   y = x(ind)
 *
 * where ind is an array of indices to x.
 *
 * Parameters:
 *   y - output vector of length k
 *   x - input vector of arbitrary length
 *   ind - array of indices into x (indices begin at 0)
 *   k - length of the array ind
 *
 **************************************************************************/
void vec_assign(double y[], double x[], mwIndex ind[], mwSize k);



/**************************************************************************
 * Matrix transpose.
 *
 * Computes Y := X'
 *
 * Parameters:
 *   X - input matrix of size n X m
 *   Y - output matrix of size m X n
 *   n, m - dimensions of X
 *
 **************************************************************************/
void transpose(double X[], double Y[], mwSize n, mwSize m);



/**************************************************************************
 * Print a matrix.
 *
 * Parameters:
 *   A - matrix of size n X m
 *   n, m - dimensions of A
 *   matname - name of matrix to display
 *
 **************************************************************************/
void printmat(double A[], int n, int m, char* matname);



/**************************************************************************
 * Print a sparse matrix.
 *
 * Parameters:
 *   A - sparse matrix of type double
 *   matname - name of matrix to display
 *
 **************************************************************************/
void printspmat(mxArray *A, char* matname);


#endif

