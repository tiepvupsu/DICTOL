/**************************************************************************
 *
 * File name: myblas.c
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Version: 1.1
 * Last updated: 13.8.2009
 *
 *************************************************************************/


#include "myblas.h"
#include <ctype.h>


/* find maximum of absolute values */

mwIndex maxabs(double c[], mwSize m)
{
  mwIndex maxid=0, k;
  double absval, maxval = SQR(*c);   /* use square which is quicker than absolute value */

  for (k=1; k<m; ++k) {
    absval = SQR(c[k]);
    if (absval > maxval) {
      maxval = absval;
      maxid = k;
    }
  }
  return maxid;
}


/* compute y := alpha*x + y */

void vec_sum(double alpha, double x[], double y[], mwSize n)
{
  mwIndex i;

  for (i=0; i<n; ++i) {
    y[i] += alpha*x[i];
  }
}


/* compute y := alpha*A*x */

void mat_vec(double alpha, double A[], double x[], double y[], mwSize n, mwSize m)
{
  mwIndex i, j, i_n;
  double *Ax;

  Ax = mxCalloc(n,sizeof(double));

  for (i=0; i<m; ++i) {
    i_n = i*n;
    for (j=0; j<n; ++j) {
      Ax[j] += A[i_n+j] * x[i];
    }
  }

  for (j=0; j<n; ++j) {
    y[j] = alpha*Ax[j];
  }

  mxFree(Ax);
}


/* compute y := alpha*A'*x */

void matT_vec(double alpha, double A[], double x[], double y[], mwSize n, mwSize m)
{
  mwIndex i, j, n_i;
  double sum0, sum1, sum2, sum3;

  for (j=0; j<m; ++j) {
    y[j] = 0;
  }

  /* use loop unrolling to accelerate computation */

  for (i=0; i<m; ++i) {
    n_i = n*i;
    sum0 = sum1 = sum2 = sum3 = 0;
    for (j=0; j+4<n; j+=4) {
      sum0 += A[n_i+j]*x[j];
      sum1 += A[n_i+j+1]*x[j+1];
      sum2 += A[n_i+j+2]*x[j+2];
      sum3 += A[n_i+j+3]*x[j+3];
    }
    y[i] += alpha * ((sum0 + sum1) + (sum2 + sum3));
    while (j<n) {
      y[i] += alpha*A[n_i+j]*x[j];
      j++;
    }
  }
}


/* compute y := alpha*A*x */

void mat_sp_vec(double alpha, double pr[], mwIndex ir[], mwIndex jc[], double x[], double y[], mwSize n, mwSize m)
{
  
  mwIndex i, j, j1, j2;

  for (i=0; i<n; ++i) {
    y[i] = 0;
  }
  
  j2 = jc[0];
  for (i=0; i<m; ++i) {
    j1 = j2; j2 = jc[i+1];
    for (j=j1; j<j2; ++j) {
      y[ir[j]] += alpha * pr[j] * x[i];
    }
  }
  
}


/* compute y := alpha*A'*x */

void matT_sp_vec(double alpha, double pr[], mwIndex ir[], mwIndex jc[], double x[], double y[], mwSize n, mwSize m)
{
  
  mwIndex i, j, j1, j2;
  
  for (i=0; i<m; ++i) {
    y[i] = 0;
  }
  
  j2 = jc[0];
  for (i=0; i<m; ++i) {
    j1 = j2; j2 = jc[i+1];
    for (j=j1; j<j2; ++j) {
      y[i] += alpha * pr[j] * x[ir[j]];
    }
  }
  
}


/* compute y := alpha*A*x */

void mat_vec_sp(double alpha, double A[], double pr[], mwIndex ir[], mwIndex jc[], double y[], mwSize n, mwSize m)
{
  
  mwIndex i, j, j_n, k, kend;
  
  for (i=0; i<n; ++i) {
    y[i] = 0;
  }
  
  kend = jc[1];
  if (kend==0) {   /* x is empty */
    return;
  }
  
  for (k=0; k<kend; ++k) {
    j = ir[k];
    j_n = j*n;
    for (i=0; i<n; ++i) {
      y[i] += alpha * A[i+j_n] * pr[k];
    }
  }

}


/* compute y := alpha*A'*x */

void matT_vec_sp(double alpha, double A[], double pr[], mwIndex ir[], mwIndex jc[], double y[], mwSize n, mwSize m)
{
  
  mwIndex i, j, j_n, k, kend;
  
  for (i=0; i<m; ++i) {
    y[i] = 0;
  }
  
  kend = jc[1];
  if (kend==0) {   /* x is empty */
    return;
  }
  
  for (j=0; j<m; ++j) {
    j_n = j*n;
    for (k=0; k<kend; ++k) {
      i = ir[k];
      y[j] += alpha * A[i+j_n] * pr[k];
    }
  }
  
}


/* compute y := alpha*A*x */

void mat_sp_vec_sp(double alpha, double pr[], mwIndex ir[], mwIndex jc[], double prx[], mwIndex irx[], mwIndex jcx[], double y[], mwSize n, mwSize m)
{
  
  mwIndex i, j, k, kend, j1, j2;

  for (i=0; i<n; ++i) {
    y[i] = 0;
  }
  
  kend = jcx[1]; 
  if (kend==0) {   /* x is empty */
    return;
  }
  
  for (k=0; k<kend; ++k) {
    i = irx[k];
    j1 = jc[i]; j2 = jc[i+1];
    for (j=j1; j<j2; ++j) {
      y[ir[j]] += alpha * pr[j] * prx[k];
    }
  }
  
}


/* compute y := alpha*A'*x */

void matT_sp_vec_sp(double alpha, double pr[], mwIndex ir[], mwIndex jc[], double prx[], mwIndex irx[], mwIndex jcx[], double y[], mwSize n, mwSize m)
{
  
  mwIndex i, j, k, jend, kend, jadd, kadd, delta;
  
  for (i=0; i<m; ++i) {
    y[i] = 0;
  }
  
  kend = jcx[1];
  if (kend==0) {   /* x is empty */
    return;
  }
  
  for (i=0; i<m; ++i) {
    j = jc[i]; 
    jend = jc[i+1];
    k = 0;
    while (j<jend && k<kend) {
      
      delta = ir[j] - irx[k];
      
      if (delta) { /* if indices differ - increment the smaller one */
        jadd = delta<0;
        kadd = 1-jadd;
        j += jadd;
        k += kadd;
      }
      
      else {    /* indices are equal - add to result and increment both */
        y[i] += alpha * pr[j] * prx[k];
        j++; k++;
      }
    }
  }
  
}


/* matrix-matrix multiplication */

void mat_mat(double alpha, double A[], double B[], double X[], mwSize n, mwSize m, mwSize k)
{
  mwIndex i1, i2, i3, iX, iA, i2_n;
  double b;
  
  for (i1=0; i1<n*k; i1++) {
    X[i1] = 0;
  }

  for (i2=0; i2<m; ++i2) {
    i2_n = i2*n;
    iX = 0;
    for (i3=0; i3<k; ++i3) {
      iA = i2_n;
      b = B[i2+i3*m];
      for (i1=0; i1<n; ++i1) {
        X[iX++] += A[iA++]*b;
      }
    }
  }
  
  for (i1=0; i1<n*k; i1++) {
    X[i1] *= alpha;
  }
}


/* matrix-transpose-matrix multiplication */

void matT_mat(double alpha, double A[], double B[], double X[], mwSize n, mwSize m, mwSize k)
{
  mwIndex i1, i2, i3, iX, iA, i2_n;
  double *x, sum0, sum1, sum2, sum3;

  for (i2=0; i2<m; ++i2) {
    for (i3=0; i3<k; ++i3) {
      sum0 = sum1 = sum2 = sum3 = 0;
      for (i1=0; i1+4<n; i1+=4) {
        sum0 += A[i1+0+i2*n]*B[i1+0+i3*n];
        sum1 += A[i1+1+i2*n]*B[i1+1+i3*n];
        sum2 += A[i1+2+i2*n]*B[i1+2+i3*n];
        sum3 += A[i1+3+i2*n]*B[i1+3+i3*n];
      }
      X[i2+i3*m] = (sum0+sum1) + (sum2+sum3);
      while(i1<n) {
        X[i2+i3*m] += A[i1+i2*n]*B[i1+i3*n];
        i1++;
      }
    }
  }
  
  for (i1=0; i1<m*k; i1++) {
    X[i1] *= alpha;
  }
}


/* tensor-matrix product */

void tens_mat(double alpha, double A[], double B[], double X[], mwSize n, mwSize m, mwSize k, mwSize l)
{
  mwIndex i1, i2, i3, i4, i2_n, nml;
  double b;
  
  nml = n*m*l;
  for (i1=0; i1<nml; ++i1) {
    X[i1] = 0;
  }

  for (i2=0; i2<m; ++i2) {
    i2_n = i2*n;
    for (i3=0; i3<k; ++i3) {
      for (i4=0; i4<l; ++i4) {
        b = B[i4+i3*l];
        for (i1=0; i1<n; ++i1) {
          X[i1 + i2_n + i4*n*m] += A[i1 + i2_n + i3*n*m] * b;
        }
      }
    }
  }
  
  for (i1=0; i1<nml; ++i1) {
    X[i1] *= alpha;
  }
}


/* tensor-matrix-transpose product */

void tens_matT(double alpha, double A[], double B[], double X[], mwSize n, mwSize m, mwSize k, mwSize l)
{
  mwIndex i1, i2, i3, i4, i2_n, nml;
  double b;
  
  nml = n*m*l;
  for (i1=0; i1<nml; ++i1) {
    X[i1] = 0;
  }

  for (i2=0; i2<m; ++i2) {
    i2_n = i2*n;
    for (i4=0; i4<l; ++i4) {
      for (i3=0; i3<k; ++i3) {
        b = B[i3+i4*k];
        for (i1=0; i1<n; ++i1) {
          X[i1 + i2_n + i4*n*m] += A[i1 + i2_n + i3*n*m] * b;
        }
      }
    }
  }
  
  for (i1=0; i1<nml; ++i1) {
    X[i1] *= alpha;
  }
}


/* dot product */

double dotprod(double a[], double b[], mwSize n)
{
  double sum = 0;
  mwIndex i;
  for (i=0; i<n; ++i)
    sum += a[i]*b[i];
  return sum;
}


/* find maximum of vector */

mwIndex maxpos(double c[], mwSize m)
{
  mwIndex maxid=0, k;
  double val, maxval = *c;

  for (k=1; k<m; ++k) {
    val = c[k];
    if (val > maxval) {
      maxval = val;
      maxid = k;
    }
  }
  return maxid;
}


/* solve L*x = b */

void backsubst_L(double L[], double b[], double x[], mwSize n, mwSize k)
{
  mwIndex i, j;
  double rhs;

  for (i=0; i<k; ++i) {
    rhs = b[i];
    for (j=0; j<i; ++j) {
      rhs -= L[j*n+i]*x[j];
    }
    x[i] = rhs/L[i*n+i];
  }
}


/* solve L'*x = b */

void backsubst_Lt(double L[], double b[], double x[], mwSize n, mwSize k)
{
  mwIndex i, j;
  double rhs;

  for (i=k; i>=1; --i) {
    rhs = b[i-1];
    for (j=i; j<k; ++j) {
      rhs -= L[(i-1)*n+j]*x[j];
    }
    x[i-1] = rhs/L[(i-1)*n+i-1];
  }
}


/* solve U*x = b */

void backsubst_U(double U[], double b[], double x[], mwSize n, mwSize k)
{
  mwIndex i, j;
  double rhs;

  for (i=k; i>=1; --i) {
    rhs = b[i-1];
    for (j=i; j<k; ++j) {
      rhs -= U[j*n+i-1]*x[j];
    }
    x[i-1] = rhs/U[(i-1)*n+i-1];
  }
}


/* solve U'*x = b */

void backsubst_Ut(double U[], double b[], double x[], mwSize n, mwSize k)
{
  mwIndex i, j;
  double rhs;

  for (i=0; i<k; ++i) {
    rhs = b[i];
    for (j=0; j<i; ++j) {
      rhs -= U[i*n+j]*x[j];
    }
    x[i] = rhs/U[i*n+i];
  }
}


/* back substitution solver */

void backsubst(char ul, double A[], double b[], double x[], mwSize n, mwSize k)
{
  if (tolower(ul) == 'u') {
    backsubst_U(A, b, x, n, k);
  }
  else if (tolower(ul) == 'l') {
    backsubst_L(A, b, x, n, k);
  }
  else {
    mexErrMsgTxt("Invalid triangular matrix type: must be ''U'' or ''L''");
  }
}


/* solve equation set using cholesky decomposition */

void cholsolve(char ul, double A[], double b[], double x[], mwSize n, mwSize k)
{
  double *tmp;

  tmp = mxMalloc(k*sizeof(double));

  if (tolower(ul) == 'l') {
    backsubst_L(A, b, tmp, n, k);
    backsubst_Lt(A, tmp, x, n, k);
  }
  else if (tolower(ul) == 'u') {
    backsubst_Ut(A, b, tmp, n, k);
    backsubst_U(A, tmp, x, n, k);
  }
  else {
    mexErrMsgTxt("Invalid triangular matrix type: must be either ''U'' or ''L''");
  }

  mxFree(tmp);
}


/* perform a permutation assignment y := x(ind(1:k)) */

void vec_assign(double y[], double x[], mwIndex ind[], mwSize k)
{
  mwIndex i;

  for (i=0; i<k; ++i)
    y[i] = x[ind[i]];
}


/* matrix transpose */

void transpose(double X[], double Y[], mwSize n, mwSize m)
{
  mwIndex i, j, i_m, j_n;
  
  if (n<m) {
    for (j=0; j<m; ++j) {
      j_n = j*n;
      for (i=0; i<n; ++i) {
        Y[j+i*m] = X[i+j_n];
      }
    }
  }
  else {
    for (i=0; i<n; ++i) {
      i_m = i*m;
      for (j=0; j<m; ++j) {
        Y[j+i_m] = X[i+j*n];
      }
    }
  }
}


/* print contents of matrix */

void printmat(double A[], int n, int m, char* matname)
{
  int i, j;
  mexPrintf("\n%s = \n\n", matname);

  if (n*m==0) {
    mexPrintf("   Empty matrix: %d-by-%d\n\n", n, m);
    return;
  }

  for (i=0; i<n; ++i) {
    for (j=0; j<m; ++j)
      mexPrintf("   %lf", A[j*n+i]);
    mexPrintf("\n");
  }
  mexPrintf("\n");
}


/* print contents of sparse matrix */

void printspmat(mxArray *a, char* matname)
{
  mwIndex *aJc = mxGetJc(a);
  mwIndex *aIr = mxGetIr(a);
  double *aPr = mxGetPr(a);

  int i;

  mexPrintf("\n%s = \n\n", matname);

  for (i=0; i<aJc[1]; ++i)
    printf("   (%d,1) = %lf\n", aIr[i]+1,aPr[i]);

  mexPrintf("\n");
}



/* matrix multiplication using Winograd's algorithm */

/*
void mat_mat2(double alpha, double A[], double B[], double X[], mwSize n, mwSize m, mwSize k)
{
  
  mwIndex i1, i2, i3, iX, iA, i2_n;
  double b, *AA, *BB;
  
  AA = mxCalloc(n,sizeof(double));
  BB = mxCalloc(k,sizeof(double));
  
  for (i1=0; i1<n*k; i1++) {
    X[i1] = 0;
  }
  
  for (i1=0; i1<n; ++i1) {
    for (i2=0; i2<m/2; ++i2) {
      AA[i1] += A[i1+2*i2*n]*A[i1+(2*i2+1)*n];
    }
  }

  for (i2=0; i2<k; ++i2) {
    for (i1=0; i1<m/2; ++i1) {
      BB[i2] += B[2*i1+i2*m]*B[2*i1+1+i2*m];
    }
  }

  for (i2=0; i2<k; ++i2) {
    for (i3=0; i3<m/2; ++i3) {
      for (i1=0; i1<n; ++i1) {
        X[i1+i2*n] += (A[i1+(2*i3)*n]+B[2*i3+1+i2*m])*(A[i1+(2*i3+1)*n]+B[2*i3+i2*m]);
      }
    }
  }
  
  if (m%2) {
    for (i2=0; i2<k; ++i2) {
      for (i1=0; i1<n; ++i1) {
        X[i1+i2*n] += A[i1+(m-1)*n]*B[m-1+i2*m];
      }
    }
  }
  
  for (i2=0; i2<k; ++i2) {
    for (i1=0; i1<n; ++i1) {
      X[i1+i2*n] -= (AA[i1] + BB[i2]);
      X[i1+i2*n] *= alpha;
    }
  }
  
  mxFree(AA);
  mxFree(BB);
}
*/




