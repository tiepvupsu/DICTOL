/**************************************************************************
 *
 * File name: omputils.h
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 18.8.2009
 *
 * Utility definitions and functions for the OMP library.
 *
 *************************************************************************/


#ifndef __OMP_UTILS_H__
#define __OMP_UTILS_H__

#include "mex.h"


/* constants for the representation mode of gamma */

extern const char FULL_GAMMA_STR[];      /* "full" */
extern const char SPARSE_GAMMA_STR[];    /* "sparse" */


#define FULL_GAMMA 1
#define SPARSE_GAMMA 2
#define INVALID_MODE 3



/**************************************************************************
 * Memory management for OMP2.
 *
 * GAMMA_INC_FACTOR:
 * The matrix GAMMA is allocated with sqrt(n)/2 coefficients per signal,
 * for a total of nzmax = L*sqrt(n)/2 nonzeros. Whenever GAMMA needs to be
 * increased, it is increased by a factor of GAMMA_INC_FACTOR.
 *
 * MAT_INC_FACTOR:
 * The matrices Lchol, Gsub and Dsub are allocated with sqrt(n)/2
 * columns each. If additional columns are needed, this number is 
 * increased by a factor of MAT_INC_FACTOR.
 **************************************************************************/

#define GAMMA_INC_FACTOR (1.4)
#define MAT_INC_FACTOR (1.6)



/**************************************************************************
 * Convert number of seconds to hour, minute and second representation.
 *
 * Parameters:
 *   sectot - total number of seconds
 *   hrs, mins, secs - output hours (whole) and minutes (whole) and seconds
 *
 **************************************************************************/
void secs2hms(double sectot, int *hrs, int *mins, double *secs);



/**************************************************************************
 * QuickSort - public-domain C implementation by Darel Rex Finley.
 *
 * Modified to sort both the array vals[] and the array data[] according 
 * to the values in the array vals[].
 *
 **************************************************************************/
void quicksort(mwIndex vals[], double data[], mwIndex n);


#endif

