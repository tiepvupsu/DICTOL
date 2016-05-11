/**************************************************************************
 *
 * File name: omputils.c
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 18.8.2009
 *
 *************************************************************************/

#include "omputils.h"
#include <math.h>


const char FULL_GAMMA_STR[] = "full";
const char SPARSE_GAMMA_STR[] = "sparse";


/* convert seconds to hours, minutes and seconds */

void secs2hms(double sectot, int *hrs, int *mins, double *secs)
{
  *hrs = (int)(floor(sectot/3600)+1e-2);
  sectot = sectot - 3600*(*hrs);
  *mins = (int)(floor(sectot/60)+1e-2);
  *secs = sectot - 60*(*mins);
}


/* quicksort, public-domain C implementation by Darel Rex Finley. */
/* modification: sorts the array data[] as well, according to the values in the array vals[] */

#define  MAX_LEVELS  300

void quicksort(mwIndex vals[], double data[], mwIndex n) {
  
  long piv, beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R, swap ;
  double datapiv;
  
  beg[0]=0;
  end[0]=n;
  
  while (i>=0) {
    
    L=beg[i]; 
    R=end[i]-1;
    
    if (L<R) {
      
      piv=vals[L];
      datapiv=data[L];
      
      while (L<R) 
      {
        while (vals[R]>=piv && L<R) 
          R--;
        if (L<R) {
          vals[L]=vals[R];
          data[L++]=data[R];
        }
        
        while (vals[L]<=piv && L<R) 
          L++;
        if (L<R) { 
          vals[R]=vals[L];
          data[R--]=data[L];
        }
      }
      
      vals[L]=piv;
      data[L]=datapiv;
      
      beg[i+1]=L+1;
      end[i+1]=end[i];
      end[i++]=L;
      
      if (end[i]-beg[i] > end[i-1]-beg[i-1]) {
        swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
        swap=end[i]; end[i]=end[i-1]; end[i-1]=swap;
      }
    }
    else {
      i--;
    }
  }
}
