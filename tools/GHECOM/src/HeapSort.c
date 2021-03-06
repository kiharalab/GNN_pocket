/*

<HeapSort.c>

 These functions are based on "Numerical Recipe in C"

 coded by Takeshi Kawabata.


=========== "ghecom" program ===========

Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under
the GNU Lesser General Public License (LGPL) version 3, see LICENSE.txt.

=========== Installation and Usage ===========

See html file under the "doc/" directory.

=========== Citing "ghecom" program ===========

Please cite:

1) Kawabata T. (2010) Detection of multi-scale pockets on protein surfaces using mathematical morphology. Proteins,78, 1195-1121.





>>> CAUTION!! <<<
  This function does not work when N = 1 !! 


*/

#include <stdio.h>
#include <stdlib.h>

void HeapSort();
void HeapSort_Indirect();
void HeapSort_Integer_NoIndex();
void BestIndex();

void HeapSort_Indirect(N,val,index,type)
 int N;
 float *val;
 int  *index;
 char type;
{ int i;
  float *X;

  X = (float *)malloc(sizeof(float)*N); 
  for (i=0;i<N;++i) { index[i] = i; X[i] = val[i];}
  if (N>1) HeapSort(N,X,index,type);
  free(X);

} /* end of HeapSort_Indirect() */







void HeapSort(N,ra,rb,Rtype)
 int N;
 float *ra;  /* array[0..N-1] to be sorted */
 int *rb;     /* array[0..N-1] sorted following ra[] */
 char Rtype;  /* 'R':reverse order */  
{
 int l,j,ir,i;
 float rra;
 int rrb;

 l = (N>>1) + 1;
 ir = N;
 
 for (;;)
 {
  if (l>1) 
   { rra = ra[--l-1];
     rrb = rb[l-1]; 
    }
 else 
  {
    rra = ra[ir-1];
    rrb = rb[ir-1];
    ra[ir-1] = ra[0];
    rb[ir-1] = rb[0];
    if (--ir ==1)
     {   
      ra[0] = rra;
      rb[0] = rrb;
      return;
      }
   }
   i = l;
   j = l << 1;
   while (j <= ir)
   {
     if (Rtype=='R') 
     { 
       if ((j<ir)&&(ra[j-1]>ra[j])) ++j;
       if (rra > ra[j-1])
        {
          ra[i-1] = ra[j-1];
          rb[i-1] = rb[j-1];
          j += (i=j);   
         }
       else j = ir+1;  
      } 
     else 
     {
      if ((j<ir)&&(ra[j-1]<ra[j])) ++j;
      if (rra < ra[j-1])
       {
         ra[i-1] = ra[j-1];
         rb[i-1] = rb[j-1];
         j += (i=j);   
        }
      else j = ir+1;  
     }
   } /* while */
 
    ra[i-1] = rra; 
    rb[i-1] = rrb; 

  } /* for (;;) */

} /* end of HeapSort() */ 


void HeapSort_Integer_NoIndex(N,ra)
 int N;
 int *ra;
{
 int l,j,ir,i,rra;

 l = (N>>1) + 1;
 ir = N;
 
 for (;;)
 {
  if (l>1) 
   { rra = ra[--l-1]; }
 else 
  {
    rra = ra[ir-1];
    ra[ir-1] = ra[0];
    if (--ir ==1)
     {   
      ra[0] = rra;
      return;
      }
   }
   i = l;
   j = l << 1;
   while (j <= ir)
   {
     if ((j<ir)&&(ra[j-1]<ra[j])) ++j;
     if (rra < ra[j-1])
      {
        ra[i-1] = ra[j-1];
        j += (i=j);   
       }
     else j = ir+1;  
    } /* while */
 
    ra[i-1] = rra; 

  } /* for (;;) */

}/* end of HeapSort_Integer_NoIndex() */





void BestIndex(N,val,best_index,type)
 int N;
 float *val;
 int  *best_index;
 char type;
{ int i;
  float BestVal;  

  /*** Return 'Best' index. Not Sorting ***/

  if (type =='R')
  {
    BestVal = val[0]; *best_index = 0;
    for (i=1;i<N;++i)
    { if (val[i]>BestVal) {BestVal = val[i]; *best_index = i;} 

      }
   } 
   else
  {
    BestVal = val[0]; *best_index = 0;
    for (i=1;i<N;++i)
    { if (val[i]<BestVal) {BestVal = val[i]; *best_index  = i;} }
   } 

} /* end of BestIndex() */



/*

>>> This function is only for testing <<<
 
void BubbleSort(N,X,index)
 int N;
 float *X;
 int  *index;
{ int i,Ncha,bi;
  float bx;

  for (i=0;i<N;++i)  index[i] = i; 
  do { Ncha = 0;
       for (i=1;i<N;++i)
       if (X[i-1]<X[i])
        { bi = index[i]; index[i] = index[i-1]; index[i-1] = bi;
          bx = X[i]; X[i] = X[i-1]; X[i-1] = bx;
          ++Ncha; }
  }while (Ncha > 0);

} 

*/
