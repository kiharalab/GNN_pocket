/*

 <VecFunc.c>

 Functions for vector manupulations

=========== "ghecom" program ===========

Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under
the GNU Lesser General Public License (LGPL) version 3, see LICENSE.txt.

=========== Installation and Usage ===========

See html file under the "doc/" directory.

=========== Citing "ghecom" program ===========

Please cite:

1) Kawabata T. (2010) Detection of multi-scale pockets on protein surfaces using mathematical morphology. Proteins,78, 1195-1121.


*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>

/*** FUNCTIONS (GLOBAL) */
void Cross_Vec();
float Dot_Prod();
void Sub_Vec(); 
void Equal_Vec();
void Min_Vec();
void Max_Vec();
float Normalize();
float Vec_Length();
void Print_Vec();
void   Make_Random_Rotation_Matrix();
void   Set_Unitary_Matrix();
void Make_RotMatrix_By_Quaternion();



void Cross_Vec(C,A,B)   /* C = A x B */
 float C[3],A[3],B[3];
{ C[0] = A[1]*B[2] - A[2]*B[1];
  C[1] = A[2]*B[0] - A[0]*B[2];
  C[2] = A[0]*B[1] - A[1]*B[0]; }


                                                                                
float Dot_Prod(A,B)  /* return A dot B */
 float A[3],B[3];
{ int i;
  float prod;
  prod = 0.0;
  for (i=0;i<3;++i){
     prod += A[i]*B[i];
  }
  return(prod);
} /* end of Dot_Prod() */



void Sub_Vec(C,A,B)    /* C = A - B */
 float C[3],A[3],B[3];
{  C[0] = A[0]-B[0];
   C[1] = A[1]-B[1];
   C[2] = A[2]-B[2];
} /* end of Sub_Vec(A,B) */

                                                                                
void Equal_Vec(A,B)  /* A := B */
 float A[3],B[3];
{  A[0] = B[0];
   A[1] = B[1];
   A[2] = B[2];
} /* end of Equal_Vec(A,B) */


void Min_Vec(C,A,B)  /* C = min[A, B] */
 float C[3],A[3],B[3];
{ 
  int i;
  for (i=0;i<3;++i){
    if (A[i]<B[i]){
      C[i] = A[i];
    }
    else{
      C[i] = B[i];
    }
  }
} /* end of Min_Vec(A,B) */


void Max_Vec(C,A,B)  /* C = max[A, B] */
 float C[3],A[3],B[3];
{ 
  int i;
  for (i=0;i<3;++i){
    if (A[i]>B[i]){
      C[i] = A[i];
    }
    else{
      C[i] = B[i];
    }
  }
} /* end of Max_Vec(A,B) */
                                                                                
                                                                                
float Normalize(A)
 float A[3];
{  float dis;
  dis = A[0]*A[0] + A[1]*A[1] + A[2]*A[2];
  if (dis>0.0)
   { dis = sqrt(dis);
     A[0] /= dis; A[1] /= dis; A[2] /= dis; }
  return(dis);

} /* end Normalize() */


float Vec_Length(A)
 float A[3];
{ float dis;
 
  dis = A[0]*A[0] + A[1]*A[1] + A[2]*A[2];
  if (dis>0.0){ dis = sqrt(dis);}
  return(dis);

} /* end Vec_Length() */


void Make_RotMatrix_By_Quaternion(R,W,X,Y,Z)
 float R[3][3];
 float W,X,Y,Z;  /* Quaternion */
{
 double len,w,x,y,z;
 /* Normalize (W,X,Y,Z) */
 len = W*W + X*X + Y*Y + Z*Z;
 if (len>0.0)
  { len = sqrt(len);
    w = W/len;
    x = X/len;
    y = Y/len;
    z = Z/len;  }
 else
 { w = W; x = X; y = Y; z = Z; }
 R[0][0] = 1.0 - 2.0*(y*y+z*z);
 R[1][0] = 2.0*(x*y+w*z);
 R[2][0] = 2.0*(x*z-w*y);
 R[0][1] = 2.0*(x*y-w*z);
 R[1][1] = 1.0 - 2.0*(x*x+z*z);
 R[2][1] = 2.0*(y*z+w*x);
 R[0][2] = 2.0*(x*z+w*y);
 R[1][2] = 2.0*(y*z-w*x);
 R[2][2] = 1.0 - 2.0*(x*x+y*y);
} /* end of Make_RotMatrix_By_Quaternion() */

void Make_Random_Rotation_Matrix(R)
 float R[3][3];
{
 float w,x,y,z;
 float len,theta2,cos_t2,sin_t2,norm;

 theta2 = 2.0*M_PI*(double)rand()/(double)RAND_MAX / 2.0;

 cos_t2 = cos(theta2);
 sin_t2 = sin(theta2);

 len = 3.0; 
 while (len>1.0){
  x = 2.0*(double)rand()/(double)RAND_MAX - 1.0;
  y = 2.0*(double)rand()/(double)RAND_MAX - 1.0;
  z = 2.0*(double)rand()/(double)RAND_MAX - 1.0;
  len = x*x+y*y+z*z;
 }

 if (len>0.0) len = sqrt(len);

 w  = cos_t2;
 norm = sin_t2/len; 
 x *= norm; 
 y *= norm; 
 z *= norm; 

/*
 printf("#quaternion : w %f x %f y %f z %f\n",w,x,y,z);
 */
  Make_RotMatrix_By_Quaternion(R,w,x,y,z);
 /*
 printf("#R0 %f %f %f\n",R[0][0],R[0][1],R[0][2]);
 printf("#R1 %f %f %f\n",R[1][0],R[1][1],R[1][2]);
 printf("#R2 %f %f %f\n",R[2][0],R[2][1],R[2][2]);
 */
} /* end of Make_Random_Rotation_Matrix() */



void Set_Unitary_Matrix(R)
 float R[3][3];
{
 int i,j;
 for (i=0;i<3;++i){
   for (j=0;j<3;++j){
     if (i==j) R[i][j] = 1.0; else R[i][j] = 0.0;
   }
  }
} /* end of Set_Unitary_Matrix() */

