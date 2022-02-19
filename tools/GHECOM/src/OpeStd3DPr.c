/*

 <OpeStd3DPr.c>
  
  Various operations for 3D-image processing 

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
#include <string.h>
#include <strings.h>
#include <math.h>
#include "globalvar.h"
#include "pdbstruct.h"
#include "PdbIO.h"
#include "Radius.h"
#include "Grid3D.h"
#include "VecFunc.h"

/*** FUNCTIONS (GLOBAL) ***/
int  Extract_Specific_3DMAP();
void Inverse_Specific_3DMAP();
void Inverse_Specific_and_Set_3DMAP();
int  SetZero_Not_Having_Specific_3DMAP();
int  Reset_Specific_3DMAP();
int  SetZero_Specific_3DMAP();
int  SetOne_BitA_for_BitB_One_3DMAP();
int  Extract_NonZero_3DMAP();
int  Count_Specific_3DMAP();
int  Count_NonZero_3DMAP();
int  Count_BitA_and_BitB_3DMAP();
int  Count_BitA_and_Not_BitB_3DMAP();
void Extract_Surface_3DMAP();
void Extract_Complement_Surface_3DMAP();
void Extract_NonZero_Surface_3DMAP();
int  Extract_A_and_B_3DMAP();
int  Extract_A_or_B_3DMAP();
int  Extract_A_and_notB_3DMAP();
int  Extract_notA_and_notB_3DMAP();
int  Extract_A_and_B_and_notC_3DMAP();
int  Extract_A_and_notB_and_notC_3DMAP();
int Increment_Map1_For_Specific_Map2_3DMAP();
void Set_NeighX_NeighY_NeighZ();


/*** FUNCTIONS (LOCAL) ***/


/*** VARIABLES (LOCAL) ***/
/*
static  int  NeighX[6],NeighY[6],NeighZ[6];  
*/


int Extract_Specific_3DMAP(Xmap,tar_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 unsigned char tar_bitmask; /* target bitmask (1 or 2 or 4 or 8 or ... or 128)*/
{
 int x,y,z;
 int  Nsurvive;
 printf("#Extract_Specific_3DMAP(Xmap,tar_bitmask %d)\n",tar_bitmask);

 Nsurvive = 0;

 for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if ((Xmap->map[x][y][z]&tar_bitmask) == tar_bitmask){ 
         Xmap->map[x][y][z] = tar_bitmask; 
         Nsurvive += 1;
       }
       else { 
         if (Xmap->map[x][y][z]>0){
           Xmap->map[x][y][z] = 0; 
         }
       }
      } 
   } 
 } 

 printf("#Extract_Specific_Bitmap Nsurvive:%d\n",Nsurvive);
 return(Nsurvive);

} /* end of Extract_Specific_3DMAP() */


void Inverse_Specific_3DMAP(Xmap,tar_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 unsigned char tar_bitmask; /* target bitmask (1 or 2 or 4 or 8 or ... or 128)*/
{
 int x,y,z;
 printf("#Inverse_Specific_3DMAP(Xmap,tar_bitmask %d)\n",tar_bitmask);

 for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if ((Xmap->map[x][y][z]&tar_bitmask) != tar_bitmask){ 
         Xmap->map[x][y][z] =  Xmap->map[x][y][z]|tar_bitmask;  /* 0 OR 1 => 1 */
       }
       else{
         Xmap->map[x][y][z] =  Xmap->map[x][y][z]^tar_bitmask; /*  1 EXOR 1 => 0 */
       }
     } 
   } 
 } 

} /* end of Inverse_Specific_3DMAP() */



void Inverse_Specific_and_Set_3DMAP(Xmap, tar_bitmask, res_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 unsigned char tar_bitmask; /* target bitmask (1 or 2 or 4 or 8 or ... or 128)*/
 unsigned char res_bitmask; /* result bitmask (1 or 2 or 4 or ...  or 128)*/
{
 int x,y,z;
 printf("#Inverse_Specific_and_Set_3DMAP(Xmap,tar_bitmask %d)\n",tar_bitmask);

 for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if ((Xmap->map[x][y][z]&tar_bitmask) != tar_bitmask){ 
         Xmap->map[x][y][z] =  Xmap->map[x][y][z]|res_bitmask;  /* 0 OR 1 => 1 */
       }
     } 
   } 
 } 

} /* end of Inverse_Specific_and_Set_3DMAP() */




int SetZero_Not_Having_Specific_3DMAP(Xmap,tar_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 unsigned char tar_bitmask; /* target bitmask (1 or 2 or 4 or 8 or ... or 128)*/
{
 int x,y,z;
 int  Nsurvive;
 printf("#SetZero_Not_Having_Specific_3DMAP(tar_bitmask %d)\n",tar_bitmask);

 Nsurvive = 0;

 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
     if ((Xmap->map[x][y][z]&tar_bitmask) != tar_bitmask) 
       Xmap->map[x][y][z] = 0; 
     else 
       ++Nsurvive;
    } /* z */
  } /* y */
 } /* x */

 printf("#SetZero_Not_Having_Specific_Bitmap Nsurvive:%d\n",Nsurvive);
 return(Nsurvive);

} /* end of SetZero_Not_Having_Specific_3DMAP() */






int Reset_Specific_3DMAP(Xmap,tar_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 unsigned char tar_bitmask; /* target bitmask (1 or 2 or 4 or 8 or ... or 128)*/
 /**** Set the bit zero for the specific tar_bitmask.  ****/
{
 int x,y,z;
 int  Nreset;
 unsigned char compl_tar_bitmask;
 compl_tar_bitmask = ~tar_bitmask;

/* 
 printf("#Reset_Specific_3DMAP(Xmap,tar_bitmask %d compl_tar_bitmask %d)\n",tar_bitmask,compl_tar_bitmask);
 */

 Nreset = 0;

 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
     if ((Xmap->map[x][y][z]&tar_bitmask) == tar_bitmask) 
      { 
       /* printf(" %d -->",Xmap->map[x][y][z]); */
       Xmap->map[x][y][z] = Xmap->map[x][y][z]&compl_tar_bitmask; 
       /* printf(" %d\n",Xmap->map[x][y][z]); */
       ++Nreset;}
     /* ^: exclusive or */
    } /* z */
  } /* y */
 } /* x */
/*
 printf("#Reset_Specific_Bitmap Nsurvive:%d\n",Nreset);
*/
  return(Nreset);
} /* end of Reset_Specific_3DMAP() */




int SetZero_Specific_3DMAP(Xmap,tar_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 unsigned char tar_bitmask; /* target bitmask (1 or 2 or 4 or 8 or ... or 128)*/
 /**** Set the map[][][]:=zero, for pixel having the specific tar_bitmask.  ****/
{
 int x,y,z;
 int  Nreset;
 printf("#SetZero_Specific_3DMAP(Xmap,tar_bitmask %d)\n",tar_bitmask);

 Nreset = 0;

 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
     if ((Xmap->map[x][y][z]&tar_bitmask) == tar_bitmask) { Xmap->map[x][y][z] = 0; Nreset += 1;}
    } /* z */
  } /* y */
 } /* x */

 printf("#SetZero_Specific_Bitmap Nsurvive:%d\n",Nreset);
 return(Nreset);
} /* end of SetZero_Specific_3DMAP() */


int SetOne_BitA_for_BitB_One_3DMAP(Xmap,A_bitmask,B_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 unsigned char A_bitmask; /* target bitmask (1 or 2 or 4 or 8 or ... or 128)*/
 unsigned char B_bitmask; /* target bitmask (1 or 2 or 4 or 8 or ... or 128)*/
 /**** Set the map[][][]:=zero, for pixel having the specific tar_bitmask.  ****/
{
 int x,y,z;
 int  Nset;
 /*
 printf("#SetOne_BitA_for_BitB_One_3DMAP(A_bitmask %d B_bitmask %d)\n",
 A_bitmask,B_bitmask);
 */

 Nset = 0;

 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
     if ((Xmap->map[x][y][z]&B_bitmask) == B_bitmask){ 
        Xmap->map[x][y][z] = Xmap->map[x][y][z]|A_bitmask ; 
        ++Nset;
      }
    } /* z */
  } /* y */
 } /* x */

 /*
 printf("#SetZero_Specific_Bitmap Nsurvive:%d\n",Nset);
 */ 
 return(Nset);
} /* end of SetZero_Specific_3DMAP() */





int Extract_NonZero_3DMAP(Xmap,res_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 unsigned char res_bitmask; /* result bitmask (1 or 2 or 4 or 8 or ... or 128)*/
{
 int x,y,z;
 int  Nsurvive;
 
 Nsurvive = 0;
 for (x=0;x<Xmap->N[0];++x){
   if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]);
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if (Xmap->map[x][y][z]!=0){ 
         Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask; 
         ++Nsurvive;
       }
    } 
  } 
 } 

 printf("#Extract_NonZero_Bitmap Nsurvive:%d\n",Nsurvive);
 return(Nsurvive);

} /* end of Extract_NonZero_3DMAP() */






int Count_Specific_3DMAP(Xmap,tar_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 unsigned char tar_bitmask; /* target bitmask (1 or 2 or 4 or 8 or ... or 128)*/
{
 int x,y,z;
 int  N;
 N = 0;
 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    if ((Xmap->map[x][y][z]&tar_bitmask) == tar_bitmask) ++N; 
    } 
  } 
 }
 return(N);
} /* end of Count_Specific_3DMAP() */



int Count_NonZero_3DMAP(Xmap)
 struct CHAR3DMAP *Xmap; /* target 3D map */
{
 int x,y,z;
 int  N;
 N = 0;
 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    if (Xmap->map[x][y][z] != 0) ++N; 
    } 
  } 
 }
 return(N);
} /* end of Count_NonZero_3DMAP() */



int Count_BitA_and_BitB_3DMAP(Xmap,A_bitmask,B_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 unsigned char A_bitmask; /* target bitmask (1 or 2 or 4 or 8 or ... or 128)*/
 unsigned char B_bitmask; /* target bitmask (1 or 2 or 4 or 8 or ... or 128)*/
{
 int x,y,z;
 int  N;
 N = 0;
 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    if (((Xmap->map[x][y][z]&A_bitmask) == A_bitmask)  &&
        ((Xmap->map[x][y][z]&B_bitmask) == B_bitmask) ) ++N; 
    } 
  } 
 }
 return(N);
} /* end of Count_BitA_and_BitB_3DMAP() */




int Count_BitA_and_Not_BitB_3DMAP(Xmap,A_bitmask,B_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 unsigned char A_bitmask; /* target bitmask (1 or 2 or 4 or 8 or ... or 128)*/
 unsigned char B_bitmask; /* target bitmask (1 or 2 or 4 or 8 or ... or 128)*/
{
 int x,y,z;
 int  N;
 N = 0;
 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    if (((Xmap->map[x][y][z]&A_bitmask) == A_bitmask)  &&
        ((Xmap->map[x][y][z]&B_bitmask) != B_bitmask) ) ++N; 
    } 
  } 
 }
 return(N);
} /* end of Count_BitA_and_Not_BitB_3DMAP() */






void Extract_Surface_3DMAP(Xmap,tar_bitmask,res_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 unsigned char tar_bitmask; /* target bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char res_bitmask; /* result bitmask (1 or 2 or 4 or ... or 128)*/
{
 int x,y,z,xx,yy,zz;
 int  Ntar_nei,n;
 printf("#Extract_Surface_3DMAP(Xmap,tar_bitmask %d)\n",res_bitmask);

 PAR.NeighborNum = 6; 
 Set_NeighX_NeighY_NeighZ(PAR.NeighborNum,PAR.NeighX,PAR.NeighY,PAR.NeighZ); 
 

 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
     if ((Xmap->map[x][y][z]&tar_bitmask) == tar_bitmask){
      Ntar_nei = 0;
      for (n=0;n< PAR.NeighborNum;++n){
          xx = x + PAR.NeighX[n];
          yy = y + PAR.NeighY[n];
          zz = z + PAR.NeighZ[n];
          if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){ 
            if ((Xmap->map[xx][yy][zz]&tar_bitmask) == tar_bitmask) ++Ntar_nei;
          }
         } /* n */

       if (Ntar_nei<PAR.NeighborNum) Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask ; 
 
       } /* Xmap == 1 */
     } /* z */
  } /* y */
 } /* x */

} /* end of Extract_Surface_3DMAP() */



void Extract_Complement_Surface_3DMAP(Xmap,tar_bitmask,res_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 unsigned char tar_bitmask; /* target bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char res_bitmask; /* result bitmask (1 or 2 or 4 or ... or 128)*/
{
 int x,y,z,xx,yy,zz;
 int  Nnontar_nei,Nnei;
 int  n;
 printf("#Extract_Complement_Surface_3DMAP(Xmap,tar_bitmask %d)\n",res_bitmask);

 PAR.NeighborNum = 6;
 Set_NeighX_NeighY_NeighZ(PAR.NeighborNum,PAR.NeighX, PAR.NeighY, PAR.NeighZ);

 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
     if ((Xmap->map[x][y][z]&tar_bitmask) != tar_bitmask){
      Nnei = Nnontar_nei =  0;
      for (n=0;n< PAR.NeighborNum;++n)
       {
          xx = x + PAR.NeighX[n];
          yy = y + PAR.NeighY[n];
          zz = z + PAR.NeighZ[n];
          if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){ 
           ++Nnei;
           if ((Xmap->map[xx][yy][zz]&tar_bitmask) != tar_bitmask) ++Nnontar_nei;
          }
         } /* n */

       if (Nnontar_nei<Nnei) Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask ; 
 
       } /* Xmap == 1 */
     } /* z */
  } /* y */
 } /* x */

} /* end of Extract_Complement_Surface_3DMAP() */










void Extract_NonZero_Surface_3DMAP(Xmap)
 struct CHAR3DMAP *Xmap; /* target 3D map */
{
 int x,y,z,xx,yy,zz;
 int  Ntar_nei,n;
 printf("#Extract_NonZero_Surface_3DMAP()\n");

 /** [1] Find completely surrounded grid to 255 **/
 PAR.NeighborNum = 6; 
 Set_NeighX_NeighY_NeighZ(PAR.NeighborNum,PAR.NeighX,PAR.NeighY,PAR.NeighZ); 


 xx = yy = zz = 0;
 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
     if (Xmap->map[x][y][z]!=0){
      Ntar_nei = 0;
      for (n=0;n<PAR.NeighborNum;++n){
          xx = x + PAR.NeighX[n];
          yy = y + PAR.NeighY[n];
          zz = z + PAR.NeighZ[n];
          if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){ 
           if (Xmap->map[xx][yy][zz] !=0) ++Ntar_nei;
          }
         } /* n */

       if (Ntar_nei == PAR.NeighborNum) Xmap->map[x][y][z] = 255;
 
       } /* Xmap == 1 */
     } /* z */
  } /* y */
 } /* x */


 /** [2] Change grids with 255 to 0 **/

 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
      if (Xmap->map[x][y][z] ==255) Xmap->map[xx][yy][zz] = 0;
     } /* z */
   } /* y */
 } /* x */

} /* end of Extract_NonZero_Surface_3DMAP() */







int Extract_A_and_B_3DMAP(Xmap,Abitmask,Bbitmask,res_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 unsigned char Abitmask; /* target bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char Bbitmask; /* substract bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char res_bitmask; /* result bitmask (1 or 2 or 4 or ...  or 128)*/
 /* tar - sub  => res */
{
 int x,y,z,Nsurvive;
 printf("#Extract_A_and_B_3DMAP(Xmap,A %d B %d res %d)\n",
 Abitmask, Bbitmask, res_bitmask);

 Nsurvive = 0;

 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    if ( ((Xmap->map[x][y][z]&Abitmask) == Abitmask) 
         &&((Xmap->map[x][y][z]&Bbitmask) == Bbitmask) ) {  
        Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask;
         ++Nsurvive;
        } 
    } 
  } 
 } 
 return(Nsurvive);
} /* end of Extract_A_and_B_3DMAP() */





int Extract_A_or_B_3DMAP(Xmap,Abitmask,Bbitmask,res_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 unsigned char Abitmask; /* target bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char Bbitmask; /* substract bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char res_bitmask; /* result bitmask (1 or 2 or 4 or ...  or 128)*/
 /* tar - sub  => res */
{
 int x,y,z,Nsurvive;
 printf("#Extract_A_or_B_3DMAP(Xmap,A %d B %d res %d)\n",
 Abitmask, Bbitmask, res_bitmask);

 Nsurvive = 0;

 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    if ( ((Xmap->map[x][y][z]&Abitmask) == Abitmask) 
         ||((Xmap->map[x][y][z]&Bbitmask) == Bbitmask) ) {  
        Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask;
         ++Nsurvive; 
      } 
    } 
  } 
 } 
 return(Nsurvive);
} /* end of Extract_A_or_B_3DMAP() */




int Extract_A_and_notB_3DMAP(Xmap,Abitmask,Bbitmask,res_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 unsigned char Abitmask; /* target bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char Bbitmask; /* substract bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char res_bitmask; /* result bitmask (1 or 2 or 4 or ...  or 128)*/
 /* tar - sub  => res */
{
 int x,y,z,Nsurvive;
 printf("#Extract_A_and_notB_3DMAP(Xmap,A %d B %d res %d)\n",
 Abitmask, Bbitmask, res_bitmask);

 Nsurvive = 0;

 for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if ( ((Xmap->map[x][y][z]&Abitmask) == Abitmask) 
          &&((Xmap->map[x][y][z]&Bbitmask) != Bbitmask) ) {  
         Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask;
         ++Nsurvive; } 
     } 
   }  
 } 
 return(Nsurvive);
} /* end of Extract_A_and_notB_3DMAP() */



int Extract_notA_and_notB_3DMAP(Xmap,Abitmask,Bbitmask,res_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 unsigned char Abitmask; /* target bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char Bbitmask; /* substract bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char res_bitmask; /* result bitmask (1 or 2 or 4 or ...  or 128)*/
 /* tar - sub  => res */
{
 int x,y,z,Nsurvive;

 printf("#Extract_notA_and_notB_3DMAP(Xmap,A %d B %d res %d)\n", Abitmask, Bbitmask, res_bitmask);

 Nsurvive = 0;

 for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if ( ((Xmap->map[x][y][z]&Abitmask) != Abitmask) 
          &&((Xmap->map[x][y][z]&Bbitmask) != Bbitmask) ) {  
         Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask;
         Nsurvive += 1; 
       } 
     } 
   }  
 } 
 return(Nsurvive);
} /* end of Extract_notA_and_notB_3DMAP() */








int Extract_A_and_B_and_notC_3DMAP(Xmap,Abitmask,Bbitmask,Cbitmask,res_bitmask)
 struct CHAR3DMAP *Xmap;    /* target 3D map */
 unsigned char Abitmask;    /* bitmask for A(1 or 2 or 4 or ... or 128)*/
 unsigned char Bbitmask;    /* bitmask for B(1 or 2 or 4 or ... or 128)*/
 unsigned char Cbitmask;    /* bitmask for C(1 or 2 or 4 or ... or 128)*/
 unsigned char res_bitmask; /* result bitmask (1 or 2 or 4 or ...  or 128)*/
 /* tar - sub  => res */
{
 int x,y,z,Nsurvive;
 
 printf("##Extract_A_and_B_and_notC_3DMAP(Xmap,A %d B %d C %d res %d)\n",
  Abitmask, Bbitmask, Cbitmask, res_bitmask);   

 Nsurvive = 0;

 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    if ( ((Xmap->map[x][y][z]&Abitmask) == Abitmask)  &&
         ((Xmap->map[x][y][z]&Bbitmask) == Bbitmask)  &&
         ((Xmap->map[x][y][z]&Cbitmask) != Cbitmask)  ){
         Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask;
         ++Nsurvive;} 
    } /* z */
  } /* y */
 } /* x */

 return(Nsurvive);

} /* end of Extract_A_and_B_and_notC_3DMAP() */



int Extract_A_and_notB_and_notC_3DMAP(Xmap,Abitmask,Bbitmask,Cbitmask,res_bitmask)
 struct CHAR3DMAP *Xmap;    /* target 3D map */
 unsigned char Abitmask;    /* bitmask for A(1 or 2 or 4 or ... or 128)*/
 unsigned char Bbitmask;    /* bitmask for B(1 or 2 or 4 or ... or 128)*/
 unsigned char Cbitmask;    /* bitmask for C(1 or 2 or 4 or ... or 128)*/
 unsigned char res_bitmask; /* result bitmask (1 or 2 or 4 or ...  or 128)*/
 /* tar - sub  => res */
{
 int x,y,z,Nsurvive;
 
 printf("##Extract_A_and_notB_and_notC_3DMAP(Xmap,A %d B %d C %d res %d)\n",
  Abitmask, Bbitmask, Cbitmask, res_bitmask);   

 Nsurvive = 0;

 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    if ( ((Xmap->map[x][y][z]&Abitmask) == Abitmask)  &&
         ((Xmap->map[x][y][z]&Bbitmask) != Bbitmask)  &&
         ((Xmap->map[x][y][z]&Cbitmask) != Cbitmask)  ){
         Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask;
         ++Nsurvive;} 
    } /* z */
  } /* y */
 } /* x */

 return(Nsurvive);

} /* end of Extract_A_and_notB_and_notC_3DMAP() */




int Increment_Map1_For_Specific_Map2_3DMAP(Map1,Map2,tar_bitmask)
 struct CHAR3DMAP *Map1; /* 3D map for incremental count*/
 struct CHAR3DMAP *Map2; /* target 3D map */
 unsigned char tar_bitmask; /* target bitmask (1 or 2 or 4 or 8 or ... or 128)*/
{
 int x,y,z,N;
 /*
 printf("#Extract_Specific_3DMAP(Xmap,tar_bitmask %d)\n",tar_bitmask);
 */
 N = 0;
 for (x=0;x<Map1->N[0];++x){
  for (y=0;y<Map1->N[1];++y){
   for (z=0;z<Map1->N[2];++z){ 
     if ((Map2->map[x][y][z]&tar_bitmask) == tar_bitmask) 
      { 
/*
        if (Map1->map[x][y][z]<255) {Map1->map[x][y][z] += 1;} 
*/ 
        Map1->map[x][y][z] += 1; 
       ++N;
       }
    } /* z */
  } /* y */
 } /* x */
 return(N);
} /* end of Increment_Xmap_For_Specific_Ymap_3DMAP() */



void Set_NeighX_NeighY_NeighZ(NeiNum,Nx,Ny,Nz)
 int NeiNum; /* Number of Neighbor. 6 or 18 or 26 */
 int Nx[26],Ny[26],Nz[26];
{
 int n;
 /*
 printf("#Set_NeighX_NeighY_NeighZ(NeiNum %d)\n",NeiNum);
 */
 
 n = 0;

 if (NeiNum>=6){
  Nx[ 0] =  1; Ny[ 0] =  0; Nz[ 0] =  0;
  Nx[ 1] = -1; Ny[ 1] =  0; Nz[ 1] =  0;
  Nx[ 2] =  0; Ny[ 2] =  1; Nz[ 2] =  0;
  Nx[ 3] =  0; Ny[ 3] = -1; Nz[ 3] =  0;
  Nx[ 4] =  0; Ny[ 4] =  0; Nz[ 4] =  1;
  Nx[ 5] =  0; Ny[ 5] =  0; Nz[ 5] = -1 ;
 }
 if (NeiNum>=18){
  Nx[ 6] =  1; Ny[ 6] =  1; Nz[ 6] =  0;
  Nx[ 7] =  0; Ny[ 7] =  1; Nz[ 7] =  1;
  Nx[ 8] =  1; Ny[ 8] =  0; Nz[ 8] =  1;
  Nx[ 9] = -1; Ny[ 9] = -1; Nz[ 9] =  0;
  Nx[10] =  0; Ny[10] = -1; Nz[10] = -1;
  Nx[11] = -1; Ny[11] =  0; Nz[11] = -1;
  Nx[12] =  1; Ny[12] = -1; Nz[12] =  0;
  Nx[13] =  0; Ny[13] =  1; Nz[13] = -1;
  Nx[14] =  1; Ny[14] =  0; Nz[14] = -1;
  Nx[15] = -1; Ny[15] =  1; Nz[15] =  0;
  Nx[16] =  0; Ny[16] = -1; Nz[16] =  1;
  Nx[17] = -1; Ny[17] =  0; Nz[17] =  1;
 }
 if (NeiNum>=26){
  Nx[18] =  1; Ny[18] =  1; Nz[18] =  1;
  Nx[19] =  1; Ny[19] =  1; Nz[19] = -1;
  Nx[20] =  1; Ny[20] = -1; Nz[20] =  1;
  Nx[21] =  1; Ny[21] = -1; Nz[21] = -1;
  Nx[22] = -1; Ny[22] =  1; Nz[22] =  1;
  Nx[23] = -1; Ny[23] =  1; Nz[23] = -1;
  Nx[24] = -1; Ny[24] = -1; Nz[24] =  1;
  Nx[25] = -1; Ny[25] = -1; Nz[25] = -1;
 }

 /*
 for (n=0;n<NeiNum;++n){
   printf("[%2d] Nx %-3d Ny %-3d Nz %-3d\n",n,Nx[n],Ny[n],Nz[n]);
 }
 */

 if ((NeiNum != 6)&&(NeiNum != 18)&&(NeiNum != 26) ){
  printf("#ERROR:Improper NeighborNum %d\n",NeiNum);
  exit(1);
 }

} /* end of Set_NeighX_NeighY_NeighZ() */
