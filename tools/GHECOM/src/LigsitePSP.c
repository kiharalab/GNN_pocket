/*

 <LigsitePSP.c>

  Calculation of LIGSITE-stylyed ray-based inaccessiblity
  
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
void Cal_LigsitePSP();

/*** FUNCTIONS (LOCAL) ***/
extern void Set_Directional_Vectors();


void Cal_LigsitePSP(Xmap,Ymap,pro_bitmask,tar_bitmask,Ndir,TermType)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 struct CHAR3DMAP *Ymap; /* result PSP value 3D map */
 unsigned char  pro_bitmask; /* protein bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char  tar_bitmask; /* target region bitmask (1 or 2 or 4 or ... 128)*/
 int  Ndir;  /* Number of direction */
 char TermType;  /* 'B'oth end 'O'ne end */
{
 int x,y,z,k,d,psp;
 int xx,yy,zz;
 int dir_vec[98][3];
 char hitPfor, hitPback,end;

 printf("#Cal_LigsitePSP(tar_bitmask %d res_bitmask %d Ndir %d TermType %c)\n", pro_bitmask,tar_bitmask,Ndir,TermType);


 /* [1] Set direction vector */
 Set_Directional_Vectors(Ndir,dir_vec);
 for (k=0;k<Ndir;++k){
  printf("[%d] %d %d %d\n",k,dir_vec[k][0],dir_vec[k][1],dir_vec[k][2]);
 } 

 /* [2] Initialize Ymap to 255 */
 for (x=0;x<Ymap->N[0];++x){
   for (y=0;y<Ymap->N[1];++y){
     for (z=0;z<Ymap->N[2];++z){
       Ymap->map[x][y][z] = 255;
     }
   }
 }


 /* [3] Scaning Xmap and write the PSP value to Ymap */
 for (x=0;x<Xmap->N[0];++x){
  /* if (PAR.OutCalProcess=='T'){if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]);} */
  if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]);
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    
    if ((Xmap->map[x][y][z]&pro_bitmask) == pro_bitmask) Ymap->map[x][y][z] = 255; 

    if ( ((Xmap->map[x][y][z]&tar_bitmask) == tar_bitmask) &&
         ((Xmap->map[x][y][z]&pro_bitmask) != pro_bitmask) ) { 
      psp = 0;
      for (k=0;k<Ndir;++k){
        hitPfor = hitPback = 0; 
        /* forward direction */
        end = 0; d = 0;
        while ((hitPfor==0)&&(end==0)){
	  xx = x + dir_vec[k][0] * d; 
          yy = y + dir_vec[k][1] * d; 
          zz = z + dir_vec[k][2] * d;
          if (  (xx<0)||(xx>=Xmap->N[0])||
                (yy<0)||(yy>=Xmap->N[1])||
                (zz<0)||(zz>=Xmap->N[2])  ) { end = 1;}
          else{
           if ((Xmap->map[xx][yy][zz]&pro_bitmask) == pro_bitmask) hitPfor = 1;
           d += 1;
          }
        } /* while */

        /* backward direction */
        end = 0; d = 0;
        while ((hitPback==0)&&(end==0)){
	  xx = x - dir_vec[k][0] * d; 
          yy = y - dir_vec[k][1] * d; 
          zz = z - dir_vec[k][2] * d;
          if (  (xx<0)||(xx>=Xmap->N[0])||
                (yy<0)||(yy>=Xmap->N[1])||
                (zz<0)||(zz>=Xmap->N[2])  ) { end = 1;}
          else{
           if ((Xmap->map[xx][yy][zz]&pro_bitmask) == pro_bitmask) hitPback = 1;
           d += 1;
          }
        } /* while */

        if (TermType=='B'){ if ((hitPfor==1) && (hitPback==1)) psp += 1; }
        if (TermType=='O'){ psp += (hitPfor + hitPback); }

      } /* k */
       Ymap->map[x][y][z] = psp;
       if (psp==0) Ymap->map[x][y][z] = 255;
     } /* Xmap == 1 */
    } /* z */
  } /* y */
 } /* x */


} /* end of Cal_LigsitePSP() */



void Set_Directional_Vectors(Ndir,dir)
 int Ndir;
 int dir[98][3];
{
 if (Ndir>=7 ){
  dir[0][0] =  1; dir[0][1] =  0; dir[0][2] =  0;
  dir[1][0] =  0; dir[1][1] =  1; dir[1][2] =  0;
  dir[2][0] =  0; dir[2][1] =  0; dir[2][2] =  1;
  dir[3][0] =  1; dir[3][1] =  1; dir[3][2] =  1;
  dir[4][0] = -1; dir[4][1] =  1; dir[4][2] =  1;
  dir[5][0] =  1; dir[5][1] = -1; dir[5][2] =  1;
  dir[6][0] =  1; dir[6][1] =  1; dir[6][2] = -1;
 }
 if (Ndir >= 13){
  dir[7][0]  =  1;  dir[7][1]  =  1;  dir[7][2]  =  0;
  dir[8][0]  =  1;  dir[8][1]  = -1;  dir[8][2]  =  0;
  dir[9][0]  =  0;  dir[9][1]  =  1;  dir[9][2]  =  1;
  dir[10][0] =  0;  dir[10][1] =  1;  dir[10][2] = -1;
  dir[11][0] =  1;  dir[11][1] =  0;  dir[11][2] =  1;
  dir[12][0] =  1;  dir[12][1] =  0;  dir[12][2] = -1;
 }

} /* end of Set_Directional_Vectors() */






