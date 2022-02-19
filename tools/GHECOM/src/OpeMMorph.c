/*

 <OpeMMorph.c>
  
  Basic operations for mathematical morphology


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
int   Make_Sphere_Probe_Map();
void  Dilation();
int   Erosion();
int   Bounded_Erosion();
int   Opening();
void  Reflection();

/*** FUNCTIONS (LOCAL) ***/


int  Make_Sphere_Probe_Map(M,Radius,grid_width)
 struct CHAR3DMAP *M;
 float  Radius;
 float  grid_width;
{
 int Ngrid,Ngrid_half,x,y,z,Nvoxel_foreground;
 float Pos[3],DD,RR; 

 Ngrid_half = (int)ceil((double)Radius/(double)grid_width); 
 
 Ngrid = 2 * Ngrid_half + 1;
 M->OrigPos[0] = 0.0; M->OrigPos[1] = 0.0; M->OrigPos[2] = 0.0;
 M->N[0] = M->N[1] = M->N[2] = Ngrid;

 Malloc_CHAR3DMAP(M,M->N[0],M->N[1],M->N[2]);
 
 RR = Radius*Radius;
 Nvoxel_foreground = 0;
 for (x=-Ngrid_half;x<=Ngrid_half;++x){
   Pos[0] = M->OrigPos[0] + grid_width * x;
   for (y=-Ngrid_half;y<=Ngrid_half;++y){ 
     Pos[1] = M->OrigPos[1] + grid_width * y;
     for (z=-Ngrid_half;z<=Ngrid_half;++z){ 
       Pos[2] = M->OrigPos[2] + grid_width * z;
       DD = Pos[0]*Pos[0] + Pos[1]*Pos[1] + Pos[2]*Pos[2];
       /* printf("x %d y %d z %d\n",x+Ngrid_half,y+Ngrid_half,z+Ngrid_half); */
       if (DD<=RR){ 
         M->map[x+Ngrid_half][y+Ngrid_half][z+Ngrid_half] = 1;  
         Nvoxel_foreground += 1;
       }
      } 
   } 
 } 
 printf("#Make_Sphere_Probe_Map(Radius %f grid_width %f) -> MapSize [%d %d %d] Nvoxel_foreground %d\n",
         Radius,grid_width,Ngrid,Ngrid,Ngrid,Nvoxel_foreground);
 
if (Nvoxel_foreground<1){
   printf("#ERROR: Sphere_Probe_Map() with Nvoxel_foreground = %d is not proper. Smaller grid size(-gw) or Large Rsmall (-rs) is recommeded.\n",Nvoxel_foreground);
   exit(1);
 }

 return(Nvoxel_foreground);

} /* end of Make_Sphere_Probe_Map() */





void Dilation(Xmap,Pmap,tar_bitmask,res_bitmask,Otype)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 struct CHAR3DMAP *Pmap; /* probe 3D map */
 unsigned char   tar_bitmask; /* target bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char   res_bitmask; /* result bitmask (1 or 2 or 4 or ... 128)*/
 char    Otype; /* 'O': output messages */
{
 int x,y,z,p,q,r;
 int xx,yy,zz;
 int Nhalf[3];

 if (Otype == 'O') printf("#Dilation(Xmap,Pmap tar_bitmask %d res_bitmask %d)\n",tar_bitmask,res_bitmask);

 Nhalf[0] = (Pmap->N[0]-1)/2;
 Nhalf[1] = (Pmap->N[1]-1)/2;
 Nhalf[2] = (Pmap->N[2]-1)/2;

 for (x=0;x<Xmap->N[0];++x){
  if (Otype=='O'){if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]);}
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    
   if ((Xmap->map[x][y][z]&tar_bitmask) == tar_bitmask){
     for (p=0;p<Pmap->N[0];++p){
      for (q=0;q<Pmap->N[1];++q){
       for (r=0;r<Pmap->N[2];++r){
        if (Pmap->map[p][q][r]==1){
          xx = x + p - Nhalf[0];
          yy = y + q - Nhalf[1];
          zz = z + r - Nhalf[2];
          /* printf("x %d p %d xx %d\n",x,p,xx); */
          if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){
           Xmap->map[xx][yy][zz] = Xmap->map[xx][yy][zz]|res_bitmask;
          }
           else {
           /*
               printf("#Dilation:out of box !! xyz %d %d %d Xmap %d xx %d yy %d zz %d\n",
               x,y,z,Xmap->map[x][y][z],xx,yy,zz);
           */
          }
         } /* Pmap == 1 */
        } /* r */
       } /* q */
      } /* p */
     } /* Xmap == 1 */
    } /* z */
  } /* y */
 } /* x */


} /* end of Dilation() */





int Erosion(Xmap,Pmap,tar_bitmask,res_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 struct CHAR3DMAP *Pmap; /* probe 3D map */
 unsigned char   tar_bitmask; /* target bitmask (1 or 2 or 4 or 8 or 16 or 32 or 64 or 128)*/
 unsigned char   res_bitmask; /* result bitmask (1 or 2 or 4 or 8 or 16 or 32 or 64 or 128)*/
{
 int x,y,z,p,q,r;
 int xx,yy,zz;
 int Nhalf[3];
 char overlap_all;
 int  Nsurvive,Ntry;
 printf("#Erosion(Xmap,Pmap tar_bitmask %d res_bitmask %d)\n",tar_bitmask,res_bitmask);

 Nhalf[0] = (Pmap->N[0]-1)/2;
 Nhalf[1] = (Pmap->N[1]-1)/2;
 Nhalf[2] = (Pmap->N[2]-1)/2;
 Nsurvive = Ntry = 0;

 for (x=0;x<Xmap->N[0];++x){
   if ((x%10)==0){printf("#x %d/%d\n",x,Xmap->N[0]);}
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if ((Xmap->map[x][y][z]&tar_bitmask) == tar_bitmask){
         ++Ntry;
         overlap_all = 1;
         p = 0; 
         while ((p<Pmap->N[0])&&(overlap_all==1)){
           q = 0;
           while ((q<Pmap->N[1])&&(overlap_all==1)){
             r = 0;
             while ((r<Pmap->N[2])&&(overlap_all==1)){
               if (Pmap->map[p][q][r]==1){
                 xx = x + p - Nhalf[0];
                 yy = y + q - Nhalf[1];
                 zz = z + r - Nhalf[2];
                 if ((xx<0)||(xx>=Xmap->N[0])||(yy<0)||(yy>=Xmap->N[1])||(zz<0)||(zz>=Xmap->N[2])){ 
                   overlap_all = 0;
            /* 
             printf("#Erosion:out of box !! xyz %d %d %d Xmap %d xx %d yy %d zz %d\n",
               x,y,z,Xmap->map[x][y][z],xx,yy,zz); 
             */
                 }
                 else{ 
                   if ((Xmap->map[xx][yy][zz]&tar_bitmask) != tar_bitmask){
                    overlap_all = 0;
                   }
                 }
               } /* Pmap == 1 */
               r += 1;
             } /* r */
           q += 1;
           } /* q */
         p += 1;
         } /* p */
   
         if (overlap_all==1) {
           Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask; 
           Nsurvive += 1;
         }
       } /* Xmap == 1 */
     } /* z */
   } /* y */
 } /* x */

 printf("#erosion Nsurvive:%d Ntry:%d\n",Nsurvive,Ntry);
 return(Nsurvive);

} /* end of Erosion() */



int Bounded_Erosion(Xmap,Pmap,tar_bitmask,res_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 struct CHAR3DMAP *Pmap; /* probe 3D map */
 unsigned char   tar_bitmask; /* target bitmask (1 or 2 or 4 or 8 or 16 or 32 or 64 or 128)*/
 unsigned char   res_bitmask; /* result bitmask (1 or 2 or 4 or 8 or 16 or 32 or 64 or 128)*/
{
 int x,y,z,p,q,r;
 int xx,yy,zz;
 int Nhalf[3];
 char overlap_all;
 int  Nsurvive,Ntry;
 printf("#Bounded_Erosion(Xmap,Pmap tar_bitmask %d res_bitmask %d)\n",tar_bitmask,res_bitmask);

 Nhalf[0] = (Pmap->N[0]-1)/2;
 Nhalf[1] = (Pmap->N[1]-1)/2;
 Nhalf[2] = (Pmap->N[2]-1)/2;
 Nsurvive = Ntry = 0;

 for (x=0;x<Xmap->N[0];++x){
   if ((x%10)==0){printf("#x %d/%d\n",x,Xmap->N[0]);}
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if ((Xmap->map[x][y][z]&tar_bitmask) == tar_bitmask){
         ++Ntry;
         overlap_all = 1;
         p = 0; 
         while ((p<Pmap->N[0])&&(overlap_all==1)){
           q = 0;
           while ((q<Pmap->N[1])&&(overlap_all==1)){
             r = 0;
             while ((r<Pmap->N[2])&&(overlap_all==1)){
               if (Pmap->map[p][q][r]==1){
                 xx = x + p - Nhalf[0];
                 yy = y + q - Nhalf[1];
                 zz = z + r - Nhalf[2];
                 if ((xx<0)||(xx>=Xmap->N[0])||(yy<0)||(yy>=Xmap->N[1])||(zz<0)||(zz>=Xmap->N[2])){ 
                   /** If (xx,yy,zz) is outside of the view, DO NOTHING !!! **/
            /* 
                   overlap_all = 0;
             printf("#Erosion:out of box !! xyz %d %d %d Xmap %d xx %d yy %d zz %d\n",
               x,y,z,Xmap->map[x][y][z],xx,yy,zz); 
             */
                 }
                 else{ 
                   if ((Xmap->map[xx][yy][zz]&tar_bitmask) != tar_bitmask){
                    overlap_all = 0;
                   }
                 }
               } /* Pmap == 1 */
               r += 1;
             } /* r */
           q += 1;
           } /* q */
         p += 1;
         } /* p */
   
         if (overlap_all==1) {
           Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask; 
           Nsurvive += 1;
         }
       } /* Xmap == 1 */
     } /* z */
   } /* y */
 } /* x */

 printf("#bounded_erosion Nsurvive:%d Ntry:%d\n",Nsurvive,Ntry);
 return(Nsurvive);

} /* end of Bounded_Erosion() */







int Opening(Xmap,Pmap,tar_bitmask,res_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 struct CHAR3DMAP *Pmap; /* probe 3D map */
 unsigned char   tar_bitmask; /* target bitmask (1 or 2 or 4 or 8 or 16 or 32 or 64 or 128)*/
 unsigned char   res_bitmask; /* result bitmask (1 or 2 or 4 or 8 or 16 or 32 or 64 or 128)*/
{
 int x,y,z,p,q,r,n;
 int xx,yy,zz;
 int Nhalf[3],NpntPmap;
 char overlap_all;
 int  Nsurvive_center;
 int  Nneighbor,*XX,*YY,*ZZ; /* array of xx, yy, and zz */
 printf("#Opening(Xmap,Pmap tar_bitmask %d res_bitmask %d)\n",tar_bitmask,res_bitmask);
 
 NpntPmap = Pmap->N[0] * Pmap->N[1] * Pmap->N[2]; 
 XX = (int *)malloc(sizeof(int)*NpntPmap);
 YY = (int *)malloc(sizeof(int)*NpntPmap);
 ZZ = (int *)malloc(sizeof(int)*NpntPmap);
 
 Nhalf[0] = (Pmap->N[0]-1)/2;
 Nhalf[1] = (Pmap->N[1]-1)/2;
 Nhalf[2] = (Pmap->N[2]-1)/2;
 Nsurvive_center = 0;

 for (x=0;x<Xmap->N[0];++x){
  if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]);
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    Nneighbor = 0;
    if ((Xmap->map[x][y][z]&tar_bitmask) == tar_bitmask){
     overlap_all = 1; 

     for (p=0;p<Pmap->N[0];++p){ 
       xx = x + p - Nhalf[0];
      for (q=0;q<Pmap->N[1];++q){ 
       yy = y + q - Nhalf[1];
       for (r=0;r<Pmap->N[2];++r){ 
       zz = z + r - Nhalf[2];
        if (Pmap->map[p][q][r]!=0){
          if ((xx<0)||(xx>=Xmap->N[0])||(yy<0)||(yy>=Xmap->N[1])||(zz<0)||(zz>=Xmap->N[2])) 
            { overlap_all = 0; p = q = r = NpntPmap;}
          else {
            if ((Xmap->map[xx][yy][zz]&tar_bitmask) == tar_bitmask){
             XX[Nneighbor] = xx; 
             YY[Nneighbor] = yy; 
             ZZ[Nneighbor] = zz; 
             ++Nneighbor;
            }
            else { overlap_all = 0; p = q = r = NpntPmap;}
          }
        } 
 
        } /* r */
       } /* q */
      } /* p */
  
     if (overlap_all==1){
      ++Nsurvive_center;
      for (n=0;n<Nneighbor;++n){ 
       xx = XX[n]; yy = YY[n]; zz = ZZ[n];
       Xmap->map[xx][yy][zz] = Xmap->map[xx][yy][zz]|res_bitmask; }
      }
     } /* Xmap == 1 */
    } /* z */
  } /* y */
 } /* x */

 free(XX);
 free(YY);
 free(ZZ);
 printf("#Opening Nsurvive_center %d\n",Nsurvive_center);
 return(Nsurvive_center);

} /* end of Opening() */






void  Reflection(R,O)
 struct CHAR3DMAP *R; /* Reflected Map */
 struct CHAR3DMAP *O; /* Original Map */
{
 int i,x,y,z,xo,yo,zo,xr,yr,zr;
 int Nhalf[3],Ngrid; 

 printf("#Reflection(R,O(%d %d %d))\n",O->N[0],O->N[1],O->N[2]);
 R->grid_width = O->grid_width;
 for (i=0;i<3;++i){ 
  R->OrigPos[i] = O->OrigPos[i];
  R->N[i] = O->N[i];
  Nhalf[i] = (O->N[i]-1)/2;
 
 }

 for (x=0;x<R->N[0];++x){
   for (y=0;y<R->N[1];++y){
     for (z=0;z<R->N[2];++z){R->map[x][y][z] = 0;}
   }
 }
 /*
 printf("#N %d %d %d Nhalf %d %d %d\n",
  O->N[0], O->N[1], O->N[2], Nhalf[0], Nhalf[1], Nhalf[2]);
 */
 Ngrid = 0; 

 for (x=-Nhalf[0];x<=Nhalf[0];++x){
   xo = Nhalf[0] + x;
   xr = Nhalf[0] - x;

   for (y=-Nhalf[1];y<=Nhalf[1];++y){
     yo = Nhalf[1] + y;
     yr = Nhalf[1] - y;
     for (z=-Nhalf[2];z<=Nhalf[2];++z){
       zo = Nhalf[2] + z;
       zr = Nhalf[2] - z;
       R->map[xr][yr][zr] = O->map[xo][yo][zo];  
   } /* z */
  } /* y */
 } /* x */

} /* end of Reflection() */
