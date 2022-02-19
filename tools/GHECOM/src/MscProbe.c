/*

 <MscProbe.c>
 
 for dealing multi scaling probes  

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
#include "HeapSort.h"


/*** FUNCTIONS (GLOBAL) ***/
int Make_MultiScale_Sphere_Probe_Map();
void Dilation_MultiScale();
void Erosion_MultiScale();
void Assign_Rinaccess_to_Atoms_From_MultiSc_Closing_Map();
void Assign_Rinaccess_to_Atom_Shells_From_MultiSc_Closing_Map();

/*** FUNCTIONS (LOCAL) ***/
static void Find_Min_Max_Of_MultiScale_Sphere_Probe();




int Make_MultiScale_Sphere_Probe_Map(M,Nradius,RadiusArray,grid_width)
 struct CHAR3DMAP *M;
 int    Nradius;      /* Number of radius          */
 float  *RadiusArray; /* RadiusArray[0..Nradius-1] (sorted by increasing order) */
 float  grid_width; 
{
 /*
   >> CAUTION << 
   CHAR3DMAP M Voxel within Radius[r] has the value (r+2). 
   This is simply due to technical reason.
   In this program, voxels with 0 corresponds to outer space, 
                    voxels with 1 corresponds to VdW volume of the molecule.
 */
 
 int Ngrid,Ngrid_half,x,y,z,number_of_grid,r;
 char hit; 
 float Pos[3],DD,RRmax,Rmax;
 float *RRarray;
 int *DotCount;

 /** [1] Malloc CHAR3DMAP **/
 Rmax = RadiusArray[Nradius-1];
 Ngrid_half = (int)(Rmax/grid_width);
 Ngrid = 2 * Ngrid_half + 1;
 M->OrigPos[0] = 0.0; M->OrigPos[1] = 0.0; M->OrigPos[2] = 0.0;

 M->N[0] = M->N[1] = M->N[2] = Ngrid;
 printf("#Make_MultiScale_Sphere_Probe_Map(Nradius %d Rmax %f grid_width %f) -> %d %d %d\n",
   Nradius, Rmax,grid_width,Ngrid,Ngrid,Ngrid);

 Malloc_CHAR3DMAP(M,M->N[0],M->N[1],M->N[2]);
                                                                                                   
 /** [2] Make RRarray **/
 RRarray  = (float *)malloc(sizeof(float)*(Nradius+2)); 
 DotCount  = (int *)malloc(sizeof(int)*(Nradius+2)); 
 for (r=0;r<Nradius;++r) 
  { RRarray[r] = RadiusArray[r]*RadiusArray[r];
    DotCount[r] = 0; }
 RRmax = Rmax*Rmax;

 /** [3] Assign radius numbet to CHAR3DMAP **/

 number_of_grid = 0;
 
 for (x=-Ngrid_half;x<=Ngrid_half;++x)
 {
  Pos[0] = M->OrigPos[0] + grid_width * x;
  for (y=-Ngrid_half;y<=Ngrid_half;++y)
  {
   Pos[1] = M->OrigPos[1] + grid_width * y;
   for (z=-Ngrid_half;z<=Ngrid_half;++z)
   {
    Pos[2] = M->OrigPos[2] + grid_width * z;
    DD = Pos[0]*Pos[0] + Pos[1]*Pos[1] + Pos[2]*Pos[2];
    /* printf("x %d y %d z %d\n",x+Ngrid_half,y+Ngrid_half,z+Ngrid_half); */
    hit = 0; r = 0;
    while ((hit==0)&&(r<Nradius)) 
    { 
     if (DD<=RRarray[r])
      { M->map[x+Ngrid_half][y+Ngrid_half][z+Ngrid_half] = r+2; 
        hit = 1; 
        ++number_of_grid; }
      else ++r; 
     }
   } /* z */
  } /* y */
 } /* x */

  for (x=0;x<M->N[0];++x)
   for (y=0;y<M->N[1];++y)
    for (z=0;z<M->N[2];++z) DotCount[M->map[x][y][z]] += 1;

 for (r=0;r<Nradius;++r) 
  printf("#radius(%3d/%3d) %7.3lf A Voxel_char_value %3d Nvoxel %d\n",
   r,Nradius,RadiusArray[r],r+2,DotCount[r+2]);
 
 free(RRarray);
 free(DotCount);
 return(number_of_grid);

} /* end of Make_MultiScale_Sphere_Probe_Map() */



void Find_Min_Max_Of_MultiScale_Sphere_Probe(M,Nradius,Min, Max)
 struct CHAR3DMAP *M;
 int    Nradius;  /* Number of radius  */
 int    *Min;     /* Min[2..Nradius+1] */
 int    *Max;     /* Max[2..Nradius+1] */
{
 int x,y,z,r;

 for (r=0;r<(Nradius+2);++r) {Min[r] = M->N[0]+1; Max[r] = 0;}
 
 for (x=0;x<M->N[0];++x)
 {
  for (y=0;y<M->N[0];++y)
  {
    for (z=0;z<M->N[0];++z)
    {
     r = M->map[x][y][z];  /* r is the radius index of Voxels */
     if (x<Min[r]) Min[r] = x;
     if (x>Max[r]) Max[r] = x;
    } /* x */
   } /* y */
 } /* z */ 

} /* end of Find_Min_Max_Of_MultiScale_Sphere_Probe() */






void Dilation_MultiScale(Xmap,Ymap,Pmap, Nradius)
 struct CHAR3DMAP *Xmap; 
   /* >> Xmap << 
    Initial state: target shape :3D map (voxels of the target shape have non-zero values) 
    Final state  : does not change from the initial state */
 struct CHAR3DMAP *Ymap; 
   /* >> Ymap << 
    Initial state: don't care. Ymap is initialized to zero in this function.
    Final state  : dilated 3D map. Voxel has minimum radius number.
                   0 -> outer space 
                   1 -> VdW volume 
                   2 -> sphere with 0-th radius.
                   3 -> sphere with 1-th radius.
                   :  
           Nradius+2 -> sphere with Nradius-th radius.
   */
  struct CHAR3DMAP *Pmap; /* probe 3D map with multi scale  */
  int    Nradius;         /* Number of Radius */
{
 int x,y,z,p,q,r;
 int xx,yy,zz;
 int Nhalf;
 int *Ncount; 
 
 printf("#Dilation_MultiScale(Xmap,Ymap,Pmap,Nradius %d)\n",Nradius);

 /** Initialize Ymap to zero **/ 
 for (x=0;x<Ymap->N[0];++x)
  for (y=0;y<Ymap->N[1];++y)
   for (z=0;z<Ymap->N[2];++z) 
   { if (Xmap->map[x][y][z] == 1) Ymap->map[x][y][z] = 1;
                             else Ymap->map[x][y][z] = 0; }

 Ncount = (int*)malloc(sizeof(int)*(Nradius+2));
 for (r=0;r<(Nradius+2);++r) Ncount[r] = 0;
 
 /** Scan Xmap and try translation by Pmap  **/ 
 Nhalf = (Pmap->N[0]-1)/2;
 
 for (x=0;x<Xmap->N[0];++x){
   if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]);
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){

       if (Xmap->map[x][y][z]>0){

         for (p=0;p<Pmap->N[0];++p){
         xx = x + p - Nhalf;
         if ((xx>=0)&&(xx<Xmap->N[0])){
           for (q=0;q<Pmap->N[1];++q){
           yy = y + q - Nhalf;
           if ((yy>=0)&&(yy<Xmap->N[1])){
             for (r=0;r<Pmap->N[2];++r){
               if (Pmap->map[p][q][r]>0){
                  zz = z + r - Nhalf;
                  /* printf("x %d p %d xx %d\n",x,p,xx); */
                  if ((zz>=0)&&(zz<Xmap->N[2])){
                    if ( (Ymap->map[xx][yy][zz]==0) ||
                         (Pmap->map[p][q][r] < Ymap->map[xx][yy][zz]) )
                        Ymap->map[xx][yy][zz] = Pmap->map[p][q][r];
                  }
                } /* Pmap > 0 */
               } /* r */
             } 
            } /* q */
          }  
        } /* p */
       } /* Xmap > 0 */
     } /* z */
   } /* y */
 } /* x */


 for (x=0;x<Ymap->N[0];++x)
  for (y=0;y<Ymap->N[1];++y)
   for (z=0;z<Ymap->N[2];++z) Ncount[Ymap->map[x][y][z]] += 1;

 for (r=0;r<(Nradius+2);++r)
  printf("Ymap val %d count %6d\n",r,Ncount[r]);
 free(Ncount);

} /* end of Dilation_MultiScale() */









void Erosion_MultiScale(Xmap,Ymap,Pmap, Nradius)
 struct CHAR3DMAP *Xmap; 
   /* >> Xmap << 
    Initial state  : dilated 3D map. Voxel has minimum radius number.
                   0 -> outer space 
                   1 -> VdW volume 
                   2 -> sphere with 0-th radius.
                   3 -> sphere with 1-th radius.
                   : :
           Nradius+2 -> sphere with Nradius-th radius.
    Final state    : does not change.
   */
 struct CHAR3DMAP *Ymap; 
  /* >> Ymap <<
    Initial state: don't care. Ymap is initialized to zero in this function.
    Final state  : closing 3D map. Voxel has minimum radius number.
                   0 -> outer space 
                   1 -> VdW volume 
                   2 -> sphere with 0-th radius.
                   3 -> sphere with 1-th radius.
                   : :
           Nradius+2 -> sphere with Nradius-th radius.
   */
 struct CHAR3DMAP *Pmap; /* probe 3D map with multi scale  */
 int    Nradius;         /* Number of Radius */
{
 int x,y,z,p,q,r,R;
 int x_Nhalf, y_Nhalf, z_Nhalf;
 int xx,yy,zz,Nround;
 int Nhalf;
 char end_Rloop,out_R,out_all; 
 unsigned char Rout;
 int  *Ncount,*MinX,*MaxX; 

 printf("#Erosion_MultiScale(Xmap,Ymap,Pmap,Nradius %d)\n",Nradius);
 
 /** Initialize Ymap to zero **/ 
 for (x=0;x<Ymap->N[0];++x)
  for (y=0;y<Ymap->N[1];++y)
   for (z=0;z<Ymap->N[2];++z) Ymap->map[x][y][z] = 0;
 
 /** [1] Malloc Temporary Variables  **/ 
 Ncount = (int*)malloc(sizeof(int)*(Nradius+2));
 MinX = (int*)malloc(sizeof(int)*(Nradius+2));
 MaxX = (int*)malloc(sizeof(int)*(Nradius+2));
 for (R=0;R<(Nradius+2);++R) Ncount[R] = 0;



 /** [2] Find Min-Max of Multiscale Sphere Probes  **/ 
 Find_Min_Max_Of_MultiScale_Sphere_Probe(Pmap,Nradius,MinX, MaxX);
 for (r=1;r<=Nradius;++r) { printf("r %d MinX %d MaxX %d\n",r,MinX[r],MaxX[r]); }

 /** [3] Scan Xmap and try translation by Pmap  **/ 
 Nhalf = (Pmap->N[0]-1)/2;
 printf("#Nhalf of Pmap : %d\n",Nhalf); 
 Rout  = 0; 
 for (x=0;x<Xmap->N[0];++x){
  if ((x%10)==0) printf("##x %d/%d\n",x,Xmap->N[0]);
  x_Nhalf = x - Nhalf;

  for (y=0;y<Xmap->N[1];++y){
   y_Nhalf = y - Nhalf;
   for (z=0;z<Xmap->N[2];++z){
         if (Xmap->map[x][y][z]==1) Ymap->map[x][y][z] = 1;
    else if (Xmap->map[x][y][z]>=2){
     z_Nhalf = z - Nhalf;
    
     R = Xmap->map[x][y][z];

     end_Rloop = 0; Nround = 0;
     
     while ((R<(Nradius+2))&&(end_Rloop==0)){
      ++Nround;
 /*
      printf(" Nround %d R %d\n",Nround,R); 
*/
      out_R  = out_all = 0; 

      for (p=MinX[R];p<=MaxX[R];++p){
       xx = x_Nhalf + p;
       if ((xx>=0)&&(xx<Xmap->N[0])){
       for (q=MinX[R];q<=MaxX[R];++q){
        yy = y_Nhalf + q;
        if ((yy>=0)&&(yy<Xmap->N[1])){
        for (r=MinX[R];r<=MaxX[R];++r){
         if ((Pmap->map[p][q][r]>0)&&(Pmap->map[p][q][r]<=R)){
           zz = z_Nhalf + r;
           if ((zz>=0)&&(zz<Xmap->N[2])){
             if (Xmap->map[xx][yy][zz]==0)
               { out_all = 1; 
                 p = Pmap->N[0]; q = Pmap->N[1]; r = Pmap->N[2];}
             if (Xmap->map[xx][yy][zz]>R) 
               { out_R = 1; 
                 Rout = Xmap->map[xx][yy][zz];
                 p = Pmap->N[0]; q = Pmap->N[1]; r = Pmap->N[2];}
            }

           } /* Pmap > 0 */
          } /* r */
         } 
         } /* q */
        }   
      } /* p */

           if ((out_R==0)&&(out_all==0)) {Ymap->map[x][y][z] = R; end_Rloop = 1;}
      else if (out_all == 1) {end_Rloop = 1;}
      else if (out_R   == 1) { R = Rout;}
 
     }  /* while */

    
    } /* Xmap > 0 */

     Ncount[Ymap->map[x][y][z]] += 1;

    } /* z */
  } /* y */
 } /* x */

 for (R=0;R<(Nradius+2);++R) printf("Ymap val %d count %6d\n",R,Ncount[R]);
 
 free(Ncount); free(MinX); free(MaxX);

} /* end of Erosion_MultiScale() */





void Assign_Rinaccess_to_Atoms_From_MultiSc_Closing_Map(AtomHead,M,Nradius,RepType)
 struct ATOM *AtomHead;
 /* Rinaccess will be written in tFactor of each atom. */
 struct CHAR3DMAP *M; 
   /* >> M << 
    Initial state  : dilated 3D map. Voxel has minimum radius number.
                   0 -> outer space 
                   1 -> VdW volume 
                   2 -> sphere with 0-th radius.
                   3 -> sphere with 1-th radius.
                   : :
           Nradius+2 -> sphere with Nradius-th radius.
    Final state    : does not change.
   */
 int Nradius;
 char RepType;  /* how to choose representative value: 'M'edian, m'I'nimum, ma'X'imum */
{
 struct ATOM *an;
 float D[3],DD[3],ddis,RR;
 int   minN[3],maxN[3];
 int   i,x,y,z,Ninside,Nmed;
 float value[512];

 printf("#Assign_Rinaccess_to_Atoms_From_MultiSc_Closing_Map(RepType %c)\n",RepType);
 an = AtomHead;
 while ( an->next != NULL){
  an = an->next;
  an->tFactor = 0.0;
  Ninside = 0; 
  
  for (i=0;i<3;++i){  
      minN[i] = (int)floor((an->Pos[i]-an->R- M->OrigPos[i])/M->grid_width) - 1;
      maxN[i] = (int)ceil((an->Pos[i]+an->R-  M->OrigPos[i])/M->grid_width) + 1; 
      if (minN[i]<0) minN[i] = 0;
      if (maxN[i]>=M->N[i]) maxN[i] = M->N[i]-1;
   }
 
  RR  = an->R*an->R; 
  for (x=minN[0];x<=maxN[0];++x){ 
   D[0] = M->OrigPos[0] + M->grid_width * x - an->Pos[0];
   DD[0] = D[0]*D[0];
   for (y=minN[1];y<=maxN[1];++y){
     D[1] = M->OrigPos[1] + M->grid_width * y - an->Pos[1];
     DD[1] = D[1]*D[1];
     for (z=minN[2];z<=maxN[2];++z){
       D[2] = M->OrigPos[2] + M->grid_width * z - an->Pos[2];
       ddis = DD[0] + DD[1] + D[2]*D[2];

       if (ddis<=RR) 
       {
        /* Inside */
        /* printf("Find skin\n");  */
        if (M->map[x][y][z]==0)  value[Ninside] = Nradius + 3;
        else value[Ninside] = M->map[x][y][z];
        ++Ninside;   
     }

     } /* z */
    } /* y */
  } /* x */

  Nmed = 0; 
  if (Ninside>1){ 
    HeapSort_Integer_NoIndex(Ninside,value);
    if ((Ninside%2)==1) Nmed = (Ninside-1)/2;
    if ((Ninside%2)==0) Nmed = Ninside/2-1;
    if (RepType == 'M') an->tFactor = value[Nmed]; 
    if (RepType == 'X') an->tFactor = value[Ninside-1]; 
    if (RepType == 'I') an->tFactor = value[0]; 
    printf("#Ninside %d Median %f Max %f Min %f\n",Ninside,value[Nmed],value[Ninside-1],value[0]);
  }
  else if (Ninside==1) an->tFactor = value[0];
  else if (Ninside==0) an->tFactor = 0.0;


} /* an */

} /* end of Assign_Rinaccess_to_Atoms_From_MultiSc_Closing_Map() */





void Assign_Rinaccess_to_Atom_Shells_From_MultiSc_Closing_Map(AtomHead,M,Nradius,Rshell,EmpType,RepType)
 struct ATOM *AtomHead;
 /* Rinaccess will be written in tFactor of each atom. */
 struct CHAR3DMAP *M; 
   /* >> M << 
    Initial state  : dilated 3D map. Voxel has minimum radius number.
                   0 -> outer space 
                   1 -> VdW volume 
                   2 -> sphere with 0-th radius.
                   3 -> sphere with 1-th radius.
                   : :
           Nradius+2 -> sphere with Nradius-th radius.
    Final state    : does not change.
   */
 int Nradius;
 float Rshell;
 char EmpType;  /* how to deal with Empty(0) value 'Z'ero, or max'X'imum radius */
 char RepType;  /* how to choose representative value: 'M'edian, m'I'nimum, ma'X'imum */
{
 struct ATOM *an;
 float D[3],DD[3],ddis,RR,RR1;
 int   minN[3],maxN[3];
 int   i,x,y,z,Nskin,Ninside,Nmed;
 float value[512];

 printf("#Assign_Rinaccess_to_Atom_Shells_From_MultiSc_Closing_Map(Rshell %f RepType %c)\n",Rshell,RepType);
 an = AtomHead;
 while ( an->next != NULL){
  an = an->next;
  an->tFactor = 0.0;
  Nskin = Ninside = 0; 
  
  for (i=0;i<3;++i){  
      minN[i] = (int)floor((an->Pos[i] - an->R - Rshell - M->OrigPos[i])/M->grid_width);
      maxN[i] = (int)ceil( (an->Pos[i] + an->R + Rshell - M->OrigPos[i])/M->grid_width); 
      if (minN[i]<0) minN[i] = 0;
      if (maxN[i]>=M->N[i]) maxN[i] = M->N[i]-1;
   }
 
  RR  = an->R*an->R; 
  RR1 = (an->R + M->grid_width)*(an->R + M->grid_width); 
  for (x=minN[0];x<=maxN[0];++x){ 
   D[0] = M->OrigPos[0] + M->grid_width * x - an->Pos[0];
   DD[0] = D[0]*D[0];
   for (y=minN[1];y<=maxN[1];++y){
     D[1] = M->OrigPos[1] + M->grid_width * y - an->Pos[1];
     DD[1] = D[1]*D[1];
     for (z=minN[2];z<=maxN[2];++z){
       D[2] = M->OrigPos[2] + M->grid_width * z - an->Pos[2];
       ddis = DD[0] + DD[1] + D[2]*D[2];

       if (ddis<=RR) ++Ninside;
       if ((ddis>RR)&&(ddis<=RR1)){
        /* Skin */
        /* printf("Find skin\n");  */
         if (M->map[x][y][z]==0){  
            if (EmpType=='X') value[Nskin] = Nradius + 3;
            if (EmpType=='Z') value[Nskin] = 0;
         }
         else  value[Nskin] = M->map[x][y][z];
         ++Nskin;
         if (Nskin>=512) {printf("#ERROR:Nskin > 512\n"); exit(1);}
       }


     } /* z */
    } /* y */
  } /* x */

  Nmed = 0;
  if (Nskin>1){ 
    HeapSort_Integer_NoIndex(Nskin,value);
    if ((Nskin%2)==1) Nmed = (Nskin-1)/2;
    if ((Nskin%2)==0) Nmed = Nskin/2-1;
  /* 
  printf(">an %s %s %s Nskin %3d min %.0f med %.0f max %.0f\n",
  an->Atom,an->Resi,an->Rnum,Nskin,value[0],value[Nmed],value[Nskin-1]); 
  */
    if (RepType == 'M') an->tFactor = value[Nmed]; 
    if (RepType == 'X') an->tFactor = value[Nskin-1]; 
    if (RepType == 'I') an->tFactor = value[0]; 
 }
 else if (Nskin==1) an->tFactor = value[0];
 else if (Nskin==0) an->tFactor = 0.0;
 
 } /* an */

} /* end of Assign_Rinaccess_to_Atom_Shells_From_MultiSc_Closing_Map() */
