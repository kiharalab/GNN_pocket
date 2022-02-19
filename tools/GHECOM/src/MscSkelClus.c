/*

 <MscSkelClus.c>
 
 for Clustering Multiscale Pocket by Skeletonization 

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
#include "MscShPrb.h"
#include "MscPockClus.h"
#include "GclsMergeSort.h"
#include "OpeStd3DPr.h"



#define CLUSTER_MASK 1
#define FIGURE_MASK  2
#define DILATE_MASK  4

/*** FUNCTIONS (GLOBAL) ***/
int Erosion_for_Non_Zero_Non_255_Voxel();
void Dilation_for_Eroded_Grid_Cluster();
void Extended_Skeleton_To_Shallower_Direction();
void Extended_Skeleton_To_Shallower_Dir_With_Volume_Restrict();

/*** FUNCTIONS (LOCAL) ***/
static void Label_Connected_Shallower_Neighbors();
static int Ngrid_for_Dilated_Skeleton_Cluster();

/*** VARIABLES (LOCAL) ***/
static  int  NeighX[26],NeighY[26],NeighZ[26];  /* for neighbor voxels */
                                                                                                                     
                                                                                                                
int Erosion_for_Non_Zero_Non_255_Voxel(Xmap,Ymap,Pmap)
 struct CHAR3DMAP *Xmap; 
 /* target 3D map with multi value, for example, 
    Initial state    : pocket 3D map. Voxel has minimum radius number.
                   0 -> outer space
                   1 -> sphere with 1-th radius.
                   2 -> sphere with 2-th radius.
                   : :
             Nradius -> sphere with Nradius-th radius.
                 255 -> VdW volume.
   
   Voxels having from 1 to 254 values are regarded as the "figure". 
 */
 struct CHAR3DMAP *Ymap; /* result eroded 3D map (to be caluclated) */ 
 struct CHAR3DMAP *Pmap; /* probe 3D map */
{
 int x,y,z,p,q,r;
 int xx,yy,zz;
 int Nhalf;
 char overlap_all;
 int  Nsurvive,Ntry;
 printf("#Erosion_for_Non_Zero_Voxel(Xmap,Pmap)\n");

 /** Initialize Ymap **/
 for (x=0;x<Xmap->N[0];++x) {
  for (y=0;y<Xmap->N[1];++y) {
   for (z=0;z<Xmap->N[2];++z) { 
   Ymap->map[x][y][z] = 0;
   }
  }
 } 


 /** Erosion **/
 Nhalf = (Pmap->N[0]-1)/2;
 Nsurvive = Ntry = 0;
 for (x=0;x<Xmap->N[0];++x){
  if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]);
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    if ((Xmap->map[x][y][z]>0)&&(Xmap->map[x][y][z]<255)){
     ++Ntry;
     overlap_all = 1;
     p = 0;
     while ((p<Pmap->N[0])&&(overlap_all==1)){
      q = 0;
      while ((q<Pmap->N[1])&&(overlap_all==1)){
       r = 0;
       while ((r<Pmap->N[2])&&(overlap_all==1)){
        if (Pmap->map[p][q][r]==1){
          xx = x + p - Nhalf;
          yy = y + q - Nhalf;
          zz = z + r - Nhalf;
          if ((xx<0)||(xx>=Xmap->N[0])||(yy<0)||(yy>=Xmap->N[1])||(zz<0)||(zz>=Xmap->N[2]))
           { overlap_all = 0;
             printf("#Erosion:out of box !! xyz %d %d %d Xmap %d xx %d yy %d zz %d\n",
               x,y,z,Xmap->map[x][y][z],xx,yy,zz); }
          else
          {
             if ((Xmap->map[xx][yy][zz]==0)||(Xmap->map[xx][yy][zz]==255)) overlap_all = 0;
          }
                                                                                                                
          } /* Pmap == 1 */
         ++r;
        } /* r */
        ++q;
       } /* q */
       ++p;
      } /* p */

     if (overlap_all==1)
         { Ymap->map[x][y][z] = Xmap->map[x][y][z]; ++Nsurvive;}
  /*
     if (overlap_all==1)
         { Ymap->map[x][y][z] = 1; ++Nsurvive;}
  */
     } /* Xmap == 1 */
    } /* z */
  } /* y */
 } /* x */
                                                                                                                
 printf("#erosion Nsurvive:%d Ntry:%d\n",Nsurvive,Ntry);
 return(Nsurvive);
                                                                                                                
} /* end of Erosion_for_Non_Zero_Non_255_Voxels() */





void Dilation_for_Eroded_Grid_Cluster(Xmap,GclsHead,Pmap)
 struct CHAR3DMAP *Xmap; 
 /* target 3D map with multi value, for example, 
    Initial state    : pocket 3D map. Voxel has minimum radius number.
                   0 -> outer space
                   1 -> sphere with 1-th radius.
                   2 -> sphere with 2-th radius.
                   : :
             Nradius -> sphere with Nradius-th radius.
                 255 -> VdW volume.
   
   Voxels having from 1 to 254 values are regarded as the "figure". 
 */
 struct GRID_CLUSTER *GclsHead;
 struct CHAR3DMAP *Pmap; /* probe 3D map */
{
 int x,y,z,p,q,r,xx,yy,zz,i,init;
 int Nhalf,Nadd,Nero,Nskel,Ndila;
 struct GRID_CLUSTER *gcls;
 struct GRID_MEMBER *gmem;
 struct CHAR3DMAP Lmap; /* for temporarily 3D map for labeling */ 
 int min[3],max[3];
 int num_index[1000];

 printf("#Dilation_for_Eroded_Grid_Cluster()\n");
 Malloc_CHAR3DMAP(&Lmap,Xmap->N[0],Xmap->N[1],Xmap->N[2]);

 init = 1;
 gcls = GclsHead;
 Nhalf = (Pmap->N[0]-1)/2;

 for (i=0;i<3;++i){
   min[i] = 0;
   max[i] = Lmap.N[i];
 } 

 while (gcls->next != NULL){
  gcls = gcls->next;
  Nadd = Nero = 0; 
  if (init==1) {Initialize_CHAR3DMAP(&Lmap); init=0;}
  else{
   for (x=min[0];x<=max[0];++x)
    for (y=min[1];y<=max[1];++y)
     for (z=min[2];z<=max[2];++z) Lmap.map[x][y][z] = 0;
  }
  for (i=0;i<3;++i) {min[i] = Xmap->N[i]; max[i] = 0;}

  Nskel = Number_Of_Grid_Member(&(gcls->HeadGridMember));
  gmem = &(gcls->HeadGridMember);
  while (gmem->next != NULL){
   gmem = gmem->next;
   x = gmem->x; y = gmem->y; z = gmem->z;
   ++Nero;
   Lmap.map[x][y][z] = Lmap.map[x][y][z]|FIGURE_MASK;
   
    for (p=0;p<Pmap->N[0];++p){
      for (q=0;q<Pmap->N[1];++q){
       for (r=0;r<Pmap->N[2];++r){
        if (Pmap->map[p][q][r]==1){
          xx = x + p - Nhalf;
          yy = y + q - Nhalf;
          zz = z + r - Nhalf;
          /* printf("x %d p %d xx %d\n",x,p,xx); */
          if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){
           Lmap.map[xx][yy][zz] = Lmap.map[xx][yy][zz]|DILATE_MASK;
           if (xx<min[0]){min[0] = xx;} if (yy<min[1]){min[1] = yy;} if (zz<min[2]){min[2] = zz;}
           if (xx>max[0]){max[0] = xx;} if (yy>max[1]){max[1] = yy;} if (zz>max[2]){max[2] = zz;}
           }
         } /* Pmap == 1 */
        } /* r */
       } /* q */
      } /* p */
   } /* gmem */

   for (x=min[0];x<=max[0];++x){
    for (y=min[1];y<=max[1];++y){
     for (z=min[2];z<=max[2];++z){  
     if (( (Lmap.map[x][y][z]&DILATE_MASK)==(DILATE_MASK) ) &&
         ( (Lmap.map[x][y][z]&FIGURE_MASK)!=(FIGURE_MASK) )  ){ 
         Add_Grid_Member(&(gcls->HeadGridMember),x,y,z,'E');
         ++Nadd;
        }
 
       Lmap.map[x][y][z] = 0;
     }
    }
  }

  Ndila = Number_Of_Grid_Member(&(gcls->HeadGridMember));
  printf("Nskel %d Nadd %d Ndila %d\n",Nskel,Nadd,Ndila);

 } /* gcls */

 Sort_GridCluster(GclsHead,num_index,'I');

 Free_CHAR3DMAP(&Lmap);

} /* end of Dilation_for_Eroded_Grid_Cluster() */





void Extended_Skeleton_To_Shallower_Direction(Xmap,GclsHead,thre_val_min,thre_val_max)
 struct CHAR3DMAP *Xmap;
  /*
    Initial state    : pocket 3D map. Voxel has minimum radius number.
                   0 -> outer space
                   1 -> sphere with 1-th radius.
                   2 -> sphere with 2-th radius.
                   : :
             Nradius -> sphere with Nradius-th radius.
                 255 -> VdW volume
      Voxels having from 1 to [thre_voxel_value] values are regarded as the "figure".
      (1 <= [thre_voxel_val] <= 254)
   */
 struct GRID_CLUSTER *GclsHead;
 int    thre_val_min;  /* Threshold Voxel Value (1 <= thre_voxel_val <= 254) */
 int    thre_val_max;  /* Threshold Voxel Value (1 <= thre_voxel_val <= 254) */
{
 struct GRID_CLUSTER *gcls;
 struct GRID_MEMBER  *gmem;
 struct CHAR3DMAP Lmap; /* for temporarily 3D map for labeling */ 
 
 printf("#Extended_Skeleton_To_Shallower_Direction(thre_val_min %d thre_val_max %d)\n"
  ,thre_val_min,thre_val_max);
 Malloc_CHAR3DMAP(&Lmap,Xmap->N[0],Xmap->N[1],Xmap->N[2]);

 Set_NeighX_NeighY_NeighZ(PAR.NeighborNum,NeighX,NeighY,NeighZ);

 gcls = GclsHead;
 while (gcls->next != NULL)
 {
  gcls = gcls->next;
  /* (i) initialize Lmap */
  gmem = &(gcls->HeadGridMember);
  while (gmem->next != NULL)
  { gmem = gmem->next;
    Lmap.map[gmem->x][gmem->y][gmem->z] = 1; }
 
  /* (ii) Labeling Neighbors for Shallower direction  */
   gmem = &(gcls->HeadGridMember);
   while (gmem->next != NULL)
   { gmem = gmem->next;
     Label_Connected_Shallower_Neighbors(
      Xmap,&Lmap,gmem->x,gmem->y,gmem->z,thre_val_max,&(gcls->HeadGridMember));
   }

  gcls->Nmember = Number_Of_Grid_Member(&(gcls->HeadGridMember));

  /* (iii) Recover Lmap */
   gmem = &(gcls->HeadGridMember);
   while (gmem->next != NULL)
   { gmem = gmem->next;
     Lmap.map[gmem->x][gmem->y][gmem->z] = 0; }

 } /* gcls */
 
 Free_CHAR3DMAP(&Lmap);

} /* end of Extended_Skeleton_To_Shallower_Direction() */




void Extended_Skeleton_To_Shallower_Dir_With_Volume_Restrict(
   Xmap,GclsHead,MSprb, Pmap,thre_val_min,VolMin, VolMax)
 struct CHAR3DMAP *Xmap;
  /*
    Initial state    : pocket 3D map. Voxel has minimum radius number.
                   0 -> outer space
                   1 -> sphere with 1-th radius.
                   2 -> sphere with 2-th radius.
                   : :
             Nradius -> sphere with Nradius-th radius.
                 255 -> VdW volume
      Voxels having from 1 to [thre_voxel_value] values are regarded as the "figure".
      (1 <= [thre_voxel_val] <= 254)
   */
 struct GRID_CLUSTER *GclsHead;
 struct MULSC_SHELL_PRB *MSprb;  
 struct CHAR3DMAP *Pmap; /* Map for small probe */
 int    thre_val_min;  /* Threshold Voxel Value (1 <= thre_voxel_val <= 254) */
 float  VolMin, VolMax; /* Minimum and Maximum Volume of pocket */ 
{
 struct GRID_CLUSTER *gcls;
 struct GRID_MEMBER  *gmem;
 struct GRID_MEMBER NewGridMember;
 struct CHAR3DMAP Lmap; /* for temporarily 3D map for labeling */ 
 int i,t,Ngrid_new, Ngrid_dila;
 char  inVrange;
 float Vcube,Vpock;

 printf("#Extended_Skeleton_To_Shallower_Dir_With_Volume_Restrict(thre_val_min %d Vmin %f Vmax %f)\n"
  ,thre_val_min,VolMin, VolMax); fflush(stdout);

 Set_NeighX_NeighY_NeighZ(PAR.NeighborNum,NeighX,NeighY,NeighZ);

 /** Malloc and Initialize Lmap **/
 Malloc_CHAR3DMAP(&Lmap,Xmap->N[0],Xmap->N[1],Xmap->N[2]);
 Lmap.grid_width = Xmap->grid_width;
 for (i=0;i<3;++i) Lmap.OrigPos[i] = Xmap->OrigPos[i]; 
 Vcube = Xmap->grid_width * Xmap->grid_width * Xmap->grid_width;

 /******************************************/
 /*** WHILE LOOP FOR GRID_CLUSTER (gcls) ***/
 /******************************************/
 
 gcls = GclsHead;
 while (gcls->next != NULL){
  gcls = gcls->next;
  gcls->flag = 0; 
  printf("#>gcls %d\n",gcls->num); fflush(stdout);

  /* (1) initialize Lmap */
  NewGridMember.next = NULL;
  Ngrid_dila = Ngrid_for_Dilated_Skeleton_Cluster(gcls,&NewGridMember,Pmap,&Lmap);
  t = thre_val_min;  
  gmem = &(gcls->HeadGridMember);
  while (gmem->next != NULL)
  { gmem = gmem->next;
    Lmap.map[gmem->x][gmem->y][gmem->z] = Lmap.map[gmem->x][gmem->y][gmem->z]|CLUSTER_MASK; }
  
  /*
  printf("g%dc t %d %f Ngrid_orig %d Ngrid_dila %d Vdila %f\n",
   gcls->num,t,MSprb->Rarray[t],Number_Of_Grid_Member(&(gcls->HeadGridMember)),Ngrid_dila,
   Ngrid_dila*Xmap->grid_width*Xmap->grid_width*Xmap->grid_width);
  */
 
 
  /* (2) Labeling Neighbors for Shallower direction  */
  /*      loop for radius_threshold_number t */ 
  

  t = MSprb->Nradius;
  inVrange = 0;
  
  while ((t>thre_val_min) && (inVrange == 0)){
    Free_Grid_Member(&NewGridMember);

   /* (i) Make extended NewGridMmber from gcls */
   gmem = &(gcls->HeadGridMember);
   while (gmem->next != NULL)
   { gmem = gmem->next;
     Label_Connected_Shallower_Neighbors(Xmap,&Lmap,gmem->x,gmem->y,gmem->z,t,&NewGridMember); 
   }

   /* (ii) Calculated Volume */
   Ngrid_new  = Number_Of_Grid_Member(&NewGridMember);
   Ngrid_dila = Ngrid_for_Dilated_Skeleton_Cluster(gcls,&NewGridMember,Pmap,&Lmap);
   Vpock = Ngrid_dila * Vcube; 
   if ((VolMin<=Vpock)&&(Vpock<=VolMax)) inVrange = 1;

   printf("g%dc t %d %f Ngrid_new %d Ngrid_dila %d Vpock %f inVrange %d\n",
    gcls->num,t,MSprb->Rarray[t],Ngrid_new,Ngrid_dila, Vpock,inVrange);
  
   /* (iii) Recover Lmap for NewGridMember */
    gmem = &NewGridMember; 
    while (gmem->next != NULL)
    { gmem = gmem->next;
      Lmap.map[gmem->x][gmem->y][gmem->z] = 0; }
 
    --t;
 
  } /* t */

  /* (3) If Volume is in the range[VolMin,VolMax], then unify NewGridMember and gcls->HeadGridMember  */
  if (inVrange==1)
   {
    gmem = &NewGridMember; 
    while (gmem->next != NULL)
    { gmem = gmem->next;
      Add_Grid_Member(&(gcls->HeadGridMember),gmem->x,gmem->y,gmem->z,'E'); }
   }
  else { gcls->flag = 1; 
        printf("#gcls %d is to deleted\n",gcls->num); 
  }
 
  Free_Grid_Member(&NewGridMember);
  
  /* (4) Recover Lmap for gcls->HeadGridMember  */
   gmem = &(gcls->HeadGridMember);
   while (gmem->next != NULL)
   { gmem = gmem->next;
     Lmap.map[gmem->x][gmem->y][gmem->z] = 0; }


 } /* gcls */


 /**************************************/
 /*** Delete Clusters  with flag = 1 ***/
 /**************************************/
 gcls = GclsHead;
 while (gcls->next != NULL)
 {
  gcls = gcls->next;
  printf("#gnum %d flag %d\n",gcls->num,gcls->flag);
  if (gcls->flag==1) 
  {
   if (gcls->prev != NULL) gcls->prev->next = gcls->next;  
   if (gcls->next != NULL) gcls->next->prev = gcls->prev;  
   }
 }

 
 Free_CHAR3DMAP(&Lmap);

} /* end of Extended_Skeleton_To_Shallower_Dir_With_Volume_Restrict() */







void Label_Connected_Shallower_Neighbors(Xmap,Lmap,x,y,z,thre_voxel_max,gmem_head)
 struct CHAR3DMAP *Xmap; /* target 3D pocket map */
 struct CHAR3DMAP *Lmap; /* cluster label 3D map */
 int    x,y,z;
 int    thre_voxel_max; /* Threshold Voxel Value (1 <= thre_voxel_val <= 254) */
 struct GRID_MEMBER *gmem_head;
{
 int  n,xx,yy,zz;
 char val0;

 val0 = Xmap->map[x][y][z]; 
 
 if ((Lmap->map[x][y][z]&CLUSTER_MASK)!=CLUSTER_MASK){ 
   Lmap->map[x][y][z] = Lmap->map[x][y][z]|CLUSTER_MASK;
   Add_Grid_Member(gmem_head,x,y,z,'E');
   }
 /*
    printf("#Label_Connected_Shallower_Neighbors(%d %d %d val %d thre_max %d val0 %d)\n",
      x,y,z,Xmap->map[x][y][z],thre_voxel_max,val0);
 */

 for (n=0;n<PAR.NeighborNum;++n){
   xx = x + NeighX[n];
   yy = y + NeighY[n];
   zz = z + NeighZ[n];
 
   if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){
    if ((Xmap->map[xx][yy][zz]>=val0)&&(Xmap->map[xx][yy][zz]<=thre_voxel_max)
        &&(Lmap->map[xx][yy][zz]==0)){
     Label_Connected_Shallower_Neighbors(Xmap,Lmap,xx,yy,zz,thre_voxel_max,gmem_head);
     }
    }
 } /* n */
} /* end of Label_Connected_Shallower_Neighbors() */






int Ngrid_for_Dilated_Skeleton_Cluster(gcls,gmem_add_head,Pmap,Lmap)
 struct GRID_CLUSTER *gcls;
 struct GRID_MEMBER  *gmem_add_head;
 struct CHAR3DMAP *Pmap; /* probe 3D map */
 struct CHAR3DMAP *Lmap; /* for temporarily 3D map for labeling */ 
{
 int p,q,r,x,y,z,i;
 int Nhalf,Ndila,Nsk_orig,Nsk_add;
 struct GRID_MEMBER *gmem;
 char not_FIGURE_MASK,not_DILATE_MASK;
 int min[3],max[3];

 Nhalf = (Pmap->N[0]-1)/2;
 
 not_FIGURE_MASK = ~FIGURE_MASK;
 not_DILATE_MASK = ~DILATE_MASK;

 /* 
printf("#Ngrid_for_Dilated_Skeleton_Cluster(Lmap %d %d %d Pmap %d %d %d) fig %d nfig %d dil %d ndil %d\n",
  Lmap->N[0],Lmap->N[1],Lmap->N[2], Pmap->N[0],Pmap->N[1],Pmap->N[2],
  FIGURE_MASK,not_FIGURE_MASK,DILATE_MASK,not_DILATE_MASK); 
 fflush(stdout);
 */
 
 for (i=0;i<3;++i) {min[i] = Lmap->N[i]; max[i] = 0;}
 Ndila = Nsk_orig = Nsk_add = 0;
 
 /** (1) Dilated for HeadGridMember of gcls **/
 gmem = &(gcls->HeadGridMember);
 
 while (gmem->next != NULL){
   gmem = gmem->next;
   ++Nsk_orig;
   Lmap->map[gmem->x][gmem->y][gmem->z] = Lmap->map[gmem->x][gmem->y][gmem->z]|FIGURE_MASK;
   Lmap->map[gmem->x][gmem->y][gmem->z] = Lmap->map[gmem->x][gmem->y][gmem->z]|DILATE_MASK;

   for (p=0;p<Pmap->N[0];++p){
     for (q=0;q<Pmap->N[1];++q){
      for (r=0;r<Pmap->N[2];++r){
       if (Pmap->map[p][q][r]==1){
          x = gmem->x + p - Nhalf; 
          y = gmem->y + q - Nhalf; 
          z = gmem->z + r - Nhalf;
          if ((x>=0)&&(x<Lmap->N[0])&&(y>=0)&&(y<Lmap->N[1])&&(z>=0)&&(z<Lmap->N[2])){
           Lmap->map[x][y][z] = Lmap->map[x][y][z]|DILATE_MASK;
           if (x<min[0]){min[0] = x;} if (y<min[1]){min[1] = y;} if (z<min[2]){min[2] = z;}
           if (x>max[0]){max[0] = x;} if (y>max[1]){max[1] = y;} if (z>max[2]){max[2] = z;}
           }
         } /* Pmap == 1 */
        } /* r */
       } /* q */
      } /* p */
   } /* gmem */
 

  /** (2) Dilated for HeadGridMember of gcls **/
  gmem = gmem_add_head;
  while (gmem->next != NULL)
  {
   gmem = gmem->next;
   ++Nsk_add;
   Lmap->map[gmem->x][gmem->y][gmem->z] = Lmap->map[gmem->x][gmem->y][gmem->z]|FIGURE_MASK;
   Lmap->map[gmem->x][gmem->y][gmem->z] = Lmap->map[gmem->x][gmem->y][gmem->z]|DILATE_MASK;
 
   for (p=0;p<Pmap->N[0];++p){
      for (q=0;q<Pmap->N[1];++q){
       for (r=0;r<Pmap->N[2];++r){
        if (Pmap->map[p][q][r]==1){
          x = gmem->x + p - Nhalf; 
          y = gmem->y + q - Nhalf; 
          z = gmem->z + r - Nhalf;
          
          if ((x>=0)&&(x<Lmap->N[0])&&(y>=0)&&(y<Lmap->N[1])&&(z>=0)&&(z<Lmap->N[2])){
           Lmap->map[x][y][z] = Lmap->map[x][y][z]|DILATE_MASK;
           if (x<min[0]){min[0] = x;} if (y<min[1]){min[1] = y;} if (z<min[2]){min[2] = z;}
           if (x>max[0]){max[0] = x;} if (y>max[1]){max[1] = y;} if (z>max[2]){max[2] = z;}
           }
         } /* Pmap == 1 */
        } /* r */
       } /* q */
      } /* p */
   } /* gmem */

   /*** (3) Counting Pocket Grids and reset Lmap to zero ***/
   for (x=min[0];x<=max[0];++x){
    for (y=min[1];y<=max[1];++y){
     for (z=min[2];z<=max[2];++z){  
     if ((Lmap->map[x][y][z]&DILATE_MASK)==(DILATE_MASK)) { ++Ndila; 
     }
     Lmap->map[x][y][z] = Lmap->map[x][y][z]&not_FIGURE_MASK;
     Lmap->map[x][y][z] = Lmap->map[x][y][z]&not_DILATE_MASK;
     }
    }
  }

 /*
 printf("#Ndila %d Nsk_orig %d Nsk_add %d\n",Ndila,Nsk_orig,Nsk_add);
 */
 
 return(Ndila);

} /* end of Ngrid_for_Dilated_Skeleton_Cluster() */
