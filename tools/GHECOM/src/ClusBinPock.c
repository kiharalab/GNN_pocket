/*

 <ClusBinPock.c>
  
  Functions for clustering of Binary Image Pocket on CHAR3DMAP 

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
#include "OpeStd3DPr.h"
#include "MapCCP4_CHAR.h"
#include "ReadOpt.h"



/*** FUNCTIONS (GLOBAL) ***/
int  Make_Clusters_for_Binary_CHAR3DMAP_Recursive();
int  Make_Clusters_for_Binary_CHAR3DMAP_Queue();
void Output_Clusters_of_CHAR3DMAP_in_PDB_format();
void Output_Clusters_of_CHAR3DMAP_in_CCP4_format();


/*** FUNCTIONS (LOCAL) ***/
static int  label_neighbor_with_nomarked_target_Recursive();
static int  label_neighbor_with_nomarked_target_Queue();
static void add_new_CHAR3DMAP_cluster();
static void Output_Cluster_CHAR3DMAP_in_PDB_format();
static int Number_of_CHAR3DMAP_Clusters();


/*** FUNCTIONS (LOCAL) FOR QUEUE_XYZ ***/
static int PUT_QUEUE_XYZ();
static int GET_QUEUE_XYZ();
static void INIT_QUEUE_XYZ();
static void EMPTY_QUEUE_XYZ();
static int STORE_NUM_QUEUE_XYZ();
static void MALLOC_QUEUE_XYZ();
static void FREE_QUEUE_XYZ();

/*** GLOBAL VARIABLES FOR QUEUE_XYZ ***/
static int *QUEUE_X, *QUEUE_Y, *QUEUE_Z;     /* [0..Q_MAX] (malloc later) */
static int Q_MAX, Q_HEAD, Q_TAIL;
static int Q_STORE_MAX;  /* maximum number of XYZs stored in QUEUE */
static int Q_COUNT;  /* number of Q-stored XYZ in total */



static int MAX_RECUR_NUM;
static int Nrecur;



/************** FUNCTIONS FOR QUEUE_XYZ (FROM HERE)  ****************/

static int PUT_QUEUE_XYZ(x,y,z)
  int x,y,z;
{
  int q_store;
   if (Q_TAIL>Q_MAX){
    printf("#ERROR:PUT_QUEUE_XYZ:Q_TAIL %d is over Q_MAX %d\n",Q_TAIL,Q_MAX);
    exit(1);
  }
  QUEUE_X[Q_TAIL] = x;
  QUEUE_Y[Q_TAIL] = y;
  QUEUE_Z[Q_TAIL] = z;
  Q_TAIL += 1; 
  if (Q_TAIL > Q_MAX){
    Q_TAIL = 0;
  }
  q_store = STORE_NUM_QUEUE_XYZ();
 /*
 *   printf("QUEUE_PUT %d %d %d HEAD %d TAIL %d N %d MAX %d\n",x,y,z,Q_HEAD,Q_TAIL,q_store,Q_MAX);  
 *     */ 
  if (q_store >Q_STORE_MAX){
    Q_STORE_MAX = q_store;
  }
  if (Q_TAIL == Q_HEAD){
   /* If one more XYZ is added to the full QUEUE (with Q_MAX XYZs),  release all the XYZs. */
    printf("#WARNING:QUEUE_XYZ IS FULL ! (HEAD %d TAIL %d N %d MAX %d) Q_STORE_MAX %d\n",Q_HEAD, Q_TAIL, STORE_NUM_QUEUE_XYZ(),Q_MAX,Q_STORE_MAX);
    exit(1);
    return(0);
  }
  return(1);
} /* end of PUT_QUEUE_XYZ() */

static int GET_QUEUE_XYZ(x, y, z)
  int *x, *y, *z;
{
  if (Q_HEAD != Q_TAIL){
    *x = QUEUE_X[Q_HEAD];
    *y = QUEUE_Y[Q_HEAD];
    *z = QUEUE_Z[Q_HEAD];
    Q_HEAD = Q_HEAD + 1;
    if (Q_HEAD>Q_MAX){
      Q_HEAD = 0;
    }
 /* 
 *    printf("QUEUE_GET %d %d %d  HEAD %d TAIL %d N %d MAX %d\n",*x,*y,*z,Q_HEAD, Q_TAIL,STORE_NUM_QUEUE_XYZ(),Q_MAX);  
 *      */
    return(1);
  }
  else{
    *x =  *y = *z = 0;
    printf("#GET --failed--\n");
    return(0);
  }
} /* end of GET_QUEUE_XYZ() */


static void INIT_QUEUE_XYZ()
{
  Q_HEAD = 0;
  Q_TAIL = 0;
} /* end of INIT_QUEUE_XYZ() */

static void EMPTY_QUEUE_XYZ()
{
  Q_HEAD = Q_TAIL;
}  /* end of EMPTY_QUEUE_XYZ() */


static int STORE_NUM_QUEUE_XYZ()
{
  if (Q_TAIL >= Q_HEAD){
    return(Q_TAIL-Q_HEAD);
  }
  else{
    return(Q_TAIL+Q_MAX - Q_HEAD);
  }
}  /* end of STORE_NUM_QUEUE_XYZ() */


static void MALLOC_QUEUE_XYZ()
{
  int i;
  printf("#MALLOC_QUEUE_XYZ(Q_MAX:%d)\n",Q_MAX);
  QUEUE_X = (int *)malloc(sizeof(int)*(Q_MAX+1));
  QUEUE_Y = (int *)malloc(sizeof(int)*(Q_MAX+1));
  QUEUE_Z = (int *)malloc(sizeof(int)*(Q_MAX+1));
  for (i=0;i<=Q_MAX;++i){
    QUEUE_X[i] = -1;
    QUEUE_Y[i] = -1;
    QUEUE_Z[i] = -1;
  }
} /* end of MALLOC_QUEUE_XYZ() */


static void FREE_QUEUE_XYZ()
{
  free(QUEUE_Z);
  free(QUEUE_Y);
  free(QUEUE_X);
} /* end of FREE_QUEUE_XYZ() */


/************** FUNCTIONS FOR QUEUE_XYZ (TO HERE)  ****************/



















int Make_Clusters_for_Binary_CHAR3DMAP_Recursive(Xmap,tar_bitmask,mark_bitmask,curr_bitmask,ClusMapHead)
 struct CHAR3DMAP *Xmap;         /* target 3D map */
 unsigned char tar_bitmask;      /* target   forground bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char mark_bitmask;     /* marked   forground bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char curr_bitmask;     /* currrent forground bitmask (1 or 2 or 4 or ... or 128)*/
 struct CHAR3DMAP *ClusMapHead; /* head of CHAR3DMAP list for clusters */
{
 int  x,y,z;
 int Ngrid1,Ngrid0;
 
 printf("#Make_Cluster_for_Binary_CHAR3DMAP_Recursive(Xmap,tar %d mark %d curr %d)\n",tar_bitmask,mark_bitmask,curr_bitmask);
 
 Set_NeighX_NeighY_NeighZ(PAR.NeighborNum,PAR.NeighX, PAR.NeighY, PAR.NeighZ);

 /*** [1] count tar_bitmask and reset mark_bitmask and curr_bitmask ***/  
 Ngrid0 = Ngrid1 = 0;
 for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if ((Xmap->map[x][y][z]&tar_bitmask)  == tar_bitmask) { ++Ngrid1;}
       if ((Xmap->map[x][y][z]&mark_bitmask )== mark_bitmask){ 
          Xmap->map[x][y][z] = Xmap->map[x][y][z]&(~mark_bitmask);
       }
       if ((Xmap->map[x][y][z]&curr_bitmask )== curr_bitmask){ 
          Xmap->map[x][y][z] = Xmap->map[x][y][z]&(~curr_bitmask);
       }
     } 
   } 
 } 
 
 /** [2] Label  Neighbors with map==1, assign ClusterNumber [2..Ncluster+2] to map **/
 /***    map 0 is for background, and map 1 is for foreground ***/

 printf("#Label Neighbors with map==1\n");
 MAX_RECUR_NUM = Xmap->N[0]*Xmap->N[1]*Xmap->N[2];

 for (x=0;x<Xmap->N[0];++x){
   if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]);
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if (((Xmap->map[x][y][z]&tar_bitmask) == tar_bitmask)&&((Xmap->map[x][y][z]&mark_bitmask) != mark_bitmask)){
         Nrecur = 0;
         label_neighbor_with_nomarked_target_Recursive(Xmap,tar_bitmask, mark_bitmask,curr_bitmask,x,y,z);
         add_new_CHAR3DMAP_cluster(Xmap,tar_bitmask, mark_bitmask, curr_bitmask, ClusMapHead);
       }
     } 
   }
 } 

 printf("#Number_of_CHAR3DMAP_Clusters: %d\n",Number_of_CHAR3DMAP_Clusters(ClusMapHead));
 return(1);
 
} /* end of Make_Cluster_for_Binary_CHAR3DMAP_Recursive() */








int label_neighbor_with_nomarked_target_Recursive(Xmap,tar_bitmask, mark_bitmask, curr_bitmask,x,y,z)
  struct CHAR3DMAP *Xmap; /* target 3D map */
  unsigned char tar_bitmask;      /* target   forground bitmask (1 or 2 or 4 or ... or 128)*/
  unsigned char mark_bitmask;     /* marked   forground bitmask (1 or 2 or 4 or ... or 128)*/
  unsigned char curr_bitmask;     /* currrent forground bitmask (1 or 2 or 4 or ... or 128)*/
  int    x,y,z;
{
 int  n,xx,yy,zz; 

 /*
 printf("#label_neighbor_with_nomarked_target(%d %d %d in %d %d %d Xmap %d)\n",x,y,z,
       Xmap->N[0], Xmap->N[1], Xmap->N[2], Xmap->map[x][y][z]); fflush(stdout);
 */
  Nrecur += 1;
  if (Nrecur>MAX_RECUR_NUM){ return(0);}
 
  Xmap->map[x][y][z] = Xmap->map[x][y][z]|curr_bitmask;
  Xmap->map[x][y][z] = Xmap->map[x][y][z]|mark_bitmask;

  for (n=0;n<PAR.NeighborNum;++n){
    xx = x + PAR.NeighX[n];  
    yy = y + PAR.NeighY[n];  
    zz = z + PAR.NeighZ[n];  

    if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){
     /*
     printf("Nrecur %d cluster_num %d n %d xx %d yy %d zz %d map %d\n",Nrecur,cluster_num,n,xx,yy,zz,Xmap->map[xx][yy][zz]); fflush(stdout); 
     */
     if (((Xmap->map[xx][yy][zz]&tar_bitmask)==tar_bitmask)&&((Xmap->map[xx][yy][zz]&mark_bitmask)!=mark_bitmask)){
        label_neighbor_with_nomarked_target_Recursive(Xmap,tar_bitmask, mark_bitmask,curr_bitmask,xx,yy,zz);
     }
    }  
  }
  return(1);

} /* end of label_neighbor_with_nomarked_target_Recursive() */




int Make_Clusters_for_Binary_CHAR3DMAP_Queue(Xmap,tar_bitmask,mark_bitmask,curr_bitmask,ClusMapHead)
 struct CHAR3DMAP *Xmap;         /* target 3D map */
 unsigned char tar_bitmask;      /* target   forground bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char mark_bitmask;     /* marked   forground bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char curr_bitmask;     /* currrent forground bitmask (1 or 2 or 4 or ... or 128)*/
 struct CHAR3DMAP *ClusMapHead; /* head of CHAR3DMAP list for clusters */
{
 int  x,y,z,xx,yy,zz;
 
 printf("#Make_Cluster_for_Binary_CHAR3DMAP_Queue(Xmap,tar %d mark %d curr %d)\n",tar_bitmask,mark_bitmask,curr_bitmask);
 
 Set_NeighX_NeighY_NeighZ(PAR.NeighborNum,PAR.NeighX, PAR.NeighY, PAR.NeighZ);

 /*  Q_MAX = pow(N[0]*N[1]*N[2],2/3) was emplically estimated. See the page 2019/10/17 in the notebook */ 
 Q_MAX = (int)pow((double)(Xmap->N[0]*Xmap->N[1]*Xmap->N[2]),2.0/3.0);   
 /* Q_MAX = Xmap->N[0]*Xmap->N[1]*Xmap->N[2]; */

  printf("#N %d %d %d N012 %d Q_MAX %d\n", Xmap->N[0],Xmap->N[1],Xmap->N[2], Xmap->N[0]*Xmap->N[1]*Xmap->N[2],Q_MAX);

  MALLOC_QUEUE_XYZ();
  Q_STORE_MAX = 0;
  INIT_QUEUE_XYZ();


 /*** [1] count tar_bitmask and reset mark_bitmask and curr_bitmask ***/  
 for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if ((Xmap->map[x][y][z]&mark_bitmask )== mark_bitmask){ 
          Xmap->map[x][y][z] = Xmap->map[x][y][z]&(~mark_bitmask);
       }
       if ((Xmap->map[x][y][z]&curr_bitmask )== curr_bitmask){ 
          Xmap->map[x][y][z] = Xmap->map[x][y][z]&(~curr_bitmask);
       }
     } 
   } 
 } 
 
 /** [2] Label  Neighbors with map==1, assign ClusterNumber [2..Ncluster+2] to map **/
 /***    map 0 is for background, and map 1 is for foreground ***/

 printf("#Label Neighbors with map==1\n");

 for (x=0;x<Xmap->N[0];++x){
   if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]);
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if (((Xmap->map[x][y][z]&tar_bitmask) == tar_bitmask)&&((Xmap->map[x][y][z]&mark_bitmask) != mark_bitmask)){
         INIT_QUEUE_XYZ();
         Xmap->map[x][y][z] = Xmap->map[x][y][z]|curr_bitmask;
         Xmap->map[x][y][z] = Xmap->map[x][y][z]|mark_bitmask;
         PUT_QUEUE_XYZ(x,y,z);

         while (Q_HEAD != Q_TAIL){
           GET_QUEUE_XYZ(&xx,&yy,&zz);
           label_neighbor_with_nomarked_target_Queue(Xmap,tar_bitmask, mark_bitmask,curr_bitmask,xx,yy,zz);
         }
         add_new_CHAR3DMAP_cluster(Xmap,tar_bitmask, mark_bitmask, curr_bitmask, ClusMapHead);
         INIT_QUEUE_XYZ();
       }
     } 
   }
 } 

 printf("#Q_STORE_MAX %d\n",Q_STORE_MAX);
  printf("#Number_of_CHAR3DMAP_Clusters: %d\n",Number_of_CHAR3DMAP_Clusters(ClusMapHead));
 return(1);
 
} /* end of Make_Cluster_for_Binary_CHAR3DMAP_Queue() */




int label_neighbor_with_nomarked_target_Queue(Xmap,tar_bitmask, mark_bitmask, curr_bitmask,x,y,z)
  struct CHAR3DMAP *Xmap; /* target 3D map */
  unsigned char tar_bitmask;      /* target   forground bitmask (1 or 2 or 4 or ... or 128)*/
  unsigned char mark_bitmask;     /* marked   forground bitmask (1 or 2 or 4 or ... or 128)*/
  unsigned char curr_bitmask;     /* currrent forground bitmask (1 or 2 or 4 or ... or 128)*/
  int    x,y,z;
{
 int  n,xx,yy,zz;

 /*
 *  printf("#label_neighbor_with_nomarked_target_Queue(%d %d %d in %d %d %d Xmap %d)\n",x,y,z,
 *         Xmap->N[0], Xmap->N[1], Xmap->N[2], Xmap->map[x][y][z]); fflush(stdout);
 *          */
 
  Xmap->map[x][y][z] = Xmap->map[x][y][z]|curr_bitmask;
  Xmap->map[x][y][z] = Xmap->map[x][y][z]|mark_bitmask;

  for (n=0;n<PAR.NeighborNum;++n){
    xx = x + PAR.NeighX[n];  
    yy = y + PAR.NeighY[n];  
    zz = z + PAR.NeighZ[n];  

    if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){
      if (((Xmap->map[xx][yy][zz]&tar_bitmask)==tar_bitmask)&&((Xmap->map[xx][yy][zz]&mark_bitmask)!=mark_bitmask)){
        Xmap->map[xx][yy][zz] = Xmap->map[xx][yy][zz]|curr_bitmask;
        Xmap->map[xx][yy][zz] = Xmap->map[xx][yy][zz]|mark_bitmask;
        PUT_QUEUE_XYZ(xx,yy,zz);
      }
    }
  }
  return(1);

} /* end of label_neighbor_with_nomarked_target_Queue() */






void add_new_CHAR3DMAP_cluster(Xmap,tar_bitmask, mark_bitmask, curr_bitmask, ClusMapHead)
  struct CHAR3DMAP *Xmap; /* target 3D map */
  unsigned char tar_bitmask;      /* target   forground bitmask (1 or 2 or 4 or ... or 128)*/
  unsigned char mark_bitmask;     /* marked   forground bitmask (1 or 2 or 4 or ... or 128)*/
  unsigned char curr_bitmask;     /* currrent forground bitmask (1 or 2 or 4 or ... or 128)*/
  struct CHAR3DMAP *ClusMapHead; /* head of CHAR3DMAP list for clusters */
{
 int  i,x,y,z,xx,yy,zz,init,Nforeground; 
 int  minN[3], maxN[3];
 struct CHAR3DMAP *cmap;

 /** [1] get minN[] and maxN[]  **/
 minN[0] = minN[1] = minN[2] = Xmap->N[0];
 maxN[0] = maxN[1] = maxN[2] = 0;
 Nforeground = 0;
 init = 1; 
 for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if ((Xmap->map[x][y][z]&curr_bitmask)==curr_bitmask){
         Nforeground += 1;
         if (init==1){
           minN[0] = maxN[0] = x;
           minN[1] = maxN[1] = y;
           minN[2] = maxN[2] = z;
           init = 0; 
         }
         else{
           if (x<minN[0]){ minN[0] = x; }
           if (x>maxN[0]){ maxN[0] = x; }
           if (y<minN[1]){ minN[1] = y; }
           if (y>maxN[1]){ maxN[1] = y; }
           if (z<minN[2]){ minN[2] = z; }
           if (z>maxN[2]){ maxN[2] = z; }
         }
       } 
    }
   }
 }

 /** [2] malloc new CHAR3DMAP  **/
 if ((minN[0]>maxN[0]) || (minN[1]>maxN[1]) || (minN[2]>maxN[2])){
   printf("#ERROR:add_new_CHAR3DMAP_cluster():Nforeground %d\n",Nforeground);
   exit(1);
 }

 cmap = ClusMapHead;
 while (cmap->next !=NULL){
   cmap = cmap->next;
 }
 cmap->next = (struct CHAR3DMAP *)malloc(sizeof(struct CHAR3DMAP));
 cmap->next->prev = cmap;
 cmap->next->next = NULL;
 cmap = cmap->next;

 cmap->grid_width = Xmap->grid_width;
 cmap->Nforeground = Nforeground;
 cmap->value = (float)cmap->Nforeground;
 
 for (i=0;i<3;++i){
   cmap->N[i] = maxN[i] - minN[i] + 1;
   cmap->OrigPos[i] = Xmap->OrigPos[i] + minN[i] * Xmap->grid_width;
   cmap->Noffset[i] = minN[i];
 }

 Malloc_CHAR3DMAP(cmap, cmap->N[0],cmap->N[1],cmap->N[2]);
 Initialize_CHAR3DMAP(cmap);
 printf("#N %d %d %d Nfore %d\n", cmap->N[0],cmap->N[1],cmap->N[2], cmap->Nforeground);

 /** [3] assign in->map  and erace curr_bitmask **/
 for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if ((Xmap->map[x][y][z]&curr_bitmask)==curr_bitmask){
          xx = x - minN[0];
          yy = y - minN[1];
          zz = z - minN[2];
          cmap->map[xx][yy][zz] = cmap->map[xx][yy][zz]|tar_bitmask;
          Xmap->map[x][y][z] = Xmap->map[x][y][z]|mark_bitmask;
          Xmap->map[x][y][z] = Xmap->map[x][y][z]&(~curr_bitmask);
       }
     }
   }
 } 
 
} /* end of add_new_CHAR3DMAP_cluster() */




void Output_Clusters_of_CHAR3DMAP_in_PDB_format(ofname,wmode,ObjMap,ClusMapHead,tar_bitmask)
  char   *ofname;
  char   wmode;
  struct CHAR3DMAP *ObjMap;     /* Original object map */
  struct CHAR3DMAP *ClusMapHead; /* head of CHAR3DMAP list for clusters */
  unsigned char tar_bitmask;
{
  struct CHAR3DMAP *cmap;
  int Ncluster;
  FILE *fpo;


 /** (0) count Ncluster **/
  Ncluster = 0;
  cmap = ClusMapHead;
  while (cmap->next != NULL){
    cmap = cmap->next;
    Ncluster += 1;
    cmap->rank = Ncluster;
  }

 /** (1) output HEADER strings **/
  if (wmode == 'a'){
    fpo = fopen(ofname,"a");
  }
  else{
    fpo = fopen(ofname,"w");
    printf("#Output_Clusters_of_CHAR3DMAP_in_PDB_format(Ncluster %d)-->'%s'\n",Ncluster,ofname);
  }
  fprintf(fpo,"REMARK Output_Clusters_of_CHAR3DMAP_in_PDB_format()\n");
  fprintf(fpo,"REMARK    COMMAND %s\n",PAR.COMMAND);
  Set_END_DATE();
  fprintf(fpo,"REMARK  DATE_START %s\n",PAR.START_DATE);
  fprintf(fpo,"REMARK  DATE_END   %s\n",PAR.END_DATE);
  fprintf(fpo,"REMARK  COMP_TIME  %lf seconds\n",PAR.COMP_TIME_SEC);
  fprintf(fpo,"REMARK Ncluster %d\n",Ncluster);
  fprintf(fpo,"REMARK  grid_size %4d %4d %4d\n",ObjMap->N[0],ObjMap->N[1],ObjMap->N[2]);
  fprintf(fpo,"REMARK  grid_width %8.3f\n",ObjMap->grid_width);
  fprintf(fpo,"REMARK  OrigPos %8.3f %8.3f %8.3f\n",ObjMap->OrigPos[0],ObjMap->OrigPos[1],ObjMap->OrigPos[2]);

  cmap = ClusMapHead;
  while (cmap->next != NULL){
    cmap = cmap->next;
    fprintf(fpo,"#REMARK CLUSTER[%d] NVOXEL  %d\n",cmap->rank,cmap->Nforeground);
  }
  fclose(fpo);

 /** (2) output each cluster **/
  cmap = ClusMapHead;
  while (cmap->next != NULL){
    cmap = cmap->next;
    /*
    Output_CHAR3DMAP_PDB_format(ofname,'a',cmap,1,"",&ComHead);
     */
    Output_Cluster_CHAR3DMAP_in_PDB_format(ofname,'a',cmap,tar_bitmask,cmap->rank);
  }


} /* end of Output_Clusters_of_CHAR3DMAP_in_PDB_format() */




void Output_Clusters_of_CHAR3DMAP_in_CCP4_format(ofname,ClusMapHead,tar_bitmask)
  char   *ofname;
  struct CHAR3DMAP *ClusMapHead; /* head of CHAR3DMAP list for clusters */
  unsigned char tar_bitmask;
{
  struct CHAR3DMAP *cmap;
  int Ncluster;
  char core[MAX_FILE_NAME],ofname_cluster[MAX_FILE_NAME];

 /** (0) count Ncluster **/
  Ncluster = 0;
  cmap = ClusMapHead;
  while (cmap->next != NULL){
    cmap = cmap->next;
    Ncluster += 1;
    cmap->rank = Ncluster;
  }

 Find_Filename_Core_with_Directory(core,ofname);

 /** (1) output each cluster **/
  cmap = ClusMapHead;
  while (cmap->next != NULL){
    cmap = cmap->next;
    sprintf(ofname_cluster,"%s_%d.map",core,cmap->rank);
    Write_MapCCP4_File_CHAR3DMAP(ofname_cluster,cmap,'D');
  }

} /* end of Output_Clusters_of_CHAR3DMAP_in_CCP4_format() */



void Output_Cluster_CHAR3DMAP_in_PDB_format(ofname,mode,M,bitmask,number)
 char   *ofname;
 char   mode;             /* 'w'rite 'a'ppend */
 struct CHAR3DMAP *M;
 unsigned char bitmask;
 int    number;  /* 1,2,3,.... */
{
 int i,j,k,Natom;
 FILE *fp;
 float Pnt[3],tFactor;
 char buff[64],core[128],AtomStr[5];
 char ABC[63];  /* changed in 2019/12/18 char ABC[62]  to ABC[63] */

 sprintf(ABC,"1234567890ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");

 printf("#Output_Cluster_CHAR3DMAP_in_PDB_format(mode %c bitmask %d) -->\"%s\" ",mode,bitmask,ofname);
 if (mode=='a'){fp = fopen(ofname,"a");}
          else {fp = fopen(ofname,"w");}
 if (fp==NULL) { printf("#ERROR:Can't write to \"%s\"\n",ofname); exit(1);}

 Find_Filename_Core(core,ofname);
 Get_Part_Of_Line(buff,ofname,0,29);
 /** Output Basic Grid Information **/
 fprintf(fp,"MODEL%9d\n",number);
 fprintf(fp,"HEADER    GRID     %-30s %s   %s\n",buff,Get_Date_String_PDB(),PAR.pdbidPro);
 fprintf(fp,"REMARK  cluster_grid_size %4d %4d %4d\n",M->N[0],M->N[1],M->N[2]);
 fprintf(fp,"REMARK  cluster_grid_width %8.3f\n",M->grid_width);
 fprintf(fp,"REMARK  cluster_Noffset %4d %4d %4d\n",M->Noffset[0], M->Noffset[1], M->Noffset[2]);
 fprintf(fp,"REMARK  cluster_OrigPos %8.3f %8.3f %8.3f\n",M->OrigPos[0],M->OrigPos[1],M->OrigPos[2]);
 fprintf(fp,"REMARK  If (a grid with value less than [min_val]), do not write the grid.\n");

 fprintf(fp,"REMARK Column of [Residue Number] corresponds to the grid's char value\n");
 fprintf(fp,"REMARK              [grid_char_value]                 [radius][grid_char_value] [i] [j] [k]\n");
 Natom = 0;
 printf("#M->N %d %d %d\n",M->N[0],M->N[1],M->N[2]);
 for (i=0;i<M->N[0];++i){
   for (j=0;j<M->N[1];++j){
     for (k=0;k<M->N[2];++k){
       /* printf("#%d %d %d val %d\n",i,j,k,M->map[i][j][k]);  */
       if ((M->map[i][j][k]&bitmask)==bitmask){
         ++Natom;
         Pnt[0] = M->OrigPos[0] + i*M->grid_width;
         Pnt[1] = M->OrigPos[1] + j*M->grid_width;
         Pnt[2] = M->OrigPos[2] + k*M->grid_width;
       /* 
         sprintf(buff,"X%d",M->map[i][j][k]);
         Get_Part_Of_Line(AtomStr,buff,0,2);
         */
         sprintf(AtomStr,"CA ");
         tFactor = (float)M->map[i][j][k];
         /*
 *          fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f %d %d %d\n",
 *                      Natom%100000,AtomStr,"GRD",' ',M->map[i][j][k],
 *                                  Pnt[0],Pnt[1],Pnt[2],M->grid_width/2.0,tFactor,i,j,k);
 *                                           */
         fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f",
            Natom%100000,AtomStr,"GRD",ABC[(number-1)%62],M->map[i][j][k],Pnt[0],Pnt[1],Pnt[2],
            M->grid_width/2.0,tFactor);

/*
         fprintf(fp," %d %d %d",i,j,k);
*/
         fprintf(fp," %d %d %d",i + M->Noffset[0],j + M->Noffset[1], k + M->Noffset[2]);
         fprintf(fp,"\n");
        }
      }
    }
  }
 fprintf(fp,"ENDMDL\n");
 fclose(fp);

} /* end of Output_Cluster_CHAR3DMAP_in_PDB_format() */



int Number_of_CHAR3DMAP_Clusters(ClusMapHead)
  struct CHAR3DMAP *ClusMapHead;
{
  int Ncluster;
  struct CHAR3DMAP *cmap;

  Ncluster = 0;
  cmap = ClusMapHead;
  while (cmap->next != NULL){
    cmap = cmap->next;
    Ncluster += 1;
    cmap->rank = Ncluster;
  }
  return(Ncluster);
} /* end of Number_of_CHAR3DMAP_Clusters() */

