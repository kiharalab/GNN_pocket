/*

 <FldFilQue.c>
  
  Functions for  "flood fill" or "clustering" using "Queue"-based Breadth First Search. 

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


#define MAX_RECUR_NUM 10000000000 

/*** FUNCTIONS (GLOBAL) ***/
void Extract_Outerspace_Region_Recursive();
void Extract_Outerspace_Region_Queue();
int  Label_Connected_Bit_Clusters_Recursive();


/*** FUNCTIONS (LOCAL) ***/
static int Recursive_Label_Neighbors();
static int  Label_Neighbor_With_Map_One();
static void ReLabel_Neighbor_With_Same_Label();
static void Bubble_Sort();
static void Put_Neighbors_in_QUEUE_XYZ();

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
  printf("QUEUE_PUT %d %d %d HEAD %d TAIL %d N %d MAX %d\n",x,y,z,Q_HEAD,Q_TAIL,q_store,Q_MAX);  
  */
  if (q_store >Q_STORE_MAX){
    Q_STORE_MAX = q_store;
  }
  if (Q_TAIL == Q_HEAD){
   /* If one more XYZ is added to the full QUEUE (with Q_MAX XYZs),  release all the XYZs. */
    printf("#WARNING:QUEUE_XYZ IS FULL ! (HEAD %d TAIL %d N %d MAX %d)\n",Q_HEAD, Q_TAIL, STORE_NUM_QUEUE_XYZ(),Q_MAX);
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
   printf("QUEUE_GET %d %d %d  HEAD %d TAIL %d N %d MAX %d\n",*x,*y,*z,Q_HEAD, Q_TAIL,STORE_NUM_QUEUE_XYZ(),Q_MAX);  
  */
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











void Extract_Outerspace_Region_Recursive(Xmap,tar_bitmask,res_bitmask)
  struct CHAR3DMAP *Xmap;     /* target 3D map */
  unsigned char tar_bitmask;  /* target bitmask (1 or 2 or 4 or ... or 128) (X closing P)   */
  unsigned char res_bitmask;  /* result bitmask (1 or 2 or 4 or ... or 128) --> 'cavity'    */
{
 int  x,y,z;
 printf("#Extract_Outerspace_Region_Recursive(Xmap,tar %d res %d)\n",tar_bitmask,res_bitmask);

 /*

  */
 /** [1] Labeling "outer space" region (not inner blank space), as out_bitmask **/
 Set_NeighX_NeighY_NeighZ(PAR.NeighborNum,PAR.NeighX,PAR.NeighY,PAR.NeighZ);
 printf("#Labeling 'outer space' region as out_bitmask.\n");
 Nrecur = 0; 
 Q_COUNT = 0;
 for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if ((x==0)||(x==(Xmap->N[0]-1))||
           (y==0)||(y==(Xmap->N[1]-1))||
           (z==0)||(z==(Xmap->N[2]-1)) ){

         if (((Xmap->map[x][y][z]&tar_bitmask)  == tar_bitmask)&&
             ((Xmap->map[x][y][z]&res_bitmask)  != res_bitmask) ){
       
             Recursive_Label_Neighbors(Xmap,x,y,z,tar_bitmask, res_bitmask);
         }
       }
   }
  }
 }
 printf("#Q_COUNT %d Nrecur %d\n",Q_COUNT,Nrecur); 
 
} /* end of Extract_Outerspace_Region_Recursive() */





int Recursive_Label_Neighbors(Xmap,x,y,z,tar_bitmask, res_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 int    x,y,z;
 unsigned char tar_bitmask;   /* target bitmask (1 or 2 or 4 or ...  or 128)*/
 unsigned char res_bitmask; /* label  bitmask (1 or 2 or 4 or ...  or 128)*/
{
 int n,xx,yy,zz; 

 Nrecur += 1;
 if (Nrecur>MAX_RECUR_NUM){ return(0);}
 /*
  >> CAUTION <<
   Recursive callings of this function "Recursive_Label_Background_Bitmask_Neighbor" 
   easily lead to the system error "stack size overflow", because there are 
   too many empty grids around the figure X.

   ( 
    To show [stack_size]
     > ulimit -s   
    To set [stack_size]
     > ulimit -s [stack_size]   
   ) 

   For avoiding the "stack size overflow", the algorithm stops when the number 
   of recursive calls (Nrecursive) is over MAX_RECUR_NUM.

 */

  /*
 printf("#Label_Bitmask_Neighbor(%8d) %d/%d %d/%d %d/%d (%d)\n",Nrecur,x,Xmap->N[0],y,Xmap->N[1],z,Xmap->N[2],Xmap->map[x][y][z]);
 fflush(stdout);
 */
 
 Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask; 
 Q_COUNT += 1;
 
 for (n=0;n<PAR.NeighborNum;++n){
   xx = x + PAR.NeighX[n];  
   yy = y + PAR.NeighY[n];  
   zz = z + PAR.NeighZ[n];  

   if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){
     if (    ((Xmap->map[xx][yy][zz]&tar_bitmask) == tar_bitmask) 
          && ((Xmap->map[xx][yy][zz]&res_bitmask) != res_bitmask))
       Recursive_Label_Neighbors(Xmap,xx,yy,zz,tar_bitmask, res_bitmask);
     }  
  }

 return(1);
} /* end of Recursive_Label_Neighbors() */




int Label_Connected_Bit_Clusters_Recursive(Xmap,tar_bitmask,MinNbit_cluster,Number_of_cluster,Nbit_cluster)
 struct CHAR3DMAP *Xmap;     /* target 3D map */
 unsigned char tar_bitmask;  /* target bitmask (1 or 2 or 4 or ... or 128)*/
 int    MinNbit_cluster;
 int    *Number_of_cluster;  /* Number of cluster detected (<MAX_UNSIGNED_CHAR)  */
 int    *Nbit_cluster;       /* Nbit_cluster[1..Number_of_cluster] */
{
 int  x,y,z,c;
 int Ncluster,Ngrid1,Ngrid0;
 char Ncluster_char;
 int Ngrid_in_cluster[MAX_UNSIGNED_CHAR],index_cluster[MAX_UNSIGNED_CHAR],old_to_new_cls[MAX_UNSIGNED_CHAR];
 
 printf("#Label_Connected_Bit_Clusters_Recursive(Xmap,tar %d MinNbit_cluster %d)\n",tar_bitmask,MinNbit_cluster);
 
 Set_NeighX_NeighY_NeighZ(PAR.NeighborNum,PAR.NeighX, PAR.NeighY, PAR.NeighZ);

 /** [1] Change Xmap to 0-1 bitmap (if tar_bitmask,  map=1, otherwise map=0) **/
 
 Ngrid0 = Ngrid1 = 0;

 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
   if((Xmap->map[x][y][z]&tar_bitmask)== tar_bitmask) {Xmap->map[x][y][z] = 1; ++Ngrid1;}
                                                 else {Xmap->map[x][y][z] = 0; ++Ngrid0;}
   } /* z */
  } /* y */
 } /* x */
 
 printf("#Ngrid1(tar_bitmask %d) %d Ngrid0 %d\n",tar_bitmask,Ngrid1,Ngrid0);
 /** [2] Label  Neighbors with map==1, assign ClusterNumber [2..Ncluster+2] to map **/
 /***    map 0 is for background, and map 1 is for foreground ***/
 /***    "Ncluster" = [number of cluster] + 1 **/

 Ncluster = 1; Ncluster_char = 1;
 printf("#Label Neighbors with map==1\n");
 for (c=0;c<MAX_UNSIGNED_CHAR;++c) Ngrid_in_cluster[c] = 0;

 for (x=0;x<Xmap->N[0];++x){
  if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]);
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    if (Xmap->map[x][y][z]==1){
      ++Ncluster;  
      
      if (Ncluster>=(MAX_UNSIGNED_CHAR-2)) { 
        printf("#ERROR:Ncluster %d is over MAX_UNSIGNED_CHAR-2 %d !!\n",Ncluster,MAX_UNSIGNED_CHAR-2); 
        Ncluster_char = (Ncluster%128)+2;
        return(0);
      }
      Nrecur = 0;
      Label_Neighbor_With_Map_One(Xmap,x,y,z,Ncluster,&Ngrid_in_cluster[Ncluster]);
      
     /** If the size of cluster is less than threshold, remove that cluster **/ 
     
    if (Ngrid_in_cluster[Ncluster] < MinNbit_cluster){
/*
        printf("#cluster Ngrid %d is smaller than %d\n", Ngrid_in_cluster[Ncluster],MinNbit_cluster);
*/
        ReLabel_Neighbor_With_Same_Label(Xmap,x,y,z,Ncluster,0);
        Ngrid_in_cluster[Ncluster] = 0;
        --Ncluster;
      }
     } 
   } /* z */
  } /* y */
 } /* x */

 /** [3] Sorting cluster_num by its grid_size **/
 printf("#Ngrid0 %d Ngrid1 %d Ncluster %d\n",Ngrid0, Ngrid1, Ncluster);

 /*
 for (c=0;c<=Ncluster;++c) { printf("cluster %d Ngrid %d\n",c,Ngrid_in_cluster[c]); }
 */

 Bubble_Sort(Ngrid_in_cluster,index_cluster,Ncluster+1);
 for (c=0;c<=Ncluster;++c){  
  /*
   printf("%d sorted_cluster %d Ngrid %d\n",c,index_cluster[c],Ngrid_in_cluster[index_cluster[c]]); 
  */
   old_to_new_cls[index_cluster[c]] = c+1;
   Nbit_cluster[c+1] = Ngrid_in_cluster[index_cluster[c]]; 
  }

 
 /** [4] rename of clusters **/
 /** Finally points in cluster n has value n (1<=n<Ncluster) */ 
 
 for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if (Xmap->map[x][y][z]>1){ 
         Xmap->map[x][y][z] = old_to_new_cls[Xmap->map[x][y][z]];
       }
     } 
   } 
 } 

 *Number_of_cluster = Ncluster-1;
 return(1);
 
} /* end of Label_Connected_Bit_Clusters_Recursive() */








int Label_Neighbor_With_Map_One(Xmap,x,y,z,cluster_num,Ngrid_in_cluster)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 int    x,y,z;
 int    cluster_num;
 int    *Ngrid_in_cluster;
{
 int  n,xx,yy,zz; 
/*
 printf("#Label_Neighbor_With_Map_One(%d %d %d map %d clus_num:%d Ngrid_in_clus %d)\n",
 x,y,z,Xmap->map[x][y][z],cluster_num,*Ngrid_in_cluster);
*/
   Nrecur += 1;
   if (Nrecur>MAX_RECUR_NUM){ return(0);}
 
   Xmap->map[x][y][z] = cluster_num;
 
  *Ngrid_in_cluster += 1; 
 
  for (n=0;n<PAR.NeighborNum;++n){
    xx = x + PAR.NeighX[n];  
    yy = y + PAR.NeighY[n];  
    zz = z + PAR.NeighZ[n];  

    if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){
     /*
     printf("Nrecur %d cluster_num %d n %d xx %d yy %d zz %d map %d\n",Nrecur,cluster_num,n,xx,yy,zz,Xmap->map[xx][yy][zz]); fflush(stdout); 
     */
     if (Xmap->map[xx][yy][zz]==1){
        Label_Neighbor_With_Map_One(Xmap,xx,yy,zz,cluster_num,Ngrid_in_cluster); 
     }
    }  
  }
  return(1);
} /* end of Label_Neighbor_With_Map_One() */





void ReLabel_Neighbor_With_Same_Label(Xmap,x,y,z,target_num, new_num)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 int    x,y,z;
 unsigned char   target_num;
 unsigned char   new_num;
{
 int  n,xx,yy,zz; 

 Xmap->map[x][y][z] = new_num;
 Set_NeighX_NeighY_NeighZ(PAR.NeighborNum,PAR.NeighX,PAR.NeighY,PAR.NeighZ); 
 
 for (n=0;n<PAR.NeighborNum;++n){
   xx = x + PAR.NeighX[n];  
   yy = y + PAR.NeighY[n];  
   zz = z + PAR.NeighZ[n];  

   if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){
    if (Xmap->map[xx][yy][zz]==target_num) 
     ReLabel_Neighbor_With_Same_Label(Xmap,xx,yy,zz,target_num,new_num); 
     }  
  }

} /* end of ReLabel_Neighbor_With_Same_Label() */















void Bubble_Sort(X,index,N)
 int  *X;
 int  *index;
 int   N; /* number of array */
{
 int i,Nchange,buff_i;
 float buff;
 int *XX; 
        
 XX = (int *)malloc(sizeof(int)*(N+1));
 for (i=0;i<N;++i) 
 { XX[i]    = X[i];
   index[i] = i; } 

 do{
 Nchange = 0;
 for (i=1;i<N;++i){
   if (XX[i-1]<XX[i]){ 
    buff   = XX[i-1]; XX[i-1] = XX[i]; XX[i] = buff;
    buff_i = index[i-1]; index[i-1] = index[i]; index[i] = buff_i;
     ++Nchange; }
  }
 }while (Nchange >0);

 free(XX);
} /* end of Bubble_Sort() */



void Extract_Outerspace_Region_Queue(Xmap,tar_bitmask,res_bitmask)
  struct CHAR3DMAP *Xmap;     /* target 3D map */
  unsigned char tar_bitmask;  /* target bitmask (1 or 2 or 4 or ... or 128) (X closing P)   */
  unsigned char res_bitmask;  /* result bitmask (1 or 2 or 4 or ... or 128) --> 'outer_space'    */
{
 int  x,y,z,xx,yy,zz;

 printf("#Extract_Outerspace_Region_Queue(Xmap,tar %d res %d)\n",tar_bitmask,res_bitmask);

 /*
  [0] T:target (X dilation P)^c

  TTTTTTTTTTT 
  TTT#TT##TTT 
  TT#######TT
  TT###T###TT
  TT###T###TT
  TT#######TT
  TTT#####TTT 
  TTTTTTTTTTT 

  [1] @:outer space and T, T:not outer space 

  @@@@@@@@@@@ 
  @@@#@@##@@@ 
  @@#######@@
  @@###T###@@
  @@###T###@@
  @@#######@@
  @@@#####@@@ 
  @@@@@@@@@@@ 
  */
  /* Q_MAX =  PAR.NeighborNum; */
/*
  Q_MAX = 2*(Xmap->N[0]*Xmap->N[1] + Xmap->N[1]*Xmap->N[2] + Xmap->N[0]*Xmap->N[2]);
  Q_MAX = 10*(Xmap->N[0]*Xmap->N[1] + Xmap->N[1]*Xmap->N[2] + Xmap->N[0]*Xmap->N[2]);
*/
 /*  Q_MAX = pow(N[0]*N[1]*N[2],2/3) was emplically estimated. See the page 2019/10/17 in the notebook */  
  Q_MAX = (int)pow((double)(Xmap->N[0]*Xmap->N[1]*Xmap->N[2]),2.0/3.0);  
 
  Q_MAX = (int)(Xmap->N[0]*Xmap->N[1]*Xmap->N[2]);


  printf("#Q_MAX %d\n",Q_MAX);
  MALLOC_QUEUE_XYZ();
  Q_STORE_MAX = 0;
  INIT_QUEUE_XYZ();
  /** [1] Labeling "outer space" region (not inner blank space), as res_bitmask **/
  Set_NeighX_NeighY_NeighZ(PAR.NeighborNum,PAR.NeighX,PAR.NeighY,PAR.NeighZ);
  Q_COUNT = 0;
  for (x=0;x<Xmap->N[0];++x){
    for (y=0;y<Xmap->N[1];++y){
      for (z=0;z<Xmap->N[2];++z){
        if ((x==0)||(x==(Xmap->N[0]-1))||
            (y==0)||(y==(Xmap->N[1]-1))||
            (z==0)||(z==(Xmap->N[2]-1)) ){
        if (((Xmap->map[x][y][z]&tar_bitmask)  == tar_bitmask)&&
           ((Xmap->map[x][y][z]&res_bitmask)  != res_bitmask) ) {

          INIT_QUEUE_XYZ();

          Put_Neighbors_in_QUEUE_XYZ(x,y,z,Xmap,tar_bitmask,  res_bitmask);

          while (Q_HEAD != Q_TAIL){
            GET_QUEUE_XYZ(&xx,&yy,&zz);
            Put_Neighbors_in_QUEUE_XYZ(xx,yy,zz,Xmap,tar_bitmask,  res_bitmask);
          }
         }
       }
     }
   }
 }

 printf("#Q_COUNT %d Q_STORE_MAX %d\n",Q_COUNT,Q_STORE_MAX);
/*
 FREE_QUEUE_XYZ();
 */
} /* end of Extract_Outerspace_Region_Queue() */


void Put_Neighbors_in_QUEUE_XYZ(x,y,z,Xmap,tar_bitmask, res_bitmask)
  int x,y,z;
  struct CHAR3DMAP *Xmap;     /* target 3D map */
  unsigned char tar_bitmask;  /* target bitmask (1 or 2 or 4 or ... or 128) (X closing P)   */
  unsigned char res_bitmask;  /* result bitmask (1 or 2 or 4 or ... or 128) --> 'cavity'    */

{
 int n,xx,yy,zz;

 Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask; 
 Q_COUNT += 1; 
 for (n=0;n<PAR.NeighborNum;++n){
   xx = x + PAR.NeighX[n];  
   yy = y + PAR.NeighY[n];  
   zz = z + PAR.NeighZ[n];  

   if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){
     if (    ((Xmap->map[xx][yy][zz]&tar_bitmask) == tar_bitmask) 
          && ((Xmap->map[xx][yy][zz]&res_bitmask) != res_bitmask)){
/*
        printf("PUT xx_yy_zz %d %d %d tar %d res %d\n",xx,yy,zz ,
             Xmap->map[xx][yy][zz]&tar_bitmask, Xmap->map[xx][yy][zz]&res_bitmask);
*/
        Xmap->map[xx][yy][zz] = Xmap->map[xx][yy][zz]|res_bitmask; 
        PUT_QUEUE_XYZ(xx,yy,zz);
     }  
  }
 }
} /* end of Put_Neighbors_in_QUEUE_XYZ() */


