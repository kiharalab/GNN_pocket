/*

 <OpeClus3DPr.c>
  
  Functions for clustering of 3D-images 

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


#define MAX_RECUR_NUM 100000 

/*** FUNCTIONS (GLOBAL) ***/
void Extract_Cavity_Region();
void Extract_Outerspace_Region();
int  Label_Connected_Bit_Clusters();
int Extract_Neighbor_Bit_Clusters_Around_PocPos();
int Mark_Bits_Around_Pos_In_Direction();
int Diffusion_From_Start_To_End_Regions();


/*** FUNCTIONS (LOCAL) ***/
static int Mark_Neighbor_Bits_Around_Pos();
static int Check_Labeled_Background_Neighbor();
static int  Recursive_Label_Background_Bitmask_Neighbor();
static int  Recursive_Label_Foreground_Bitmask_Neighbor();
static int  Label_Neighbor_With_Map_One();
static void ReLabel_Neighbor_With_Same_Label();
static void Bubble_Sort();


static int Nrecur;

void Extract_Cavity_Region(Xmap,tar_bitmask,out_bitmask,res_bitmask)
  struct CHAR3DMAP *Xmap;     /* target 3D map */
  unsigned char tar_bitmask;  /* target bitmask (1 or 2 or 4 or ... or 128) (X closing P)   */
  unsigned char out_bitmask;  /* outer space  bitmask (1 or 2 or 4 or ... or 128) --> 'outer space'  */
  unsigned char res_bitmask;  /* result bitmask (1 or 2 or 4 or ... or 128) --> 'cavity'    */
{
 int  x,y,z;
 int Ncavity;
 int outsp_point, no_tar_outsp; 
 printf("#Extract_Cavity_Region(Xmap,tar %d res %d)\n",tar_bitmask,res_bitmask);

 /*

  [0] X:target (P-closed X), #:'background'  
  ########### 
  ##XXXXXXX##
  #XXX###XXX#
  #XXX###XXX#
  ##XXXXXXX##
  ###########
  
  [1] o:outer space
  ooooooooooo 
  ooXXXXXXXoo
  oXXX###XXXo
  oXXX###XXXo
  ooXXXXXXXoo
  ooooooooooo 

  [2] @:result (cavity region) 
  ooooooooooo 
  ooXXXXXXXoo
  oXXX@@@XXXo
  oXXX@@@XXXo
  ooXXXXXXXoo
  ooooooooooo 

  */

 /** [1] Labeling "outer space" region (not inner blank space), as out_bitmask **/
 Set_NeighX_NeighY_NeighZ(PAR.NeighborNum,PAR.NeighX,PAR.NeighY,PAR.NeighZ);
 printf("#Labeling 'outer space' region as out_bitmask.\n");
 Nrecur = 0; 
 for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       outsp_point = no_tar_outsp = 0;

       if ((x==0)||(x==(Xmap->N[0]-1))||
           (y==0)||(y==(Xmap->N[1]-1))||
           (z==0)||(z==(Xmap->N[2]-1)) ) {outsp_point = 1;}

       if (((Xmap->map[x][y][z]&tar_bitmask)  != tar_bitmask)&&
           ((Xmap->map[x][y][z]&out_bitmask)  != out_bitmask) ) {no_tar_outsp = 1;}
   
       if ((no_tar_outsp==1) && (outsp_point==0)){
           outsp_point = Check_Labeled_Background_Neighbor(Xmap,x,y,z,tar_bitmask,out_bitmask);
       } 
 
       if ((no_tar_outsp==1) && (outsp_point==1)){ 
          Nrecur = 0;
          Recursive_Label_Background_Bitmask_Neighbor(Xmap,x,y,z,tar_bitmask, out_bitmask);
       }
   }
  }
 }
 
 /** (2) Labeling "cavity" region **/
 Ncavity = 0;
 printf("#Labeling 'cavity' region\n");
 for (x=0;x<Xmap->N[0];++x){
   printf("#x %d/%d\n",x,Xmap->N[0]);
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if (     ((Xmap->map[x][y][z]&tar_bitmask)  != tar_bitmask)
            && ((Xmap->map[x][y][z]&out_bitmask)  != out_bitmask) ){
         Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask; 
        ++Ncavity; 
       } 
     } 
   } 
 } 
 printf("#Nrecur %d\n",Nrecur);
 printf("#Ncavity %d\n",Ncavity);
 
} /* end of Extract_Cavity_Region() */



void Extract_Outerspace_Region(Xmap,tar_bitmask,out_bitmask,res_bitmask)
  struct CHAR3DMAP *Xmap;     /* target 3D map */
  unsigned char tar_bitmask;  /* target bitmask (1 or 2 or 4 or ... or 128) (X closing P)   */
  unsigned char out_bitmask;  /* outer space  bitmask (1 or 2 or 4 or ... or 128) --> 'outer space'  */
  unsigned char res_bitmask;  /* result bitmask (1 or 2 or 4 or ... or 128) --> 'cavity'    */
{
 int  x,y,z;
 int Noutsp_target;
 int outsp_point, tar_no_outsp; 
 printf("#Extract_Outerspace_Region(Xmap,tar %d res %d)\n",tar_bitmask,res_bitmask);

 /*

  [0] X:protein, P:probe, #:'background'  
   P
  PPP
   P 

  ########### 
  ########### 
  ###X##XX###
  ###X###X###
  ###X###X###
  ###XXXXX###
  ########### 
  ########### 

  [1] 1:(X dilation P), #:'background'  

  ########### 
  ###1##11### 
  ##1111111##
  ##111#111##
  ##111#111##
  ##1111111##
  ###11111### 
  ########### 


  [2] T:target (X dilation P)^c

  TTTTTTTTTTT 
  TTT#TT##TTT 
  TT#######TT
  TT###T###TT
  TT###T###TT
  TT#######TT
  TTT#####TTT 
  TTTTTTTTTTT 

 [3] @:outer space and T, ?:not outer space and T

  @@@@@@@@@@@ 
  @@@#@@##@@@ 
  @@#######@@
  @@###?###@@
  @@###?###@@
  @@#######@@
  @@@#####@@@ 
  @@@@@@@@@@@ 

  */

 /** [1] Labeling "outer space" region (not inner blank space), as out_bitmask **/
 Set_NeighX_NeighY_NeighZ(PAR.NeighborNum,PAR.NeighX,PAR.NeighY,PAR.NeighZ);
 printf("#Labeling 'outer space' region as out_bitmask.\n");
 Nrecur = 0; 
 for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       outsp_point = tar_no_outsp = 0;

       if ((x==0)||(x==(Xmap->N[0]-1))||
           (y==0)||(y==(Xmap->N[1]-1))||
           (z==0)||(z==(Xmap->N[2]-1)) ) {outsp_point = 1;}

       if (((Xmap->map[x][y][z]&tar_bitmask)  == tar_bitmask)&&
           ((Xmap->map[x][y][z]&out_bitmask)  != out_bitmask) ) {tar_no_outsp = 1;}
   
       if ((tar_no_outsp==1) && (outsp_point==0)){
           outsp_point = Check_Labeled_Background_Neighbor(Xmap,x,y,z,tar_bitmask,out_bitmask);
       } 
 
       if ((tar_no_outsp==1) && (outsp_point==1)){ 
          Nrecur = 0;
          Recursive_Label_Foreground_Bitmask_Neighbor(Xmap,x,y,z,tar_bitmask, out_bitmask);
       }
   }
  }
 }
 
 /** (2) Labeling "outerspace" target region **/
 Noutsp_target = 0;
 printf("#Labeling outerspace region\n");
 for (x=0;x<Xmap->N[0];++x){
   /* printf("#x %d/%d\n",x,Xmap->N[0]); */
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if (     ((Xmap->map[x][y][z]&tar_bitmask)  == tar_bitmask)
             && ((Xmap->map[x][y][z]&out_bitmask)  == out_bitmask) ){
         Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask; 
         Noutsp_target += 1; 
       } 
     } 
   } 
 } 
 printf("#Nrecur %d\n",Nrecur);
 printf("#Noutsp_target %d\n",Noutsp_target);
 
} /* end of Extract_Outerspace_Region() */





























int Check_Labeled_Background_Neighbor(Xmap,x,y,z,tar_bitmask, label_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 int    x,y,z;
 unsigned char tar_bitmask;   /* target bitmask (1 or 2 or 4 or ...  or 128)*/
 unsigned char label_bitmask; /* label  bitmask (1 or 2 or 4 or ...  or 128)*/
{
  int n,xx,yy,zz;

  for (n=0;n<PAR.NeighborNum;++n){
     xx = x + PAR.NeighX[n];  
     yy = y + PAR.NeighY[n];  
     zz = z + PAR.NeighZ[n];  
        if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){
          if (    ((Xmap->map[xx][yy][zz]&tar_bitmask)  !=  tar_bitmask) 
               && ((Xmap->map[xx][yy][zz]&label_bitmask) == label_bitmask) ){ return(1);}
         }
  }

  return(0);
} /* end of Check_Labeled_Background_Neighbor() */



int Recursive_Label_Background_Bitmask_Neighbor(Xmap,x,y,z,tar_bitmask, label_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 int    x,y,z;
 unsigned char tar_bitmask;   /* target bitmask (1 or 2 or 4 or ...  or 128)*/
 unsigned char label_bitmask; /* label  bitmask (1 or 2 or 4 or ...  or 128)*/
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

   For example,

[takawaba@crambin ~]$ uname -a
Linux crambin 2.6.32-573.el6.x86_64 #1 SMP Thu Jul 23 15:44:03 UTC 2015 x86_64 x86_64 x86_64 GNU/Linux
[takawaba@crambin ~]$ ulimit -s
10240
[kawabata@empiar ~]$ uname -a
Linux empiar.pdbj.org 3.10.0-693.21.1.el7.x86_64 #1 SMP Wed Mar 7 19:03:37 UTC 2018 x86_64 x86_64 x86_64 GNU/Linux
[kawabata@empiar ~]$ ulimit -s
8192

   For avoiding the "stack size overflow", the algorithm stops when the number 
   of recursive calls (Nrecursive) is over MAX_RECUR_NUM.

 */


  /*
 printf("#Label_Bitmask_Neighbor(%8d) %d/%d %d/%d %d/%d (%d)\n",Nrecur,x,Xmap->N[0],y,Xmap->N[1],z,Xmap->N[2],Xmap->map[x][y][z]);
 fflush(stdout);
 */
 
 Xmap->map[x][y][z] = Xmap->map[x][y][z]|label_bitmask; 
 
 for (n=0;n<PAR.NeighborNum;++n){
   xx = x + PAR.NeighX[n];  
   yy = y + PAR.NeighY[n];  
   zz = z + PAR.NeighZ[n];  

   if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){
     if (    ((Xmap->map[xx][yy][zz]&tar_bitmask)  != tar_bitmask) 
          && ((Xmap->map[xx][yy][zz]&label_bitmask)!= label_bitmask))
       Recursive_Label_Background_Bitmask_Neighbor(Xmap,xx,yy,zz,tar_bitmask, label_bitmask);
     }  
  }

 return(1);
} /* end of Recursive_Label_Background_Bitmask_Neighbor() */



int Recursive_Label_Foreground_Bitmask_Neighbor(Xmap,x,y,z,tar_bitmask, label_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 int    x,y,z;
 unsigned char tar_bitmask;   /* target bitmask (1 or 2 or 4 or ...  or 128)*/
 unsigned char label_bitmask; /* label  bitmask (1 or 2 or 4 or ...  or 128)*/
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
 
 Xmap->map[x][y][z] = Xmap->map[x][y][z]|label_bitmask; 
 
 for (n=0;n<PAR.NeighborNum;++n){
   xx = x + PAR.NeighX[n];  
   yy = y + PAR.NeighY[n];  
   zz = z + PAR.NeighZ[n];  

   if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){
     if (    ((Xmap->map[xx][yy][zz]&tar_bitmask)   == tar_bitmask) 
          && ((Xmap->map[xx][yy][zz]&label_bitmask) != label_bitmask))
       Recursive_Label_Foreground_Bitmask_Neighbor(Xmap,xx,yy,zz,tar_bitmask, label_bitmask);
     }  
  }

 return(1);
} /* end of Recursive_Label_Foreground_Bitmask_Neighbor() */




int Label_Connected_Bit_Clusters(Xmap,tar_bitmask,MinNbit_cluster,Number_of_cluster,Nbit_cluster)
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
 
 printf("#Label_Connected_Bit_Clusters(Xmap,tar %d MinNbit_cluster %d)\n",tar_bitmask,MinNbit_cluster);
 
 Set_NeighX_NeighY_NeighZ(PAR.NeighborNum,PAR.NeighX, PAR.NeighY, PAR.NeighZ);

 /** [1] Change Xmap to 0-1 bitmap (if tar_bitmask,  map=1, otherwise map=0) **/
 
 Ngrid0 = Ngrid1 = 0;

 for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
         if((Xmap->map[x][y][z]&tar_bitmask)== tar_bitmask) {Xmap->map[x][y][z] = 1; ++Ngrid1;}
     } 
   } 
 } 
 
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
 
} /* end of Label_Connected_Bit_Clusters() */








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



int Extract_Neighbor_Bit_Clusters_Around_PocPos(Xmap,PocPos,Daround,Number_of_cluster,Nbit_cluster)
 struct CHAR3DMAP *Xmap;     /* target 3D map */
 float  PocPos[3];           /* Forcusing Point */
 float  Daround;             /* Distance around PocPos[] */
 int    *Number_of_cluster;  /* Number of cluster detected (<MAX_UNSIGNED_CHAR)  */
 int    *Nbit_cluster;       /* Nbit_cluster[1..Number_of_cluster] */
{
 int  x,y,z,c,i,Nnew_cluster,Nbit_cluster_old[MAX_UNSIGNED_CHAR];
 int  MaxPoc[3],MinPoc[3],IntPoc[3]; 
 float DifPos[3],DDaround,DD;
 char mark[MAX_UNSIGNED_CHAR],Nnew_from_old[MAX_UNSIGNED_CHAR]; 

 printf("#Extract_Neighbor_Bit_Clusters_Around_PocPos(%f %f %f D %f)\n",
  PocPos[0],PocPos[1],PocPos[2],Daround);
 DDaround = Daround*Daround;

 for (c=0;c<MAX_UNSIGNED_CHAR;++c){ 
    mark[c] = 0; 
    Nbit_cluster_old[c] = Nbit_cluster[c];
    Nbit_cluster[c] = 0;
    Nnew_from_old[c] = 0;
 }

 /** (1) Decide Rough Region for Focus **/
 for (i=0;i<3;++i){
   IntPoc[i] = (int)((PocPos[i] + Daround  - Xmap->OrigPos[i])/Xmap->grid_width);
   MaxPoc[i] = (int)ceil((PocPos[i] + Daround  - Xmap->OrigPos[i])/Xmap->grid_width)  + 1;
   MinPoc[i] = (int)floor((PocPos[i]- Daround  - Xmap->OrigPos[i])/Xmap->grid_width) - 1;
/*
   printf("i %d PocPos %f OrigPos %f MaxPos %d MinPos %d\n",i,PocPos[i],Xmap->OrigPos[i],MaxPoc[i],MinPoc[i]);
*/
   if (MaxPoc[i]<0)           MaxPoc[i] = 0;
   if (MaxPoc[i]>=Xmap->N[i]) MaxPoc[i] = Xmap->N[i]-1; 
   if (MinPoc[i]<0)           MinPoc[i] = 0;
   if (MinPoc[i]>=Xmap->N[i]) MinPoc[i] = Xmap->N[i]-1; 
 }

/*
 printf("MinPoc %d %d %d MaxPoc %d %d %d\n",
 MinPoc[0], MinPoc[1], MinPoc[2], MaxPoc[0], MaxPoc[1], MaxPoc[2]);
*/
 
 if  ((IntPoc[0]<0)||(IntPoc[0]>=Xmap->N[0])
   ||(IntPoc[1]<0)||(IntPoc[1]>=Xmap->N[1])
   ||(IntPoc[2]<0)||(IntPoc[2]>=Xmap->N[2]))
 {
  printf("#ERROR:PocPos(%f %f %f)->(%d %d %d) is out of range.\n",
   PocPos[0], PocPos[1], PocPos[2], IntPoc[0], IntPoc[1], IntPoc[2]);
  exit(1); 
 }

 /** (2) Find bit around PocPos  **/
 for (x=MinPoc[0];x<=MaxPoc[0];++x){
   DifPos[0] = Xmap->OrigPos[0] + Xmap->grid_width * x - PocPos[0];
  for (y=MinPoc[1];y<=MaxPoc[1];++y){
    DifPos[1] = Xmap->OrigPos[1] + Xmap->grid_width * y - PocPos[1];
   for (z=MinPoc[2];z<=MaxPoc[2];++z){
    if (Xmap->map[x][y][z]>0){
     DifPos[2] = Xmap->OrigPos[2] + Xmap->grid_width * z - PocPos[2];
     DD = DifPos[0]*DifPos[0] + DifPos[1]*DifPos[1] + DifPos[2]*DifPos[2];
     if (DD<=DDaround) mark[Xmap->map[x][y][z]] = 1;
    }
   } /* z */
  } /* y */
 } /* x */
 
 /** (3) Change to zero for non_marked bit  **/
 Nnew_cluster = 0;
 for (c=0;c<256;++c)
 {
  if (mark[c]>0) {
   printf("c %d\n",c);
   ++Nnew_cluster; 
   Nnew_from_old[c] = Nnew_cluster; 
   Nbit_cluster[(int)Nnew_from_old[c]] = Nbit_cluster_old[(int)c];
  } 
 }
 printf("#Nnew_cluster %d\n",Nnew_cluster);
 if (Nnew_cluster==0)
 {
  printf("#WARNING:No tunnel cluster is found around (%f %f %f)\n",
   PocPos[0],PocPos[1],PocPos[2]);
 }

 for (x=0;x<Xmap->N[0];++x){
 for (y=0;y<Xmap->N[1];++y){
 for (z=0;z<Xmap->N[2];++z){
  if (Xmap->map[x][y][z]>0) Xmap->map[x][y][z] = Nnew_from_old[Xmap->map[x][y][z]];
   } /* z */
  } /* y */
 } /* x */

 *Number_of_cluster = Nnew_cluster;
 return(Nnew_cluster);


} /* end of Extract_Neighbor_Bit_Clusters_Around_PocPos() */





int Mark_Neighbor_Bits_Around_Pos(Xmap,Pos,Dmin,Dmax,tar_bitmask, res_bitmask)
 struct CHAR3DMAP *Xmap;    /* target 3D map */
 float  Pos[3];             /* Forcusing Point */
 float  Dmin;               /* Min Distance around PocPos[] */
 float  Dmax;               /* Max Distance around PocPos[] */
 unsigned char tar_bitmask; /* target bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char res_bitmask; /* result bitmask (1 or 2 or 4 or ... or 128)*/
{
 int  x,y,z,i;
 int  MaxPoc[3],MinPoc[3],IntPoc[3],Nmark_bit; 
 float DifPos[3],DDmin,DDmax,DD,minDD;
 char  min_ok,max_ok,bit_ok; 

 printf("#Mark_Neighbor_Bits_Around_Pos(%f %f %f Dmin %f Dmax %f tar_bit %d res_bit %d)\n",
  Pos[0],Pos[1],Pos[2],Dmin,Dmax,tar_bitmask,res_bitmask);

 DDmin = Dmin*Dmin;
 DDmax = Dmax*Dmax;

 /** (1) Decide Rough Region for Focus **/
 if (Dmax > 0.0){
  for (i=0;i<3;++i){
    IntPoc[i] = (int)((Pos[i] + Dmax  - Xmap->OrigPos[i])/Xmap->grid_width);
    MaxPoc[i] = (int)ceil((Pos[i] + Dmax  - Xmap->OrigPos[i])/Xmap->grid_width)  + 1;
    MinPoc[i] = (int)floor((Pos[i]- Dmax  - Xmap->OrigPos[i])/Xmap->grid_width) - 1;
    if (MaxPoc[i]<0)           MaxPoc[i] = 0;
    if (MaxPoc[i]>=Xmap->N[i]) MaxPoc[i] = Xmap->N[i]-1; 
    if (MinPoc[i]<0)           MinPoc[i] = 0;
    if (MinPoc[i]>=Xmap->N[i]) MinPoc[i] = Xmap->N[i]-1; 
  /*
   printf("i %d PocPos %f OrigPos %f MaxPos %d MinPos %d\n",
    i,PocPos[i],Xmap->OrigPos[i],MaxPoc[i],MinPoc[i]);
  */
  }
 }
 else
  for (i=0;i<3;++i)
  {
   MinPoc[i] = 0;
   MaxPoc[i] = Xmap->N[i]-1;
  }

 /** (2) Find bit around PocPos  **/
 Nmark_bit = 0; minDD = -1;
 for (x=MinPoc[0];x<=MaxPoc[0];++x){
   DifPos[0] = Xmap->OrigPos[0] + Xmap->grid_width * x - Pos[0];
  for (y=MinPoc[1];y<=MaxPoc[1];++y){
    DifPos[1] = Xmap->OrigPos[1] + Xmap->grid_width * y - Pos[1];
   for (z=MinPoc[2];z<=MaxPoc[2];++z){
    if (Xmap->map[x][y][z]>0){
     DifPos[2] = Xmap->OrigPos[2] + Xmap->grid_width * z - Pos[2];
     DD = DifPos[0]*DifPos[0] + DifPos[1]*DifPos[1] + DifPos[2]*DifPos[2];
     min_ok = max_ok = bit_ok = 0; 
     if ((minDD==-1.0)||(DD<minDD)) minDD = DD;
     if ((Dmin<0) || (DD>=DDmin)) min_ok = 1;
     if ((Dmax<0) || (DD<=DDmax)) max_ok = 1;
     if ((tar_bitmask==0)||((Xmap->map[x][y][z]&tar_bitmask) == tar_bitmask)) bit_ok = 1;
     if (min_ok*max_ok*bit_ok==1) { Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask ; ++Nmark_bit;}
    }
   } /* z */
  } /* y */
 } /* x */

 if (minDD>0.0) printf("#minD %f\n",sqrt(minDD)); 
 printf("#Nmark_bit %d\n",Nmark_bit);
 return(Nmark_bit);

} /* end of Mark_Neighbor_Bits_Around_Pos() */




int Mark_Bits_Around_Pos_In_Direction(Xmap,Pos,Dir,AngMax,Dmin,Dmax,tar_bitmask, res_bitmask)
 struct CHAR3DMAP *Xmap;    /* target 3D map */
 float  Pos[3];             /* Forcusing Point */
 float  Dir[3];             /* Direction */
 float  AngMax;             /* Max Angle (degree) between Dir and Pos->X  */
 float  Dmin;               /* minimum  Distance */
 float  Dmax;               /* maximum Distance */
 unsigned char tar_bitmask; /* target bitmask (1 or 2 or 4 or ... or 128)*/
 unsigned char res_bitmask; /* result bitmask (1 or 2 or 4 or ... or 128)*/
{
 int  x,y,z;
 int  Nmark_bit,Nreject; 
 float XP[3],DotProd,dir[3],cos_thre,D;
 char ang_ok,dmin_ok,dmax_ok;
 
 dir[0] = Dir[0]; dir[1] = Dir[1]; dir[2] = Dir[2];
 Normalize(dir);
 
 cos_thre = cos(AngMax*M_PI/180.0);

 printf("#Mark_Bits_Around_Pos_In_Direction(Pos %f %f %f Dir %f %f %f AngMax %f Dmin %f Dmax %f tar_bit %d res_bit %d)\n",
  Pos[0],Pos[1],Pos[2],Dir[0],Dir[1],Dir[2],AngMax,Dmin,Dmax,tar_bitmask,res_bitmask);

 /* printf("#AngMax %f cos_thre %f\n",AngMax,cos_thre); */

 Nmark_bit = Nreject = 0; 
 for (x=0;x<Xmap->N[0];++x){
   XP[0] = Xmap->OrigPos[0] + Xmap->grid_width * x - Pos[0];
  for (y=0;y<Xmap->N[1];++y){
    XP[1] = Xmap->OrigPos[1] + Xmap->grid_width * y - Pos[1];
   for (z=0;z<Xmap->N[2];++z){
     if ((tar_bitmask==0)||((Xmap->map[x][y][z]&tar_bitmask) == tar_bitmask)){ 
      XP[2] = Xmap->OrigPos[2] + Xmap->grid_width * z - Pos[2];
      D  = Vec_Length(XP);
      if (D>0.0) DotProd = (XP[0]*dir[0] + XP[1]*dir[1] + XP[2]*dir[2])/D; else DotProd = 0.0;

      ang_ok = dmin_ok = dmax_ok = 0;
       if (DotProd>=cos_thre)    ang_ok = 1; 
      if ((Dmin<0.0)||(D>=Dmin)) dmin_ok = 1; 
      if ((Dmax<0.0)||(D<=Dmax)) dmax_ok = 1; 
 
     if ((ang_ok==1)&&(dmin_ok==1)&&(dmax_ok==1))
       { Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask ; 
 /*
      printf("%f %f %f D %f Dprod %f ang %d dmin %d dmax %d\n",
        Xmap->OrigPos[0]+Xmap->grid_width * x,
        Xmap->OrigPos[1]+Xmap->grid_width * y,
        Xmap->OrigPos[2]+Xmap->grid_width * z,
        D,DotProd,ang_ok,dmin_ok,dmax_ok);
 */
         ++Nmark_bit;}
      else {++Nreject;}  
      }
   } /* z */
  } /* y */
 } /* x */

 printf("#Nmark_bit %d Nreject %d\n",Nmark_bit,Nreject);
 return(Nmark_bit);

} /* Mark_Bits_Around_Pos_In_Direction() */



int Diffusion_From_Start_To_End_Regions(Xmap,rgn_bit,sta_bit,end_bit,mrk_bit,tmp_bit)
 struct CHAR3DMAP *Xmap;    /* target 3D map */
 unsigned char rgn_bit; /* region points bitmask (1 or 2 or ... or 128) */
 unsigned char sta_bit; /* start points bitmask (1 or 2 or ... or 128) */
 unsigned char end_bit; /* end   points bitmask (1 or 2 or ... or 128) */
 unsigned char mrk_bit; /* mark  points bitmask (1 or 2 or ... or 128) */
 unsigned char tmp_bit; /* tmp   points bitmask (1 or 2 or ... or 128) */
{
 int  x,y,z,n,xx,yy,zz;
 int  Nround,Nmark_total,Nmark_round;
 char arrive_end,compl_tmp_bit,compl_mrk_bit;

 printf("#Diffusion_From_Start_To_End_Regions(sta_bitmask %d end_bitmask %d mrk_bitmask %d)\n",
  sta_bit,end_bit,mrk_bit);
 
 Set_NeighX_NeighY_NeighZ(PAR.NeighborNum,PAR.NeighX, PAR.NeighY, PAR.NeighZ);
 compl_tmp_bit = ~tmp_bit;
 compl_mrk_bit = ~mrk_bit;

 /** Initialize **/
 for (x=0;x<Xmap->N[0];++x)
  for (y=0;y<Xmap->N[1];++y)
   for (z=0;z<Xmap->N[2];++z) 
   if ((Xmap->map[x][y][z]&rgn_bit) == rgn_bit)   { 
    Xmap->map[x][y][z] = Xmap->map[x][y][z]&compl_tmp_bit;
    Xmap->map[x][y][z] = Xmap->map[x][y][z]&compl_mrk_bit;
    if((Xmap->map[x][y][z]&sta_bit) == sta_bit)  
     Xmap->map[x][y][z] = Xmap->map[x][y][z]|mrk_bit;
   }

 /** Scanning Repeat **/
 arrive_end = 0; 
 Nmark_total = 0; Nround = 0;

 do{ 
 ++Nround;
 Nmark_round = 0; 

  /** (i) if map[xyz] = sta or mrk, then scan neighbors **/
  for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
    for (z=0;z<Xmap->N[2];++z){
      if ((Xmap->map[x][y][z]&mrk_bit)==mrk_bit){
       for (n=0;n<PAR.NeighborNum;++n){
        xx = x + PAR.NeighX[n];  
        yy = y + PAR.NeighY[n];  
        zz = z + PAR.NeighZ[n];  
        if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2]) &&
            ((Xmap->map[xx][yy][zz]&rgn_bit)==rgn_bit)&&((Xmap->map[xx][yy][zz]&mrk_bit)!=mrk_bit))
             Xmap->map[xx][yy][zz] = Xmap->map[xx][yy][zz]|tmp_bit; 
       }/* n */
     } /* if sta or mrk */
    } /* z */
   } /* y */
  } /* x */


 /** (ii) If tmp_bitmask, then mark mrk_bitmask and erase tmp_bitmask.  **/
  for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
    for (z=0;z<Xmap->N[2];++z){
      if ( (Xmap->map[x][y][z]&tmp_bit) == tmp_bit)
       {
          Xmap->map[x][y][z] = Xmap->map[x][y][z]|mrk_bit; 
          Xmap->map[x][y][z] = Xmap->map[x][y][z]&compl_tmp_bit;
          if ((Xmap->map[x][y][z]&end_bit) == end_bit) {arrive_end = 1;}  
          ++Nmark_round;
       }
    } /* z */
   } /* y */
  } /* x */
 
 printf("#Round %d Nmark %d\n",Nround,Nmark_round); 
 Nmark_total += Nmark_round;
 
 if (arrive_end==1) { printf("#Min_Step %d\n",Nround); return(Nround);}

 } while ((Nmark_round>0)&&(arrive_end==0));

 printf("#Min_Step -1. No Pathway from the start to the end\n");
 return(-1);

} /* end of Diffusion_From_Start_To_End_Regions() */



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


