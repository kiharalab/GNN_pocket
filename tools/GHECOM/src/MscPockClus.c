/*

 <MscPockClus.c>
 
 for Clustering Multiscale Pocket

 20140326: The 'Recurrent' algorithm for Label_Connected_Neighbors() 
            was replaced into the "Stack" algorithm.


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

/*** FUNCTIONS (GLOBAL) ***/
void Cluster_MultiScale_Pocket_Single_Linkage();
void Sort_GridCluster();
void Output_CHAR3DMAP_PDB_format_By_Grid_Cluster();
void Set_Map_to_Zero_for_Low_Ranked_Grid_Cluster();
void Add_Grid_Member();
void Free_Grid_Member();
int Check_Overlap_of_Grid_Member();
int Number_Of_Grid_Member();

/*** FUNCTIONS (LOCAL) ***/
static void Label_Connected_Neighbors_Recurrent();
static void Label_Connected_Neighbors_Stack();
static void Set_Grid_Cluster_Properties();


/*** FUNCTIONS (LOCAL) for the function Label_Connected_Neighbors_Stack() ***/
static void malloc_stack_arrayXYZ();
static void free_stack_arrayXYZ();
static void pushXYZ();
static void popXYZ();
static void stack_initXYZ();
static int  stack_emptyXYZ();
static int not_over_max();

/*** VARIABLES (LOCAL) for the function Label_Connected_Neighbors_Stack() ***/
int Nbyte_malloc;
static short *stack_arrayX, *stack_arrayY, *stack_arrayZ;
static int   top_stackX,top_stackY,top_stackZ;



void malloc_stack_arrayXYZ(N)
  int N;
{
  stack_arrayX = (short *)malloc(sizeof(short)*N);
  stack_arrayY = (short *)malloc(sizeof(short)*N);
  stack_arrayZ = (short *)malloc(sizeof(short)*N);
} 

void free_stack_arrayXYZ()
{
  free(stack_arrayX);
  free(stack_arrayY);
  free(stack_arrayZ);
}

void pushXYZ(x,y,z)
 int x,y,z;
{
  stack_arrayX[top_stackX++] = x;
  stack_arrayY[top_stackY++] = y;
  stack_arrayZ[top_stackZ++] = z;
}

void popXYZ(x,y,z)
 int *x,*y,*z;
{
  *x = stack_arrayX[--top_stackX];
  *y = stack_arrayY[--top_stackY];
  *z = stack_arrayZ[--top_stackZ];
}

void stack_initXYZ()
{
  top_stackX = 0;
  top_stackY = 0;
  top_stackZ = 0;
}

int stack_emptyXYZ()
{
  return(!top_stackX);
}








void Cluster_MultiScale_Pocket_Single_Linkage(Xmap,Zmap,GclsHead,thre_voxel_val,thre_type,MS,SortType)
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
 struct CHAR3DMAP *Zmap;  
  /* map for clustering (should be initialized to 0) 
  >> After calculation << 
  =0     --> not labeled.
  =255   --> in the stack. 
  =1.254 --> cluster label.
 */
 struct GRID_CLUSTER *GclsHead;
 int  thre_voxel_val; /* Threshold Voxel Value (1 <= thre_voxel_val <= 254) */  
 char thre_type;      /*" m'I'n:map[][][]>=thre, ma'X':map[][][]<=thre */ 
 struct MULSC_SHELL_PRB *MS;  /* this is necessary, only for SortType == 'I' */
 char   SortType;        /* sorted by 'N'member, 'V'olume, 'I':invRvolume */
{
 int x,y,z,Ncluster,Nnonzero;
 struct GRID_CLUSTER *gcls,*gcls0;
 struct GRID_MEMBER *gmem;
 int *num_index;
 char   AlgoType;        /* 'R'ecurrent, 'S'tack' */
 
 Set_NeighX_NeighY_NeighZ(PAR.NeighborNum,PAR.NeighX,PAR.NeighY,PAR.NeighZ);
 AlgoType = 'S'; 

 printf("#Cluster_MultiScale_Pocket_Single_Linkage(thre_voxel_val %d thre_type %c AlgoType %c)\n",thre_voxel_val,thre_type,AlgoType);

 if (AlgoType=='S'){
   Nnonzero = Count_NonZero_3DMAP(Xmap);
   malloc_stack_arrayXYZ(Nnonzero);
   stack_initXYZ();
 }

 /*** (1) Cluster Xmap and write cluster_num into Zmap ***/
 Initialize_CHAR3DMAP(Zmap);
 Ncluster = 0;
 Nbyte_malloc = 0; 
 GclsHead->next = NULL;
 gcls = GclsHead; 

 for (x=0;x<Xmap->N[0];++x){
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
   
   if ((Xmap->map[x][y][z]>0)&&(Xmap->map[x][y][z]<255)&&(Zmap->map[x][y][z]==0)){
    if ( ((thre_type=='X')&&(Xmap->map[x][y][z]<= thre_voxel_val)) ||
         ((thre_type=='I')&&(Xmap->map[x][y][z]>= thre_voxel_val))  ){
        gcls0 = GclsHead->next;
        GclsHead->next = (struct GRID_CLUSTER*)malloc(sizeof(struct GRID_CLUSTER));
        gcls = GclsHead->next;
        gcls->next = gcls0;
        gcls->prev = GclsHead; 
        if (gcls0 != NULL) gcls0->prev = gcls;
        gcls->HeadGridMember.next = NULL;
        gcls->Nmember = 0;
        ++Ncluster;  
        gcls->num = Ncluster;
        if (AlgoType=='R'){
          Label_Connected_Neighbors_Recurrent(Xmap,Zmap,x,y,z,thre_voxel_val,thre_type,gcls);
        }
        else if (AlgoType=='S'){
          Label_Connected_Neighbors_Stack(Xmap,Zmap,x,y,z,thre_voxel_val,thre_type,gcls);
        }
/*
        printf("Ncluster %d gcls->Nmember %d malloc %d byte\n",Ncluster,gcls->Nmember,Nbyte_malloc); fflush(stdout); 
*/
      }
     }
   } /* z */ 
  } /* y */
 } /* x */

 printf("#Ncluster %d\n",Ncluster);
 
 /*** (2) Set Grid Grid Cluster Properties ***/
  gcls = GclsHead;
  while (gcls->next != NULL){
   gcls = gcls->next;
   gcls->Nmember = 0;
   gmem = &(gcls->HeadGridMember);
   while (gmem->next !=NULL){ 
     gmem = gmem->next;
     ++(gcls->Nmember); 
   } 
   if (SortType=='I') Set_Grid_Cluster_Properties(gcls,Xmap,MS);
  }

 /*** (3) Sorting cluster  ***/
 /*
 gcls = GclsHead;
 while (gcls->next!=NULL){
  gcls = gcls->next;
  printf("#cluster [%d] Nmember %d\n",gcls->num,gcls->Nmember);
 }
 */

 num_index = (int*)malloc(sizeof(int)*(Ncluster+1));
 Sort_GridCluster(GclsHead,num_index,SortType);

 /*
 for (x=1;x<=Ncluster;++x){ printf("#[old] %d -> [new] %d\n",x,num_index[x]); }
 */

 /* 
 gcls = GclsHead;
 while (gcls->next!=NULL){
  gcls = gcls->next;
  printf("#cluster [%d] Nmember %d\n",gcls->num,gcls->Nmember);
 }
 */

 /*** (3) Renumbering Zmap cluster number ***/
 for (x=0;x<Zmap->N[0];++x){
  for (y=0;y<Zmap->N[1];++y){
   for (z=0;z<Zmap->N[2];++z){
     if (Zmap->map[x][y][z]>0) Zmap->map[x][y][z] = num_index[Zmap->map[x][y][z]];  
   }
  }
 }
 free(num_index);

 if (AlgoType=='S'){ free_stack_arrayXYZ();}

 /*** test for Free_Grid_Member() ***/
 /*
 gcls = GclsHead;
 while (gcls->next!=NULL){
   gcls = gcls->next;
   Free_Grid_Member(&(gcls->HeadGridMember)); 

 }
 */
 /*
 gcls = GclsHead;
 while (gcls->next!=NULL){
   gcls = gcls->next;
   printf("#>gcls %d Nmember %d\n",gcls->num,gcls->Nmember);
   gmem = &(gcls->HeadGridMember);
   while (gmem->next != NULL){
     gmem = gmem->next;
     printf(" %d %d %d\n",gmem->x,gmem->y,gmem->z);
   }
 }
 exit(1);
 */
} /* end of Cluster_MultiScale_Pocket_Single_Linkage() */







void Sort_GridCluster(Gcls,num_index,SortType)
 struct GRID_CLUSTER *Gcls;
 int    *num_index;  /* num_index[old_num] -> new_num */ 
 char   SortType;    /* sorted by 'N'member, 'V'olume, 'I':invRvolume */
{
 struct GRID_CLUSTER *gcls;
 int Ncluster;
 
 /* (1) Recount Nmember  */
 gcls = Gcls; 
 while (gcls->next != NULL){
  gcls = gcls->next;
  if (SortType=='N') gcls->value = (float)(gcls->Nmember); 
  if (SortType=='V') gcls->value = gcls->Volume; 
  if (SortType=='I') gcls->value = gcls->invRvolume; 
 }

 /* (2) Sort */
 Merge_Sort_Double_Linked_List_GRID_CLUSTER(Gcls,'R');
 
 /* (3) Renumbering gcls->num */
 Ncluster = 0;
 gcls = Gcls;
 while (gcls->next != NULL){
  gcls = gcls->next;
  ++Ncluster;
  num_index[gcls->num] = Ncluster;
  gcls->num = Ncluster;
 }
} /* end of Sort_GridCluster() */








void Label_Connected_Neighbors_Recurrent(Xmap,Ymap,x,y,z,thre_voxel_val,thre_type,gclus)
 struct CHAR3DMAP *Xmap; /* target 3D pocket map */
 struct CHAR3DMAP *Ymap; /* cluster label 3D map */
 int    x,y,z;
 unsigned char  thre_voxel_val; /* Threshold Voxel Value (1 <= thre_voxel_val <= 254) */  
 char   thre_type;      /*" m'I'n:map[][][]>=thre, ma'X':map[][][]<=thre */ 
 struct GRID_CLUSTER *gclus;
{
 int  n,xx,yy,zz;

 /*
 printf("#Label_Connected_Neighbors_Recurrent(%d %d %d val %d gnum %d Nmember %d thre %d)\n",
   x,y,z,Xmap->map[x][y][z],gclus->num,gclus->Nmember,thre_voxel_val); fflush(stdout);
 */
 
 Ymap->map[x][y][z] = gclus->num;
 Add_Grid_Member(&(gclus->HeadGridMember),x,y,z);
 gclus->Nmember += 1;
 
 for (n=0;n<PAR.NeighborNum;++n){
   xx = x + PAR.NeighX[n];
   yy = y + PAR.NeighY[n];
   zz = z + PAR.NeighZ[n];
   if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){
    if ((Xmap->map[xx][yy][zz]>0)&&(Xmap->map[xx][yy][zz]<255)&&(Ymap->map[xx][yy][zz]==0)){
      if ( ((thre_type=='X')&&(Xmap->map[xx][yy][zz] <= thre_voxel_val)) ||
           ((thre_type=='I')&&(Xmap->map[xx][yy][zz] >= thre_voxel_val))  ){
          Label_Connected_Neighbors_Recurrent(Xmap,Ymap,xx,yy,zz,thre_voxel_val,thre_type,gclus);
      } 
     }
    }
 } /* n */

} /* end of Label_Connected_Neighbors_Recurrent() */


void Label_Connected_Neighbors_Stack(Xmap,Ymap,x0,y0,z0,thre_voxel_val,thre_type,gclus)
 struct CHAR3DMAP *Xmap; /* target 3D pocket map */
 struct CHAR3DMAP *Ymap; /* cluster label 3D map. 
   Ymap[x][y][z]=0     --> not labeled.
   Ymap[x][y][z]=255   --> in the stack. 
   Ymap[x][y][z]=1.254 --> cluster label.
 Therefore, the number of cluster should be 1<=cluster_num<=254.
*/
 int    x0,y0,z0;
 unsigned char  thre_voxel_val; /* Threshold Voxel Value (1 <= thre_voxel_val <= 254) */  
 char   thre_type;      /*" m'I'n:map[][][]>=thre, ma'X':map[][][]<=thre */ 
 struct GRID_CLUSTER *gclus;
{
 int  n,x,y,z,xx,yy,zz;

/*
 printf("#Label_Connected_Neighbors_Stack(%d %d %d val %d gnum %d Nmember %d thre %d)\n",
   x0,y0,z0,Xmap->map[x][y][z],gclus->num,gclus->Nmember,thre_voxel_val); fflush(stdout);
 */

 pushXYZ(x0,y0,z0);
 /* printf("#push %d %d %d\n",x0,y0,z0); fflush(stdout); */

 while (!stack_emptyXYZ()){
   popXYZ(&x,&y,&z);
   /* printf("#pop %d %d %d\n",x,y,z); fflush(stdout); */
   Ymap->map[x][y][z] = gclus->num;
   gclus->Nmember += 1;
   Add_Grid_Member(&(gclus->HeadGridMember),x,y,z);

   for (n=0;n<PAR.NeighborNum;++n){
     xx = x + PAR.NeighX[n];
     yy = y + PAR.NeighY[n];
     zz = z + PAR.NeighZ[n];
     if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){
      if ((Xmap->map[xx][yy][zz]>0)&&(Xmap->map[xx][yy][zz]<255)&&(Ymap->map[xx][yy][zz]==0)){
        if ( ((thre_type=='X')&&(Xmap->map[xx][yy][zz] <= thre_voxel_val)) ||
             ((thre_type=='I')&&(Xmap->map[xx][yy][zz] >= thre_voxel_val))  ){
            pushXYZ(xx,yy,zz);
            Ymap->map[xx][yy][zz] = 255; /* If Ymap[x][y][z]=255, (x,y,z) is in the stack. */
        } 
       }
      }
   } 
 }

} /* end of Label_Connected_Neighbors_Stack() */









void Add_Grid_Member(gmem_head,x,y,z)
 struct GRID_MEMBER *gmem_head;
 int x,y,z;
{
 struct GRID_MEMBER *gm_pre;

 /*
 printf("Add_Grid_Member(%d %d %d Nbyte_malloc %d) \n",x,y,z,Nbyte_malloc);  fflush(stdout); 
 */

 gm_pre = gmem_head->next;  
 gmem_head->next = (struct GRID_MEMBER*)malloc(sizeof(struct GRID_MEMBER));
 Nbyte_malloc += sizeof(struct GRID_MEMBER);
 if (Nbyte_malloc >= PAR.MAX_MEMORY_MEGA * 1024 * 1024){
   printf("#ERROR: memory over. Nbyte_malloc %d byte MAX_MEMORY_MEGA %f\n",Nbyte_malloc,PAR.MAX_MEMORY_MEGA);
   exit(1);
 }
 if (gmem_head->next == NULL){
   printf("#ERROR: malloc failed. (Add_Grid_Member())\n");
   exit(1);
 }
 gmem_head->next->next = gm_pre;
 /*
 gmem_head->next->prev = gmem_head;
 if (gm_pre != NULL) gm_pre->prev = gmem_head->next;
  */
 /*
 gmem_head->next->x = x;
 gmem_head->next->y = y;
 gmem_head->next->z = z;
 */
 gmem_head->next->x = (short)x;
 gmem_head->next->y = (short)y;
 gmem_head->next->z = (short)z;
 /*
 gmem_head->next->value = value;
 gmem_head->next->type = type;
 */

} /* end of Add_Grid_Member() */



void Free_Grid_Member(gmem_head)
 struct GRID_MEMBER *gmem_head;
{
 /*
 struct GRID_MEMBER *gm,*gm_prev;
 */
 struct GRID_MEMBER *gm,*gm_next;

 /* 
 printf("#Free_Grid_Member()\n"); fflush(stdout);
 */

/*
 gm = gmem_head;
 while (gm->next != NULL) { gm = gm->next; }
 
 while (gm != gmem_head){ 
   gm_prev = gm->prev; 
   free(gm);
   gm = gm_prev;
  }
 gmem_head->next = NULL;
*/
 gm = gmem_head;
 while ((gm!=NULL) && (gm->next != NULL)){ 
   gm = gm->next; 
   gm_next = gm->next;
   free(gm);
   gm = gm_next;
 } 
 gmem_head->next = NULL;

} /* end of Free_Grid_Member() */


int Check_Overlap_of_Grid_Member(gmem_head,x,y,z)
 struct GRID_MEMBER *gmem_head;
 int x,y,z;
{
 struct GRID_MEMBER *gm;
 
 gm = gmem_head;
 while (gm->next != NULL) 
 { gm = gm->next; 
   if ((gm->x == x) && (gm->y == y) && (gm->z == z)) return(1);
  }
 return(0); 

} /* end of Check_Overlap_of_Grid_Member() */



int Number_Of_Grid_Member(gmem_head)
 struct GRID_MEMBER *gmem_head;
{
 struct GRID_MEMBER *gm;
 int Ngmem;

 Ngmem = 0;
 gm = gmem_head;
 while (gm->next != NULL) { gm = gm->next; ++Ngmem;}
 return(Ngmem);

} /* end of Number_Of_Grid_Member() */



void Set_Grid_Cluster_Properties(gcls,M,MS)
 struct GRID_CLUSTER *gcls;
 struct CHAR3DMAP *M;
 struct MULSC_SHELL_PRB *MS;
{
 struct GRID_MEMBER *gmem;
 int x,y,z,v,init,n,k;
 float Rinacc,Vvox;
 float sum_invRinacc,ave_invRinacc;

 Vvox = M->grid_width * M->grid_width * M->grid_width;
 sum_invRinacc = 0.0;
 init = 1;
 for (v=0;v<256;++v) gcls->NRinacc[v] = 0;
 n = 0;
 for (k=0;k<3;++k) gcls->Pos[k] = 0.0;
 gcls->Nmember = 0; 
 
 gmem = &(gcls->HeadGridMember);
 while (gmem->next != NULL){
   gmem = gmem->next;
   x = gmem->x;
   y = gmem->y;
   z = gmem->z;
   if ((x>=0)&&(x<M->N[0])&&(y>=0)&&(y<M->N[1])&&(z>=0)&&(z<M->N[2])){
     v      = M->map[x][y][z];
     Rinacc = MS->Rarray[v];
     gcls->NRinacc[v] += 1;
     sum_invRinacc += 1.0/Rinacc;
     gcls->Pos[0] += M->OrigPos[0] + M->grid_width * x; 
     gcls->Pos[1] += M->OrigPos[1] + M->grid_width * y; 
     gcls->Pos[2] += M->OrigPos[2] + M->grid_width * z; 
     (gcls->Nmember) += 1; 
    /* printf("[%d-%d/%d]Rinacc %f\n",gcls->num,n,gcls->Nmember,Rinacc); */
     n += 1; 
     if (init==1) {gcls->minRinacc = gcls->maxRinacc = Rinacc; init = 0;}
     else{
       if (Rinacc < gcls->minRinacc) gcls->minRinacc = Rinacc;
       if (Rinacc > gcls->maxRinacc) gcls->maxRinacc = Rinacc;
     } 
   }
 } /* gmem */

   for (k=0;k<3;++k){ gcls->Pos[k] /= gcls->Nmember;}

   gcls->invRvolume = sum_invRinacc * Vvox;
   ave_invRinacc   = sum_invRinacc/(float)gcls->Nmember;
   gcls->aveRinacc = 1.0/ave_invRinacc;
   gcls->Volume = (float)gcls->Nmember * Vvox;
 /*
  printf("sum_invRinacc %f ave_invvRinacc %f Nmember %d gcls->aveRinacc %f\n",
   sum_invRinacc,ave_invRinacc,gcls->Nmember,gcls->aveRinacc);
 */

} /* end of Set_Grid_Cluster_Properties() */



void Output_CHAR3DMAP_PDB_format_By_Grid_Cluster(ofname,mode,M,Gclus,MS,minNgridClus,Nrank_keep,title,CommentHead)
 char   *ofname;
 char   mode;  /* 'w'rite, 'a'ppend */
 struct CHAR3DMAP *M;
 struct GRID_CLUSTER *Gclus; 
 struct MULSC_SHELL_PRB *MS;
 int    minNgridClus;     /* Minimum Number of Grids in One Cluster */ 
 int    Nrank_keep;       /* if (Nrank_keep < -1), output everything */
 char   *title;
 struct STRING_NODE *CommentHead;
{
 int x,y,z,r,Natom,Nclus;
 FILE *fp;
 float Pnt[3],R,invR;
 char buff[128],AtomStr[5],core[128];
 struct STRING_NODE *sn;
 struct GRID_CLUSTER *gcls;
 struct GRID_MEMBER  *gmem;
 char  chain_str[256];
 int   Ngrid_total, len_chain_str;
 float Volume_total;
 

 sprintf(chain_str,"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
 len_chain_str = strlen(chain_str);
 
 printf("#Output_CHAR3DMAP_PDB_format_By_Grid_Cluster() -->\"%s\"\n",ofname);
 if (mode=='a') fp = fopen(ofname,"a");
           else fp = fopen(ofname,"w");
 if (fp==NULL) { printf("#ERROR:Can't write to \"%s\"\n",ofname); exit(1);}
 Find_Filename_Core(core,ofname);
 Get_Part_Of_Line(buff,ofname,0,29);
 fprintf(fp,"HEADER    GRID     %-30s %s   %s\n",buff,Get_Date_String_PDB(),PAR.pdbidPro);
 if (title[0]!='\0') fprintf(fp,"TITLE   %s\n",title);
 fprintf(fp,"REMARK  minNgridClus %d Nrank_keep %d\n",minNgridClus,Nrank_keep);
 fprintf(fp,"REMARK  COMMAND %s\n",PAR.COMMAND);
 fprintf(fp,"REMARK  DATE    %s\n",Get_Date_String());
 fprintf(fp,"REMARK  grid_size %4d %4d %4d\n",M->N[0],M->N[1],M->N[2]);
 fprintf(fp,"REMARK  grid_width %8.3f\n",M->grid_width);
 fprintf(fp,"REMARK  OrigPos %8.3f%8.3f %8.3f\n",M->OrigPos[0],M->OrigPos[1],M->OrigPos[2]);
 
 sn = CommentHead;
 while (sn->next != NULL){ 
   sn = sn->next;
   fprintf(fp,"REMARK  COMMENT %s\n",sn->line); }
 
 /*** Write MultiScale Probe information ***/
 fprintf(fp,"REMARK  MULSC_PROBE: NRADIUS %d\n",MS->Nradius);
 for (r=1;r<=MS->Nradius;++r){ 
    fprintf(fp,"REMARK  MULSC_PROBE: RADIUS %3d th %6.3f Ang %6.3f 1/Ang \n",
    r,MS->Rarray[r],1.0/MS->Rarray[r]); }
 
 fprintf(fp,"REMARK  MIN_NUM_OF_GRID_IN_CLUSTER %d\n",minNgridClus);
 /*** Set Grid Grid Cluster Properties ***/
 gcls = Gclus;
 while (gcls->next != NULL){
  gcls = gcls->next;
  Set_Grid_Cluster_Properties(gcls,M,MS);
 }

 
 /*** Write Grid Cluster Information ***/
 gcls = Gclus;
 Ngrid_total = 0;
 Volume_total = 0.0;
 while (gcls->next != NULL){
   gcls = gcls->next;
   if (gcls->Nmember >= minNgridClus)
   fprintf(fp,"REMARK  CLUSTER_PROPERTY %3d Ngrid %6d Vol(AAA) %7.1f Rinac(A) av %5.2f mi %5.2f invRvolume(AA) %7.1f \n",
     gcls->num,gcls->Nmember,gcls->Volume,gcls->aveRinacc,gcls->minRinacc,gcls->invRvolume);

   Ngrid_total  += gcls->Nmember;
   Volume_total += gcls->Volume;

 } /* gcls */


 fprintf(fp,"REMARK   Rinac(A) av     : harmonic average value of Rinaccess in the cluster (angstrom)\n");
 fprintf(fp,"REMARK   Rinac(A) mi     : minimum value of Rinaccess in the cluster (angstrom)\n");
 fprintf(fp,"REMARK   invRvolume(AA)  : sum of [1/Rinaccc] * Vol_of_one_voxel for the grids in the cluster (angstrom * angstrom)\n");
 fprintf(fp,"REMARK   Atoms with name ' CK ' correspond to the skeleton or core atoms.\n");
 fprintf(fp,"REMARK   Atoms with name ' CE ' correspond to the extended atoms from the skeleton or core.\n");
 fprintf(fp,"REMARK   Ngrid_total %d Volume_total(AAA) %lf\n",Ngrid_total, Volume_total);
 fprintf(fp,"TER\n");

 /*** Write Grid Position for Each Cluster ***/                                                                                                          
 Natom = Nclus = 0;
 gcls = Gclus;
 while (gcls->next != NULL){
  gcls = gcls->next;
  ++Nclus;
  if ((gcls->Nmember >= minNgridClus) && ((Nclus<=Nrank_keep)||(Nrank_keep<0))){
  fprintf(fp,"REMARK  CLUSTER %3d\n",Nclus);
  fprintf(fp,"REMARK  CLUSTER_PROPERTY %3d Ngrid %6d Vol(AAA) %7.1f Rinac(A) av %5.2f mi %5.2f invRvolume(AA) %7.1f \n",
      gcls->num,gcls->Nmember,gcls->Volume,gcls->aveRinacc,gcls->minRinacc,gcls->invRvolume);
  for (r=0;r<256;++r)
   if (gcls->NRinacc[r]>0){
    fprintf(fp,"REMARK  CLUSTER_RINACC_HISTOGRAM %3d Rinacc %5.2f Ngrid %6d %6.2f %%\n",
     gcls->num,MS->Rarray[r],gcls->NRinacc[r],100.0*gcls->NRinacc[r]/gcls->Nmember);
    }  

  fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
         Natom%100000," CA ","CEN",chain_str[not_over_max(Nclus,len_chain_str)],
            255, gcls->Pos[0],gcls->Pos[1],gcls->Pos[2],M->grid_width/2.0,gcls->aveRinacc);
  Natom += 1;
  fprintf(fp,"REMARK  [ChainID:ClusNum],[ResNum:char value]       [Radius][Rinacc][x][y][z]\n");
 
  gmem = &(gcls->HeadGridMember);
  while (gmem->next != NULL){
    gmem = gmem->next;
    x = gmem->x;
    y = gmem->y;
    z = gmem->z;
    if ((x>=0)&&(x<M->N[0])&&(y>=0)&&(y<M->N[1])&&(z>=0)&&(z<M->N[2])){
      ++Natom;
      /*
      if (gmem->type == 'K') sprintf(AtomStr,"CK ");
                        else sprintf(AtomStr,"CE ");
      */
      sprintf(AtomStr,"CE ");

       Pnt[0] = M->OrigPos[0] + x*M->grid_width;
       Pnt[1] = M->OrigPos[1] + y*M->grid_width;
       Pnt[2] = M->OrigPos[2] + z*M->grid_width;
       R = MS->Rarray[M->map[x][y][z]];
       if (R>0.0) invR = 1.0/R; else invR = 0.0;

       fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f %d %d %d",
           Natom%100000,AtomStr,"GRD",chain_str[not_over_max(Nclus,len_chain_str)],M->map[x][y][z],
           Pnt[0],Pnt[1],Pnt[2],M->grid_width/2.0,R,x,y,z);
       fprintf(fp,"\n");

   }
  } /* gmem */
  fprintf(fp,"TER\n");
  }
 } /* gclus */ 
 fprintf(fp,"END\n");
 fclose(fp);
                                                                                                           
} /* end of Output_CHAR3DMAP_PDB_format_By_Grid_Cluster() */



void Set_Map_to_Zero_for_Low_Ranked_Grid_Cluster(M,Gclus,Nrank_keep)
 struct CHAR3DMAP *M;
 struct GRID_CLUSTER *Gclus; 
 int   Nrank_keep;
{
 int Nrank,Nerase;
 struct GRID_CLUSTER *gcls;
 struct GRID_MEMBER  *gmem;

 printf("#Set_Map_to_Zero_for_Low_Ranked_Grid_Cluster(Nrank_keep:%d)\n",Nrank_keep);
 Nrank = 0; 
 Nerase = 0;
 gcls = Gclus;
 while (gcls->next != NULL){
  gcls = gcls->next;
  ++Nrank;
  if (Nrank > Nrank_keep){
    gmem = &(gcls->HeadGridMember);
    while (gmem->next != NULL){
       gmem = gmem->next;
       M->map[gmem->x][gmem->y][gmem->z] = 0;
       ++Nerase;
     } 
   }
 }
 printf("#Nerase %d\n",Nerase);

} /* end of Set_Map_to_Zero_for_Low_Ranked_Grid_Cluster() */


int not_over_max(i, max)
  int i,max;
{
  if (i>=max){ return(max-1); }
  return(i);
}
