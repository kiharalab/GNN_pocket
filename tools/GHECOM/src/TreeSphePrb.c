/*

 <TreeSphePrb.c>

 Functions for Hierarchical Tree Clustering for Spherical Probes

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
#include "globalvar.h"
#include "pdbstruct.h"
#include "PdbIO.h"
#include "GenSphePrb.h"
#include "MscPockClus.h"
#include "Grid3D.h"
#include "GclsMergeSort.h"

/** struct TREENODE : this is for constructing tree structure */

struct TREENODE { 
            int    num;                /* Node Number (0..NNODE-1) */
            int    dnum;               /* Node Number for distance matrix. (0..NLEAF-1)*/
                                       /* For leaf nodes, num = dnum, whereas dnum of ancestor is dnum of child1 (ch1.dnum).*/
            float  height;             /* Height (for UPGMA). heights for leaves are 0.*/
            struct TREENODE *ch1,*ch2;     /* Child Pointer  (for TREE graph )*/
            struct TREENODE *par;          /* Parent Pointer (for TREE graph ) */
            float  len1,len2;          /* Length from parent to child */
            struct TREENODE *prev, *next;  /* Back and next Pointer (for LINEAR graph) */
            struct TREENODE *cprev,*cnext; /* Back and next Pointer (for CLUSTER ) */
            char   type;               /* type for nodes        :'L'eaf, 'A'ncestor  */
            int    nmem;               /* Number of leaves under thie node  */
            float  Pos[3];             /* 3D-coordinate of the data */
            float  weight;             /* weight for data */
            float  score;              /* score for sorting */
            struct ATOM *atom;         /* pointer to ATOM */
           };



/*** Functions (global) ***/
int Make_Weighted_Hierarichical_Cluster_Of_Atoms();
int Project_Sphe_Probe_Cluster_into_Grid();

/*** Functions (local) ***/
static void Renumbering_Gclus_by_Sum_of_Weights();
static int   Make_TREENODEs_from_ATOMs();
static void Calculate_Euclid_Dmap_From_3Dpoints();
static int  Weighted_Centroid_Hierarchical_Clustering();
static void Add_TREENODE_Llist();
static void Delete_TREENODE_Llist();
static void Copy_TREENODE();
static int  Number_TREENODE_Llist();
static void Make_Cluster();

/**** LOCAL VARIABLES ****/
/** Array-style data **/
static float **DCEN;   /* distance matrix [Nnode][Nnode] : (malloc later) : for centroid */
static float **DSLC;   /* distance matrix [Nnode][Nnode] : (malloc later) : for single linkage clustering */
static int NLEAF;      /* Number of LEAF */
static int NNODE;      /* Number of NODE (in Llist or Tree ) */

static char  TREE_TYPE;     /* 'C'entroid, 'H'ybrid with single_linkage */
static float SIG_THRE_SLC;  /* sigmoid threshold       for single linkage clustereing */
static float SIG_ALPH_SLC;  /* sigmoid alpha parameter for single linkage clustereing */
static float SIG_BETA_SLC;  /* sigmoid beta  parameter for single linkage clustereing */


int Make_Weighted_Hierarichical_Cluster_Of_Atoms(HeadGclus,HeadAtomNode,Dthre,tree_type)
  struct GRID_CLUSTER *HeadGclus;
  struct ATOM         *HeadAtomNode;
  float  Dthre;       /* Threshold for distance */
  char   tree_type;   /* 'C'entroid, 'H'ybrid with single_linkage */
{
  struct TREENODE HeadLinearTreeNode,RootTreeNode;
  struct GRID_CLUSTER *curr_gclus;
  int Ncluster;

  TREE_TYPE = tree_type;

  HeadGclus->next = NULL;
  Ncluster = 0;
  
  Make_TREENODEs_from_ATOMs(&HeadLinearTreeNode,HeadAtomNode);
  printf("#NLEAF %d\n",NLEAF);
  if (NLEAF==0){return(0);}

  SIG_THRE_SLC = HeadAtomNode->next->R * 4.0;
  printf("#SIG_THRE_SLC %f\n",SIG_THRE_SLC);
  SIG_ALPH_SLC = 1.0;
  SIG_BETA_SLC = 1.0;

  Calculate_Euclid_Dmap_From_3Dpoints(&HeadLinearTreeNode);

  Weighted_Centroid_Hierarchical_Clustering(&HeadLinearTreeNode,&RootTreeNode);
  curr_gclus = NULL;
  Make_Cluster(&RootTreeNode,Dthre,0,&Ncluster,curr_gclus,HeadGclus);
  printf("#Ncluster %d\n",Ncluster);
       
  /** renumbering gclus by the sum of weights **/
  Merge_Sort_Double_Linked_List_GRID_CLUSTER(HeadGclus,'R');
  Renumbering_Gclus_by_Sum_of_Weights(Ncluster,HeadGclus,HeadAtomNode);
  return(Ncluster); 

} /* end of Make_Weighted_Hierarichical_Cluster_Of_Atoms() */




int Project_Sphe_Probe_Cluster_into_Grid(HeadGclus,HeadAtomNode,Xmap,Zmap)
 struct GRID_CLUSTER *HeadGclus;
 struct ATOM         *HeadAtomNode;
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
 struct CHAR3DMAP *Zmap;  /* map for clustering (should be initialized to 0) */
{
 int x,y,z,i,Ncluster,min_cluster_num;
 struct ATOM *an;
 struct GRID_CLUSTER *gcls, **gcls_pointer;
 float min,max,DD,DifPos[3],RR,minDD,threDD;
 int Nmin[3],Nmax[3];
 
 printf("#Project_Sphe_Probe_Cluster_into_Grid(HeadGclus,HeadAtomNode)\n");

 if (HeadAtomNode->next == NULL) return(0);
 if (HeadGclus->next    == NULL) return(0);

 /** [1] Initializing Zmap **/
 for (x=0;x<Zmap->N[0];++x){
   for (y=0;y<Zmap->N[1];++y){
     for (z=0;z<Zmap->N[2];++z){
       Zmap->map[x][y][z] = 0;
       /* printf("init %d %d %d Xmap %d Zmap %d\n",x,y,z,Xmap->map[x][y][z],Zmap->map[x][y][z]); */
     }
   }
 }    

 /*** [2] Count Ncluster and make gcls_pointer[] **/
 Ncluster = 0;
 gcls = HeadGclus;
 while (gcls->next != NULL){
   gcls = gcls->next;
   Ncluster += 1;
 }
 gcls_pointer = (struct GRID_CLUSTER**)malloc(sizeof(struct GRID_CLUSTER*)*(Ncluster+1));

 gcls = HeadGclus;
 while (gcls->next != NULL){
   gcls = gcls->next;
   gcls_pointer[gcls->num] = gcls; 
   /* printf("gcls->num %d Ncluster %d\n",gcls->num,Ncluster);  */
 }

 /** [2] Project Sphe Probe Cluster into Zmzp **/
 an = HeadAtomNode;
 while (an->next != NULL){
   an = an->next; 
   RR = an->R * an->R;
   for (i=0;i<3;++i){
     min = an->Pos[i] - an->R;
     max = an->Pos[i] + an->R;
     Nmin[i] = (int)floor((min - Zmap->OrigPos[i])/Zmap->grid_width);
     Nmax[i] = (int)ceil( (max - Zmap->OrigPos[i])/Zmap->grid_width);
     if (Nmin[i]<0)        Nmin[i] = 0;
     if (Nmax[i]>=Zmap->N[i]) Nmax[i] = Zmap->N[i]-1;
   }
 
   /* printf("Nmin %d %d %d Nmax %d %d %d\n",Nmin[0],Nmin[1],Nmin[2],Nmax[0],Nmax[1],Nmax[2]); */

   for (x=Nmin[0];x<=Nmax[0];++x){
     for (y=Nmin[1];y<=Nmax[1];++y){
       for (z=Nmin[2];z<=Nmax[2];++z){
         DifPos[0] = Zmap->OrigPos[0] + x*Zmap->grid_width - an->Pos[0];
         DifPos[1] = Zmap->OrigPos[1] + y*Zmap->grid_width - an->Pos[1];
         DifPos[2] = Zmap->OrigPos[2] + z*Zmap->grid_width - an->Pos[2];
         DD = DifPos[0]*DifPos[0] + DifPos[1]*DifPos[1] + DifPos[2]*DifPos[2];
         if ((Xmap->map[x][y][z]!=0)&&(Xmap->map[x][y][z]!=255)&&(DD <= RR)){
           /* if overlap with more than two probes, set Zmap==255 */
            if ((Zmap->map[x][y][z]!=0)&&(Zmap->map[x][y][z]!=an->cluster_num)) Zmap->map[x][y][z] = 255;  
                                                                          else  Zmap->map[x][y][z] = an->cluster_num;
         }
      } 
     } 
   }
 
 }

 for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if ((Zmap->map[x][y][z]>0)&&(Zmap->map[x][y][z]<255)){
           Add_Grid_Member(&(gcls_pointer[Zmap->map[x][y][z]]->HeadGridMember),x,y,z,Zmap->map[x][y][z],'-');
        }
        if (Zmap->map[x][y][z]==255) Zmap->map[x][y][z] = 0;
     }
   }
 }

 /** [3] Find closest sphe-probe for the residual pocket grids **/
 for (x=0;x<Xmap->N[0];++x){
   for (y=0;y<Xmap->N[1];++y){
     for (z=0;z<Xmap->N[2];++z){
       if ((Zmap->map[x][y][z]==0)&&(Xmap->map[x][y][z]!=0)&&(Xmap->map[x][y][z]!=255)){
         minDD = -1.0; min_cluster_num = 0;
         an = HeadAtomNode;
         while (an->next != NULL){
           an = an->next;
           DifPos[0] = Xmap->OrigPos[0] + x*Xmap->grid_width - an->Pos[0];
           DifPos[1] = Xmap->OrigPos[1] + y*Xmap->grid_width - an->Pos[1];
           DifPos[2] = Xmap->OrigPos[2] + z*Xmap->grid_width - an->Pos[2];
           DD = DifPos[0]*DifPos[0] + DifPos[1]*DifPos[1] + DifPos[2]*DifPos[2];
           threDD = (an->R*2.0)*(an->R*2.0); 
           if ((DD < threDD) && ((minDD<0.0) || ((DD<minDD)))) {minDD = DD; min_cluster_num = an->cluster_num;}
         }   
         if (minDD>0.0){
           Zmap->map[x][y][z] = min_cluster_num; 
           /* printf("add %d minD %f (%d %d %d)\n",min_cluster_num,sqrt(minDD),x,y,z); */
           Add_Grid_Member(&(gcls_pointer[min_cluster_num]->HeadGridMember),x,y,z,Xmap->map[x][y][z],'-');
         }
       }
     }
   }
 }

 free(gcls_pointer);
 return(Ncluster);

} /* end of Project_Sphe_Probe_Cluster_into_Grid(HeadGclus,HeadAtomNode,Xmap) */



void Renumbering_Gclus_by_Sum_of_Weights(Ncluster,HeadGclus,HeadAtomNode)
  int   Ncluster;
  struct GRID_CLUSTER *HeadGclus;
  struct ATOM *HeadAtomNode; 
{
  int    n,k,*new_from_old;  
  struct GRID_CLUSTER *gn;
  struct ATOM *an;
  
  new_from_old = (int *)malloc(sizeof(int)*(Ncluster+1));
  gn = HeadGclus;
  n = 0;
  while (gn->next != NULL){
    gn = gn->next;
    n += 1;
    new_from_old[gn->num] = n; 
    gn->num = n; 
    for (k=0;k<3;++k) gn->Pos[k] /= gn->value;  /* calculate gcls->Pos */
  } 

  an = HeadAtomNode;
  while (an->next != NULL){
    an = an->next;
    an->cluster_num = new_from_old[an->cluster_num];
  }

  free(new_from_old);

} /* end of Renumbering_Gclus_by_Sum_of_Weights() */




int Make_TREENODEs_from_ATOMs(HeadLinearTreeNode,HeadAtomNode)
  struct TREENODE *HeadLinearTreeNode;
  struct ATOM     *HeadAtomNode;
{
 struct ATOM *an;
 struct TREENODE *newtn;

 NLEAF = 0;
 an = HeadAtomNode;
 HeadLinearTreeNode->next = NULL;
 while (an->next != NULL){
   an = an->next;
  
   newtn = (struct TREENODE*)malloc(sizeof(struct TREENODE)); 
  
   newtn->atom = an;
   newtn->Pos[0] = an->Pos[0];
   newtn->Pos[1] = an->Pos[1];
   newtn->Pos[2] = an->Pos[2];
   newtn->num = newtn->dnum = NLEAF;
   newtn->ch1 = newtn->ch2 = newtn->par = NULL;
   newtn->weight = 1.0/an->Rinacc/an->Rinacc;
   newtn->type = 'L';
   newtn->height = 0.0;  
   /* printf("weight %f Rinacc %f\n",newtn->weight,an->Rinacc); */
   if (HeadLinearTreeNode->next != NULL) HeadLinearTreeNode->next->prev = newtn;
   newtn->next = HeadLinearTreeNode->next;
   newtn->prev = HeadLinearTreeNode;

   HeadLinearTreeNode->next = newtn; 
   NLEAF += 1;
 }
  printf("#NLEAF %d\n",NLEAF);
  NNODE = NLEAF;
  return(NLEAF);
} /* end of Make_TREENODEs_from_ATOMs() */






void Calculate_Euclid_Dmap_From_3Dpoints(HeadLinearTreeNode)
  struct TREENODE *HeadLinearTreeNode;
{
 int i;
 struct TREENODE *an,*bn;
 float dx,dy,dz,memoryGB;

 memoryGB = (float)sizeof(float)*(float)(NLEAF * NLEAF)/1024.0/1024.0/1024.0;
 printf("#Calculate_Euclid_Dmap_From_3Dpoints(NLEAF:%d, %.3f GB required)\n",NLEAF,memoryGB); fflush(stdout);

 DCEN = (float **)malloc(sizeof(float *)*NLEAF);
 for (i=0;i<NLEAF;++i)  DCEN[i] = (float *)malloc(sizeof(float )*NLEAF);
 if (TREE_TYPE=='H'){
   DSLC = (float **)malloc(sizeof(float *)*NLEAF);
   for (i=0;i<NLEAF;++i)  DSLC[i] = (float *)malloc(sizeof(float )*NLEAF);
 } 

 an = HeadLinearTreeNode;
 while (an->next != NULL){
   an = an->next;
   DCEN[an->dnum][an->dnum] = 0.0;
   bn = an;
   while ((bn!=NULL) && (bn->next != NULL)){
     bn = bn->next;
     dx = an->Pos[0] - bn->Pos[0];
     dy = an->Pos[1] - bn->Pos[1];
     dz = an->Pos[2] - bn->Pos[2];
     DCEN[an->dnum][bn->dnum] = sqrt(dx*dx + dy*dy + dz*dz);
     
     if (TREE_TYPE=='H'){ 
       DSLC[an->dnum][bn->dnum]  = DCEN[an->dnum][bn->dnum];
       DSLC[bn->dnum][an->dnum]  = DSLC[an->dnum][bn->dnum];
       DCEN[an->dnum][bn->dnum]  =  
          DCEN[an->dnum][bn->dnum] * (1.0 + SIG_BETA_SLC/(1+exp(-SIG_ALPH_SLC*(DSLC[an->dnum][bn->dnum]-SIG_THRE_SLC))));
     } 

     DCEN[bn->dnum][an->dnum]  =  DCEN[an->dnum][bn->dnum];
  }
 }
} /* end of Calculate_Euclid_Dmap_From_3Dpoints() */




int Weighted_Centroid_Hierarchical_Clustering(lstanode,trootnode)
 struct TREENODE *lstanode,*trootnode;
{
 int end,NLnode,k;
 float minDij;
 struct TREENODE *I,*J,*minI,*minJ,*newN,*M;
  
 printf("#Weighted_Hierarchical_Clustering(lstanode,trootnode)\n");
 NLnode  = Number_TREENODE_Llist(lstanode);
 I = NULL; J = NULL;
 end = 0;
 while (end == 0){

  /** (1) Finding the closest node pairs (minI, minJ) **/
   minI = NULL; minJ = NULL;
   minDij = 0.0;
   /** Find Mimimum Dij **/
   I = lstanode;
   while (I->next !=NULL){
     I = I->next;
     J = I;
     while (J->next != NULL) {
       J = J->next;
       if ( ((minI==NULL)&&(minJ==NULL))||
            (DCEN[I->dnum][J->dnum] < minDij)) {minDij = DCEN[I->dnum][J->dnum]; minI = I; minJ = J;}
      }

    }

  /*
  printf("min I %d J %d minDij %f\n",minI->dnum,minJ->dnum,minDij);
  if (minI==NULL) {printf("#woops!! minI is NULL (NLnode %d)\n",NLnode); exit(1);}
  if (minJ==NULL) {printf("#woops!! minJ is NULL (NLnode %d)\n",NLnode); exit(1);}
  */

 
  /** (2) Making New Node **/
   newN = (struct TREENODE *) malloc(sizeof(struct TREENODE));
   newN->num = NNODE;
   newN->dnum = minI->dnum;  /* newN uses 'dnum' of minI */
   newN->type = 'A';
   newN->nmem   = minI->nmem   + minJ->nmem;
   newN->weight = minI->weight + minJ->weight;

   if (newN->weight == 0){
     printf("woops! newN %d minI %d minJ %d weight newN %f minI %f minJ %f\n",
      newN->num,minI->num,minJ->num,newN->weight,minI->weight,minJ->weight);
     exit(1);
   }

   for (k=0;k<3;++k) newN->Pos[k] = (minI->weight * minI->Pos[k]  + minJ->weight * minJ->Pos[k])/(float)(newN->weight);
   newN->ch1 = minI;
   newN->ch2 = minJ;
   minI->par = newN;
   minJ->par = newN;
  
   newN->height = DCEN[minI->dnum][minJ->dnum];

   newN->len1 = newN->height - minI->height;
   newN->len2 = newN->height - minJ->height;

   ++NNODE;

   /** (3) Set new DCEN for Linear List **/
   M = lstanode;
   while (M->next !=NULL){
     M = M->next;
    /* printf("I %d J %d NLnode %d newN %d M %d /NNODE %d NLEAF %d nmax_node %d minDij %f\n", 
      minI->num,minJ->num,NLnode, newN->num,M->num,NNODE,NLEAF, 2*NLEAF+1,minDij); fflush(stdout); */
     /*** Calculate DCEN ***/
     DCEN[newN->dnum][M->dnum]  =  0.0;
     for (k=0;k<3;++k)
       DCEN[newN->dnum][M->dnum]  +=  (newN->Pos[k] - M->Pos[k]) * (newN->Pos[k] - M->Pos[k]);
     DCEN[newN->dnum][M->dnum]  =  sqrt(DCEN[newN->dnum][M->dnum]);

     if (TREE_TYPE=='H'){
       /*** Calculate DSLC ***/
       if (DSLC[I->dnum][M->dnum]<DSLC[I->dnum][M->dnum]) DSLC[newN->dnum][M->dnum] = DSLC[I->dnum][M->dnum];
                                                     else DSLC[newN->dnum][M->dnum] = DSLC[J->dnum][M->dnum];
       DSLC[M->dnum][newN->dnum] = DSLC[newN->dnum][M->dnum];
     
       /*** Combine DCEN and DSLC, and store as DCEN ***/
       DCEN[newN->dnum][M->dnum]  =  
          DCEN[newN->dnum][M->dnum] * (1.0 + SIG_BETA_SLC/(1+exp(-SIG_ALPH_SLC*(DSLC[newN->dnum][M->dnum]-SIG_THRE_SLC))));
     } 

     DCEN[M->dnum][newN->dnum] = DCEN[newN->dnum][M->dnum];
   } 

   /** (4) Add newN and delete minI and minJ **/
   Add_TREENODE_Llist(lstanode,newN);
   Delete_TREENODE_Llist(minI);
   Delete_TREENODE_Llist(minJ);

   NLnode  = Number_TREENODE_Llist(lstanode);

   /** Terminal Operation **/
   if (NLnode <=1 )
    { end = 1;
      I = lstanode->next;
      Copy_TREENODE(trootnode,I);
      trootnode->type = 'A';
     }

  } /* while (end==0) */

  return(1);
} /* end of Weighted_Centroid_Hierarchical_Clustering() */




void Add_TREENODE_Llist(stanode,newnode)
 struct TREENODE *stanode,*newnode;
{
 struct TREENODE *bn;

 bn = stanode;
 while (bn->next != NULL) bn = bn->next;
 bn->next = newnode;
 newnode->prev = bn;
 newnode->next = NULL;

} /* end of Add_TREENODE_Llist() */




void Delete_TREENODE_Llist(delnode)
 struct TREENODE *delnode;
{
 if (delnode->prev != NULL)
   delnode->prev->next = delnode->next;
 if (delnode->next != NULL)
   delnode->next->prev = delnode->prev;
} /* end of Delete_TREENODE_Llist() */




void Copy_TREENODE(A,B)
 struct TREENODE *A,*B;
{
 A->type  = B->type;
 A->num   = B->num;
 A->dnum  = B->dnum;
 A->par  = B->par;
 A->height = B->height;
 A->nmem = B->nmem;
 A->Pos[0] = B->Pos[0];
 A->Pos[1] = B->Pos[1];
 A->Pos[2] = B->Pos[2];
 if (A->type == 'A') {
    A->ch1 = B->ch1;
    A->ch2 = B->ch2;
    B->ch1->par = A;
    B->ch2->par = A;
    A->len1 = B->len1;
    A->len2 = B->len2; }
} /* end of Copy_TREENODE() */





int Number_TREENODE_Llist(stanode)
 struct TREENODE *stanode;
{
 struct TREENODE *bn;
 int n;
 n = 0;
 bn = stanode;
 while (bn->next != NULL) { bn = bn->next; ++n;}
 return(n);
} /* end of Number_TREENODE_Llist() */




void Make_Cluster(node,height_thre,out_cluster,Ncluster,curr_gclus,head_gclus)
 struct TREENODE *node;
 float height_thre;
 char  out_cluster;
 int    *Ncluster;    /* number of cluster */
 struct GRID_CLUSTER *curr_gclus;  /* pointer of current grid_cluster */
 struct GRID_CLUSTER *head_gclus;  /* pointer of head of llist of grid_clusters */
{
  char new_out_cluster;
  int k;

 /*
  printf("#Make_Cluster(cstanode,node %c,height %f height_thre %fout_cluster %d Ncluster %d)\n",
     node->type,node->height,height_thre,out_cluster,*Ncluster);
 */

  new_out_cluster = out_cluster;

  if ((out_cluster==0) && (node->height <= height_thre)){
     /** making a new cluster **/
     curr_gclus = (struct GRID_CLUSTER*)malloc(sizeof(struct GRID_CLUSTER));
     if (head_gclus->next != NULL) head_gclus->next->prev = curr_gclus;
     curr_gclus->prev = head_gclus;
     curr_gclus->next = head_gclus->next;
     head_gclus->next = curr_gclus;      
     curr_gclus->Nmember = 0; 
     curr_gclus->value   = 0.0; 
     curr_gclus->HeadGridMember.next = NULL;
     curr_gclus->Pos[0] = curr_gclus->Pos[1] = curr_gclus->Pos[2] =  0.0;
     (*Ncluster) += 1;
     curr_gclus->num   = *Ncluster;   /* gclus->num ranges [1..Ncluster] */
     new_out_cluster = 1;
  }

  if (node->type == 'A'){
    Make_Cluster(node->ch1,height_thre,new_out_cluster,Ncluster,curr_gclus,head_gclus);
    Make_Cluster(node->ch2,height_thre,new_out_cluster,Ncluster,curr_gclus,head_gclus);
   }

  if (node->type == 'L'){ 
      curr_gclus->Nmember += 1;
      curr_gclus->value += node->weight; 
      for (k=0;k<3;++k) curr_gclus->Pos[k] += node->weight * node->Pos[k];
      node->atom->cluster_num = curr_gclus->num;
 /*
      printf("curr_gclus->num %d node->num %d node->atom->cluster_num  %d\n",
      curr_gclus->num,node->num,node->atom->cluster_num); 
*/
     }

} /* end of Make_Cluster() */
