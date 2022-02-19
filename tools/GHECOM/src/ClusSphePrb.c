/*
 
 <ClusSphePrb.c>

 for clustering spherical probes using Single Linkage

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


/*** FUNCTIONS (GLOBAL) ***/

void Make_Single_Linkage_Cluster_Of_Atoms();
void Extract_Rep_Atoms_in_Cluster_Center();

/*** FUNCTIONS (LOCAl) ***/
void Mark_Neighbor();

void Make_Single_Linkage_Cluster_Of_Atoms(Ahead,Rhead,Dthre)
 struct ATOM    *Ahead;
 struct RESIDUE *Rhead;
 float  Dthre; 
{
 struct FLOAT2DMAP Dmat;
 struct ATOM    *an;
 struct RESIDUE *rn;
 int Natom,Ncluster;
printf("#Make_Single_Linkage_Cluster_Of_Atoms(Dthre %f)\n",Dthre); fflush(stdout);
 /*** [1] Inirialize rnum  and calculating distance map ***/
 Rhead->next = NULL;
 an = Ahead;
 while (an->next != NULL) { an = an->next; an->rnum = 0; an->res = NULL;}
 Renumber_Atom_Number(Ahead);
 Natom = Number_Of_ATOMs(Ahead);
 Malloc_FLOAT2DMAP(&Dmat,Natom);
 Cal_Distance_FLOAT2DMAP(Ahead,&Dmat);
 printf("#Natom %d\n",Natom);
 /*** [2] Make Single Linkage Clustering ***/
 Ncluster = 0;
 an = Ahead;  rn = Rhead;
 while (an->next != NULL){
  an = an->next;
  if (an->res == NULL){
     /* (1) Making new RESIDUE */
     ++Ncluster;
     rn->next = (struct RESIDUE *)malloc(sizeof(struct RESIDUE));
     rn->next->prev = rn;
     rn->next->next = NULL;
     rn = rn->next;
     an->rnum =  rn->num  = Ncluster;
     rn->Natom = 0;
     sprintf(rn->Resi,"PRB");
     sprintf(rn->Rnum,"%4d ",rn->num%10000);
     rn->Ahead.rnext = NULL;
     an->res = rn;
     /* printf(">Ncluster %d %d\n",Ncluster,rn->num); */
     Add_Atom_To_Residue_Tail(rn,an);

     /* (2) Marking Neighbor Atoms */
     Mark_Neighbor(Ahead,&Dmat,an,rn,Dthre);
   }

 } /* an */

 printf("#Ncluster %d\n",Ncluster);

 Free_FLOAT2DMAP(&Dmat);

} /* end of void Make_Single_Linkage_Cluster_Of_Atoms() */




void Extract_Rep_Atoms_in_Cluster_Center(Ahead,Rhead)
 struct ATOM    *Ahead;
 struct RESIDUE *Rhead;
{
 struct ATOM *an;
 struct RESIDUE *rn;
 int k;
 float DD,minDD;
 struct ATOM *mn;

 printf("#Extract_Rep_Atoms_in_Cluster_Center(Ahead,Rhead)\n");

 /* (1) Calculate center of cluster */
 rn = Rhead;
 while (rn->next != NULL){
   rn = rn->next;
   for (k=0;k<3;++k) rn->Pos[k] = 0.0;
   an = &(rn->Ahead);
   while (an->rnext != NULL){
    an = an->rnext;
    an->mark = 0;
    /* printf("atom [%d]  %f %f %f\n",an->num,an->Pos[0],an->Pos[1],an->Pos[2]);  */
    for (k=0;k<3;++k) rn->Pos[k] += an->Pos[k]; 
   }
  for (k=0;k<3;++k) rn->Pos[k] /= rn->Natom;
}

 /* (2) Mark center-closest atom */
 rn = Rhead;
 mn = NULL;
 while (rn->next != NULL){
   rn = rn->next;
   minDD = -1.0;
   an = &(rn->Ahead);
   while (an->rnext != NULL){
    an = an->rnext;
    DD = 0.0;
    for (k=0;k<3;++k) DD += (rn->Pos[k] - an->Pos[k]) * (rn->Pos[k] - an->Pos[k]);
    if ((minDD < 0.0) || (DD < minDD)) {minDD = DD; mn = an;}
   }

  /* printf("[%d] Natom %d %f %f %f minDD %f\n",rn->num,rn->Natom,rn->Pos[0],rn->Pos[1],rn->Pos[2],minDD); */
  mn->mark = 1; 
}


 /* (3) Remove non-representative (mark==0) atom  */
 an = Ahead;
 while (an->next != NULL){
   an = an->next; 
   if (an->mark==0){
    if (an->prev != NULL) an->prev->next = an->next;
    if (an->next != NULL) an->next->prev = an->prev;
   }
 }

} /* end of Extract_Rep_Atoms_in_Cluster_Center() */



void Mark_Neighbor(Ahead,Dmat,an,rn,Dthre)
 struct ATOM    *Ahead;
 struct FLOAT2DMAP  *Dmat;
 struct ATOM    *an;
 struct RESIDUE *rn;
 float  Dthre;
{
 struct ATOM *bn;

/* printf("Mark_Neighbor(an %d,rn %d)\n",an->num,rn->num); */

 bn = Ahead;
 while (bn->next != NULL){
   bn = bn->next;
   if ((bn->res==NULL)&&(Dmat->m[an->num][bn->num]<= Dthre)){ 
      Add_Atom_To_Residue_Tail(rn,bn);
      Mark_Neighbor(Ahead,Dmat,bn,rn,Dthre); }
 }

} /* end of Mark_Neighbor() */

 
