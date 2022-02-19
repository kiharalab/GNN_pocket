/*

 <MscLabelRes.c>

 for labeling MolVol and Pocket values into focused atoms (residues). 


=========== "ghecom" program ===========

Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under
the GNU Lesser General Public License (LGPL) version 3, see LICENSE.txt.

=========== Installation and Usage ===========

See html file under the "doc/" directory.

=========== Citing "ghecom" program ===========

Please cite:

1) Kawabata T. (2010) Detection of multi-scale pockets on protein surfaces using mathematical morphology. Proteins,78, 1195-1121.


-------------------------------------------------------------------------------
 Surrounding shells are placed around the focused atoms as follows:

   X : points of the focused atoms
   P : points of the protein 
   s : points of the shell 
   p : points of the protein in the shell
   e : empty points in the shell

              s         e 
    X        sXs       eXe
   XXX  ::  sXXXs  :: eXXXe  
    XPP      sXs       eXp 
   PPPP       s         p
   PPPP      
   
   NPinside        : [# of 'X'] Npoints inside of the focused atoms
   NPshell         : [# of 's'] Npoints outside atom_sphere and inside shell_sphere
   NPshell_prot    : [# of 'p'] Npoints of NPshell overlap target protein.
   NPshell_access  : [# of 'e'] Npoints of NPshell not overlap target protein (empty points)


   *ShellAcc (both for protein and ligand atoms) is defined as follows :
      ShellAcc = [#e]/[#s] = NPshell_access/NPshell
   

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


/*** FUNCTIONS (GLOBAL) ***/
void Label_Multiscale_MolVol_or_Pocket_to_Protein_Residue_Shells();
void Label_Multiscale_Each_Pocket_Cluster_to_Protein_Residue_Shells();
void Write_Residue_File_Assigned_Multiscale_MolVol_and_Pocket();



void Label_Multiscale_MolVol_or_Pocket_to_Protein_Residue_Shells(ProResHead,X,MS,Wshell,MolPockType)
 struct RESIDUE *ProResHead;
 /* Rinaccess will be written in tFactor of each atom. */
 struct CHAR3DMAP *X; 
   /* >> X << 
    Initial state  : multiscale closing 3D map. Voxel has minimum radius number.
                   0 -> outer space
                   1 -> sphere with 1-th radius.
                   2 -> sphere with 2-th radius.
                   : :
             Nradius -> sphere with Nradius-th radius.
                 255 -> VdW volume
   */
 struct MULSC_SHELL_PRB *MS;
 float  Wshell;  /* Width of the shell (be added to protein atoms (shell radius)) */
 char   MolPockType;  /* 'M'olVol, 'P'ocket */
{
 struct ATOM *an;
 struct RESIDUE *rn;
 float D[3],DD[3],ddis,RR,RRshell;
 int   minN[3],maxN[3],minNres[3],maxNres[3];
 int   i,x,y,z,k;
 int   NPshell,NPshell_access,NPinside;
 float invRsum,Rsum;
 int   Ncount[256];
 float Rmax, Rmax_sigma, invRmax_sigma,Rmin;
 struct CHAR3DMAP M;
 
 printf("#Label_Multiscale_MolVol_or_Pocket_to_Protein_Residue_Shells(MolPockType %c Wshell %f)\n",MolPockType,Wshell);
 Rmin = MS->Rarray[1];
 Rmax = MS->Rarray[MS->Nradius];
 Rmax_sigma = MS->Rarray[MS->Nradius] + (MS->Rarray[MS->Nradius] - MS->Rarray[MS->Nradius-1]);
 invRmax_sigma = 1.0/Rmax_sigma;
 printf("#Rmax %f Rmax_sigma %f\n",Rmax,Rmax_sigma);
 Malloc_CHAR3DMAP(&M,X->N[0],X->N[1],X->N[2]);
 
 rn = ProResHead;
 while (rn->next != NULL){
  rn = rn->next; 
  for (k=0;k<256;++k) Ncount[k] = 0;
  NPshell = NPshell_access = 0;
  NPinside = 0;
 
  for (i=0;i<3;++i) {maxNres[i] = 0; minNres[i] = X->N[i];} 
   
  an = &(rn->Ahead);

  while (an->rnext != NULL){
   an = an->rnext;
   an->tFactor = 0.0;
   NPshell = NPshell_access = 0; 

 /** (1) Make MAP 'M' 
        mask 1 for inside atom_sphere;
        mask 2 for inside shell_sphere;
     -->
       M[x][y][z] == 0     : void (outer space) 
       M[x][y][z] == 1 or 3: inside_residue 
       M[x][y][z] == 2     : accessible_shell 
  **/


   for (i=0;i<3;++i)
   {  minN[i] = (int)floor((an->Pos[i] - an->R - Wshell - X->OrigPos[i])/X->grid_width) - 1;
      if (minN[i]<0) minN[i] = 0;
      maxN[i] = (int)ceil ((an->Pos[i] + an->R + Wshell - X->OrigPos[i])/X->grid_width) + 1;
      if (maxN[i]>=X->N[i]) maxN[i] = X->N[i]-1;
      
      if (minN[i]<minNres[i]) minNres[i] = minN[i]; 
      if (maxN[i]>maxNres[i]) maxNres[i] = maxN[i]; 
   } /* i */

   RR       = an->R*an->R; 
   RRshell  = (an->R+Wshell)*(an->R+Wshell);
   for (x=minN[0];x<=maxN[0];++x){ 
    D[0] = X->OrigPos[0] + X->grid_width * x - an->Pos[0];
    DD[0] = D[0]*D[0];
    for (y=minN[1];y<=maxN[1];++y){
     D[1] = X->OrigPos[1] + X->grid_width * y - an->Pos[1];
     DD[1] = D[1]*D[1];
     for (z=minN[2];z<=maxN[2];++z){
        D[2] = X->OrigPos[2] + X->grid_width * z - an->Pos[2];
        ddis = DD[0] + DD[1] + D[2]*D[2];
             if (ddis<=RR)      M.map[x][y][z] = M.map[x][y][z]|1;
        else if (ddis<=RRshell) M.map[x][y][z] = M.map[x][y][z]|2;
     } /* z */
    } /* y */
  } /* x */

 } /* an */

 /** (2) Calculate averaged 1/R for "accessible shell" (M[x][y][z]==1) **/
  for (x=minNres[0];x<=maxNres[0];++x){ 
   for (y=minNres[1];y<=maxNres[1];++y){ 
    for (z=minNres[2];z<=maxNres[2];++z){ 
      if ((M.map[x][y][z]&1)==1) ++NPinside; 
      if (((M.map[x][y][z]&1)!=1)&&((M.map[x][y][z]&2)==2)){
         ++NPshell;
         Ncount[X->map[x][y][z]] += 1;  
         if (X->map[x][y][z]!=255) ++NPshell_access;
         }
     } /* z */
   }  /* y */
  } /* x */
 


    if (MolPockType == 'M'){
     rn->shellAcc =  rn->Rinacc =  0.0;

     if (NPshell>0){
       rn->shellAcc = (float)NPshell_access/(float)NPshell;
       Rsum = 0.0;
       for (k=1;k<255;++k)  {if (Ncount[k]>0)    Rsum += (Ncount[k]*MS->Rarray[k]);}
       Rsum    += (Ncount[0]*Rmax_sigma);
       rn->Rinacc = Rsum/NPshell;
       rn->rw_shellAcc = Rsum/(NPshell*Rmax_sigma);
     }

   }

   if (MolPockType == 'P'){
    rn->pocketness = 0.0;
    invRsum = 0.0;
    for (k=1;k<255;++k)  {if (Ncount[k]>0) invRsum += (Ncount[k]*1.0/MS->Rarray[k]);}
    if (NPshell>0){
     rn->pocketness = invRsum * MS->Rarray[1]/NPshell;
    }
   }
 

 /** (3) Re-initialize map M  **/
  for (x=minNres[0];x<=maxNres[0];++x){ 
   for (y=minNres[1];y<=maxNres[1];++y){ 
    for (z=minNres[2];z<=maxNres[2];++z){ M.map[x][y][z] = 0; } 
   }  /* y */
  } /* x */

} /* rn */

 Free_CHAR3DMAP(&M);
} /* end of Label_Multiscale_MolVol_or_Pocket_to_Protein_Residue_Shells() */







void Label_Multiscale_Each_Pocket_Cluster_to_Protein_Residue_Shells(ProResHead,X,Z,M,MS,Wshell)
 struct RESIDUE *ProResHead;
 /* Rinaccess will be written in tFactor of each atom. */
 struct CHAR3DMAP *X; 
   /* >> X << 
    Initial state  : multiscale closing 3D map. Voxel has minimum radius number.
                   0 -> outer space
                   1 -> sphere with 1-th radius.
                   2 -> sphere with 2-th radius.
                   : :
             Nradius -> sphere with Nradius-th radius.
                 255 -> VdW volume
   */
 struct CHAR3DMAP *Z;  /* Map containing the cluster number */
 struct CHAR3DMAP *M;  /* Map tempolariry used for labeling  */
 struct MULSC_SHELL_PRB *MS;
 float  Wshell;  /* Width of the shell (be added to protein atoms (shell radius)) */
{
 struct ATOM *an;
 struct RESIDUE *rn;
 float D[3],DD[3],ddis,RR,RRshell;
 int   minN[3],maxN[3],minNres[3],maxNres[3];
 int   i,x,y,z,k,c;
 int   NPshell,NPshell_access,NPinside;
 float invRsum;
 int   Ncount[MAX_CLUSTER_NUM][256];
 float Rmax, Rmax_sigma, invRmax_sigma,Rmin;
 
 printf("#Label_Multiscale_MolVol_or_Pocket_to_Protein_Residue_Shells(Wshell %f Nmax_cluster %d)\n",Wshell,PAR.Nmax_cluster);
 Rmin = MS->Rarray[1];
 Rmax = MS->Rarray[MS->Nradius];
 Rmax_sigma = MS->Rarray[MS->Nradius] + (MS->Rarray[MS->Nradius] - MS->Rarray[MS->Nradius-1]);
 invRmax_sigma = 1.0/Rmax_sigma;
 /*
 Malloc_CHAR3DMAP(&M,X->N[0],X->N[1],X->N[2]);
 */
 Initialize_CHAR3DMAP(M);
 
 rn = ProResHead;

 while (rn->next != NULL){
  rn = rn->next; 
  for (c=0;c<=PAR.Nmax_cluster;++c){
    for (k=0;k<256;++k) Ncount[c][k] = 0;
  }

  NPshell = NPshell_access = 0;
  NPinside = 0;
 
  for (i=0;i<3;++i) {maxNres[i] = 0; minNres[i] = X->N[i];} 
   
  an = &(rn->Ahead);

  while (an->rnext != NULL){
   an = an->rnext;
   an->tFactor = 0.0;
   NPshell = NPshell_access = 0; 

 /** (1) Make MAP 'M' 
        mask 1 for inside atom_sphere;
        mask 2 for inside shell_sphere;
     -->
       M[x][y][z] == 0     : void (outer space) 
       M[x][y][z] == 1 or 3: inside_residue 
       M[x][y][z] == 2     : accessible_shell 
  **/


   for (i=0;i<3;++i)
   {  minN[i] = (int)floor((an->Pos[i] - an->R - Wshell - X->OrigPos[i])/X->grid_width) - 1;
      if (minN[i]<0) minN[i] = 0;
      maxN[i] = (int)ceil ((an->Pos[i] + an->R + Wshell - X->OrigPos[i])/X->grid_width) + 1;
      if (maxN[i]>=X->N[i]) maxN[i] = X->N[i]-1;
      
      if (minN[i]<minNres[i]) minNres[i] = minN[i]; 
      if (maxN[i]>maxNres[i]) maxNres[i] = maxN[i]; 
   } /* i */

   RR       = an->R*an->R; 
   RRshell  = (an->R+Wshell)*(an->R+Wshell);
   for (x=minN[0];x<=maxN[0];++x){ 
    D[0] = X->OrigPos[0] + X->grid_width * x - an->Pos[0];
    DD[0] = D[0]*D[0];
    for (y=minN[1];y<=maxN[1];++y){
     D[1] = X->OrigPos[1] + X->grid_width * y - an->Pos[1];
     DD[1] = D[1]*D[1];
     for (z=minN[2];z<=maxN[2];++z){
        D[2] = X->OrigPos[2] + X->grid_width * z - an->Pos[2];
        ddis = DD[0] + DD[1] + D[2]*D[2];
             if (ddis<=RR)      M->map[x][y][z] = M->map[x][y][z]|1;
        else if (ddis<=RRshell) M->map[x][y][z] = M->map[x][y][z]|2;
     } /* z */
    } /* y */
  } /* x */

 } /* an */

 /** (2) Calculate averaged 1/R for "accessible shell" (M[x][y][z]==1) **/
  for (x=minNres[0];x<=maxNres[0];++x){ 
   for (y=minNres[1];y<=maxNres[1];++y){ 
    for (z=minNres[2];z<=maxNres[2];++z){ 
      if ((M->map[x][y][z]&1)==1) ++NPinside; 
      if (((M->map[x][y][z]&1)!=1)&&((M->map[x][y][z]&2)==2)){
         ++NPshell;
         /* if (Z->map[x][y][z]>0) printf("Z %d\n",Z->map[x][y][z]); */
         /* if (Z->map[x][y][z]==ClusterNum) Ncount[X->map[x][y][z]] += 1; */  
         if ((Z->map[x][y][z]>0)&&(Z->map[x][y][z]<=PAR.Nmax_cluster)) Ncount[Z->map[x][y][z]][X->map[x][y][z]] += 1;  

         if (X->map[x][y][z]!=255) ++NPshell_access;
         }
     } /* z */
   }  /* y */
  } /* x */
 

    for (c=1;c<=PAR.Nmax_cluster;++c){
      rn->pocketness_clus[c] = 0.0;
      invRsum = 0.0;
      for (k=1;k<255;++k)  {if (Ncount[c][k]>0) invRsum += (Ncount[c][k]*1.0/MS->Rarray[k]);}
      if (NPshell>0){
       rn->pocketness_clus[c] = invRsum * MS->Rarray[1]/NPshell;
      }
    } 

 /** (3) Re-initialize map M  **/
  for (x=minNres[0];x<=maxNres[0];++x){ 
   for (y=minNres[1];y<=maxNres[1];++y){ 
    for (z=minNres[2];z<=maxNres[2];++z){ M->map[x][y][z] = 0; } 
   }  /* y */
  } /* x */

} /* rn */

 /*
 Free_CHAR3DMAP(&M);
 */
} /* end of Label_Multiscale_Each_Pocket_Cluster_to_Protein_Residue_Shells() */








void Write_Residue_File_Assigned_Multiscale_MolVol_and_Pocket(fname,Rhead,MS,Out_pocketness_clus,title,CommentHead)
 char *fname;
 struct RESIDUE *Rhead;
 struct MULSC_SHELL_PRB *MS;
 char   Out_pocketness_clus; /* 'T','S':output pocketness_clus, '-' or 'F': not output pocket_clus */
 char   *title;
 struct STRING_NODE *CommentHead;
{
 FILE *fp;
 struct RESIDUE *rn;
 int r,c;
 char core[128],buff[512],buff2[512],chain;
 struct STRING_NODE *sn;
 float Rmax,Rmax_sigma,invRmax_sigma;
                                                                                                                          
 printf("#Write_Residue_File_Assigned_Multiscale_Value() -> \"%s\"\n",fname);
 Rmax = MS->Rarray[MS->Nradius];
 Rmax_sigma = MS->Rarray[MS->Nradius] + (MS->Rarray[MS->Nradius] - MS->Rarray[MS->Nradius-1]);
 invRmax_sigma = 1.0/Rmax_sigma;
 
 if (fname[0]=='-') fp = stdout;
 else{
  fp = fopen(fname,"w");
  if (fp==NULL) { printf("#ERROR:Can't write to \"%s\"\n",fname); exit(1);}
 }
 
 Find_Filename_Core(core,fname);
 Get_Part_Of_Line(buff,fname,0,40);
 Get_Part_Of_Line(buff2,core,0,9);
 if (title[0]!='\0') fprintf(fp,"#TITLE   %s\n",title);
 fprintf(fp,"#OUTPUTFILENAME    %s\n",fname);
 fprintf(fp,"#DATE    %s\n",Get_Date_String());
 if (PAR.COMMAND[0] != '\0') fprintf(fp,"#COMMAND %s\n",PAR.COMMAND);
sn = CommentHead;
 while (sn->next != NULL){ 
   sn = sn->next;
   fprintf(fp,"#COMMENT %s\n",sn->line); }
 /*** Write MultiScale Probe information ***/
 fprintf(fp,"#MULSC_PROBE: NRADIUS %d\n",MS->Nradius);
 for (r=1;r<=MS->Nradius;++r){ 
    fprintf(fp,"#MULSC_PROBE: RADIUS %3d th %6.3f Ang %6.3f 1/Ang \n",
      r,MS->Rarray[r],1.0/MS->Rarray[r]); }
       fprintf(fp,"#OUTSIDE_OF_MAX_RADIUS      %6.3f Ang %6.3f 1/Ang \n",Rmax_sigma,invRmax_sigma);
   
 fprintf(fp,"#COLUMN  1|RNUM              |Residue Number\n");
 fprintf(fp,"#COLUMN  2|CHAIN             |Chain Identifier\n");
 fprintf(fp,"#COLUMN  3|RES               |Three-letter residue name\n");
 fprintf(fp,"#COLUMN  4|shellAcc          |shell accessibility (%%)\n");
 fprintf(fp,"#COLUMN  5|Rinacc            |averaged Rinaccess (A)\n");
 fprintf(fp,"#COLUMN  6|Natom             |Number of atoms\n");
 fprintf(fp,"#COLUMN  7|Natom_contact     |Number of contacting atoms with ligand\n");
 fprintf(fp,"#COLUMN  8|pocketness        |sum of 1/[Rpocket] /(1/[Rmin]*[vol of shell]) (%%)\n");
 if ((Out_pocketness_clus!='F') && (Out_pocketness_clus!='-')){ 
   for (c=1;c<=PAR.Nmax_cluster;++c) 
   fprintf(fp,"#COLUMN %2d|pocketness_clus %d |sum of 1/[Rpocket_for_pocketcluster %d] /(1/[Rmin]*[vol of shell]) (%%)\n",c+8,c,c);
 }

 /*** Write ATOM/HETATM ***/
 rn = Rhead;
 while (rn->next!=NULL){
   rn = rn->next;
   chain = rn->Chain;
   if (chain == ' ') chain = '-';
   fprintf(fp,"%s %c %s %6.2f %6.3f %2d %2d %6.2f",
   rn->Rnum,chain,rn->Resi,100.0*rn->shellAcc,rn->Rinacc,rn->Natom,rn->Natom_contact,100.0*rn->pocketness);
   if ((Out_pocketness_clus!='F') && (Out_pocketness_clus!='-')){ 
     for (c=1;c<=PAR.Nmax_cluster;++c) fprintf(fp," %6.2f",100.0*rn->pocketness_clus[c]);  
   }
   fprintf(fp,"\n");
 }
 
 if (fp!=stdout) fclose(fp);
                                                                                                                          
} /* end of Write_Residue_File_Assigned_Multiscale_MolVol_and_Pocket() */

