/*

 <MscLabelAtm.c>
 
 for labeling MolVol and Pocket values into focused atoms. 

=========== "ghecom" program ===========

Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under
the GNU Lesser General Public License (LGPL) version 3, see LICENSE.txt.

=========== Installation and Usage ===========

See html file under the "doc/" directory.

=========== Citing "ghecom" program ===========

Please cite:

1) Kawabata T. (2010) Detection of multi-scale pockets on protein surfaces using mathematical morphology. Proteins,78, 1195-1121.



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
  
   *Rinacess in the protein region ('X') is defined as 0.
 
   *Rinaccess for protein atoms is defined as : mean of Rinacess for the NPshell
     
   *Rinaccess for ligand atoms is defined as : mean of Rinacess in the focused atom points (X).

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
void Label_Multiscale_MolVol_or_Pocket_to_Ligand_Atoms();
void Label_Multiscale_Each_Pocket_Cluster_to_Ligand_Atoms();
void Label_Multiscale_MolVol_or_Pocket_to_Protein_Atom_Shells();
void Label_Multiscale_Each_Pocket_Cluster_to_Protein_Atom_Shells();
void Cal_Number_of_Contacting_Atoms_with_Ligand();
void Write_PDB_File_with_Multiscale_MolVol_and_Pocket();


void Label_Multiscale_MolVol_or_Pocket_to_Ligand_Atoms(LigAtmHead,X,MS,Wshell,MolPockType)
 struct ATOM *LigAtmHead; /* (input /to be calculated) Rinaccess will be written in tFactor of each atom. */
 struct CHAR3DMAP *X;     /* (input) MolVol or Pocket */
 struct MULSC_SHELL_PRB *MS;
 float  Wshell;       /* Width of the shellto be added to ligand atoms (only for calculating shell acc)*/
 char   MolPockType;  /* 'M'olVol, 'P'ocket */

   /* >> X for MolPockType = 'M'olVol << 

    Initial state  : multiscale closing 3D map (molecular volume). Voxel has minimum radius number.
                   0 -> outer space
                   1 -> sphere with 1-th radius.
                   2 -> sphere with 2-th radius.
                   : :
             Nradius -> sphere with Nradius-th radius.
                 255 -> VdW volume

    For the hamonic mean calculationg, we assume that 
    the region with 0 (outer space) has Rinaccess = Rlarge_max + bin_of_Rlarge. 
  
  >> X for MolPockType for 'P'ocket  << 
    Initial state  : multiscale closing 3D map. Voxel has minimum radius number.
                   0 -> outer space or inside of the protein 
                   1 -> sphere with 1-th radius.
                   2 -> sphere with 2-th radius.
                   : :
             Nradius -> sphere with Nradius-th radius.
 */

{
 struct ATOM *an;
 float D[3],DD[3],ddis,RR,RRshell;
 int   minN[3],maxN[3];
 int   i,x,y,z,k;
 int   NPinside, NPshell, NPshell_access;
 float Rsum,invRsum;
 int   Ncount[256];
 float Rmax, Rmax_sigma,invRmax_sigma;
 printf("#Label_Multiscale_MolVol_to_Ligand_Atoms(Wshell %f MolPockType '%c')\n",Wshell,MolPockType);

 
 Rmax = MS->Rarray[MS->Nradius];
 Rmax_sigma = MS->Rarray[MS->Nradius] + (MS->Rarray[MS->Nradius] - MS->Rarray[MS->Nradius-1]);
 invRmax_sigma = 1.0/Rmax_sigma;
 

 an = LigAtmHead;
 /** BIG LOOP for each atom as 'an' **/
 while ( an->next != NULL){
   an = an->next;
   an->tFactor = 0.0;
  
   NPinside = NPshell = NPshell_access = 0; 
   /* (1) set up minN[] and maxN[] for an */ 
   for (i=0;i<3;++i) {  
     minN[i] = (int)floor((an->Pos[i] - an->R - Wshell - X->OrigPos[i])/X->grid_width) - 1;
     maxN[i] = (int)ceil((an->Pos[i]  + an->R + Wshell - X->OrigPos[i])/X->grid_width) + 1; 
     if (minN[i]<0){ minN[i] = 0;}
     if (maxN[i]>=X->N[i]){ maxN[i] = X->N[i]-1;}
   }

   for (k=0;k<256;++k){
     Ncount[k] = 0;
   }

   RR  = an->R*an->R; 
   RRshell  = (an->R + Wshell)*(an->R + Wshell);
   
   /* (2) checking neighboring voxels around an using minN[] and maxN[] */ 
   
   for (x=minN[0];x<=maxN[0];++x){ 
     D[0] = X->OrigPos[0] + X->grid_width * x - an->Pos[0];
     DD[0] = D[0]*D[0];
     for (y=minN[1];y<=maxN[1];++y){
       D[1] = X->OrigPos[1] + X->grid_width * y - an->Pos[1];
       DD[1] = D[1]*D[1];
       for (z=minN[2];z<=maxN[2];++z){
         D[2] = X->OrigPos[2] + X->grid_width * z - an->Pos[2];
         ddis = DD[0] + DD[1] + D[2]*D[2];
 
         if (ddis<=RR){
           if (X->map[x][y][z]<255) { 
             Ncount[X->map[x][y][z]] += 1; 
	     NPinside += 1;
	   }  
         }
         else if (ddis<=RRshell){ 
           NPshell += 1;
           if (X->map[x][y][z]<255){
             NPshell_access += 1;  
	   }
         }
      } /* z */
     } /* y */
   } /* x */
 
   /* (3M) Calculting value for MolVol */ 
    if (MolPockType == 'M'){
      an->shellAcc =  an->Rinacc =  0.0;  
 
      if (NPshell>0){
	 an->shellAcc = (float)NPshell_access/(float)NPshell;
      }
      
      if (NPinside>0){
         Rsum = 0.0;
         invRsum = 0.0;
         for (k=1;k<255;++k){
            if (Ncount[k]>0){ 
              Rsum    += (Ncount[k]*MS->Rarray[k]);
              invRsum += (Ncount[k]*1.0/MS->Rarray[k]);
	    }
	 }
         Rsum     += (Ncount[0]*Rmax_sigma);
         invRsum  += (Ncount[0]*1.0/Rmax_sigma);
         /* an->Rinacc = Rsum/NPinside;  */
         an->Rinacc = 1.0/(invRsum/NPinside);  
         an->rw_shellAcc = Rsum/(NPinside*Rmax_sigma);
      }
      else if (NPinside==0){
        an->shellAcc = 1.0;
        an->Rinacc = Rmax_sigma;
        an->rw_shellAcc = 1.0;
      } 
    }
 
   /* (3P) Calculting value for Pocket */ 
    if (MolPockType == 'P'){
      an->pocketness = 0.0;
      invRsum = 0.0;
      for (k=1;k<255;++k){
        if (Ncount[k]>0){ 
          invRsum += (Ncount[k]*1.0/MS->Rarray[k]);
 	}
      }
      if (NPinside>0){
        an->pocketness = invRsum * MS->Rarray[1]/NPinside;
      }
    }

 } /* an */

} /* end of Label_Multiscale_MolVol_or_Pocket_to_Ligand_Atoms() */




void Label_Multiscale_Each_Pocket_Cluster_to_Ligand_Atoms(LigAtmHead,X,Z,MS)
 struct ATOM *LigAtmHead; /* Rinaccess will be written in tFactor of each atom. */
 struct CHAR3DMAP *X; 
 struct CHAR3DMAP *Z;  /* Map containing the cluster number */
 struct MULSC_SHELL_PRB *MS;

   /* >> X for MolPockType = 'M'olVol << 

    Initial state  : multiscale closing 3D map (molecular volume). Voxel has minimum radius number.
                   0 -> outer space
                   1 -> sphere with 1-th radius.
                   2 -> sphere with 2-th radius.
                   : :
             Nradius -> sphere with Nradius-th radius.
                 255 -> VdW volume

    For the hamonic mean calculationg, we assume that 
    the region with 0 (outer space) has Rinaccess = Rlarge_max + bin_of_Rlarge. 
  
  >> X for MolPockType for 'P'ocket  << 
    Initial state  : multiscale closing 3D map. Voxel has minimum radius number.
                   0 -> outer space or inside of the protein 
                   1 -> sphere with 1-th radius.
                   2 -> sphere with 2-th radius.
                   : :
             Nradius -> sphere with Nradius-th radius.
 */

{
 struct ATOM *an;
 float D[3],DD[3],ddis,RR;
 int   minN[3],maxN[3];
 int   i,x,y,z,k,c;
 int   NPinside,Nmax_cluster;
 float invRsum;
 int   Ncount[256][256];  /* count for index of radius */
 float Rmax, Rmax_sigma,invRmax_sigma;
 float pocketness_clus[256],max_pocketness;
 printf("#Label_Multiscale_Each_Pocket_Cluster_to_Ligand_Atoms()\n");
 
 
 Nmax_cluster = 254; 
 Rmax = MS->Rarray[MS->Nradius];
 Rmax_sigma = MS->Rarray[MS->Nradius] + (MS->Rarray[MS->Nradius] - MS->Rarray[MS->Nradius-1]);
 invRmax_sigma = 1.0/Rmax_sigma;

 an = LigAtmHead;
 while ( an->next != NULL){
   an = an->next;
   an->tFactor = 0.0;
  
   NPinside =  0; 
   
   for (i=0;i<3;++i) {  
       minN[i] = (int)floor((an->Pos[i] - an->R - X->OrigPos[i])/X->grid_width) - 1;
       maxN[i] = (int)ceil((an->Pos[i]  + an->R - X->OrigPos[i])/X->grid_width) + 1; 
       if (minN[i]<0) minN[i] = 0;
       if (maxN[i]>=X->N[i]) maxN[i] = X->N[i]-1;
    }
 
   for (c=0;c<=Nmax_cluster;++c){
     for (k=0;k<256;++k) Ncount[c][k] = 0;
   } 
  
   RR  = an->R*an->R; 
   
   for (x=minN[0];x<=maxN[0];++x){ 
     D[0] = X->OrigPos[0] + X->grid_width * x - an->Pos[0];
     DD[0] = D[0]*D[0];
     for (y=minN[1];y<=maxN[1];++y){
       D[1] = X->OrigPos[1] + X->grid_width * y - an->Pos[1];
       DD[1] = D[1]*D[1];
       for (z=minN[2];z<=maxN[2];++z){
         D[2] = X->OrigPos[2] + X->grid_width * z - an->Pos[2];
         ddis = DD[0] + DD[1] + D[2]*D[2];
 
         if (ddis<=RR) {
           if ((Z->map[x][y][z]>0)&&(Z->map[x][y][z]<=Nmax_cluster)){
             Ncount[Z->map[x][y][z]][X->map[x][y][z]] += 1;
           }
           NPinside += 1;
         }
      } /* z */
     } /* y */
   } /* x */
 
    /** Calculate Values **/
     max_pocketness = 0.0;
     for (c=1;c<=Nmax_cluster;++c){
       pocketness_clus[c] = 0.0;
       invRsum = 0.0;
       for (k=1;k<255;++k){
         if (Ncount[c][k]>0){ 
 	  invRsum += (Ncount[c][k]*1.0/MS->Rarray[k]);
 	 }
       }
       if (NPinside>0){ 
         pocketness_clus[c] = invRsum * MS->Rarray[1]/NPinside; 
       }
       if (c<=PAR.Nmax_cluster){ 
         an->pocketness_clus[c] = pocketness_clus[c]; 
       }
       if (pocketness_clus[c]>max_pocketness) {
 	  an->cluster_num = c; max_pocketness = pocketness_clus[c];
       } 
    }
 
 } /* an */

} /* end of Label_Multiscale_Each_Pocket_Cluster_to_Ligand_Atoms() */













void Label_Multiscale_MolVol_or_Pocket_to_Protein_Atom_Shells(ProAtmHead,X,MS,Wshell,MolPockType)
 struct ATOM *ProAtmHead;
 struct CHAR3DMAP *X; 
 struct MULSC_SHELL_PRB *MS;
 float  Wshell; /* Width of the shellto be added to protein atoms */
 char   MolPockType;  /* 'M'olVol, 'P'ocket */

  /* >> X for MolPockType = 'M'olVol << 

    Initial state  : multiscale closing 3D map (molecular volume). Voxel has minimum radius number.
                   0 -> outer space
                   1 -> sphere with 1-th radius.
                   2 -> sphere with 2-th radius.
                   : :
             Nradius -> sphere with Nradius-th radius.
                 255 -> VdW volume

    For the hamonic mean calculationg, we assume that 
    the region with 0 (outer space) has Rinaccess = Rlarge_max + bin_of_Rlarge. 
 
    Calculate "shellAcc", "Rinacc"

 
  >> X for MolPockType for 'P'ocket  << 
    Initial state  : multiscale closing 3D map. Voxel has minimum radius number.
                   0 -> outer space or inside of the protein 
                   1 -> sphere with 1-th radius.
                   2 -> sphere with 2-th radius.
                   : :
             Nradius -> sphere with Nradius-th radius.
 
    Calculate only "pocketness".
 */ 
{
 struct ATOM *an;
 float D[3],DD[3],ddis,RR,RRshell;
 int   minN[3],maxN[3];
 int   i,x,y,z,k;
 int   NPshell,NPshell_access;
 int   Ncount[256];  /* count for index of radius */
 float Rmax, Rmax_sigma, invRmax_sigma, Rmin, Rmin_sigma, invRmin_sigma;
 float invRsum,Rsum;
 printf("#Label_Multiscale_MolVol_to_Protein_Atom_Shells(MolPockType %c Wshell %f)\n",MolPockType,Wshell);
 Rmin = MS->Rarray[1];
 Rmax = MS->Rarray[MS->Nradius];

 Rmax_sigma = MS->Rarray[MS->Nradius] + (MS->Rarray[MS->Nradius] - MS->Rarray[MS->Nradius-1]);
 Rmin_sigma = MS->Rarray[1] - (MS->Rarray[MS->Nradius] - MS->Rarray[MS->Nradius-1]);

 invRmax_sigma = 1.0/Rmax_sigma;
 invRmin_sigma = 1.0/Rmin_sigma;
 printf("#Rmax %f Rmax_sigma %f Rmin %f Rmin_sigma %f\n",Rmax,Rmax_sigma,Rmin,Rmin_sigma);

 an = ProAtmHead;
 while ( an->next != NULL){
   an = an->next;
   an->tFactor = 0.0;
   NPshell = NPshell_access = 0; 
   
   for (i=0;i<3;++i) 
    {  minN[i] = (int)floor((an->Pos[i] - an->R - Wshell - X->OrigPos[i])/X->grid_width) - 1;
       maxN[i] = (int)ceil ((an->Pos[i] + an->R + Wshell - X->OrigPos[i])/X->grid_width) + 1; 
       if (minN[i]<0) minN[i] = 0;
       if (maxN[i]>=X->N[i]) maxN[i] = X->N[i]-1;
    }
   for (k=0;k<256;++k) Ncount[k] = 0;
  
   RR       = an->R*an->R; 
   RRshell  = (an->R + Wshell)*(an->R + Wshell);
   for (x=minN[0];x<=maxN[0];++x){ 
     D[0] = X->OrigPos[0] + X->grid_width * x - an->Pos[0];
     DD[0] = D[0]*D[0];
     for (y=minN[1];y<=maxN[1];++y){
       D[1] = X->OrigPos[1] + X->grid_width * y - an->Pos[1];
       DD[1] = D[1]*D[1];
       for (z=minN[2];z<=maxN[2];++z){
         D[2] = X->OrigPos[2] + X->grid_width * z - an->Pos[2];
         ddis = DD[0] + DD[1] + D[2]*D[2];
         if ((RR<ddis)&&(ddis<=RRshell)){
           NPshell += 1;
           Ncount[X->map[x][y][z]] += 1;
         if (X->map[x][y][z]<255){  
            NPshell_access += 1; 
	 } 
        }
      } /* z */
     } /* y */
   } /* x */
 
 
   /** Calculate Values **/
 
    if (MolPockType == 'M'){
      an->shellAcc =  an->Rinacc =  an->rw_shellAcc = 0.0;  
 
      if (NPshell>0){
        an->shellAcc = (float)NPshell_access/(float)NPshell;
        invRsum = Rsum = 0.0;
        for (k=1;k<255;++k){
	  if (Ncount[k]>0){ 
	    invRsum += (Ncount[k]*1.0/MS->Rarray[k]);
	  }
	}
        for (k=1;k<255;++k){
	  if (Ncount[k]>0){
	    Rsum    += Ncount[k]*MS->Rarray[k];
	  }
	}
        Rsum      += Ncount[0]*Rmax_sigma;
        invRsum   += Ncount[0]*1.0/Rmax_sigma;
        an->Rinacc = Rsum/NPshell;
        an->rw_shellAcc   = Rsum/(NPshell*Rmax_sigma);
      }
    }
 
    if (MolPockType == 'P'){
      an->pocketness = 0.0;
      invRsum = 0.0;
      for (k=1;k<255;++k){
         if (Ncount[k]>0){ 
           invRsum += (Ncount[k]*1.0/MS->Rarray[k]);
 	}
      }
      if (NPshell>0){
        an->pocketness = invRsum * MS->Rarray[1]/NPshell;
      }
    }
 
 } /* an */

} /* end of Label_Multiscale_MolVol_or_Pocket_to_Protein_Atom_Shells() */



void Label_Multiscale_Each_Pocket_Cluster_to_Protein_Atom_Shells(ProAtmHead,X,Z,MS,Wshell)
 struct ATOM *ProAtmHead;
 struct CHAR3DMAP *X; 
 struct CHAR3DMAP *Z;  /* Map containing the cluster number */
 struct MULSC_SHELL_PRB *MS;
 float  Wshell; /* Width of the shellto be added to protein atoms */
 /* 
  >> X for MolPockType for 'P'ocket  << 
    Initial state  : multiscale closing 3D map. Voxel has minimum radius number.
                   0 -> outer space or inside of the protein 
                   1 -> sphere with 1-th radius.
                   2 -> sphere with 2-th radius.
                   : :
             Nradius -> sphere with Nradius-th radius.
 
    Calculate only "pocketness".
 */ 
{
 struct ATOM *an;
 float D[3],DD[3],ddis,RR,RRshell;
 int   minN[3],maxN[3];
 int   i,x,y,z,k,c;
 int   NPshell,NPshell_access;
 int   Ncount[MAX_CLUSTER_NUM][256];  /* count for index of radius */
 float Rmax, Rmax_sigma, invRmax_sigma, Rmin, Rmin_sigma, invRmin_sigma;
 float invRsum;
 printf("#Label_Multiscale_Each_Pocket_Cluster_to_Protein_Atom_Shells(Wshell %f Nmax_cluster %d)\n",Wshell,PAR.Nmax_cluster);
 Rmin = MS->Rarray[1];
 Rmax = MS->Rarray[MS->Nradius];

 Rmax_sigma = MS->Rarray[MS->Nradius] + (MS->Rarray[MS->Nradius] - MS->Rarray[MS->Nradius-1]);
 Rmin_sigma = MS->Rarray[1] - (MS->Rarray[MS->Nradius] - MS->Rarray[MS->Nradius-1]);

 invRmax_sigma = 1.0/Rmax_sigma;
 invRmin_sigma = 1.0/Rmin_sigma;

 an = ProAtmHead;
 while ( an->next != NULL){
   an = an->next;
   an->tFactor = 0.0;
   NPshell = NPshell_access = 0; 
   
   for (i=0;i<3;++i){  
     minN[i] = (int)floor((an->Pos[i] - an->R - Wshell - X->OrigPos[i])/X->grid_width) - 1;
     maxN[i] = (int)ceil ((an->Pos[i] + an->R + Wshell - X->OrigPos[i])/X->grid_width) + 1; 
     if (minN[i]<0) minN[i] = 0;
     if (maxN[i]>=X->N[i]) maxN[i] = X->N[i]-1;
    }
 
   for (c=0;c<=PAR.Nmax_cluster;++c){
     for (k=0;k<256;++k) Ncount[c][k] = 0;
   } 
 
   RR       = an->R*an->R; 
   RRshell  = (an->R + Wshell)*(an->R + Wshell);
   for (x=minN[0];x<=maxN[0];++x){ 
     D[0] = X->OrigPos[0] + X->grid_width * x - an->Pos[0];
     DD[0] = D[0]*D[0];
     for (y=minN[1];y<=maxN[1];++y){
       D[1] = X->OrigPos[1] + X->grid_width * y - an->Pos[1];
       DD[1] = D[1]*D[1];
       for (z=minN[2];z<=maxN[2];++z){
         D[2] = X->OrigPos[2] + X->grid_width * z - an->Pos[2];
         ddis = DD[0] + DD[1] + D[2]*D[2];
         if ((RR<ddis)&&(ddis<=RRshell)){
           NPshell += 1;
           if ((Z->map[x][y][z]>0)&&(Z->map[x][y][z]<=PAR.Nmax_cluster)){
		   Ncount[Z->map[x][y][z]][X->map[x][y][z]] += 1;
	   }
          /* Ncount[X->map[x][y][z]] += 1; */
           if (X->map[x][y][z]<255){ 
	     NPshell_access += 1;  
          }
        }
      } /* z */
     } /* y */
   } /* x */
 
 
   /** Calculate Values **/
       
     for (c=1;c<=PAR.Nmax_cluster;++c){
       an->pocketness_clus[c] = 0.0;
       invRsum = 0.0;
       for (k=1;k<255;++k) {
         if (Ncount[c][k]>0){ 
           invRsum += (Ncount[c][k]*1.0/MS->Rarray[k]);
	 }
       }	
       if (NPshell>0){
        an->pocketness_clus[c] = invRsum * MS->Rarray[1]/NPshell;
       }
     }
 
 } /* an */

} /* end of Label_Multiscale_Each_Pocket_Cluster_to_Protein_Atom_Shells() */




void Cal_Number_of_Contacting_Atoms_with_Ligand(ProAtmHead,LigAtmHead)
   struct ATOM *ProAtmHead, *LigAtmHead;
{
 struct ATOM *an,*ln;
 float D,Dthre,dx,dy,dz;

 Dthre = 1.4;
 /* (1) initialize */
 an = ProAtmHead;
 while (an->next != NULL){
   an = an->next; 
   an->res->Natom_contact = 0;
 }
 
 ln = LigAtmHead;
 while (ln->next != NULL){
  ln = ln->next;
  ln->contact = 1;
 }  
 

 /* (2) Cal distance between an and ln  */
 an = ProAtmHead;
 while (an->next != NULL){
   an = an->next;
   an->contact = 0;
   ln = LigAtmHead;

   while (ln->next != NULL){
    ln = ln->next;
    dx = an->Pos[0] - ln->Pos[0]; 
    dy = an->Pos[1] - ln->Pos[1]; 
    dz = an->Pos[2] - ln->Pos[2];
    D  = sqrt(dx*dx + dy*dy + dz*dz);
    if ( D < (an->R + ln->R + Dthre)) {an->contact = 1;  ln->contact = 1; }
   } 
  
   if (an->contact==1) an->res->Natom_contact += 1;
 }


} /* end Cal_Number_of_Contacting_Atoms_with_Ligand() */




void Write_PDB_File_with_Multiscale_MolVol_and_Pocket(fname,mode,Ahead,Type,MS,tFacType,title,CommentHead,UseResRinacc,Out_pocketness_clus)
 char *fname;
 char   mode; /* 'w'rite or 'a'ppend */
 struct ATOM *Ahead;
 char   Type;  /* 'A': output with "ATOM" (receptor), 'H':output with "HETATM" (ligand), 'P:output with "HETATM"(probe)
               otherwise:output with the original AtmHet */
 struct MULSC_SHELL_PRB *MS;
 char   tFacType;    /* Value of tFactor 'R'inaccess, 'I'nvRinacc, 'P'ocketness, 'S'hell accessiblity */
 char   *title;
 struct STRING_NODE *CommentHead;
 char   UseResRinacc;   /* 'R': use Rinacc assiged for Residue */
 char   Out_pocketness_clus; /* 'T','S':output pocketness_clus, '-' or 'F': not output pocket_clus */
{
 FILE *fp;
 struct ATOM *an;
 int r,c;
 char core[128],buff[512],tFacStr[32];
 struct STRING_NODE *sn;
 float Rmax,Rmax_sigma,invRmax_sigma;
 float Rinacc,shellAcc,pocketness,tFactor,rw_shellAcc;
  
  printf("#Write_PDB_File_with_Multiscale_MolVol_and_Pocket(tFac '%c' useRes '%c' Out_pocket_clus '%c') -> \"%s\"\n",
    tFacType,UseResRinacc,Out_pocketness_clus,fname);
  if (tFacType == 'R') sprintf(tFacStr,"Riacc");
  if (tFacType == 'I') sprintf(tFacStr,"iRiac");
  if (tFacType == 'S') sprintf(tFacStr,"shAcc");
  if (tFacType == 'P') sprintf(tFacStr,"pocke");
  if (tFacType == 'Z') sprintf(tFacStr,"zero");
  Rmax = MS->Rarray[MS->Nradius];
  Rmax_sigma = MS->Rarray[MS->Nradius] + (MS->Rarray[MS->Nradius] - MS->Rarray[MS->Nradius-1]);
  invRmax_sigma = 1.0/Rmax_sigma;
 
  if (fname[0]=='-') fp = stdout;
  else{
        if (mode=='w') fp = fopen(fname,"w");
   else if (mode=='a') fp = fopen(fname,"a");
   else                fp = fopen(fname,"w");
   if (fp==NULL) { printf("#ERROR:Can't write to \"%s\"\n",fname); exit(1);}
  }
                                                                                                                   
  Find_Filename_Core(core,fname);
  Get_Part_Of_Line(buff,fname,0,29);
  if (Type=='A'){
    fprintf(fp,"HEADER    RECEPTOR %-30s %s   %s\n",buff,Get_Date_String_PDB(),PAR.pdbidPro);
  }
  else if (Type=='H'){
    fprintf(fp,"HEADER    LIGAND   %-30s %s   %s\n",buff,Get_Date_String_PDB(),PAR.pdbidPro);
  }
  else if (Type=='P'){ 
    fprintf(fp,"HEADER    PROBE    %-30s %s   %s\n",buff,Get_Date_String_PDB(),PAR.pdbidPro);
  }
  else{ 
    fprintf(fp,"HEADER             %-30s %s   %s\n",buff,Get_Date_String_PDB(),PAR.pdbidPro);
  }

  if (title[0]!='\0') fprintf(fp,"TITLE   %s\n",title);
  fprintf(fp,"REMARK  DATE    %s\n",Get_Date_String());
  if (PAR.COMMAND[0] != '\0') fprintf(fp,"REMARK  COMMAND %s\n",PAR.COMMAND);
  sn = CommentHead;
  while (sn->next != NULL){ 
    sn = sn->next;
    fprintf(fp,"REMARK  COMMENT %s\n",sn->line); 
  }
  /*** Write MultiScale Probe information ***/
  fprintf(fp,"REMARK  MULSC_PROBE: NRADIUS %d\n",MS->Nradius);
  for (r=1;r<=MS->Nradius;++r){ 
     fprintf(fp,"REMARK  MULSC_PROBE: RADIUS %3d th %6.3f Ang %6.3f 1/Ang \n", r,MS->Rarray[r],1.0/MS->Rarray[r]); 
  }

  fprintf(fp,"REMARK  OUTSIDE_OF_MAX_RADIUS      %6.3f Ang %6.3f 1/Ang \n", Rmax_sigma,invRmax_sigma); 

  /*** Write Explanations for ATOM/HETATM lines ***/
  fprintf(fp,"REMARK  UseResidueBasedRinacc:%c\n",UseResRinacc);
  fprintf(fp,"REMARK   Occupancy[55-60]:Radius of the atom [A]\n");
  fprintf(fp,"REMARK  For protein atoms,\n");
  fprintf(fp,"REMARK  [Ratom](10: 55 - 60):Radius of atom[A]\n");
  fprintf(fp,"REMARK  [tFact](11: 61 - 66):[%s] is assigned as tFactor\n",tFacStr);
  fprintf(fp,"REMARK  [shAcc](12: 67 - 72):Shell accessibility. (Ratio of noVdW grids in the shell)[%%]\n");
  fprintf(fp,"REMARK  [Riacc](13: 73 - 78):averaged Rinaccess in the 'shell'[A]\n");
  fprintf(fp,"REMARK  [co]   (14: 79 - 81):Contact with other molecules (0 or 1)\n"); 
  fprintf(fp,"REMARK  [pocke](15: 82 - 87):pocketness. sum of 1/[Rpocket] /(1/[Rmin]*[vol of shell]) (%%)\n");
  if ((Out_pocketness_clus!='F') && (Out_pocketness_clus!='-')){
    for (c=1;c<=PAR.Nmax_cluster;++c){
      fprintf(fp,"REMARK  [pock%d](%d:%3d -%3d):pocketness_clus[%d] sum of 1/[Rpocket_for_pocke_tcluster%d] /(1/[Rmin]*[vol of shell]) (%%)\n",
     c,15+c,82+6*c,87+6*c,c,c);
    }
  }

  an = Ahead->next;
  if ((an!=NULL) && (an->conatm1!=NULL) && (an->conatm2!=NULL) && (an->conatm3 != NULL)){
    fprintf(fp,"REMARK  [cAt1] (%2d:%3d -%3d):Atom number of contact protein atom 1\n",
     16+PAR.Nmax_cluster, 88+6*PAR.Nmax_cluster,92+6*PAR.Nmax_cluster);
    fprintf(fp,"REMARK  [cAt2] (%2d:%3d -%3d):Atom number of contact protein atom 2\n",
     17+PAR.Nmax_cluster, 93+6*PAR.Nmax_cluster,97+6*PAR.Nmax_cluster);
    fprintf(fp,"REMARK  [cAt3] (%2d:%3d -%3d):Atom number of contact protein atom 3\n",
     18+PAR.Nmax_cluster, 98+6*PAR.Nmax_cluster,103+6*PAR.Nmax_cluster);
  }

  /*** Write One-line Explanation for ATOM/HETATM lines ***/
  fprintf(fp,"REMARK                                                 Ratom|%s|shAcc|Riacc|co|pocke|",tFacStr);
 
  if ((Out_pocketness_clus!='F') && (Out_pocketness_clus!='-')){
   for (c=1;c<=PAR.Nmax_cluster;++c)  fprintf(fp,"pock%d|",c);
  }

  an = Ahead->next;
  if ((an!=NULL) && (an->conatm1!=NULL) && (an->conatm2!=NULL) && (an->conatm3 != NULL))
      fprintf(fp,"cAt1|cAt2|cAt3|");
 
  fprintf(fp,"\n");

  /*** Write ATOM/HETATM ***/
  an = Ahead;

  while (an->next != NULL){
   an = an->next;
   if (an->AHtype == 'T'){
    if ((an->prev != Ahead)&&(an->prev!=NULL)&&(an->prev->AHtype != 'T')) fprintf(fp,"TER\n");
   }
   else{
   if ((an->prev != Ahead)&&(an->prev!=NULL)&&(an->AHtype != an->prev->AHtype)&&(an->prev->AHtype !='T')) 
    fprintf(fp,"TER\n");
   if (UseResRinacc== 'R'){
      Rinacc = an->res->Rinacc;
      shellAcc = an->res->shellAcc; 
      rw_shellAcc =  an->res->rw_shellAcc; 
      pocketness = an->res->pocketness; 
    }
    else{
      Rinacc = an->Rinacc;
      shellAcc = an->shellAcc; 
      rw_shellAcc =  an->rw_shellAcc; 
      pocketness = an->pocketness; 
    }
 
         if (tFacType == 'R'){ tFactor = Rinacc;}
    else if (tFacType == 'S'){ tFactor = 100.0*shellAcc;}
    else if (tFacType == 'P'){ tFactor = 100.0*pocketness;}
    else if (tFacType == 'r'){ tFactor = 100.0*rw_shellAcc;}
    else if (tFacType == 'Z'){ tFactor = 0.0;}
    else{ tFactor = 0.0;}

         if (Type=='A'){ fprintf(fp,"ATOM  ");}
    else if (Type=='H'){ fprintf(fp,"HETATM");}
    else if (Type=='P'){ fprintf(fp,"HETATM");}
    else {
           if (an->AHtype=='A'){ fprintf(fp,"ATOM  ");}
      else if (an->AHtype=='H'){ fprintf(fp,"HETATM");}
      else { fprintf(fp,"ATOM? ");}
    }
 
    fprintf(fp,"%5s %4s %3s %c%5s   %8.3f%8.3f%8.3f%5.2f %6.2f",
         an->Anum,an->Atom,an->Resi,an->Chain,an->Rnum,
         an->Pos[0],an->Pos[1],an->Pos[2],an->R,tFactor);
    
    fprintf(fp,"%6.1f%6.2f  %1d%6.2f",100.0*shellAcc,Rinacc,an->contact,100.0*pocketness);

    if ((Out_pocketness_clus!='F') && (Out_pocketness_clus!='-')){
      for (c=1;c<=PAR.Nmax_cluster;++c) fprintf(fp,"%6.2f",100.0*an->pocketness_clus[c]);
     }

    if ((an->conatm1!=NULL) && (an->conatm2!=NULL) && (an->conatm3 != NULL))
      fprintf(fp,"%5s%5s%5s",an->conatm1->Anum,an->conatm2->Anum,an->conatm3->Anum); 
    fprintf(fp,"\n");
   } 
  }
  fprintf(fp,"TER\n");
  if (fp!=stdout) fclose(fp);

} /* end of Write_PDB_File_with_Multiscale_MolVol_and_Pocket() */
