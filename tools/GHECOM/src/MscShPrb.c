/*

 <MscShPrb.c>
 
 for dealing multi scaling "shell" sphere probes  

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


/*** FUNCTIONS (GLOBAL) ***/
void Set_Rarray_By_Min_Max_Bin();
void Set_Rarray_By_Nradius_Min_Max();
void Make_MultiScale_Sphere_Shell_Probe_Map();
void Write_PDB_File_MultiSc_Shell_Spheres();
void Dilation_MultiScaleShell();
void Erosion_MultiScaleShell();
void Opening_MultiScaleShell();

/*** FUNCTIONS (LOCAL) ***/
static int Find_Index_Number_for_Radius();
static void Add_POS_INT();



void Set_Rarray_By_Min_Max_Bin(Rarray,Nradius,Rmin,Rmax,Rbin)
 float *Rarray;
 int   *Nradius;
 float Rmin,Rmax,Rbin;
{
    int k;
    Rarray[0] = 0.0;
    Rarray[1] = Rmin;
    k = 1;
    do
    {++k;
     Rarray[k] = Rarray[k-1] + Rbin; }
     while (Rarray[k]<Rmax);
    if (Rarray[k]>Rmax) Rarray[k] = Rmax;
    *Nradius = k;
    printf("#Nradius %d\n",*Nradius);
    for (k=1;k<=(*Nradius);++k){ printf("#%d %f\n",k,Rarray[k]);}

} /* end of Set_Rarray_By_Min_Max_Bin() */


void Set_Rarray_By_Nradius_Min_Max(Rarray,Nradius,Rmin,Rmax,Rbin)
 float *Rarray;
 int   Nradius;
 float Rmin,Rmax,Rbin;
{
    int k;
    for (k=0;k<MAX_RADIUS_TYPE;++k) Rarray[k] = 0.0;
    for (k=1;k<=Nradius;++k)
       Rarray[k] = Rmin + (k-1)*(Rmax - Rmin)/(Nradius-1);
    printf("#Nradius %d\n",Nradius);
    for (k=1;k<=Nradius;++k) printf("%d %f\n",k,Rarray[k]);
} /* end of Set_Rarray_By_Nradius_Min_Max() */



void Make_MultiScale_Sphere_Shell_Probe_Map(S)
 struct MULSC_SHELL_PRB  *S;
/*
  S->Rarray[] is already malloced and inputted.
 */
{
 float Pos[3],OrigPos[3],DD,*RRarray,Rmax,RRmax;
 int Ngrid,Ngrid_half,i,x,y,z,r;
 
 /** [1] Decide Grid Number **/
 Rmax = S->Rarray[S->Nradius];
 Ngrid_half = (int)(Rmax/S->grid_width);
 Ngrid = 2 * Ngrid_half + 1;
 for (i=0;i<3;++i) 
  { S->Ngrid[i]      = Ngrid;
    S->Ngrid_half[i] = Ngrid_half; }
 OrigPos[0] = OrigPos[1] = OrigPos[2] = 0.0;
 printf("#Make_MultiScale_Sphere_Shell_Probe_Map(Nradius %d N %d %d %d)\n",
  S->Nradius,S->Ngrid[0],S->Ngrid[1],S->Ngrid[2]);
 
 /** [2] Make RRarray **/
 RRarray  = (float *)malloc(sizeof(float)*(S->Nradius+2));
 for (r=1;r<=S->Nradius;++r)
  { S->HeadPosInt[r].next = NULL; 
    S->HeadPosInt[r].num   = -1;
    /* pn = &(S->HeadPosInt[r]);  */
    RRarray[r] = S->Rarray[r] * S->Rarray[r];
    S->Npoint[r] = 0; }
 RRmax = Rmax*Rmax;


 /** [3] Make Linked list of POS_INT **/
 for (x=-S->Ngrid_half[0];x<= S->Ngrid_half[0];++x){
  Pos[0] = OrigPos[0] + S->grid_width * x;
  for (y=-S->Ngrid_half[1]; y<= S->Ngrid_half[1];++y){
   Pos[1] = OrigPos[1] + S->grid_width * y;
   for (z=-S->Ngrid_half[2];z<= S->Ngrid_half[2];++z){
    Pos[2] = OrigPos[2] + S->grid_width * z;
    DD = Pos[0]*Pos[0] + Pos[1]*Pos[1] + Pos[2]*Pos[2];
    r = Find_Index_Number_for_Radius(DD,S->Nradius,RRarray);
    /* printf("DD %lf D %lf r %d\n",DD,sqrt(DD),r); */
    if (r>=0)
    { /* printf("x %d y %d z %d DD %f D %f r %d\n",x,y,z,DD,sqrt(DD),r);  */
      Add_POS_INT(&(S->HeadPosInt[r]),x,y,z,r);
      S->Npoint[r] += 1; }
   }
  } 
 }

 S->Rarray[S->Nradius+1]  = S->Rarray[S->Nradius] + 1.0;
 free(RRarray);

} /* end of Make_MultiScale_Sphere_Shell_Probe_Map() */




int Find_Index_Number_for_Radius(DD,Nradius,RRarray)
 float DD;
 int   Nradius;
 float *RRarray;
{
 int r;
 
 r = 1; 
 do{
  if (DD<=RRarray[r]) return(r);
  else ++r; 
 } while (r<=Nradius);

 return(-1);

} /* end of Find_Index_Number_for_Radius() */




void Add_POS_INT(HeadPos,x,y,z,n)
 struct POS_INT *HeadPos;
 int x,y,z,n;
{
 struct POS_INT *pn;

 /* 
 printf("#Add_POS_INT(HeadPos %d %d %d)\n",x,y,z);
 */

 pn = HeadPos;
 while (pn->next != NULL) pn = pn->next; 

 pn->next = (struct POS_INT*)malloc(sizeof(struct POS_INT));
 pn->next->prev = pn;
 pn->next->iPos[0] = x; 
 pn->next->iPos[1] = y; 
 pn->next->iPos[2] = z;
 pn->next->num     = n;
 pn->next->next = NULL;
} /* end of Add_POS_INT() */






void Write_PDB_File_MultiSc_Shell_Spheres(fname,S,comment,CommentHead)
 char *fname;
 struct MULSC_SHELL_PRB  *S;
 char *comment;
 struct STRING_NODE *CommentHead;
{
 FILE *fp;
 int i,r,Natom;
 char core[128],buff[512],buff2[512],AtomStr[32];
 struct STRING_NODE *sn;
 struct POS_INT *pn; 
 float fPos[3];
 
 printf("#Write_PDB_File_MultiSc_Shell_Spheres() -> \"%s\"\n",fname);
 fp = fopen(fname,"w");
 if (fp==NULL) { printf("#ERROR:Can't write to \"%s\"\n",fname); exit(1);}
 Find_Filename_Core(core,fname);
 Get_Part_Of_Line(buff,fname,0,40);
 Get_Part_Of_Line(buff2,core,0,9);
 fprintf(fp,"HEADER    %-40s%s   %-10s0XXX   1\n",buff,Get_Date_String_PDB(),buff2);
 fprintf(fp,"REMARK    DATE %s\n",Get_Date_String());
 fprintf(fp,"REMARK    Nradius %d\n",S->Nradius);
 for (r=1;r<=S->Nradius;++r) fprintf(fp,"REMARK    %2d -th radius %f\n",r,S->Rarray[r]);
 if (PAR.COMMAND[0] != '\0') fprintf(fp,"REMARK  COMMAND %s\n",PAR.COMMAND);
 if (comment[0]!='\0') fprintf(fp,"REMARK  COMMENT %s\n",comment);
                                                                                                                   
 sn = CommentHead;
 while (sn->next != NULL)
 { sn = sn->next;
   fprintf(fp,"REMARK  COMMENT %s\n",sn->line); }

 Natom = 0;
 
 for (r=1;r<=S->Nradius;++r){
  pn =  &(S->HeadPosInt[r]);
 
  while (pn->next != NULL){
   pn = pn->next;
   ++Natom;
   sprintf(buff,"C%d",r);
   Get_Part_Of_Line(AtomStr,buff,0,2);
   for (i=0;i<3;++i) fPos[i] = S->grid_width * pn->iPos[i]; 
   fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f %2d %2d %2d\n",
        Natom%100000,AtomStr,"GRD",' ',r,
        fPos[0],fPos[1],fPos[2],S->grid_width/2.0,S->Rarray[r],pn->iPos[0],pn->iPos[1],pn->iPos[2]);
  } /* pn */

 } /* r */


} /* end of Write_PDB_File_MultiSc_Shell_Spheres() **/



                                                                                             
void Dilation_MultiScaleShell(Xmap,Ymap,MS)
 struct CHAR3DMAP *Xmap;
   /* >> Xmap <<
    Initial state: target shape :3D map (voxel in shape has non-zero value)
                   0           -> outer space
                   more than 0 -> within the target shape 
    Final state  : just change to [more_than_0] to 255
                   0           -> outer space
                   255         -> within the target shape 
   */ 
 struct CHAR3DMAP *Ymap;
   /* >> Ymap <<
    Initial state: don't care. Ymap is initialized to zero in this function.

    Final state  : closed 3D map. Voxel has minimum radius number.
                   0 -> outer space
                   1 -> sphere with 1-th radius
                   2 -> sphere with 2-th radius.
                   :
             Nradius -> sphere with Nradius-th radius.
   */
  struct MULSC_SHELL_PRB *MS;
{
 int x,y,z,r;
 int xx,yy,zz;
 struct POS_INT *pn;
 
 printf("#Dilation_MultiScale(Xmap,Ymap,MS,MS->Nradius %d)\n",MS->Nradius);
 
 /** Initialize Ymap to zero **/
 for (x=0;x<Ymap->N[0];++x)
  for (y=0;y<Ymap->N[1];++y)
   for (z=0;z<Ymap->N[2];++z) 
    { if (Xmap->map[x][y][z]>0) Xmap->map[x][y][z] = 255; 
      Ymap->map[x][y][z] = 0; }
 
 /** Scan Xmap and try translation by Pmap  **/
 for (x=0;x<Xmap->N[0];++x){
  if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]);
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){

    if (Xmap->map[x][y][z]>0){

     for (r=1;r<=MS->Nradius;++r){  
      pn = &(MS->HeadPosInt[r]);
      while (pn->next != NULL){
       pn = pn->next;
       xx = x + pn->iPos[0]; 
       if ((xx>=0)&&(xx<Xmap->N[0])){
        yy = y + pn->iPos[1]; 
        if ((yy>=0)&&(yy<Xmap->N[1])){
          zz = z + pn->iPos[2]; 
          if ((zz>=0)&&(zz<Xmap->N[2])){ 
           if ((Ymap->map[xx][yy][zz]==0)||(r < Ymap->map[xx][yy][zz])) 
               Ymap->map[xx][yy][zz] = r;
          } /* zz */
         } /* yy */
        } /* xx */

       } /* while(pn) */

      } /* r */
     } /* Xmap > 0 */
    } /* z */
  } /* y */
 } /* x */

} /* end of Dilation_MultiScaleShell() */





void Erosion_MultiScaleShell(Xmap,Ymap,MS)
 struct CHAR3DMAP *Xmap;
   /* >> Xmap <<
    Initial state  : dilated 3D map. Voxel has minimum radius number.
                   0 -> outer space
                   1 -> sphere with 1-th radius
                   2 -> sphere with 2-th radius.
                   :
             Nradius -> sphere with Nradius-th radius.
    Final state    : does not change.
   */
 struct CHAR3DMAP *Ymap;
  /* >> Ymap <<
    Initial state: 0           -> outer space
                   255         -> within the target shape 

    Final state  : opening 3D map. Voxel has minimum radius number.
                   0 -> outer space
                   1 -> sphere with 1-th radius.
                   2 -> sphere with 2-th radius.
                   : :
             Nradius -> sphere with Nradius-th radius.
                 255 -> VdW volume
   */
 struct MULSC_SHELL_PRB *MS;
{
 int x,y,z,xx,yy,zz,r;
 struct POS_INT *pn;
 char in_shell;    /* flag for inside of the shell */
 
 printf("#Erosion_MultiScaleShell(Xmap,Ymap,MS->Nradius %d)\n",MS->Nradius);
 
 for (x=0;x<Xmap->N[0];++x){
  if ((x%10)==0) printf("##x %d/%d\n",x,Xmap->N[0]);
  
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    if ((Xmap->map[x][y][z]>0)&&(Ymap->map[x][y][z]==0)){
     r = 1; 
      
     while (r<=MS->Nradius){
       in_shell = 1; 
       pn = &(MS->HeadPosInt[r]);

       while (pn->next != NULL){
         pn = pn->next;
         xx = x + pn->iPos[0]; 
         yy = y + pn->iPos[1]; 
         zz = z + pn->iPos[2]; 
        
         if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){ 
              if (Xmap->map[xx][yy][zz]> r) {in_shell = 0;  break;}
              if (Xmap->map[xx][yy][zz]==0) {in_shell = 0; r = MS->Nradius + 1; break;}
         }
         else {in_shell = 0; r = MS->Nradius + 1;  break;}
 
        } /* pn */ 

       if (in_shell==1) { Ymap->map[x][y][z] = r; r = MS->Nradius + 1;}
       else  r += 1;

     } /* while r */

    } /* Xmap > 0 */
    
    } /* z */
  } /* y */
 } /* x */


} /* end of Erosion_MultiScaleShell() */
                                                                                               




void Opening_MultiScaleShell(Xmap,Ymap,Pmap)
 struct CHAR3DMAP *Xmap;
  /* >> Xmap <<
    Initial state  : closing 3D map. Voxel has minimum radius number.
                   0 -> outer space
                   1 -> sphere with 1-th radius.
                   2 -> sphere with 2-th radius.
                   : :
             Nradius -> sphere with Nradius-th radius.
                 255 -> VdW volume
   */
 struct CHAR3DMAP *Ymap;
  /* 
    Initial state  : Don't care.
    Final state    : pocket 3D map. Voxel has minimum radius number.
                   0 -> outer space
                   1 -> sphere with 1-th radius.
                   2 -> sphere with 2-th radius.
                   : :
             Nradius -> sphere with Nradius-th radius.
                 255 -> VdW volume
   */
 struct CHAR3DMAP *Pmap; /* Map for small probe spheres */
{
 int x,y,z,xx,yy,zz,p,q,r;
 int Nhalf,kmax;
 char stop;
 int Nzero_out, Nzero_emp,Nzero_in, Npock;
 
 printf("#Opening_MultiScaleShell(Xmap,Ymap,Pmap)\n");
 Nhalf = (Pmap->N[0]-1)/2;                                                                      
 
 /** [1] Initialize Ymap to zero **/
 for (x=0;x<Ymap->N[0];++x)
  for (y=0;y<Ymap->N[1];++y)
   for (z=0;z<Ymap->N[2];++z) 
    { Ymap->map[x][y][z] = 0;
      if (Xmap->map[x][y][z]==255) Ymap->map[x][y][z] = 255; }
 
 /** [2] Scan Xmap and try translation by Pmap  **/
 Nzero_emp = Nzero_out = Nzero_in = Npock = 0; 
 for (x=0;x<Xmap->N[0];++x){
  if ((x%10)==0) printf("##x %d/%d\n",x,Xmap->N[0]);
  
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    if ((Xmap->map[x][y][z]>0)&&(Xmap->map[x][y][z]<255)){
     kmax = 0; stop = 0;
     /** (1) Scan by probe Pmap, and find maximum value(radius) of molecular surface **/
     for (p=0;p<Pmap->N[0];++p)
      for (q=0;q<Pmap->N[1];++q)
       for (r=0;r<Pmap->N[2];++r)
       if (Pmap->map[p][q][r]==1){ 
        xx = x + p - Nhalf;
        yy = y + q - Nhalf;
        zz = z + r - Nhalf;
        if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){
               if (Xmap->map[xx][yy][zz]==255) { stop = 1;  kmax = 0; ++Nzero_in;}
          else if (Xmap->map[xx][yy][zz]==0)   { stop = 1;  kmax = 0; ++Nzero_out;}
          else if (Xmap->map[xx][yy][zz]>kmax) { kmax = Xmap->map[xx][yy][zz];}
         }
        else { stop = 1; kmax=0; ++Nzero_emp;}
      
       if (stop==1) {p = Pmap->N[0]; q = Pmap->N[1]; r = Pmap->N[2];}
      } /* p,q,r */

     /** (2) Fill the I(x) by a probe, and take minimum **/
      if (kmax>0){
        ++Npock; 
        /* printf("xyz %d %d %d kmax %d\n",x,y,z,kmax); */  
        for (p=0;p<Pmap->N[0];++p)
         for (q=0;q<Pmap->N[1];++q)
          for (r=0;r<Pmap->N[2];++r)
           if (Pmap->map[p][q][r]==1){
            xx = x + p - Nhalf;
            yy = y + q - Nhalf;
            zz = z + r - Nhalf;
            if ((xx>=0)&&(xx<Ymap->N[0])&&(yy>=0)&&(yy<Ymap->N[1])&&(zz>=0)&&(zz<Ymap->N[2])){
             if ((Ymap->map[xx][yy][zz]==0)||(kmax < Ymap->map[xx][yy][zz]))
               Ymap->map[xx][yy][zz] = kmax;
             }
           } /* p,q,r */
       } /* kmax > 0 */
      } /* Xmap >=2 */
    } /* z */
  } /* y */
 } /* x */

 printf("#Npock %d Nzero_emp %d Nzero_out %d Nzero_in %d\n",Npock,Nzero_emp,Nzero_out,Nzero_in);

} /* end of Opening_MultiScaleShell() */
