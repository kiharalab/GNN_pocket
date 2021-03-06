/*

 <OpnRotPrb.c>
 
  Functions for opening for 
  multiple probes generated by rotations

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
#include "OpnRotPrb.h"
#include "GrPtMergeSort.h"
#include "OpeStd3DPr.h"
#include "OpeMMorph.h"


/*** FUNCTIONS (GLOBAL) ***/
void Make_CHAR4DMAP_from_ProbeAtoms_and_RmatArray();
void Make_GRIDPNT3D_List_from_CHAR4DMAP_Multiple_Probes();
void Dilation_CHAR4DMAP();
int Opening_Using_CHAR4DMAP_and_GRIDPNT3D();
int Opening_Using_CHAR4DMAP_and_GRIDPNT3D_with_Fitness();
void Mark_Minimum_Filtering_GRIDPNT3Ds();
void Make_Another_GRID3DPNT_List_from_Marked_List();
int NonAllZero_Region_GRID3DPNT_Listed_Probes();
int Opening_for_Focused_Region();
int Malloc_CHAR4DMAP();
void Initialize_CHAR4DMAP();
void Make_Rod_Probe();
void Make_Disc_Probe();
void Make_CHAR4DMAP_for_Rotated_RodDisc_by_RmatArray();


/*** FUNCTIONS (LOCAL) ***/
static struct GRIDPNT3D* Add_GRIDPNT3D_Stack();



void Make_CHAR4DMAP_from_ProbeAtoms_and_RmatArray(MapPS4D,MapPS,Ahead,Nrotate,RmatArray,Rscale,bitmask)
 struct CHAR4DMAP *MapPS4D;
 struct CHAR3DMAP *MapPS;
 struct ATOM *Ahead;
 int    Nrotate;
 struct ROT_MATRIX *RmatArray;
 float  Rscale;
 char   bitmask;
{
 int r,x,y,z;
 printf("#Make_CHAR4DMAP_from_ProbeAtoms_and_RmatArray()\n");

 for (r=0;r<Nrotate;++r){
  Rotate_Molecule(Ahead,RmatArray[r].R);     
  Set_Bit_3DMAP_Inside_Atoms(MapPS,Ahead,bitmask,Rscale);
  for (x=0;x<MapPS->N[0];++x){
   for (y=0;y<MapPS->N[1];++y){
    for (z=0;z<MapPS->N[2];++z){
     if ((MapPS->map[x][y][z]&bitmask)==bitmask) 
       MapPS4D->map[x][y][z][r] = MapPS4D->map[x][y][z][r]|bitmask;
    } /* z */
   } /* y */
  } /* x */
  InvRotate_Molecule(Ahead,RmatArray[r].R);
  Reset_Specific_3DMAP(MapPS,bitmask);
 } /* r */

} /* end of Make_CHAR4DMAP_from_ProbeAtoms_and_RmatArray() */




void Make_GRIDPNT3D_List_from_CHAR4DMAP_Multiple_Probes(MapPS4D,Nrotate,PntHead,bitmask)
 struct CHAR4DMAP *MapPS4D;
 int Nrotate;
 struct GRIDPNT3D *PntHead;
 char bitmask;
{
 int r,x,y,z,N;
 struct GRIDPNT3D *pnt; 

 printf("#Make_GRIDPNT3D_List_from_CHAR4DMAP_Multiple_Probes()\n");

  PntHead->next = NULL; 
  for (x=0;x<MapPS4D->N[0];++x){
   for (y=0;y<MapPS4D->N[1];++y){
    for (z=0;z<MapPS4D->N[2];++z){
     N = 0; 
     for (r=0;r<Nrotate;++r) 
      {if ((MapPS4D->map[x][y][z][r]&bitmask)==bitmask) ++N;}
     if (N>0){
      pnt = Add_GRIDPNT3D_Stack(PntHead,x,y,z);
      pnt->Nprobe = N;
      pnt->value = pnt->Nprobe; 
     } 
    }
   }
  }
 
 Merge_Sort_Double_Linked_List_GRIDPNT3D(PntHead,'R');

} /* end of Make_GRIDPNT3D_List_from_CHAR4DMAP_Multiple_Probes() */



void Dilation_CHAR4DMAP(MapPS4D,MapPS,MapPSD,Nrotate,tar_bit,res_bit)
 struct CHAR4DMAP *MapPS4D;
 struct CHAR3DMAP *MapPS;
 struct CHAR3DMAP *MapPSD;
 int    Nrotate;
 char tar_bit;
 char res_bit;
{
 int r,x,y,z;
 printf("#Dilation_CHAR4DMAP()\n");

 for (r=0;r<Nrotate;++r){
  
 for (x=0;x<MapPS->N[0];++x){
   for (y=0;y<MapPS->N[1];++y){
    for (z=0;z<MapPS->N[2];++z){
      MapPS->map[x][y][z]= MapPS4D->map[x][y][z][r];
    } /* z */
   } /* y */
  } /* x */
 
  Dilation(MapPS,MapPSD,tar_bit,res_bit);
 
  for (x=0;x<MapPS->N[0];++x){
   for (y=0;y<MapPS->N[1];++y){
    for (z=0;z<MapPS->N[2];++z){
      MapPS4D->map[x][y][z][r] = MapPS->map[x][y][z];
    } /* z */
   } /* y */
  } /* x */

} /* r */

} /* end of Dilation_CHAR4DMAP() */










int Opening_Using_CHAR4DMAP_and_GRIDPNT3D(Xmap,Ymap,MapPS4D,PntHeadPS,tar_bit,foc_bit,res_bit,prb_bit)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 struct CHAR3DMAP *Ymap; /* 3D map for counting */
 struct CHAR4DMAP *MapPS4D;
 struct GRIDPNT3D *PntHeadPS;
 unsigned char   tar_bit; /* target bitmask (1,2,4,8,16,32,64,128)*/
 unsigned char   foc_bit; /* focused region bitmask (1,2,4,8,16,32,64,128). if 0, for every points*/
 unsigned char   res_bit; /* result bitmask (1,2,4,8,16,32,64,128)*/
 unsigned char   prb_bit; /* probe bitmask (1,2,4,8,16,32,64,128)*/
{
 int x,y,z,p,q,r,k,xx,yy,zz;
 int Nhalf[3],Nprobe,Nreject_probe,Naccept_probe,Nsurvive_center,Ntry;
 char end,out,hit; 
 struct GRIDPNT3D *pnt; 
 char *Nprobe_mark; 
 printf("#Opening_Using_CHAR4DMAP_and_GRIDPNT3D(tar %d foc %d res %d prb %d)\n",
     tar_bit,foc_bit,res_bit,prb_bit);

 /*** Initialize ***/ 
 for (k=0;k<3;++k) Nhalf[k] = (MapPS4D->N[k]-1)/2;
 Nsurvive_center = 0;
 Nprobe = MapPS4D->N[3];
 Nprobe_mark = (char *)malloc(sizeof(char)*Nprobe); 


 /*** Scan x,y,z,p,q,r ***/
 for (x=0;x<Xmap->N[0];++x){
  if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]);
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    if ((foc_bit==0)||((Xmap->map[x][y][z]&foc_bit) == foc_bit)){

    /** [1] Scanning **/
    Nreject_probe = 0; Ntry = 0;
    for (k=0;k<Nprobe;++k) Nprobe_mark[k] = 1;  
    pnt = PntHeadPS; end = 0;

    while ((pnt->next!=NULL)&&(end==0)){
     pnt = pnt->next;
     ++Ntry;
     xx = x + pnt->x - Nhalf[0];
     yy = y + pnt->y - Nhalf[1];
     zz = z + pnt->z - Nhalf[2];
     out = 0;
        if ((xx<0)||(xx>=Xmap->N[0])||(yy<0)||(yy>=Xmap->N[1])||(zz<0)||(zz>=Xmap->N[2])) 
     out = 1;
     else if ((Xmap->map[xx][yy][zz]&tar_bit) != tar_bit) out = 1;
      
     if (out==1){
       if (pnt->Nprobe==Nprobe) {Nreject_probe = Nprobe;}
       else{
        for (k=0;k<Nprobe;++k){
         if ((MapPS4D->map[pnt->x][pnt->y][pnt->z][k]&prb_bit)==prb_bit){ 
           if (Nprobe_mark[k] == 1) ++Nreject_probe; 
          Nprobe_mark[k] = 0;
         }   
        }  
      } 
     }
     
     if (Nreject_probe >= Nprobe) end = 1;
     
    } /* pnt */

   /*
      printf("xyz %d %d %d Ntry %d Nreject_probe %d  Nprobe %d\n",x,y,z,Ntry,Nreject_probe,Nprobe); 
    printf("#Ntry %d Nreject_probe %d Nprobe %d\n",Ntry,Nreject_probe,Nprobe); 
   */ 
     /** [2] Marking for surviving point **/
    Naccept_probe = Nprobe - Nreject_probe;

     if (Nreject_probe<Nprobe){  
     
      ++Nsurvive_center;

      /*
       if (Naccept_probe>255) Naccept_probe = 255;
       Ymap->map[x][y][z] = Naccept_probe;
      */ 
       for (p=0;p<MapPS4D->N[0];++p){
        for (q=0;q<MapPS4D->N[1];++q){
         for (r=0;r<MapPS4D->N[2];++r){
           k = 0; hit = 0;  
           while ((k<Nprobe)&&(hit==0)){     
             if ((Nprobe_mark[k]==1)&&((MapPS4D->map[p][q][r][k]&prb_bit)==prb_bit)) hit = 1;
             ++k;
           }

          if (hit==1){   
           xx = x + p - Nhalf[0];
           yy = y + q - Nhalf[1];
           zz = z + r - Nhalf[2];
           if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])){
            Xmap->map[xx][yy][zz] = Xmap->map[xx][yy][zz]|res_bit; 
            if (Ymap->map[xx][yy][zz]<255) Ymap->map[xx][yy][zz] += 1;  
           }  
         }

        } /* r */
       } /* q */
      } /* p */
     }  

    } /* foc_bit */

   } /* z */
  } /* y */
 } /* x */

 free(Nprobe_mark);
 printf("#Nsurvive_center %d\n",Nsurvive_center);
 return(Nsurvive_center);

} /* end of Opening_Using_CHAR4DMAP_and_GRIDPNT3D() */





int Opening_Using_CHAR4DMAP_and_GRIDPNT3D_with_Fitness(
   Xmap,Ymap,MapPS4D,PntHeadPS,
   tar_bit,foc_bit,res_bit,prt_bit,prb_bit,fit_bit,fit_ratio_thre,Ave_fit_ratio,Max_fit_ratio)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 struct CHAR3DMAP *Ymap; /* 3D map for counting */
 struct CHAR4DMAP *MapPS4D;
 struct GRIDPNT3D *PntHeadPS;
 unsigned char   tar_bit; /* target bitmask (1,2,4,8,16,32,64,128)*/
 unsigned char   foc_bit; /* focused region bitmask (1,2,4,8,16,32,64,128). if 0, for every points*/
 unsigned char   res_bit; /* result bitmask (1,2,4,8,16,32,64,128)*/
 unsigned char   prt_bit; /* protein bitmask (1,2,4,8,16,32,64,128)*/
 unsigned char   prb_bit; /* probe bitmask (1,2,4,8,16,32,64,128)*/
 unsigned char   fit_bit; /* fitshell bitmask (1,2,4,8,16,32,64,128)*/
 float  fit_ratio_thre;
 float  *Ave_fit_ratio,*Max_fit_ratio;
{
 int x,y,z,p,q,r,k,xx,yy,zz;
 int Nhalf[3],Nprobe,Nreject_probe,Naccept_probe,Ntry,Nprot_shell,Nnonprot_shell,Nshell;
 int Nsurvive_center,Nall_center; 
 char end,out,char_fit_ratio; 
 struct GRIDPNT3D *pnt; 
 char *Nprobe_mark; 
 float Nthre_nonprot,Nthre_prot,fit_ratio,fit_ratio_max,fit_ratio_ave;
 printf("#Opening_Using_CHAR4DMAP_and_GRIDPNT3D_with_Fitness(tar %d foc %d res %d prt %d prb %d fit %d\n",
     tar_bit,foc_bit,res_bit,prt_bit,prb_bit,fit_bit);

 /*** [I] Initialize ***/ 
 for (k=0;k<3;++k) Nhalf[k] = (MapPS4D->N[k]-1)/2;
 Nsurvive_center = Nall_center = 0;
 Nprobe = MapPS4D->N[3];
 Nprobe_mark = (char *)malloc(sizeof(char)*Nprobe); 

 Nshell = 0;
 for (p=0;p<MapPS4D->N[0];++p){
  for (q=0;q<MapPS4D->N[1];++q){
   for (r=0;r<MapPS4D->N[2];++r){
    if (((MapPS4D->map[p][q][r][0]&fit_bit)==fit_bit)
        &&((MapPS4D->map[p][q][r][0]&prb_bit)!=prb_bit)) ++Nshell;  
   }
  }
}
 Nthre_nonprot = (1.0-fit_ratio_thre)*Nshell;
 Nthre_prot = fit_ratio_thre*Nshell;
 printf("#Nshell %d thre noprot %f prot %f\n",Nshell,Nthre_nonprot,Nthre_prot);
 
 /*** [II] Scan x,y,z,p,q,r ***/
 fit_ratio_max = fit_ratio_ave = 0.0;

 for (x=0;x<Xmap->N[0];++x){
  if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]);
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    if ((foc_bit==0)||((Xmap->map[x][y][z]&foc_bit) == foc_bit)){

    /** (1) Scanning **/
    Nreject_probe = 0; Ntry = 0;
    for (k=0;k<Nprobe;++k) Nprobe_mark[k] = 1;  
    end = 0;
   
    pnt = PntHeadPS; 
    while ((pnt->next!=NULL)&&(end==0)){
     pnt = pnt->next;
     ++Ntry;
     xx = x + pnt->x - Nhalf[0];
     yy = y + pnt->y - Nhalf[1];
     zz = z + pnt->z - Nhalf[2];
     out = 0;
        if ((xx<0)||(xx>=Xmap->N[0])||(yy<0)||(yy>=Xmap->N[1])||(zz<0)||(zz>=Xmap->N[2])) 
     out = 1;
     else if ((Xmap->map[xx][yy][zz]&tar_bit) != tar_bit) out = 1;
      
     if (out==1){
       if (pnt->Nprobe==Nprobe) {Nreject_probe = Nprobe;}
       else{
        for (k=0;k<Nprobe;++k){
         if ((MapPS4D->map[pnt->x][pnt->y][pnt->z][k]&prb_bit)==prb_bit){ 
           if (Nprobe_mark[k] == 1) ++Nreject_probe; 
          Nprobe_mark[k] = 0;
         }   
        }  
      } 
     }
     
     if (Nreject_probe >= Nprobe) end = 1;
     
    } /* pnt */

   /*
    printf("#Ntry %d Nreject_probe %d Nprobe %d\n",Ntry,Nreject_probe,Nprobe); 
      printf("xyz %d %d %d Ntry %d Nreject_probe %d  Nprobe %d\n",x,y,z,Ntry,Nreject_probe,Nprobe); 
   */ 
   
    /** (2) Marking for surviving point **/
    Naccept_probe = Nprobe - Nreject_probe;

    if (Nreject_probe<Nprobe){  
      for (k=0;k<Nprobe;++k){    
       if (Nprobe_mark[k]==1){
        Nprot_shell = Nnonprot_shell = 0;
        ++Nsurvive_center;
       for (p=0;p<MapPS4D->N[0];++p){
           xx = x + p - Nhalf[0];
        for (q=0;q<MapPS4D->N[1];++q){
           yy = y + q - Nhalf[1];
         for (r=0;r<MapPS4D->N[2];++r){
           zz = z + r - Nhalf[2];
           if((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2]) &&
              ((MapPS4D->map[p][q][r][k]&fit_bit)==fit_bit)&&((MapPS4D->map[p][q][r][k]&prb_bit)!=prb_bit)
              && ((Xmap->map[xx][yy][zz]&prt_bit)==prt_bit)) ++Nprot_shell; else ++Nnonprot_shell;
         /* 
         if (Nnonprot_shell > Nthre_nonprot) p = q = r = MapPS4D->N[0]*MapPS4D->N[1]*MapPS4D->N[2];
          */
        } /* r */
       } /* q */
      } /* p */


       fit_ratio = (float)Nprot_shell/(float)Nshell;
       if (fit_ratio>fit_ratio_max) fit_ratio_max = fit_ratio;
       fit_ratio_ave += fit_ratio;
       ++Nall_center;
 
       char_fit_ratio = (int)floor(255.0*fit_ratio); 
       
       if (Nprot_shell  >= Nthre_prot){

       printf("Accept xyz %d %d %d k %d Nprot_shell %d Nshell %d ratio %f\n",
        x,y,z,k,Nprot_shell,Nshell,(float)Nprot_shell/(float)Nshell); 
       for (p=0;p<MapPS4D->N[0];++p){
           xx = x + p - Nhalf[0];
        for (q=0;q<MapPS4D->N[1];++q){
           yy = y + q - Nhalf[1];
         for (r=0;r<MapPS4D->N[2];++r){
           zz = z + r - Nhalf[2];
           if((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2])
               &&((MapPS4D->map[p][q][r][k]&prb_bit)==prb_bit) )
            { Xmap->map[xx][yy][zz] = Xmap->map[xx][yy][zz]|res_bit; 
              if (char_fit_ratio>Ymap->map[xx][yy][zz]) Ymap->map[xx][yy][zz] = char_fit_ratio;
            }
         } /* r */
        } /* q */
       } /* p */

       } /* if > fit_ratio */
      }
     }  /* k */ 
     }

    } /* foc_bitmask */

   } /* z */
  } /* y */
 } /* x */


 if (Nall_center>0) fit_ratio_ave /= Nall_center;

 printf("#Nsurvive_center %d\n",Nsurvive_center);
 printf("#Nall_center %d fit_ratio_ave %f fit_ratio_ave %f\n",Nall_center,fit_ratio_ave,fit_ratio_max);

 *Ave_fit_ratio = fit_ratio_ave; 
 *Max_fit_ratio = fit_ratio_max; 
 free(Nprobe_mark);
 return(Nsurvive_center);

} /* end of Opening_Using_CHAR4DMAP_and_GRIDPNT3D_with_Fitness() */








void Mark_Minimum_Filtering_GRIDPNT3Ds(MapPS4D,PntHead,Nrotate,bitmask)
 struct CHAR4DMAP *MapPS4D;
 struct GRIDPNT3D *PntHead;
 int    Nrotate;
 char   bitmask;
{
 int k,Nmark,Nmark_new;
 struct GRIDPNT3D *pnt; 
 char  *probe_mark,end;

 /** [1] Initialize **/
 probe_mark = (char *)malloc(sizeof(char)*(Nrotate+1));
 for (k=0;k<Nrotate;++k) probe_mark[k] = 0;
 pnt = PntHead;
 while (pnt->next != NULL){
  pnt = pnt->next;
  pnt->mark = 0;
 }

 /** [2] Mark GRIDPNT3D **/
 pnt = PntHead; Nmark = 0; end = 0;
 while ((pnt->next != NULL)&&(end==0)){
  pnt = pnt->next;
  Nmark_new = 0;
  for (k=0;k<Nrotate;++k){
   if ((MapPS4D->map[pnt->x][pnt->y][pnt->z][k]&bitmask)==bitmask){
    if (probe_mark[k]==0) {++Nmark; ++Nmark_new; probe_mark[k] = 1;}
   } 
  }
  if (Nmark_new>0) pnt->mark = 1; 

  if (Nmark >= Nrotate) end = 1;
 
 } /* pnt */

 /** [3] Remove Non-marked GRIDPNT3D **/
 /*
 pnt = PntHead;
 while (pnt->next!=NULL){
  pnt = pnt->next;
  if (pnt->mark==0){
    if (pnt->prev !=NULL) pnt->prev->next = pnt->next; 
    if (pnt->next !=NULL) pnt->next->prev = pnt->prev; 
   }

 }
 */

 /** [4] Output  **/
 pnt = PntHead;
 while (pnt->next!=NULL){
  pnt = pnt->next;
  if (pnt->mark==1) 
  printf("#xyz %2d %2d %2d Nprobe %3d mark %d\n",
   pnt->x,
   pnt->y,
   pnt->z,
   pnt->Nprobe,
   pnt->mark);
 }

 free(probe_mark);

} /* end of Mark_Minimum_Filtering_GRIDPNT3Ds() */

void Make_Another_GRID3DPNT_List_from_Marked_List(NewListHead,MarkListHead)
 struct GRIDPNT3D *NewListHead;
 struct GRIDPNT3D *MarkListHead;
{
 struct GRIDPNT3D *pnt,*new;

 NewListHead->next = NULL;
 pnt = MarkListHead;

 while (pnt->next !=NULL){
  pnt = pnt->next;
  if (pnt->mark == 1){
    new = Add_GRIDPNT3D_Stack(NewListHead,pnt->x,pnt->y,pnt->z);
    new->Nprobe = pnt->Nprobe;
    new->value  = pnt->value;
    new->mark   = pnt->mark;
   }
 }


} /* end of Make_Another_GRID3DPNT_List_from_Marked_List() */



int NonAllZero_Region_GRID3DPNT_Listed_Probes(Xmap,Phead,PmapN,tar_bitmask,res_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 struct GRIDPNT3D *Phead; /* list of probe 3D point */
 int    PmapN[3];  
 unsigned char   tar_bitmask; /* target bitmask (1 or 2 or 4 or 8 or 16 or 32 or 64 or 128)*/
 unsigned char   res_bitmask; /* result bitmask (1 or 2 or 4 or 8 or 16 or 32 or 64 or 128)*/
{
 int x,y,z,xx,yy,zz;
 int Nhalf[3];
 char allzero;
 int  Nsurvive;
 struct GRIDPNT3D *pnt;
 printf("#NonAllZero_GRID3DPNT_Listed_Probes(Xmap,Pmap tar_bitmask %d res_bitmask %d)\n",
  tar_bitmask,res_bitmask);
 
 Nhalf[0]  = (PmapN[0]-1)/2;
 Nhalf[1]  = (PmapN[1]-1)/2;
 Nhalf[2]  = (PmapN[2]-1)/2;
 printf("#Nhalf %d %d %d\n",Nhalf[0],Nhalf[1],Nhalf[2]);
 
 Nsurvive = 0;
 for (x=0;x<Xmap->N[0];++x){
  if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]);
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    allzero = 1;
    pnt = Phead; 
    while (pnt->next != NULL){
     pnt = pnt->next;
     xx = x + pnt->x - Nhalf[0];
     yy = y + pnt->y - Nhalf[1];
     zz = z + pnt->z - Nhalf[2];
     if ((xx>=0)&&(xx<Xmap->N[0])&&(yy>=0)&&(yy<Xmap->N[1])&&(zz>=0)&&(zz<Xmap->N[2]))
     { 
      if ((Xmap->map[xx][yy][zz]&tar_bitmask) == tar_bitmask) allzero = 0;
     }
    } /* pnt */

    if (allzero==0)
     { Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask;
       ++Nsurvive;}

   } /* z */
  } /* y */
 } /* x */
                                                                                                                  
 printf("#Nsurvive:%d\n",Nsurvive);
 return(Nsurvive);

} /* end of NonAllZeroRegion_GRID3DPNT_Listed_Probes() */




int Opening_for_Focused_Region(Xmap,Pmap,tar_bitmask,foc_bitmask,res_bitmask)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 struct CHAR3DMAP *Pmap; /* probe 3D map  */
 unsigned char   tar_bitmask; /* target bitmask [pocket] (1,2,4,8,16,32,64,128)*/
 unsigned char   foc_bitmask; /* focused region bitmask  (1,2,4,8,16,32,64,128)*/
 unsigned char   res_bitmask; /* result bitmask (1,2,4,8,16,32,64,128)*/
{
 int x,y,z,p,q,r,n,i,xx,yy,zz,Nmark;
 char allin;
 int Nhalf[3];
 int  Nsurvive_center,NpntPmap;
 int  *XX,*YY,*ZZ; /* array of xx, yy, and zz */
 printf("#Opening_for_Focused_Region(tar %d res %d)", tar_bitmask,res_bitmask);

 /*** Initialize ***/ 

 for (i=0;i<3;++i) Nhalf[i] = (Pmap->N[i]-1)/2;
 Nsurvive_center = 0;

 NpntPmap = Pmap->N[0]*Pmap->N[1]*Pmap->N[2];
 XX = (int *)malloc(sizeof(int)*NpntPmap);
 YY = (int *)malloc(sizeof(int)*NpntPmap);
 ZZ = (int *)malloc(sizeof(int)*NpntPmap);

 /* 
  printf("#Npnt %d ThreIn %.1f\n",NpntPmap,NpntPmapThreIn);
 */

 /*** Scan x,y,z,p,q,r ***/
 for (x=0;x<Xmap->N[0];++x){
 /* if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]); */
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    if ((Xmap->map[x][y][z]&foc_bitmask) == foc_bitmask){ 
     allin = 1;  Nmark = 0;
     for (p=0;p<Pmap->N[0];++p){ 
       xx = x + p - Nhalf[0];
      for (q=0;q<Pmap->N[1];++q){ 
       yy = y + q - Nhalf[1];
       for (r=0;r<Pmap->N[2];++r){ 
       zz = z + r - Nhalf[2];
        if (Pmap->map[p][q][r]!=0){
         if ((xx<0)||(xx>=Xmap->N[0])||(yy<0)||(yy>=Xmap->N[1])||(zz<0)||(zz>=Xmap->N[2])) 
           {  
            p = q = r = NpntPmap;  allin = 0;
            }
	  else if ((Xmap->map[xx][yy][zz]&tar_bitmask) == tar_bitmask)
	  { 
            XX[Nmark] = xx; YY[Nmark] = yy; ZZ[Nmark] = zz; ++Nmark; 
          }
          else {
            p = q = r = NpntPmap; allin = 0; 
          }
       	 } 
        } /* r */
       } /* q */
      } /* p */
     
    if (allin==1){ 
      ++Nsurvive_center;
      for (n=0;n<Nmark;++n){ 
       xx = XX[n]; yy = YY[n]; zz = ZZ[n];
       Xmap->map[xx][yy][zz] = Xmap->map[xx][yy][zz]|res_bitmask; }
      }
     } /* Xmap == 1 */
    } /* z */
  } /* y */
 } /* x */

 free(XX);
 free(YY);
 free(ZZ);
 printf("# Nsurvive_center %d\n",Nsurvive_center);
 return(Nsurvive_center);

} /* end of Opening_for_Focused_Region() */






struct GRIDPNT3D* Add_GRIDPNT3D_Stack(PntHead,x,y,z)
 struct GRIDPNT3D *PntHead;
 int x,y,z;
{
 struct GRIDPNT3D *pre;
 
 pre = PntHead->next;
 PntHead->next = (struct GRIDPNT3D *)malloc(sizeof(struct GRIDPNT3D));
 PntHead->next->prev = PntHead;
 PntHead->next->next = pre;
 if (pre != NULL)
  { pre->prev = PntHead->next;}
 PntHead->next->x = x; PntHead->next->y = y; PntHead->next->z = z;
 return(PntHead->next);

} /* end of Add_GRIDPNT3D_POINTER() */




int Malloc_CHAR4DMAP(M,x,y,z,w)
 struct CHAR4DMAP *M;
 int x,y,z,w;
{
 int i,j,k;
 double byte,Mbyte;
 byte = (double)x*y*z*w*(double)sizeof(char);
 Mbyte = byte/1024/1024;
 printf("#Malloc_CHAR4DMAP(%d %d %d %d) %.2lf byte -> %.2lf Mbyte\n",x,y,z,w,byte,Mbyte);
 if (Mbyte > 512) { printf("#ERROR:Memory size is too large\n"); exit(1); }
 M->N[0] = x; M->N[1] = y; M->N[2] = z; M->N[3] = w;

 M->map = (unsigned char ****)malloc(sizeof(unsigned char ***) * x);
 for (i=0;i<x;++i)
 {
  M->map[i] = (unsigned char ***)malloc(sizeof(unsigned char **) * y);
  for (j=0;j<y;++j){
   M->map[i][j] = (unsigned char **)malloc(sizeof(unsigned char *) * z);
   for (k=0;k<z;++k){
    M->map[i][j][k] = (unsigned char *)malloc(sizeof(unsigned char) * w);
   }
  } 
 }


  Initialize_CHAR4DMAP(M);
  return(1);

} /* end of Malloc_CHAR4DMAP() */


void Initialize_CHAR4DMAP(M)
 struct CHAR4DMAP *M;
{
 int i,j,k,m;
 for (i=0;i<M->N[0];++i)
  for (j=0;j<M->N[1];++j)
   for (k=0;k<M->N[2];++k)
    for (m=0;m<M->N[3];++m) M->map[i][j][k][m] = 0;

} /* end of Initialize_CHAR4DMAP() */



void Make_Rod_Probe(Pmap,Rdire,Rwidth,dVec,bitmask)
 struct CHAR3DMAP *Pmap; /* target 3D map */
 float  Rdire;   /* Radius for rod direcion */
 float  Rwidth;  /* Radius for width */
 float  dVec[3]; /* Rod direction vector */
 char   bitmask;
{
 int i,x,y,z;
 float P[3],O[3],PO[3],PY[3],dvec[3],dVecLen;
 float Ddire,Dnorm,Rdire_width,RRwidth;

 Rdire_width = Rdire + Rwidth;
 RRwidth = Rwidth * Rwidth;
 dVecLen = Vec_Length(dVec);

 for (i=0;i<3;++i){
  O[i] = Pmap->OrigPos[i] + Pmap->grid_width * (Pmap->N[i]-1)/2;
  dvec[i] = dVec[i]/dVecLen;
 }
 /*
 printf("O %f %f %f gw %f\n",O[0],O[1],O[2],Pmap->grid_width);
 printf("dVec %f %f %f dvec %f %f %f\n",dVec[0],dVec[1],dVec[2],dvec[0],dvec[1],dvec[2]);
 */ 
 for (x=0;x<Pmap->N[0];++x){
  P[0] = Pmap->OrigPos[0] + Pmap->grid_width * x;
  for (y=0;y<Pmap->N[1];++y){
   P[1] = Pmap->OrigPos[1] + Pmap->grid_width * y;
    for (z=0;z<Pmap->N[2];++z){
     P[2] = Pmap->OrigPos[2] + Pmap->grid_width * z;
     Sub_Vec(PO,P,O);
     Ddire = Dot_Prod(PO,dvec);
     for (i=0;i<3;++i) PY[i] = PO[i] - Ddire*dvec[i];
     Dnorm  = Vec_Length(PY);
     if (Ddire<0.0) Ddire = -Ddire;
     /* printf("P %f %f %f Ddire %f Dnorm %f\n",P[0],P[1],P[2],Ddire,Dnorm); */
          if ((Ddire<=Rdire)&&(Dnorm<=Rwidth)) Pmap->map[x][y][z] = 1;
     else if ((Ddire<=Rdire_width)&&
              (((Ddire-Rdire)*(Ddire-Rdire)+Dnorm*Dnorm)<=RRwidth)) 
               Pmap->map[x][y][z] = bitmask;
    else Pmap->map[x][y][z] = 0;

   } /* z */
  } /* y */
 } /* x */

} /* end of Make_Rod_Probe() */



void Make_Disc_Probe(Pmap,Rdire,Rwidth,nVec,bitmask)
 struct CHAR3DMAP *Pmap; /* target 3D map */
 float  Rdire;   /* Radius for disc  */
 float  Rwidth;  /* Radius for width */
 float  nVec[3]; /* Normal vector for disc */
 char bitmask;
{
 int i,x,y,z;
 float P[3],O[3],PO[3],PY[3],nvec[3],nVecLen;
 float Ddire,Dnorm,Rdire_width,RRwidth;

 Rdire_width = Rdire + Rwidth;
 RRwidth = Rwidth * Rwidth;
 nVecLen = Vec_Length(nVec);

 for (i=0;i<3;++i){
  O[i] = Pmap->OrigPos[i] + Pmap->grid_width * (Pmap->N[i]-1)/2;
  nvec[i] = nVec[i]/nVecLen;
 }
 /*
 printf("O %f %f %f gw %f\n",O[0],O[1],O[2],Pmap->grid_width);
 printf("nVec %f %f %f nvec %f %f %f\n",nVec[0],nVec[1],nVec[2],nvec[0],nvec[1],nvec[2]);
 */
 
 for (x=0;x<Pmap->N[0];++x){
  P[0] = Pmap->OrigPos[0] + Pmap->grid_width * x;
  for (y=0;y<Pmap->N[1];++y){
   P[1] = Pmap->OrigPos[1] + Pmap->grid_width * y;
    for (z=0;z<Pmap->N[2];++z){
     P[2] = Pmap->OrigPos[2] + Pmap->grid_width * z;
     Sub_Vec(PO,P,O);
     Dnorm = Dot_Prod(PO,nvec);
     for (i=0;i<3;++i) PY[i] = PO[i] - Dnorm*nvec[i];
     Ddire  = Vec_Length(PY);
     if (Dnorm<0.0) Dnorm = -Dnorm;
/*      printf("P %f %f %f Ddire %f Dnorm %f\n",P[0],P[1],P[2],Ddire,Dnorm);  */
  
         if ((Ddire<=Rdire)&&(Dnorm<=Rwidth)) Pmap->map[x][y][z] = 1;
     else if ((Ddire<=Rdire_width)&&
              (((Ddire-Rdire)*(Ddire-Rdire)+Dnorm*Dnorm)<=RRwidth)) Pmap->map[x][y][z] = bitmask;
     else Pmap->map[x][y][z] = 0;
   } /* z */
  } /* y */
 } /* x */
} /* end of Make_Disc_Probe() */








void Make_CHAR4DMAP_for_Rotated_RodDisc_by_RmatArray(MapPS4D,MapPS,RodDiscType,Rdire,Rwidth,Nrotate,RmatArray,bitmask)
 struct CHAR4DMAP *MapPS4D;
 struct CHAR3DMAP *MapPS;
 char   RodDiscType; /* 'R'od or 'D'isc */
 float  Rdire;   /* Radius for rod direcion */
 float  Rwidth;  /* Radius for width */
 int    Nrotate;
 struct ROT_MATRIX *RmatArray;
 char bitmask;
{
 int r,x,y,z;
 printf("#Make_CHAR4DMAP_for_Rotated_Rods_by_RmatArray(RodDisc '%c' Rdire %f Rwidth %f Nrotate %d)\n",
  RodDiscType,Rdire,Rwidth,Nrotate);

 for (r=0;r<Nrotate;++r){
       if (RodDiscType=='R') Make_Rod_Probe(MapPS,Rdire,Rwidth,RmatArray[r].R[0],bitmask);
  else if (RodDiscType=='D') Make_Disc_Probe(MapPS,Rdire,Rwidth,RmatArray[r].R[1],bitmask);
  else {printf("#ERROR:RodDiscType %c is improper\n",RodDiscType); exit(1);} 
  for (x=0;x<MapPS->N[0];++x){
   for (y=0;y<MapPS->N[1];++y){
    for (z=0;z<MapPS->N[2];++z){
     if ((MapPS->map[x][y][z]&bitmask)==bitmask) 
      MapPS4D->map[x][y][z][r] = MapPS4D->map[x][y][z][r]|bitmask;
    } /* z */
   } /* y */
  } /* x */
  Reset_Specific_3DMAP(MapPS,bitmask);
 } /* r */

} /* end of Make_CHAR4DMAP_for_Rotated_RodDisc_by_RmatArray() */
