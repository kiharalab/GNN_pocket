/*

 <TauOpnRotPrb.c>
 
  Functions for opening for probes with many 
  rotational variants 

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
#include "TauOpnRotPrb.h"
#include "GrPtMergeSort.h"

/*** FUNCTIONS (GLOBAL) ***/
int Tau_Erosion();
int Tau_Opening();
int Tau_Epsilon_Opening();
int Tau_Epsilon_Opening_for_Focused_Region();
int Number_of_Nonzero_Grids();




int Tau_Erosion(Xmap,Pmap,tar_bitmask,res_bitmask,tau)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 struct CHAR3DMAP *Pmap; /* probe 3D map */
 unsigned char   tar_bitmask; /* target bitmask (1 or 2 or 4 or 8 or 16 or 32 or 64 or 128)*/
 unsigned char   res_bitmask; /* result bitmask (1 or 2 or 4 or 8 or 16 or 32 or 64 or 128)*/
 float    tau;
{
 int x,y,z,p,q,r,NpntPmap;
 int xx,yy,zz;
 int Nhalf[3];
 char end;
 int  Nsurvive,Ntry,Nin,Nout;
 float Nin_thre,Nout_thre;
 printf("#Tau_Erosion(Xmap,Pmap tar_bitmask %d res_bitmask %d tau %f)\n",tar_bitmask,res_bitmask,tau);
                                                                                                                  
 Nhalf[0]  = (Pmap->N[0]-1)/2;
 Nhalf[1]  = (Pmap->N[1]-1)/2;
 Nhalf[2]  = (Pmap->N[2]-1)/2;
 NpntPmap  = Number_of_Nonzero_Grids(Pmap);
 Nin_thre  = tau * NpntPmap;
 Nout_thre = (1.0-tau) * NpntPmap;
 
 Nsurvive = Ntry = 0;
 for (x=0;x<Xmap->N[0];++x){
  if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]);
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    /* if ((Xmap->map[x][y][z]&tar_bitmask) == tar_bitmask) */
    {
     ++Ntry;
     end = 0; Nin = Nout = 0;
     p = 0; 
     while ((p<Pmap->N[0])&&(end==0)){
      q = 0;
      while ((q<Pmap->N[1])&&(end==0)){
       r = 0;
       while ((r<Pmap->N[2])&&(end==0)){
        if (Pmap->map[p][q][r]==1){
          xx = x + p - Nhalf[0];
          yy = y + q - Nhalf[1];
          zz = z + r - Nhalf[2];
          if ((xx<0)||(xx>=Xmap->N[0])||(yy<0)||(yy>=Xmap->N[1])||(zz<0)||(zz>=Xmap->N[2])){
             ++Nout; 
             /* printf("#Erosion:out of box !! xyz %d %d %d Xmap %d xx %d yy %d zz %d\n",
               x,y,z,Xmap->map[x][y][z],xx,yy,zz); 
              */
             }
          else{
           if ((Xmap->map[xx][yy][zz]&tar_bitmask) == tar_bitmask) ++Nin; else ++Nout;
           if (Nout>Nout_thre) end = 1; 
          }
                                                                                                                  
          } /* Pmap == 1 */
         ++r;
        } /* r */
        ++q;
       } /* q */
       ++p;
      } /* p */
                                                                                                                  
     if ((end==0)&&(Nin>=Nin_thre))
         { Xmap->map[x][y][z] = Xmap->map[x][y][z]|res_bitmask;
           ++Nsurvive;}
     } /* Xmap == 1 */
    } /* z */
  } /* y */
 } /* x */
                                                                                                                  
 printf("#tau_erosion Nsurvive:%d Ntry:%d\n",Nsurvive,Ntry);
 return(Nsurvive);
} /* end of Tau_Erosion() */



int Tau_Opening(Xmap,Pmap,tar_bitmask,res_bitmask,tau)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 struct CHAR3DMAP *Pmap; /* probe 3D map */
 unsigned char   tar_bitmask; /* target bitmask (1 or 2 or 4 or 8 or 16 or 32 or 64 or 128)*/
 unsigned char   res_bitmask; /* result bitmask (1 or 2 or 4 or 8 or 16 or 32 or 64 or 128)*/
 float tau;   /* |X and P| >= (1.0-tau)|P| */
{
 int x,y,z,p,q,r,n,i;
 int xx,yy,zz;
 int Nhalf[3],NpntPmap,Nin,Nout;
 float NpntPmapThreIn,NpntPmapThreOut;
 int  Nsurvive_center;
 int  *XX,*YY,*ZZ; /* array of xx, yy, and zz */
 printf("#Tau_Opening(Xmap,Pmap tar_bitmask %d res_bitmask %d tau %f)\n",tar_bitmask,res_bitmask,tau);

 /*** Initialize ***/ 

 for (i=0;i<3;++i) Nhalf[i] = (Pmap->N[i]-1)/2;
 Nsurvive_center = 0;
 
 NpntPmap = Number_of_Nonzero_Grids(Pmap);

 XX = (int *)malloc(sizeof(int)*NpntPmap);
 YY = (int *)malloc(sizeof(int)*NpntPmap);
 ZZ = (int *)malloc(sizeof(int)*NpntPmap);

 NpntPmapThreIn  = (1.0-tau) * NpntPmap;
 NpntPmapThreOut = tau * NpntPmap;
 printf("#NpntPmap %d NpntPmapThre %f\n",NpntPmap,NpntPmapThreIn);

 /*** Scan x,y,z,p,q,r ***/
 for (x=0;x<Xmap->N[0];++x){
  if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]);
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    Nin = Nout = 0;
    /* if ((Xmap->map[x][y][z]&tar_bitmask) == tar_bitmask)*/
    {
     for (p=0;p<Pmap->N[0];++p){ 
       xx = x + p - Nhalf[0];
      for (q=0;q<Pmap->N[1];++q){ 
       yy = y + q - Nhalf[1];
       for (r=0;r<Pmap->N[2];++r){ 
       zz = z + r - Nhalf[2];
        if (Pmap->map[p][q][r]!=0){
          if ((xx<0)||(xx>=Xmap->N[0])||(yy<0)||(yy>=Xmap->N[1])||(zz<0)||(zz>=Xmap->N[2])) 
            {  p = q = r = NpntPmap; }
          else {
            if ((Xmap->map[xx][yy][zz]&tar_bitmask) == tar_bitmask){
             XX[Nin] = xx; 
             YY[Nin] = yy; 
             ZZ[Nin] = zz; 
             ++Nin;
            }
            else  { 
             ++Nout;  
             if (Nout>NpntPmapThreOut) p = q = r = NpntPmap;
            }
          }
        } 
 
        } /* r */
       } /* q */
      } /* p */
   /* 
     printf("#Nneighbor %d NpntPmapThre %.1f\n",Nneighbor,NpntPmapThre); 
   */ 
     if (Nin>=NpntPmapThreIn){
      ++Nsurvive_center;
      for (n=0;n<Nin;++n){ 
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
 printf("#Tau_Opening Nsurvive_center %d\n",Nsurvive_center);
 return(Nsurvive_center);

} /* end of Tau_Opening() */



int Tau_Epsilon_Opening(Xmap,Pmap,pro_bitmask,tar_bitmask,res_bitmask,tau,eps)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 struct CHAR3DMAP *Pmap; /* probe 3D map  */
 unsigned char   pro_bitmask; /* protein bitmask (1,2,4,8,16,32,64,128)*/
 unsigned char   tar_bitmask; /* target bitmask [pocket] (1,2,4,8,16,32,64,128)*/
 unsigned char   res_bitmask; /* result bitmask (1,2,4,8,16,32,64,128)*/
 float tau;       /* |X              and P| < tau|P| */
 float eps;   /* |NonPocketSpace and P| < eps|P| */
{
 int x,y,z,p,q,r,n,i;
 int xx,yy,zz;
 int Nhalf[3],NpntPmap,Nin,Nout,Nwithin;
 float NpntPmapThreIn,NpntPmapThreWithin,NpntPmapThreOut;
 int  Nsurvive_center;
 int  *XX,*YY,*ZZ; /* array of xx, yy, and zz */
 printf("#Tau_Epsilon_Opening(tar %d res %d tau %f eps %f)",
  tar_bitmask,res_bitmask,tau,eps);

 /*** Initialize ***/ 

 for (i=0;i<3;++i) Nhalf[i] = (Pmap->N[i]-1)/2;
 Nsurvive_center = 0;

 NpntPmap = Number_of_Nonzero_Grids(Pmap);

 XX = (int *)malloc(sizeof(int)*NpntPmap);
 YY = (int *)malloc(sizeof(int)*NpntPmap);
 ZZ = (int *)malloc(sizeof(int)*NpntPmap);

 NpntPmapThreIn     = (1.0-tau-eps) * NpntPmap;
 NpntPmapThreWithin = tau * NpntPmap;
 NpntPmapThreOut    = eps * NpntPmap;
 printf("#NpntPmap %d Thre %.1f\n",NpntPmap,NpntPmapThreIn);

 /*** Scan x,y,z,p,q,r ***/
 for (x=0;x<Xmap->N[0];++x){
  if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]);
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    Nin = Nwithin = Nout = 0;

/*   if ((Xmap->map[x][y][z]&tar_bitmask) == tar_bitmask){ */
  if (1){ 

     for (p=0;p<Pmap->N[0];++p){ 
       xx = x + p - Nhalf[0];
      for (q=0;q<Pmap->N[1];++q){ 
       yy = y + q - Nhalf[1];
       for (r=0;r<Pmap->N[2];++r){ 
       zz = z + r - Nhalf[2];
        if (Pmap->map[p][q][r]!=0){
          if ((xx<0)||(xx>=Xmap->N[0])||(yy<0)||(yy>=Xmap->N[1])||(zz<0)||(zz>=Xmap->N[2])) 
            {  p = q = r = NpntPmap; }
          else if ((Xmap->map[xx][yy][zz]&tar_bitmask) == tar_bitmask){
             XX[Nin] = xx; 
             YY[Nin] = yy; 
             ZZ[Nin] = zz; 
             ++Nin;
          }
          else if ((Xmap->map[xx][yy][zz]&pro_bitmask) == pro_bitmask){
             ++Nwithin;  
              if (Nwithin>NpntPmapThreWithin) p = q = r = NpntPmap; 
          }
          else {
             ++Nout;  
              if (Nout>NpntPmapThreOut) p = q = r = NpntPmap; 
          }
        } 
        } /* r */
       } /* q */
      } /* p */
     
    
   if ((Nin>=NpntPmapThreIn)&&(Nwithin<=NpntPmapThreWithin)&&(Nout<=NpntPmapThreOut)){
   /* 
    printf("#Nin %d Nwithin %d Nout %d Ntotal %d ",Nin,Nwithin,Nout,Nin+Nwithin+Nout);
    if ((Nin>=NpntPmapThreIn)&&(Nwithin<=NpntPmapThreWithin)&&(Nout<=NpntPmapThreOut)) printf("HIT!!");
    printf("\n");
   */

      ++Nsurvive_center;
      for (n=0;n<Nin;++n){ 
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
 printf("#Tau_Epsilon_Opening Nsurvive_center %d\n",Nsurvive_center);
 return(Nsurvive_center);

} /* end of Tau_Epsilon_Opening() */




int Tau_Epsilon_Opening_for_Focused_Region(Xmap,Pmap,pro_bitmask,tar_bitmask,foc_bitmask,res_bitmask,tau,eps)
 struct CHAR3DMAP *Xmap; /* target 3D map */
 struct CHAR3DMAP *Pmap; /* probe 3D map  */
 unsigned char   pro_bitmask; /* protein bitmask (1,2,4,8,16,32,64,128)*/
 unsigned char   tar_bitmask; /* target bitmask [pocket] (1,2,4,8,16,32,64,128)*/
 unsigned char   foc_bitmask; /* focused region bitmask  (1,2,4,8,16,32,64,128)*/
 unsigned char   res_bitmask; /* result bitmask (1,2,4,8,16,32,64,128)*/
 float tau;       /* |X              and P| < tau|P| */
 float eps;   /* |NonPocketSpace and P| < eps|P| */
{
 int x,y,z,p,q,r,n,i;
 int xx,yy,zz;
 int Nhalf[3],NpntPmap,Nin,Nout,Nwithin,Nmark;
 float NpntPmapThreIn,NpntPmapThreWithin,NpntPmapThreOut;
 int  Nsurvive_center;
 int  *XX,*YY,*ZZ; /* array of xx, yy, and zz */
 printf("#Tau_Epsilon_Opening_for_Focused_Region(tar %d res %d tau %.2f eps %.2f)",
  tar_bitmask,res_bitmask,tau,eps);

 /*** Initialize ***/ 

 for (i=0;i<3;++i) Nhalf[i] = (Pmap->N[i]-1)/2;
 Nsurvive_center = 0;

 NpntPmap = Number_of_Nonzero_Grids(Pmap);

 XX = (int *)malloc(sizeof(int)*NpntPmap);
 YY = (int *)malloc(sizeof(int)*NpntPmap);
 ZZ = (int *)malloc(sizeof(int)*NpntPmap);

 NpntPmapThreIn     = (1.0-tau-eps) * NpntPmap;
 NpntPmapThreWithin = tau * NpntPmap;
 NpntPmapThreOut    = eps * NpntPmap;
 /* 
  printf("#Npnt %d ThreIn %.1f\n",NpntPmap,NpntPmapThreIn);
 */

 /*** Scan x,y,z,p,q,r ***/
 for (x=0;x<Xmap->N[0];++x){
 /* if ((x%10)==0) printf("#x %d/%d\n",x,Xmap->N[0]); */
  for (y=0;y<Xmap->N[1];++y){
   for (z=0;z<Xmap->N[2];++z){
    Nin = Nwithin = Nout = Nmark = 0;
    if ((Xmap->map[x][y][z]&foc_bitmask) == foc_bitmask){ 
     for (p=0;p<Pmap->N[0];++p){ 
       xx = x + p - Nhalf[0];
      for (q=0;q<Pmap->N[1];++q){ 
       yy = y + q - Nhalf[1];
       for (r=0;r<Pmap->N[2];++r){ 
       zz = z + r - Nhalf[2];
        if (Pmap->map[p][q][r]!=0){
         if ((xx<0)||(xx>=Xmap->N[0])||(yy<0)||(yy>=Xmap->N[1])||(zz<0)||(zz>=Xmap->N[2])) 
           {  
            p = q = r = NpntPmap; 
            }
	  else if ((Xmap->map[xx][yy][zz]&tar_bitmask) == tar_bitmask)
	   { 
            XX[Nmark] = xx; YY[Nmark] = yy; ZZ[Nmark] = zz; ++Nmark; 
	    ++Nin; 
	    }
          else if ((Xmap->map[xx][yy][zz]&pro_bitmask) == pro_bitmask){
             /* XX[Nmark] = xx; YY[Nmark] = yy; ZZ[Nmark] = zz; ++Nmark; */
             ++Nwithin;  
              if (Nwithin>NpntPmapThreWithin) p = q = r = NpntPmap; 
          }
          else {
             XX[Nmark] = xx; YY[Nmark] = yy; ZZ[Nmark] = zz; ++Nmark; 
             ++Nout;  
              if (Nout>NpntPmapThreOut) p = q = r = NpntPmap; 
          }
       	 } 
        } /* r */
       } /* q */
      } /* p */
     
    
   if ((Nin>=NpntPmapThreIn)&&(Nwithin<=NpntPmapThreWithin)&&(Nout<=NpntPmapThreOut)){
   /* 
    printf("#Nin %d Nwithin %d Nout %d Ntotal %d ",Nin,Nwithin,Nout,Nin+Nwithin+Nout);
    if ((Nin>=NpntPmapThreIn)&&(Nwithin<=NpntPmapThreWithin)&&(Nout<=NpntPmapThreOut)) printf("HIT!!");
    printf("\n");
   */

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

} /* end of Tau_Epsilon_Opening_for_Focused_Region() */








int Number_of_Nonzero_Grids(Pmap)
 struct CHAR3DMAP *Pmap; /* probe 3D map  */
{
 int p,q,r,Npntmap;

 Npntmap = 0; 
 for (p=0;p<Pmap->N[0];++p){
  for (q=0;q<Pmap->N[1];++q){
   for (r=0;r<Pmap->N[2];++r){
    if (Pmap->map[p][q][r]!=0) {++Npntmap;}
   }
  }
 }
 return(Npntmap);
} /* end of Number_of_Nonzero_Grids() */


