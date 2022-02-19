/*

 <Grid3D.c>

  functions for dealing with  struct CHAR3DMAP.


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
#include "Grid3D.h" 
#include "VecFunc.h" 
#include "PdbIO.h" 
#include "MscShPrb.h" 
#include "ReadOpt.h" 


/*** FUNCTIONS (GLOBAL) ***/
int  Malloc_CHAR3DMAP();
int  Malloc_INT3DMAP();
int  Malloc_VOXEL();
void Initialize_CHAR3DMAP();
void Initialize_INT3DMAP();
void Initialize_VOXEL();
void Free_CHAR3DMAP();
void Free_INT3DMAP();
void Free_VOXEL();
void Decide_Num_Grid_3D();
int  Set_Bit_3DMAP_Inside_Atoms();
void Set_Bit_By_Another_Map_in_CHAR3DMAP();
void Output_CHAR3DMAP_PDB_format();
void Output_CHAR3DMAP_with_tar_bitmask_PDB_format(); 
void Output_CHAR3DMAP_PDB_format_Min_Max();
void Output_CHAR3DMAP_PDB_format_with_MULSC_info();
void Output_CHAR3DMAP_PDB_format_Grouped_By_Char_Value();
void Output_Two_CHAR3DMAPs_in_PDB_format();
void MinMaxXYZ_Of_Atoms();
void Count_CHAR3DMAP_Insides_Overlap();
int  Assign_Map_Value_to_tFactor_of_Atoms();
void Read_CHAR3DMAP_PDB_format();
float Corr_Coeff_of_Two_CHAR3DMAPs();
void Output_Histogram_of_Two_CHAR3DMAPs();
void Make_Difference_CHAR3DMAP();
double Point_Correlation_Coefficient();
void Make_Difference_CHAR3DMAP_MultiScale();
void Resize_CHAR3DMAP_by_MinMax();
void Copy_CHAR3DMAP();
void Copy_CHAR3DMAP_with_Nshift();
void Compare_Two_CHAR3DMAP();
void Compare_Two_Bit_in_CHAR3DMAP();
void assign_CHAR3DMAP_from_VOXEL();
void Foreground_MinMax_CHAR3DMAPL();
void resize_CHAR3DMAP_from_MinMax_and_MarginWidth();
void copy_CHAR3DMAP_Same_GRID_and_Diff_OrigPos();


/*** FUNCTIONS (LOCAL) ***/
static int model_num_str_match();



int Malloc_CHAR3DMAP(M,x,y,z)
 struct CHAR3DMAP *M;
 int x,y,z;
{
 int i,j;
 double byte,Mbyte;
 
 byte = (double)x*y*z*(double)sizeof(char);
 Mbyte = byte/1024/1024;
 printf("#Malloc_CHAR3DMAP(%d %d %d) %.2lf byte -> %.2lf Mbyte\n",x,y,z,byte,Mbyte);
 if (Mbyte > PAR.MAX_MEMORY_MEGA) { printf("#ERROR:Memory size (%lf)is larger than MAX_MEMORY_MEGA(%f)\n",Mbyte, PAR.MAX_MEMORY_MEGA); exit(1); }
 
 M->N[0] = x; M->N[1] = y; M->N[2] = z; 
                                                                               
 M->map = (unsigned char ***)malloc(sizeof(unsigned char **) * x);
 for (i=0;i<x;++i){
   M->map[i] = (unsigned char **)malloc(sizeof(unsigned char *) * y);
   for (j=0;j<y;++j){
     M->map[i][j] = (unsigned char *)malloc(sizeof(unsigned char) * z);
   }
 }

 Initialize_CHAR3DMAP(M);
 M->malloced = 1; 
  return(1);
} /* end of Malloc_CHAR3DMAP() */




                                                                                                              
int Malloc_INT3DMAP(M,x,y,z)
 struct INT3DMAP *M;
 int x,y,z;
{
 int i,j;
 double byte,Mbyte;
 byte = (double)x*y*z*(double)sizeof(int);
 Mbyte = byte/1024/1024;
 printf("#Malloc_INT3DMAP(%d %d %d) %.2lf byte -> %.2lf Mbyte\n",x,y,z,byte,Mbyte);
 if (Mbyte > PAR.MAX_MEMORY_MEGA) { printf("#ERROR:Memory size (%lf)is larger than MAX_MEMORY_MEGA(%f)\n",Mbyte, PAR.MAX_MEMORY_MEGA); exit(1); }
 M->N[0] = x; M->N[1] = y; M->N[2] = z;
                                                                                                              
 M->map = (int ***)malloc(sizeof(int **) * x);
 for (i=0;i<x;++i){ 
   M->map[i] = (int **)malloc(sizeof(int *) * y);
   for (j=0;j<y;++j) M->map[i][j] = (int *)malloc(sizeof(int) * z);
 }
 M->malloced = 1;
  return(1);
} /* end of Malloc_INT3DMAP() */





int Malloc_VOXEL(vox,x,y,z)
 struct VOXEL *vox;
 int x,y,z;
{
 int i,j;
 double byte,Mbyte;

 byte = (double)x*y*z*(double)sizeof(float);
 Mbyte = byte/1024/1024;
/*
 *  printf("#Malloc_VOXEL(%d %d %d) %.2lf byte %.2lf Mbyte\n",x,y,z,byte,Mbyte);
 *  */
 if (Mbyte > PAR.MAX_MEMORY_MEGA){
   printf("#ERROR(Malloc_VOXEL):Memory size (%d %d %d) is too large (%.2lf Mbyte > MAX_MEMORY_MEGA %f Mbyte).\n",x,y,z,Mbyte,PAR.MAX_MEMORY_MEGA);
   return(0);
 }

 vox->dat = (float ***)malloc(sizeof(float **) * x);
 for (i=0;i<x;++i){
   vox->dat[i] = (float **)malloc(sizeof(float *) * y);
   for (j=0;j<y;++j){
     vox->dat[i][j] = (float *)malloc(sizeof(float) * z);
   }
 }

 vox->mal = 1;
 vox->N[0] = x; vox->N[1] = y; vox->N[2] = z;
 return(1);
} /* end of Malloc_VOXEL() */




void Initialize_CHAR3DMAP(M)
 struct CHAR3DMAP *M;
{
 int i,j,k;

 for (i=0;i<M->N[0];++i){
   for (j=0;j<M->N[1];++j){
     for (k=0;k<M->N[2];++k){ 
       M->map[i][j][k] = 0;
     }
   }
 }
} /* end of Initialize_CHAR3DMAP() */


void Initialize_INT3DMAP(M)
 struct INT3DMAP *M;
{
 int i,j,k;

 for (i=0;i<M->N[0];++i){
   for (j=0;j<M->N[1];++j){
     for (k=0;k<M->N[2];++k){ 
       M->map[i][j][k] = 0;
     }
   }
 }
} /* end of Initialize_INT3DMAP() */



void Initialize_VOXEL(M)
 struct VOXEL *M;
{
 int i,j,k;

 for (i=0;i<M->N[0];++i){
   for (j=0;j<M->N[1];++j){
     for (k=0;k<M->N[2];++k){ 
       M->dat[i][j][k] = 0;
     }
   }
 }
} /* end of Initialize_VOXEL() */





void Free_CHAR3DMAP(M)
 struct CHAR3DMAP *M;
{
 int i,j;
 printf("#Free_CHAR3DMAP(%d %d %d)\n",M->N[0],M->N[1],M->N[2]);
 for (i=M->N[0]-1;i>=0;--i){
  for (j=M->N[1]-1;j>=0;--j) free(M->map[i][j]);
  free(M->map[i]);
 }
 free(M->map);
} /* end of Free_CHAR3DMAP() */


void Free_INT3DMAP(M)
 struct CHAR3DMAP *M;
{
 int i,j;
 printf("#Free_INT3DMAP(%d %d %d)\n",M->N[0],M->N[1],M->N[2]);
 for (i=M->N[0]-1;i>=0;--i){
  for (j=M->N[1]-1;j>=0;--j) free(M->map[i][j]);
  free(M->map[i]);
 }
 free(M->map);
} /* end of Free_INT3DMAP() */



void Free_VOXEL(vox)
 struct VOXEL *vox;
{
 int i,j;

 for (i=0;i<vox->N[0];++i){
   for (j=0;j<vox->N[1];++j){free(vox->dat[i][j]);}
   free(vox->dat[i]);
 }

 free(vox->dat);
 vox->mal = 0;

} /* end of Free_VOXEL() */





void Decide_Num_Grid_3D(Ngrid,OrigPos,grid_width,margin_width,Min,Max,OddsEvenNgrid)
 int   Ngrid[3];        /* (to be calculated) */
 float OrigPos[3];      /* (to be calculated) */
 float grid_width;      /* (given) */
 float margin_width;    /* (given) */ 
 float Min[3],Max[3];   /* (given) */
 char  OddsEvenNgrid;   /* 'O':Ngrid should be odds number */
{
 int i;
 float Width[3];

 printf("#Decide_Num_Grid_3D(OddEven %c)\n",OddsEvenNgrid);
 if (OddsEvenNgrid == 'O'){
   for (i=0;i<3;++i){
     OrigPos[i] = Min[i] - margin_width;
     Width[i] = (Max[i] - Min[i])/2.0 + margin_width;
     Ngrid[i] = (int)ceil(Width[i]/grid_width) * 2 + 1;
    }
 }else{
   for (i=0;i<3;++i){
     OrigPos[i] = Min[i] - margin_width;
     Width[i] = Max[i]-Min[i] + 2*margin_width;
     Ngrid[i] = (int)ceil(Width[i]/grid_width);
    }
 }

 printf("#gwidth %f margin %f Width %f %f %f Ngrid %d %d %d OrigPos %f %f %f\n",
 grid_width, margin_width,Width[0],Width[1],Width[2],Ngrid[0],Ngrid[1],Ngrid[2], OrigPos[0],OrigPos[1],OrigPos[2]);
} /* end of Decide_Num_Grid_3D() */









int Set_Bit_3DMAP_Inside_Atoms(M,Ahead,res_bitmask,Rscale)
 struct CHAR3DMAP *M;
 struct ATOM *Ahead;
 unsigned char res_bitmask; /* result bitmask (1 or 2 or 4 or 8 or ... or 128)*/
 float  Rscale; /* scale for radius of atoms. The normal value is 1.0. */
{
 struct ATOM *an; 
 int i,j,k,Ninside,Natom;
 float min,max,DD,DifPos[3],R,RR;
 int Nmin[3],Nmax[3];

 /*
 printf("#Set_Bit_3DMAP_Inside_Atoms(M,Ahead,res_bitmask:%d Rscale:%f)\n",res_bitmask,Rscale);
 */
 for (i=0;i<3;++i){ 
   M->minN[i] = M->N[i];
   M->maxN[i] = 0;
 }

 Ninside = 0;
 an = Ahead; Natom = 0;
 while (an->next != NULL){
   an = an->next;
   Natom += 1;
   R = an->R * Rscale; 
   RR = R * R;

   for (i=0;i<3;++i){ 
     min = an->Pos[i] - R;
     max = an->Pos[i] + R; 
     Nmin[i] = (int)floor((min - M->OrigPos[i])/M->grid_width); 
     Nmax[i] = (int)ceil( (max - M->OrigPos[i])/M->grid_width); 
     if (Nmin[i]<0){        Nmin[i] = 0;} 
     if (Nmax[i]>=M->N[i]){ Nmax[i] = M->N[i]-1;}
   }
 /*
   printf("#R %f min %d %d %d max %d %d %d\n", an->R,Nmin[0], Nmin[1], Nmin[2], Nmax[0], Nmax[1], Nmax[2]);
 */
   for (i=Nmin[0];i<=Nmax[0];++i){
     DifPos[0] = M->OrigPos[0] + i*M->grid_width - an->Pos[0];
     for (j=Nmin[1];j<=Nmax[1];++j){ 
       DifPos[1] = M->OrigPos[1] + j*M->grid_width - an->Pos[1];
       for (k=Nmin[2];k<=Nmax[2];++k){
         DifPos[2] = M->OrigPos[2] + k*M->grid_width - an->Pos[2];
         if ((M->map[i][j][k]&res_bitmask) != res_bitmask){
           DD = DifPos[0]*DifPos[0] + DifPos[1]*DifPos[1] + DifPos[2]*DifPos[2];
           if (DD <= RR){
             Ninside += 1; 
             M->map[i][j][k] = M->map[i][j][k]|res_bitmask;
             if (i<M->minN[0]){ M->minN[0] = i;}
             if (i>M->maxN[0]){ M->maxN[0] = i;}
             if (j<M->minN[1]){ M->minN[1] = j;}
             if (j>M->maxN[1]){ M->maxN[1] = j;}
             if (k<M->minN[2]){ M->minN[2] = k;}
             if (k>M->maxN[2]){ M->maxN[2] = k;}
           }
         } 
      }  
     } 
   } 

 } 
 
 printf("#Set_Bit_3DMAP_Inside_Atoms(res_bitmask %d Rscale %f Natom %d) -->",res_bitmask,Rscale,Natom);
 printf("#Ngrid_inside %d\n",Ninside);
 return(Ninside);

} /* end of Set_Bit_3DMAP_Inside_Atoms() */



void Set_Bit_By_Another_Map_in_CHAR3DMAP(M,T,res_bitmask,threTval)
 struct CHAR3DMAP *M;   /* to be set bit */
 struct CHAR3DMAP *T;   /* another map. It has the same grid_width */
 unsigned char res_bitmask; /* result bitmask (1 or 2 or 4 or 8 or ... or 128)*/
 int    threTval;        /* if T->map[][][] > threTval, then it is foreground */
{
 int i, offset[3],x,y,z,xx,yy,zz;

 printf("#Set_Bit_By_Another_Map_in_CHAR3DMAP(res_bitmask %d threTval %d)\n",res_bitmask, threTval);
 if (M->grid_width != T->grid_width){
   printf("#Set_Bit_By_Another_Map_in_CHAR3DMAP(): Different grid_width (%f and %f)\n",M->grid_width, T->grid_width);
   exit(1); 
 }

 for (i=0;i<3;++i){
   offset[i] = (int)((T->OrigPos[i]-M->OrigPos[i])/T->grid_width);
 }
 printf("#M->OrigPos %f %f %f T->OrigPos %f %f %f\n",
   M->OrigPos[0], M->OrigPos[1], M->OrigPos[2], 
   T->OrigPos[0], T->OrigPos[1], T->OrigPos[2]);


 printf("#offset %d %d %d\n",offset[0], offset[1], offset[2]);


 for (x=0;x<T->N[0];++x){
   xx = x + offset[0];
   if ((xx>=0)&&(xx<M->N[0])){
     for (y=0;y<T->N[1];++y){
       yy = y + offset[1];
       if ((yy>=0)&&(yy<M->N[1])){
         for (z=0;z<T->N[2];++z){
           zz = z + offset[2];
           if ((zz>=0)&&(zz<M->N[2])){
             if (T->map[x][y][z] > threTval){
               M->map[xx][yy][zz] = M->map[xx][yy][zz]|res_bitmask;
             }
           }     
         }
       } 
     }
   }
 }



} /* end of Set_Bit_By_Another_Map_in_CHAR3DMAP() */



void Output_CHAR3DMAP_PDB_format(ofname,mode,M,min_val,title,CommentHead)
 char   *ofname;
 char   mode;  /* 'w'rite 'a'ppend */
 struct CHAR3DMAP *M;
 unsigned char   min_val;   /* if map[][][] >=min_value, then write  */
 char   *title; 
 struct STRING_NODE *CommentHead;  
{
 int i,j,k,v,Natom;
 FILE *fp;
 float Pnt[3],tFactor;
 char buff[64],core[128],AtomStr[5];
 struct STRING_NODE *sn;
 int Ngrid[256],Nsum;
 double unit_vol;
 
 printf("#Output_CHAR3DMAP_PDB_format(mode %c min_val %d) -->\"%s\" ",mode,min_val,ofname);
 if (mode=='a'){fp = fopen(ofname,"a");}
            else{ fp = fopen(ofname,"w");}
 if (fp==NULL) { printf("#ERROR:Can't write to \"%s\"\n",ofname); exit(1);}
 
 Find_Filename_Core(core,ofname);
 Get_Part_Of_Line(buff,ofname,0,29);

  Set_END_DATE();

 /** Output Basic Grid Information **/

 fprintf(fp,"HEADER    GRID     %-30s %s   %s\n",buff,Get_Date_String_PDB(),PAR.pdbidPro);
 if (title[0]!='\0') 
   fprintf(fp,"TITLE   %s\n",title);
 fprintf(fp,"REMARK  COMMAND %s\n",PAR.COMMAND);
 fprintf(fp,"REMARK  HOSTNAME %s\n",PAR.HOSTNAME);
 fprintf(fp,"REMARK  DATE    %s\n",Get_Date_String());
 fprintf(fp,"REMARK  DATE_START %s\n",PAR.START_DATE);
 fprintf(fp,"REMARK  DATE_END   %s\n",PAR.END_DATE);
 fprintf(fp,"REMARK  COMP_TIME  %lf seconds\n",PAR.COMP_TIME_SEC);

 fprintf(fp,"REMARK  grid_size %4d %4d %4d\n",M->N[0],M->N[1],M->N[2]);
 fprintf(fp,"REMARK  grid_width %8.3f\n",M->grid_width);
 fprintf(fp,"REMARK  OrigPos %8.3f %8.3f %8.3f\n",M->OrigPos[0],M->OrigPos[1],M->OrigPos[2]);
 unit_vol = M->grid_width * M->grid_width * M->grid_width;
 fprintf(fp,"REMARK  min_val %d\n",min_val);
 fprintf(fp,"REMARK  If (a grid with value less than [min_val]), do not write the grid.\n");
 
 /** Output Comment **/
 sn = CommentHead;
 while (sn->next != NULL){ 
   sn = sn->next;
   fprintf(fp,"REMARK  COMMENT %s\n",sn->line);
  }
 

 /** Count and Output Number_of_Grid for each map values **/
 for (v=0;v<256;++v){ Ngrid[v] = 0;}
 for (i=0;i<M->N[0];++i){
   for (j=0;j<M->N[1];++j){
     for (k=0;k<M->N[2];++k){ 
       if (M->map[i][j][k]>=min_val){ Ngrid[M->map[i][j][k]] += 1; }
     }
   }
 }

 Nsum = 0;
 fprintf(fp,  "REMARK  [CNgrd] Cumulative Number of Grids\n");
 fprintf(fp,  "REMARK  [CV]    Cumulative Volume\n");
 for (v=1;v<256;++v){ 
  if (Ngrid[v]>0){ 
   Nsum += Ngrid[v];
   fprintf(fp,  "REMARK  HISTOGRAM_NGRID VAL %3d Ngrd %8d V %10.3lf CNgrd %8d CV %10.3lf\n",
    v,Ngrid[v],Ngrid[v]*unit_vol,Nsum,Nsum*unit_vol);
  }
 }

 /** Output Grid Values **/ 
 fprintf(fp,"REMARK Column of [Residue Number] corresponds to the grid's char value\n");
 fprintf(fp,"REMARK              [grid_char_value]                 [radius][grid_char_value] [i] [j] [k]\n");

 Natom = 0;
 printf("#M->N %d %d %d\n",M->N[0],M->N[1],M->N[2]);
 for (i=0;i<M->N[0];++i){
   for (j=0;j<M->N[1];++j){
     for (k=0;k<M->N[2];++k){
       /* printf("#%d %d %d val %d\n",i,j,k,M->map[i][j][k]);  */
       if (M->map[i][j][k]>=min_val){  
         ++Natom; 
         Pnt[0] = M->OrigPos[0] + i*M->grid_width;
         Pnt[1] = M->OrigPos[1] + j*M->grid_width;
         Pnt[2] = M->OrigPos[2] + k*M->grid_width;
         /* sprintf(buff,"C%d",M->map[i][j][k]); */
         /* sprintf(buff,"X%d",M->map[i][j][k]);  */
         /* Get_Part_Of_Line(AtomStr,buff,0,2); */
         sprintf(AtomStr,"CA ");
          tFactor = (float)M->map[i][j][k];
         /*
         fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f %d %d %d\n",
            Natom%100000,AtomStr,"GRD",' ',M->map[i][j][k],
            Pnt[0],Pnt[1],Pnt[2],M->grid_width/2.0,tFactor,i,j,k);
         */
         fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f",
            Natom%100000,AtomStr,"GRD",' ',M->map[i][j][k],Pnt[0],Pnt[1],Pnt[2],
            M->grid_width/2.0,tFactor);

         fprintf(fp," %d %d %d",i,j,k);

         fprintf(fp,"\n");
      }
    
      } 
    } 
  } 

 fprintf(fp,"END\n");
 printf("# Natom %d\n",Natom);
 fclose(fp);

} /* end of Output_CHAR3DMAP_PDB_format() */





void Output_CHAR3DMAP_with_tar_bitmask_PDB_format(ofname,M,tar_bitmask)
 char   *ofname;
 struct CHAR3DMAP *M;
 unsigned char   tar_bitmask; /* target bitmask (1 or 2 or 4 or ... or 128)*/
{
 int i,j,k,Natom;
 FILE *fp;
 float Pnt[3],tFactor;
 char buff[64],AtomStr[5];
 
 printf("#Output_CHAR3DMAP_tar_bitmask_PDB_format(tar_bitmask %d) -->\"%s\" ",tar_bitmask,ofname);
 fp = fopen(ofname,"w");
 if (fp==NULL) { printf("#ERROR:Can't write to \"%s\"\n",ofname); exit(1);}
 
 Set_END_DATE();
 /** Output Basic Grid Information **/

 fprintf(fp,"HEADER    GRID     %-30s %s\n",buff,Get_Date_String_PDB());
 fprintf(fp,"REMARK  COMMAND %s\n",PAR.COMMAND);
 fprintf(fp,"REMARK  HOSTNAME %s\n",PAR.HOSTNAME);
 fprintf(fp,"REMARK  DATE    %s\n",Get_Date_String());
 fprintf(fp,"REMARK  DATE_START %s\n",PAR.START_DATE);
 fprintf(fp,"REMARK  DATE_END   %s\n",PAR.END_DATE);
 fprintf(fp,"REMARK  COMP_TIME  %lf seconds\n",PAR.COMP_TIME_SEC);

 fprintf(fp,"REMARK  grid_size %4d %4d %4d\n",M->N[0],M->N[1],M->N[2]);
 fprintf(fp,"REMARK  grid_width %8.3f\n",M->grid_width);
 fprintf(fp,"REMARK  OrigPos %8.3f %8.3f %8.3f\n",M->OrigPos[0],M->OrigPos[1],M->OrigPos[2]);
 fprintf(fp,"REMARK  tar_bitmask %d\n",tar_bitmask);
 fprintf(fp,"REMARK  If (a grid with value less than [min_val]), do not write the grid.\n");
 
 Natom = 0;
 for (i=0;i<M->N[0];++i){
   for (j=0;j<M->N[1];++j){
     for (k=0;k<M->N[2];++k){
       if ((M->map[i][j][k]&tar_bitmask) == tar_bitmask){
         ++Natom; 
         Pnt[0] = M->OrigPos[0] + i*M->grid_width;
         Pnt[1] = M->OrigPos[1] + j*M->grid_width;
         Pnt[2] = M->OrigPos[2] + k*M->grid_width;
         /* sprintf(buff,"C%d",M->map[i][j][k]); */
         sprintf(buff,"X%d",M->map[i][j][k]); 
         Get_Part_Of_Line(AtomStr,buff,0,2);
         tFactor = (float)M->map[i][j][k];
         /*
         fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f %d %d %d\n",
            Natom%100000,AtomStr,"GRD",' ',M->map[i][j][k],
            Pnt[0],Pnt[1],Pnt[2],M->grid_width/2.0,tFactor,i,j,k);
         */
         fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f",
            Natom%100000,AtomStr,"GRD",' ',M->map[i][j][k],Pnt[0],Pnt[1],Pnt[2],
            M->grid_width/2.0,tFactor);

         fprintf(fp," %d %d %d",i,j,k);

         fprintf(fp,"\n");
      }
    
      } 
    } 
  } 

 fprintf(fp,"END\n");
 printf("# Natom %d\n",Natom);
 fclose(fp);

} /* end of Output_CHAR3DMAP_with_tar_bitmask_PDB_format() */































void Output_CHAR3DMAP_PDB_format_Min_Max(ofname,mode,M,min_val,max_val,title,CommentHead)
 char   *ofname;
 char   mode;  /* 'w'rite 'a'ppend */
 struct CHAR3DMAP *M;
 unsigned char   min_val;/* if map[][][] >=min_value, then write                    */
 unsigned char   max_val;/* if map[][][] <=max_value, then write                    */
 char   *title; 
 struct STRING_NODE *CommentHead;  
{
 int i,j,k,v,Natom;
 FILE *fp;
 float Pnt[3],tFactor;
 char buff[64],core[128],AtomStr[5];
 struct STRING_NODE *sn;
 int Ngrid[256],Nsum;
 double unit_vol;
 
 printf("#Output_CHAR3DMAP_PDB_format_Min_Max(min %d max %d) -->\"%s\" ",min_val,max_val,ofname);
 if (mode=='a'){  fp = fopen(ofname,"a");}
            else{ fp = fopen(ofname,"w");}
 if (fp==NULL) { printf("#ERROR:Can't write to \"%s\"\n",ofname); exit(1);}
 
 Find_Filename_Core(core,ofname);
 Get_Part_Of_Line(buff,ofname,0,29);
 /** Output Basic Grid Information **/
 fprintf(fp,"HEADER    GRID     %-30s %s   %s\n",buff,Get_Date_String_PDB(),PAR.pdbidPro);
 if (title[0]!='\0') fprintf(fp,"TITLE   %s\n",title);
 fprintf(fp,"REMARK  COMMAND %s\n",PAR.COMMAND);
 fprintf(fp,"REMARK  DATE    %s\n",Get_Date_String());
 fprintf(fp,"REMARK  grid_size %4d %4d %4d\n",M->N[0],M->N[1],M->N[2]);
 fprintf(fp,"REMARK  grid_width %8.3f\n",M->grid_width);
 fprintf(fp,"REMARK  OrigPos %8.3f%8.3f %8.3f\n",M->OrigPos[0],M->OrigPos[1],M->OrigPos[2]);
 unit_vol = M->grid_width * M->grid_width * M->grid_width;
 fprintf(fp,"REMARK  min_val %d\n",min_val);
 fprintf(fp,"REMARK  If (a grid with value less than [min_val]), do not write the grid.\n");
 
 /** Output Comment **/
 sn = CommentHead;
 while (sn->next != NULL){ 
   sn = sn->next;
   fprintf(fp,"REMARK  COMMENT %s\n",sn->line); 
 }
 

 /** Count and Output Number_of_Grid for each map values **/
 for (v=0;v<256;++v) Ngrid[v] = 0;
 for (i=0;i<M->N[0];++i){
   for (j=0;j<M->N[1];++j){
     for (k=0;k<M->N[2];++k){ 
       if ((M->map[i][j][k]>=min_val)&&(M->map[i][j][k]<=max_val))
  	    Ngrid[M->map[i][j][k]] += 1; 
    }
   } 
 }
 Nsum = 0;
 fprintf(fp,  "REMARK  [CNgrd] Cumulative Number of Grids\n");
 fprintf(fp,  "REMARK  [CV]    Cumulative Volume\n");
 for (v=1;v<256;++v){ 
  if (Ngrid[v]>0){ 
   Nsum += Ngrid[v];
   fprintf(fp,  "REMARK  HISTOGRAM_NGRID VAL %3d Ngrd %8d V %10.3lf CNgrd %8d CV %10.3lf\n",
    v,Ngrid[v],Ngrid[v]*unit_vol,Nsum,Nsum*unit_vol);
  }
 }

 /** Output Grid Values **/ 
 fprintf(fp,"REMARK Column of [Residue Number] corresponds to the grid's char value\n");
 fprintf(fp,"REMARK              [grid_char_value]                 [radius][grid_char_value] [i] [j] [k]\n");

 Natom = 0;
 printf("#M->N %d %d %d\n",M->N[0],M->N[1],M->N[2]);

 for (i=0;i<M->N[0];++i){
   for (j=0;j<M->N[1];++j){
     for (k=0;k<M->N[2];++k){
       /* printf("#%d %d %d val %d\n",i,j,k,M->map[i][j][k]);  */
   if ((M->map[i][j][k]>=min_val)&&(M->map[i][j][k]<=max_val)){  
     ++Natom; 
     Pnt[0] = M->OrigPos[0] + i*M->grid_width;
     Pnt[1] = M->OrigPos[1] + j*M->grid_width;
     Pnt[2] = M->OrigPos[2] + k*M->grid_width;
     /* sprintf(buff,"C%d",M->map[i][j][k]); */
     sprintf(buff,"X%d",M->map[i][j][k]); 
     Get_Part_Of_Line(AtomStr,buff,0,2);
     tFactor = (float)M->map[i][j][k];
     fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f %d %d %d\n",
        Natom%100000,AtomStr,"GRD",' ',M->map[i][j][k],
        Pnt[0],Pnt[1],Pnt[2],M->grid_width/2.0,tFactor,i,j,k);
    }
    
     } 
   } 
 } 

 fprintf(fp,"END\n");
 printf("#Natom %d\n",Natom);
 fclose(fp);

} /* end of Output_CHAR3DMAP_PDB_format_Min_Max() */





void Output_CHAR3DMAP_PDB_format_with_MULSC_info(ofname,mode,M,min_val,max_val,title,MS,CommentHead)
 char   *ofname;
 char   mode;  /* 'w'rite 'a'ppend */
 struct CHAR3DMAP *M;
 unsigned char   min_val;   /* if map[][][] >=min_value, then write */
 unsigned char   max_val;   /* if map[][][] <=max_value, then write */
 char   *title; 
 struct  MULSC_SHELL_PRB *MS; 
 struct STRING_NODE *CommentHead;  
{
 int i,j,k,v,Natom;
 FILE *fp;
 float Pnt[3],tFactor;
 char buff[64],core[128],AtomStr[5];
 struct STRING_NODE *sn;
 int Ngrid[256],Nsum;
 double unit_vol;
 
 printf("#Output_CHAR3DMAP_PDB_format_with_MULSC_info(min %d max %d) -->\"%s\" ",min_val,max_val,ofname);
 if (mode=='a')  fp = fopen(ofname,"a");
            else fp = fopen(ofname,"w");
 if (fp==NULL) { printf("#ERROR:Can't write to \"%s\"\n",ofname); exit(1);}
 
 Find_Filename_Core(core,ofname);
 Get_Part_Of_Line(buff,ofname,0,29);
 /** Output Basic Grid Information **/
 fprintf(fp,"HEADER    GRID     %-30s %s   %s\n",buff,Get_Date_String_PDB(),PAR.pdbidPro);
 if (title[0]!='\0') fprintf(fp,"TITLE   %s\n",title);
 fprintf(fp,"REMARK  COMMAND %s\n",PAR.COMMAND);
 fprintf(fp,"REMARK  DATE    %s\n",Get_Date_String());
 fprintf(fp,"REMARK  grid_size %4d %4d %4d\n",M->N[0],M->N[1],M->N[2]);
 fprintf(fp,"REMARK  grid_width %8.3f\n",M->grid_width);
 fprintf(fp,"REMARK  OrigPos %8.3f%8.3f %8.3f\n",M->OrigPos[0],M->OrigPos[1],M->OrigPos[2]);
 unit_vol = M->grid_width * M->grid_width * M->grid_width;
 fprintf(fp,"REMARK  min_val %d max_val %d\n",min_val,max_val);
 fprintf(fp,"REMARK  If (a grid with value in [min_val<=value<=max_val]), do not write the grid.\n");
 
 /** Output Comment **/
 sn = CommentHead;
 while (sn->next != NULL)
 { sn = sn->next;
   fprintf(fp,"REMARK  COMMENT %s\n",sn->line); }
 
 /** Count and Output Number_of_Grid for each map values **/
 for (v=0;v<256;++v){ Ngrid[v] = 0;}
 for (i=0;i<M->N[0];++i){
  for (j=0;j<M->N[1];++j){
   for (k=0;k<M->N[2];++k){ 
    if (M->map[i][j][k]>=min_val){ Ngrid[M->map[i][j][k]] += 1; }
   }
  }
 }
 Nsum = 0;
 fprintf(fp,  "REMARK  [CNgrd] Cumulative Number of Grids\n");
 fprintf(fp,  "REMARK  [CV]    Cumulative Volume\n");
 for (v=1;v<256;++v){ 
  if (Ngrid[v]>0){ 
   Nsum += Ngrid[v];
   if (v<=MS->Nradius) tFactor = MS->Rarray[v]; else tFactor = 0.0; 
   fprintf(fp,  "REMARK  HISTOGRAM_NGRID VAL %3d R %6.2f Ngrd %8d V %10.3lf CNgrd %8d CV %10.3lf\n",
    v,tFactor,Ngrid[v],Ngrid[v]*unit_vol,Nsum,Nsum*unit_vol);
  }
 }

 /** Output Grid Values **/ 
 fprintf(fp,"REMARK Column of [Residue Number] corresponds to the grid's char value\n");
 fprintf(fp,"REMARK              [grid_char_value]                 [radius][R(grid)] [i] [j] [k]\n");
 Natom = 0;

 for (i=0;i<M->N[0];++i){
  for (j=0;j<M->N[1];++j){
   for (k=0;k<M->N[2];++k){
   
   if ((M->map[i][j][k]>=min_val) && (M->map[i][j][k]<=max_val)){  
     ++Natom; 
     Pnt[0] = M->OrigPos[0] + i*M->grid_width;
     Pnt[1] = M->OrigPos[1] + j*M->grid_width;
     Pnt[2] = M->OrigPos[2] + k*M->grid_width;
     /* sprintf(buff,"C%d",M->map[i][j][k]); */
     sprintf(buff,"X%d",M->map[i][j][k]); 
     Get_Part_Of_Line(AtomStr,buff,0,2);
     if ((M->map[i][j][k] >= 1) 
         && (M->map[i][j][k] <=MS->Nradius) ){
       tFactor = MS->Rarray[M->map[i][j][k]];
     }
     else tFactor = 0.0;
     /*
     if (M->map[i][j][k]>0) tFactor = 1.0/M->map[i][j][k];
                      else  tFactor = 0.0;
     tFactor = M->map[i][j][k];
     */

     /* This is mainly for dealing with Voxels of Protein VdW volume having 255 values */
     fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f",
        Natom%100000,AtomStr,"GRD",' ',M->map[i][j][k],
        Pnt[0],Pnt[1],Pnt[2],M->grid_width/2.0,tFactor);
     fprintf(fp," %d %d %d",i,j,k);
     fprintf(fp,"\n");
    }
    
    } 
   } 
  } 

 fprintf(fp,"END\n");
 printf("# Natom %d\n",Natom);
 fclose(fp);

} /* end of Output_CHAR3DMAP_PDB_format_with_MULSC_info() */












void Output_CHAR3DMAP_PDB_format_Grouped_By_Char_Value(ofname,mode,M,Ncluster, Nbit_cluster,title,CommentHead)
 char   *ofname;
 char   mode;  /* 'w'rite 'a'ppend */
 struct CHAR3DMAP *M;
 int    Ncluster;
 int    Nbit_cluster[MAX_UNSIGNED_CHAR]; /* Nbit_cluster[1..Ncluster] */
 char   *title;
 struct STRING_NODE *CommentHead;  
{
 int i,j,k,c,Natom;
 FILE *fp;
 float Pnt[3],Vvox,tFactor;
 char buff[128],AtomStr[5],core[128];
 struct STRING_NODE *sn;

 printf("#Output_CHAR3DMAP_PDB_format_Grouped_By_Char_Value() -->\"%s\" ",ofname);
 if (mode=='a')  fp = fopen(ofname,"a");
            else fp = fopen(ofname,"w");
 if (fp==NULL) { printf("#ERROR:Can't write to \"%s\"\n",ofname); exit(1);}

 Vvox = M->grid_width * M->grid_width * M->grid_width;
 
 Find_Filename_Core(core,ofname);
 Get_Part_Of_Line(buff,ofname,0,29);
 fprintf(fp,"HEADER    GRID     %-30s %s   %s\n",buff,Get_Date_String_PDB(),PAR.pdbidPro);
 if (title[0]!='\0') fprintf(fp,"TITLE   %s\n",title);
 fprintf(fp,"REMARK  COMMAND %s\n",PAR.COMMAND);
 fprintf(fp,"REMARK  DATE    %s\n",Get_Date_String());
 fprintf(fp,"REMARK  grid_size %4d %4d %4d\n",M->N[0],M->N[1],M->N[2]);
 fprintf(fp,"REMARK  grid_width %8.3f\n",M->grid_width);
 fprintf(fp,"REMARK  OrigPos %8.3f%8.3f %8.3f\n",M->OrigPos[0],M->OrigPos[1],M->OrigPos[2]);

 sn = CommentHead;
 while (sn->next != NULL){ 
   sn = sn->next;
   fprintf(fp,"REMARK  COMMENT %s\n",sn->line);
 }

 fprintf(fp,"REMARK  Ncluster %d\n",Ncluster);
 for (c=1;c<=Ncluster;++c){
  if (c<MAX_UNSIGNED_CHAR){
    fprintf(fp,"REMARK  CLUSTER %3d Ngrid %5d Volume %8.3f AAA\n",
      c,Nbit_cluster[c],Nbit_cluster[c]*Vvox);
    }
 }

 fprintf(fp,"REMARK                                              [Radius][MapVal] [i] [j] [k]\n");
 Natom = 0;

 for (c=1;c<=Ncluster;++c){

   tFactor = 10.0 - (float)c + 1.0; 
   if (tFactor < 0.0) tFactor = 0;

   for (i=0;i<M->N[0];++i){
     for (j=0;j<M->N[1];++j){
       for (k=0;k<M->N[2];++k){
         if (M->map[i][j][k]==c){  
           Natom += 1; 
           Pnt[0] = M->OrigPos[0] + i*M->grid_width;
           Pnt[1] = M->OrigPos[1] + j*M->grid_width;
           Pnt[2] = M->OrigPos[2] + k*M->grid_width;
           sprintf(buff," C%d",M->map[i][j][k]);
           Get_Part_Of_Line(AtomStr,buff,0,2);
           fprintf(fp,"HETATM%5d %-4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f",
              Natom%100000,AtomStr,"GRD",' ',M->map[i][j][k],
              Pnt[0],Pnt[1],Pnt[2],M->grid_width/2.0,tFactor);
           fprintf(fp," %d %d %d",i,j,k); 
           fprintf(fp,"\n");
         } 
       } 
     }
   } 
   fprintf(fp,"TER\n");
  } /* c */
 
 fprintf(fp,"END\n");
 printf("#Natom %d\n",Natom);
 fclose(fp);

} /* end of Output_CHAR3DMAP_PDB_format_Grouped_By_Char_Value() */



void Output_Two_CHAR3DMAPs_in_PDB_format(ofname,A,B,min_val,max_val,title)
 char   *ofname;
 struct CHAR3DMAP *A;
 struct CHAR3DMAP *B;
 unsigned char   min_val;   /* if map[][][] >=min_value, then write    */
 unsigned char   max_val;   /* if map[][][] <=max_value, then write    */
 char   *title; 
{
 int i,j,k,Natom;
 FILE *fp;
 float Pnt[3],tFactor;
 char buff[64],core[128],AtomStr[5],chain;
 unsigned char a,b;
 char chain_str[512];
 
 sprintf(chain_str,"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
 
 printf("#Output_Two_CHAR3DMAPs_PDB_format(min %d) -->\"%s\" ",min_val,ofname);
 if ((A->N[0]!=B->N[0])||(A->N[1]!=B->N[1])||(A->N[2]!=B->N[2])){
   printf("#ERROR:grid sizes are different.(%d %d %d) (%d %d %d)\n",
        A->N[0],A->N[1],A->N[2], B->N[0],B->N[1],B->N[2]);
   exit(1);
 }


 fp = fopen(ofname,"w");
 if (fp==NULL) { printf("#ERROR:Can't write to \"%s\"\n",ofname); exit(1);}
 
 Find_Filename_Core(core,ofname);
 Get_Part_Of_Line(buff,ofname,0,29);
 /** Output Basic Grid Information **/
 fprintf(fp,"HEADER    GRID     %-30s %s   %s\n",buff,Get_Date_String_PDB(),PAR.pdbidPro);
 if (title[0]!='\0') fprintf(fp,"TITLE   %s\n",title);
 fprintf(fp,"REMARK  COMMAND %s\n",PAR.COMMAND);
 fprintf(fp,"REMARK  DATE    %s\n",Get_Date_String());
 fprintf(fp,"REMARK  grid_size %4d %4d %4d\n",A->N[0],A->N[1],A->N[2]);
 fprintf(fp,"REMARK  grid_width %8.3f\n",A->grid_width);
 fprintf(fp,"REMARK  OrigPos_A %8.3f %8.3f %8.3f\n",A->OrigPos[0],A->OrigPos[1],A->OrigPos[2]);
 fprintf(fp,"REMARK  OrigPos_B %8.3f %8.3f %8.3f\n",B->OrigPos[0],B->OrigPos[1],B->OrigPos[2]);
 fprintf(fp,"REMARK  min_val %d\n",min_val);
 fprintf(fp,"REMARK  If (a grid with value less than [min_val]), do not write the grid.\n");
 fprintf(fp,"REMARK                                                            [x][y][z][valA][valB]\n");
 
 Natom = 0;
 for (i=0;i<A->N[0];++i){
  for (j=0;j<A->N[1];++j){
   for (k=0;k<A->N[2];++k){
     a = A->map[i][j][k];
     b = B->map[i][j][k];
/* 
  if ((a>=min_val)&&(a<=max_val)&&(b>=min_val)&&(b<=max_val) ){
 */
  if (((a==0)&&(b>0))||((b==0)&&(a>0))){
     printf("%d %d %d a %d b %d min_val %d max_val %d\n",i,j,k,A->map[i][j][k],B->map[i][j][k],min_val,max_val);  
     ++Natom; 
     Pnt[0] = A->OrigPos[0] + i*A->grid_width;
     Pnt[1] = A->OrigPos[1] + j*A->grid_width;
     Pnt[2] = A->OrigPos[2] + k*A->grid_width;
     /* sprintf(buff,"C%d",a); */
     sprintf(buff,"X%d",a); 
     Get_Part_Of_Line(AtomStr,buff,0,2);
     chain = chain_str[b%62];
     tFactor = (float)a - (float)b;
     fprintf(fp,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f %d %d %d %d %d\n",
        Natom%100000,AtomStr,"GRD",chain,a,
        Pnt[0],Pnt[1],Pnt[2],A->grid_width/2.0,tFactor,i,j,k,a,b);
    }
    
   }
  }
 }
 fprintf(fp,"END\n");
 printf("#Natom %d\n",Natom);
 fclose(fp);

} /* end of Output_Two_CHAR3DMAPs_in_PDB_format() */












void MinMaxXYZ_Of_Atoms(Ahead,Min,Max)
 struct ATOM *Ahead;
 float Min[3],Max[3];
{
 struct ATOM *an;
 char init;
 int i;
 float minpos,maxpos; 

 an = Ahead;
 init = 1;
 
 while (an->next != NULL){
  an = an->next;
  for (i=0;i<3;++i){ 
    minpos = an->Pos[i] - an->R;
    maxpos = an->Pos[i] + an->R; 
    if ((init==1)||(minpos<Min[i])) Min[i] = minpos;
    if ((init==1)||(maxpos>Max[i])) Max[i] = maxpos;
  }
  if (init==1) init = 0;

 } /* an */

} /* end of MinMaxXYZ_Of_Atoms(Ahead,Min,Max) */



void Count_CHAR3DMAP_Insides_Overlap(M,inside1,inside2,overlap,Nin1,Nin2,Nover,Vin1,Vin2,Vover)
 struct CHAR3DMAP *M;
 int inside1,inside2,overlap;   /* number of M->map[][][]   */
 int   *Nin1, *Nin2, *Nover;    /* Counted number in M-Map */
 float *Vin1, *Vin2, *Vover;    /* Volume in M-Map */
{
 int i,j,k; 
 float Vvox;

 *Nin1 = *Nin2 = *Nover = 0;

 for (i=0;i<M->N[0];++i)
  for (j=0;j<M->N[1];++j)
   for (k=0;k<M->N[2];++k){
     if (M->map[i][j][k] == inside1) *Nin1  += 1;
     if (M->map[i][j][k] == inside2) *Nin2  += 1;
     if (M->map[i][j][k] == overlap) *Nover += 1;
   }  

 Vvox = M->grid_width * M->grid_width * M->grid_width;

 *Vin1  = Vvox * (*Nin1); 
 *Vin2  = Vvox * (*Nin2); 
 *Vover = Vvox * (*Nover); 

} /* end of Count_CHAR3DMAP_Insides_Overlap() */



int Assign_Map_Value_to_tFactor_of_Atoms(AtomHead,M,Type,tar_bitmask)
 struct ATOM *AtomHead;
 /* Map char value will be written in tFactor of each atom. */
 struct CHAR3DMAP *M;
 char   Type;   /* 'B'it specified, 'C'luster_num of pockets */
 unsigned char   tar_bitmask;
 /*  
   case 'B' : [ Bit specfied ] 
     if a skin lattice around the atom has a tar_bitmask, tFactor = 1; otherwise tFactor = 0;  

   case 'C' : [ cluster number of pockets ] 
     if a skin lattice around the atom has a non zero value, 
     find the smallest number. ( cluster num is sorted by incresing order for its size.)
     If atoms contact with more than one cluster, larger cluster has a priority. 
 */ 
{
 struct ATOM *an;
 float D[3],DD[3],ddis,RR,RR1;
 int   minN[3],maxN[3];
 int   i,x,y,z,Nskin,Npocket,Ninside,Nmargrid, Natom_marked;
 
 printf("#Assign_Map_Value_to_tFactor_of_Atoms()\n");
 Nmargrid = 1;


 Natom_marked = 0;
 an = AtomHead;
 while ( an->next != NULL){
  an = an->next;
  if (Type == 'B') an->tFactor = 0.0;
  if (Type == 'C') an->tFactor = 255.0;
  
  Nskin = Ninside = Npocket = 0;
                                                                                                          
  for (i=0;i<3;++i){  
      minN[i] = (int)floor((an->Pos[i] - an->R - M->OrigPos[i])/M->grid_width) - Nmargrid;
      maxN[i] = (int)ceil( (an->Pos[i] + an->R - M->OrigPos[i])/M->grid_width) + Nmargrid;
      if (minN[i]<0) minN[i] = 0;
      if (maxN[i]>=M->N[i]) maxN[i] = M->N[i]-1;
  }
                                                                                                          
  RR  = an->R*an->R;
  RR1 = (an->R + Nmargrid * M->grid_width)*(an->R + Nmargrid * M->grid_width);
  for (x=minN[0];x<=maxN[0];++x){
   D[0] = M->OrigPos[0] + M->grid_width * x - an->Pos[0];
   DD[0] = D[0]*D[0];
   for (y=minN[1];y<=maxN[1];++y){
     D[1] = M->OrigPos[1] + M->grid_width * y - an->Pos[1];
     DD[1] = D[1]*D[1];
     for (z=minN[2];z<=maxN[2];++z){
       D[2] = M->OrigPos[2] + M->grid_width * z - an->Pos[2];
       ddis = DD[0] + DD[1] + D[2]*D[2];
                                                                                                          
       if (ddis<=RR) ++Ninside;
       if ((ddis>RR)&&(ddis<=RR1)){
         /* Skin */
         ++Nskin;
        /*
         printf("%d %d %d map %d bitmask %d mask %d\n",x,y,z,M->map[x][y][z],tar_bitmask,M->map[x][y][z]&tar_bitmask); 
       */

         if (Type == 'B'){
           if ((M->map[x][y][z]&tar_bitmask) == tar_bitmask) 
            {an->tFactor = 1.0; ++Npocket;}
          }
         if (Type == 'C'){
           if ((M->map[x][y][z]>0) && (M->map[x][y][z]< an->tFactor)) 
              { an->tFactor = M->map[x][y][z]; ++Npocket;}
          }
        }
     } /* z */
    } /* y */
  } /* x */
 
  if ((Type == 'C')&&(Npocket==0)) an->tFactor = 0.0;

  if (Npocket>0) Natom_marked += 1;
  /*
  printf(">an num %s atom %s resi %s R %f Nskin %d Ninside %d\n",an->Anum,an->Atom,an->Resi,an->R,Nskin,Ninside);
  printf(" minN %d %d %d maxN %d %d %d RR %f RR1 %f tFactor %f\n",
  minN[0],minN[1],minN[2],maxN[0],maxN[1],maxN[2],RR,RR1,an->tFactor);
 */
                                                                                                          
  } /* an */
  
  return(Natom_marked);                                                                                                       
} /* end of Assign_Map_Value_to_tFactor_of_Atoms() */




void Read_CHAR3DMAP_PDB_format(ifname,M, ChainSelect, model_num_str)
 char   *ifname;
 struct CHAR3DMAP *M;
 char   ChainSelect; /* '-': accept everything, '1':only accept 1, '2':accpet 1 and 2, '3': accept 1,2,3 */
 char   *model_num_str;   /* model number. if -1, not selecting  '1':model_num=1 ,'<=2' model_num <=2 */
{
 int i,j,k,L,Ngrid_read,model_num;
 FILE *fp;
 char buff[128];
 char line[256],subline[256];
 int Nword,Wsta[100],Wend[100];
 char malloc_ok,accept,ch;
 float val;
/*
>> FILE FORMAT EXAMPLE <<
HEADER    mclos.gdb                               26-JUL-06   mclos     0XXX   1
REMARK  COMMAND MathMorph 1qjaA.pdb -M M -rs 4 -rl 12 -br 4
REMARK  DATE    Jul 26,2006 12:51:3
REMARK  grid_size  102  101  115
REMARK  grid_width    1.000
REMARK  OrigPos  -36.173 -58.631  -11.357
REMARK  COMMENT Closing Shape with MultiScale Probe Sphere
REMARK  COMMENT pdbfile 1qjaA.pdb Natom 1745 VdW volume: 20610 grids, 20610.00 AAA
REMARK                                                [Radius][MapVal] [i] [j] [k]
HETATM    1   C1 GRD     1      -9.173  -4.631  43.643  0.50  1.00 27 54 55
HETATM    2   C1 GRD     1      -9.173  -4.631  48.643  0.50  1.00 27 54 60
HETATM    3   C1 GRD     1      -9.173  -3.631  42.643  0.50  1.00 27 55 54
HETATM    4   C1 GRD     1      -9.173  -3.631  43.643  0.50  1.00 27 55 55
:
HETATM36881  C11 GRD   116     223.178 138.356 107.227  0.50116.00 258 148 121
HETATM36882  C11 GRD   116     223.178 138.356 108.227  0.50116.00 258 148 122
HETATM36883  C11 GRD   116     223.178 139.356 107.227  0.50116.00 258 149 121
TER
END
*/

 printf("#Read_CHAR3DMAP_PDB_format(\"%s\" ChainSelect:'%c' select_model_num:%s)\n",ifname,ChainSelect,model_num_str);
 fp = fopen(ifname,"r");
 if (fp==NULL) { printf("#ERROR:Can't read pdbfile \"%s\"\n",ifname); exit(1);}
 malloc_ok = 0;
 Ngrid_read = 0;
 model_num = 0;

 while (feof(fp)==0){
   fgets(line,255,fp);
   L = strlen(line);
   if ((L>0)&&(line[L-1]=='\n')) {line[L-1] = '\0'; L = strlen(line);}
   /* printf("%s\n",line); */
   
   /** Read grid_size and malloc and initialize **/
   if (strncmp(line,"REMARK  grid_size",17) == 0){
    /* REMARK  grid_size  102  101  115 */
     Get_Part_Of_Line(subline,line,17,L-1);
     Split_to_Words(subline,' ',&Nword,Wsta,Wend,100);
     Get_Part_Of_Line(buff,subline,Wsta[0],Wend[0]); M->N[0] = atoi(buff);
     Get_Part_Of_Line(buff,subline,Wsta[1],Wend[1]); M->N[1] = atoi(buff);
     Get_Part_Of_Line(buff,subline,Wsta[2],Wend[2]); M->N[2] = atoi(buff);
     printf("#gird_size %d %d %d\n",M->N[0],M->N[1],M->N[2]);
     malloc_ok = Malloc_CHAR3DMAP(M,M->N[0],M->N[1],M->N[2]);
     for (k=0;k<3;++k){ M->minN[k] = M->N[k]; M->maxN[k] = 0;}
   }
 
   if (strncmp(line,"REMARK  grid_width",18) == 0){
    /* REMARK  grid_width    1.000*/
    Get_Part_Of_Line(subline,line,18,L-1);
    M->grid_width = atof(subline); 
    printf("#grid_width %f\n",M->grid_width);
   } 
   if (strncmp(line,"REMARK  OrigPos",15) == 0){
    /* REMARK  OrigPos  -36.173 -58.631  -11.357*/
    Get_Part_Of_Line(subline,line,15,L-1);
    Split_to_Words(subline,' ',&Nword,Wsta,Wend,100);
    Get_Part_Of_Line(buff,subline,Wsta[0],Wend[0]); M->OrigPos[0] = atof(buff);
    Get_Part_Of_Line(buff,subline,Wsta[1],Wend[1]); M->OrigPos[1] = atof(buff);
    Get_Part_Of_Line(buff,subline,Wsta[2],Wend[2]); M->OrigPos[2] = atof(buff);
    printf("#OrigPos %f %f %f\n",M->OrigPos[0],M->OrigPos[1],M->OrigPos[2]);
     
   }
 
      if (strncmp(line,"MODEL",5)==0){
/*
            1
            012345678901234
            MODEL        1
            MODEL        2
            MODEL       15
            */
      Get_Part_Of_Line(buff,line,5,13);
      model_num = atoi(buff);
   }





 
   if ((strncmp(line,"HETATM",6) == 0) && (malloc_ok==1)){
   /*
             10        20        30        40        50        60        70
   01234567890123456789012345678901234567890123456789012345678901234567890123456789
   HETATM    1   C1 GRD     1      -9.173  -4.631  43.643  0.50  1.00 27 54 55
   HETATM    2   C1 GRD     1      -9.173  -4.631  48.643  0.50  1.00 27 54 60
   HETATM    3   C1 GRD     1      -9.173  -3.631  42.643  0.50  1.00 27 55 54
   :
   HETATM36881  C11 GRD   116     223.178 138.356 107.227  0.50116.00 258 148 121
   :
   REMARK  [ChainID:ClusNum],[ResNum:char value]       [Radius][Rinacc][x][y][z]
             10        20        30        40        50        60        70
   01234567890123456789012345678901234567890123456789012345678901234567890123456789
   HETATM    2  CE  GRD 1   3      16.549  55.162  53.883  0.40  3.00 63 65 62
   HETATM    3  CE  GRD 1   3      15.749  58.362  56.283  0.40  3.00 62 69 65
   HETATM    4  CE  GRD 1   3      15.749  57.562  57.083  0.40  3.00 62 68 66
   HETATM    5  CE  GRD 1   3      14.949  59.162  56.283  0.40  3.00 61 70 65
   HETATM  371  CE  GRD 2   9      16.549  61.562  65.883  0.40  6.00 63 73 77
   HETATM  372  CE  GRD 2   5      16.549  60.762  65.883  0.40  4.00 63 72 77
   HETATM  373  CE  GRD 2   5      16.549  59.162  65.083  0.40  4.00 63 70 76
   HETATM  374  CE  GRD 2   8      16.549  58.362  65.083  0.40  5.50 63 69 76
   HETATM  375  CE  GRD 2   9      15.749  61.562  66.683  0.40  6.00 62 73 78
   HETATM  376  CE  GRD 2   9      15.749  61.562  65.883  0.40  6.00 62 73 77
   :
   */
 
    Get_Part_Of_Line(buff,line,17,19);
    if (strcmp(buff,"GRD")==0){
      accept = 0;
      ch = line[21]; 

      if (ChainSelect=='-'){ accept = 1;}
      if ((ChainSelect=='1')&&(ch=='1')){ accept = 1;}
      if ((ChainSelect=='2')&&((ch=='1')||(ch=='2'))){ accept = 1;}
      if ((ChainSelect=='3')&&((ch=='1')||(ch=='2')||(ch=='3'))){ accept = 1;}
      if ((ChainSelect=='4')&&((ch=='1')||(ch=='2')||(ch=='3')||(ch=='4'))){ accept = 1;}
      if ((ChainSelect=='5')&&((ch=='1')||(ch=='2')||(ch=='3')||(ch=='4')||(ch=='5'))){ accept = 1;}

      if (model_num_str_match(model_num_str, model_num)==0){
         accept = 0;
      }


      if (accept==1){
        Get_Part_Of_Line(subline,line,66,L-1);
        Split_to_Words(subline,' ',&Nword,Wsta,Wend,100);
        Get_Part_Of_Line(buff,subline,Wsta[0],Wend[0]); i = atoi(buff);
        Get_Part_Of_Line(buff,subline,Wsta[1],Wend[1]); j = atoi(buff);
        Get_Part_Of_Line(buff,subline,Wsta[2],Wend[2]); k = atoi(buff);
        Get_Part_Of_Line(buff,line,22,25); val = atof(buff);
        M->map[i][j][k] = (int)val;
    /*   
        printf("#%s:%d %d %d val %f M->map %d\n",line,i,j,k,val,M->map[i][j][k]);  fflush(stdout); 
    */
        Ngrid_read += 1;
        if (i>M->maxN[0]){ M->maxN[0] = i;}
        if (j>M->maxN[1]){ M->maxN[1] = j;}
        if (k>M->maxN[2]){ M->maxN[2] = k;}
        if (i<M->minN[0]){ M->minN[0] = i;}
        if (j<M->minN[1]){ M->minN[1] = j;}
        if (k<M->minN[2]){ M->minN[2] = k;}
      }
     }
   }
 
  
 } /* while */

 printf("#Ngrid_read %d\n",Ngrid_read);
 fclose(fp);

} /* end of Read_CHAR3DMAP_PDB_format() */


int model_num_str_match(model_num_str, model_num)
  char *model_num_str;
  int  model_num;
{
  char buff[MAX_FILE_NAME]; 
  int L,m;
  
  L = strlen(model_num_str);

  if (model_num_str[0]=='-'){
     return(1);
  }

  else if ((model_num_str[0]=='l')&&(model_num_str[1]=='e')){
     Get_Part_Of_Line(buff,model_num_str,2,L-1);
     m = atoi(buff);
     if (model_num <= m){
        return(1);
     } 
     else{
      return(0);
    }
  }
  
  else{
   m = atoi(model_num_str);
   if (model_num == m){
      return(1);
   } 
   else{
    return(0);
  }

  } 

} /* end of model_num_str_match() */



float Corr_Coeff_of_Two_CHAR3DMAPs(A,B)
 struct CHAR3DMAP *A,*B;
{
 int x,y,z,a,b,N;
 float Sa,Sb,Saa,Sbb,Sab;
 float Ma,Mb,Va,Vb,Cov,CC; 
 N = 0;
 Sa = Sb =  0.0;
 Saa = Sbb = Sab = 0.0;
 for (x=0;x<A->N[0];++x){
    for (y=0;y<A->N[1];++y){
      for (z=0;z<A->N[2];++z){
        a = A->map[x][y][z]; 
        b = B->map[x][y][z];
        if ((a>0)&&(a<255)&&(b>0)&&(b<255)){
          N += 1;
          Sa  += (float)a;
          Sb  += (float)b;
          Saa  += (float)(a*a);
          Sbb  += (float)(b*b);
          Sab += (float)(a*b);
        }
      }
    }
  }

 Ma = Sa/N; 
 Mb = Sb/N; 
 Va = Saa/N - Ma*Ma;
 Vb = Sbb/N - Mb*Mb;
 Cov = Sab/N - Ma * Mb;

 if ((Va>0.0) && (Vb>0.0)) CC = Cov/sqrt(Va)/sqrt(Vb);  else CC = 0.0;
 printf("#Ma %f Mb %f Va %f Vb %f Cov %f CC %f\n",Ma,Mb,Va,Vb,Cov,CC);
 return(CC);
} /* end of Corr_Coeff_of_Two_CHAR3DMAPs() */


void Output_Histogram_of_Two_CHAR3DMAPs(ofname,A,B)
 char   *ofname;
 struct CHAR3DMAP *A,*B;
{
  FILE *fpo;
  int Na[256],Nb[256],Nab[256][256];
  int a,b,x,y,z;
  printf("#Output_Histogram_of_Two_CHAR3DMAPs()-->'%s'\n",ofname);
  if ((A->N[0]!=B->N[0])||(A->N[1]!=B->N[1])||(A->N[2]!=B->N[2])){
    printf("#ERROR:Grid sizes are different. (%d %d %d)-(%d %d %d)\n",
    A->N[0], A->N[1], A->N[2], B->N[0], B->N[1], B->N[2]);
    exit(1);
  }
	  
  /** [1] Initialize **/
  for (a=0;a<255;++a){
    Na[a] = 0;
    for (b=0;b<255;++b){
      Nb[b] = 0;
      Nab[a][b] = 0;
    }
  }
  /** [2] Count **/
  for (x=0;x<A->N[0];++x){
    for (y=0;y<A->N[1];++y){
      for (z=0;z<A->N[2];++z){
        a = A->map[x][y][z]; 
        b = B->map[x][y][z];
        if ((a>=0)&&(a<255)&&(b>=0)&&(b<255)){
          Na[a] += 1;
          Nb[b] += 1;
          Nab[a][b] += 1;
        }
      }
    }
  }
  /** [3] Output to file **/
  fpo = fopen(ofname,"w");
  fprintf(fpo,"#HISTOGRAM FOR mapA\n");
  fprintf(fpo,"#COMMAND %s\n",PAR.COMMAND);
  fprintf(fpo,"#DATE    %s\n",Get_Date_String());
  for (a=0;a<255;++a){
    if (Na[a]>0) fprintf(fpo,"#A %d count %d\n",a,Na[a]);
  }
  fprintf(fpo,"#HISTOGRAM FOR mapB\n");
  for (b=0;b<255;++b){
    if (Nb[b]>0) fprintf(fpo,"#B %d count %d\n",b,Nb[b]);
  }

  for (a=0;a<255;++a){
    for (b=0;b<255;++b){
      if ((Na[a]>0)&&(Nb[b]>0)) fprintf(fpo,"%d %d %d\n",a,b,Nab[a][b]);
    }
    fprintf(fpo,"\n");
  }
  fclose(fpo);
} /* end of Output_Histogram_of_Two_CHAR3DMAPs() */
                                                                                                             
     

                                                                                                        
void Make_Difference_CHAR3DMAP(Amap,Bmap,Nab)
 struct CHAR3DMAP *Amap; /* target 3D map */
 struct CHAR3DMAP *Bmap; /* reference 3D map */
 double Nab[2][2]; 
{
 int x,y,z,vA,vB,a,b;
 double CC;

 printf("#Make_Difference_CHAR3DMAP(Xmap,Pmap)\n");
 printf("#Amap %d %d %d\n",Amap->N[0],Amap->N[1],Amap->N[2]);
 printf("#Bmap %d %d %d\n",Bmap->N[0],Bmap->N[1],Bmap->N[2]);
 for (a=0;a<2;++a)
  for (b=0;b<2;++b) Nab[a][b] = 0;

 if ((Amap->N[0]!=Bmap->N[0])||(Amap->N[1]!=Bmap->N[1])||(Amap->N[2]!=Bmap->N[2])){
  printf("#ERROR:Grid sizes for Xmap and Ymap are different.\n");
  exit(1);
 }
 
 for (x=0;x<Amap->N[0];++x){
  /* if ((x%10)==0) printf("#x %d/%d\n",x,Amap->N[0]); */
  for (y=0;y<Amap->N[1];++y){
   for (z=0;z<Amap->N[2];++z){
    vA = Amap->map[x][y][z];
    vB = Bmap->map[x][y][z];
    if (vA!=0) vA = 1;  
    if (vB!=0) vB = 1;  
         if ((vA==1)&&(vB==1)) Amap->map[x][y][z] = 3; 
    else if ((vA==0)&&(vB==1)) Amap->map[x][y][z] = 2; 
    else if ((vA==1)&&(vB==0)) Amap->map[x][y][z] = 1; 
    else Amap->map[x][y][z] = 0; 
    Nab[vA][vB] += 1.0;
   }
  }
 } 
 
 CC = Point_Correlation_Coefficient(Nab);
 printf("#Nab 00 %8.0lf 01 %8.0lf\n",Nab[0][0],Nab[0][1]);
 printf("#Nab 10 %8.0lf 11 %8.0lf\n",Nab[1][0],Nab[1][1]);
 printf("#CC  %lf\n",CC);
} /* end of Make_Difference_CHAR3DMAP() */





double Point_Correlation_Coefficient(Nab)
 double Nab[2][2];
{
 double bunshi,bunbo,CC;
 
 bunshi = (double)(Nab[0][0]*Nab[1][1]) - (double)(Nab[0][1]*Nab[1][0]);
 bunbo  = (double)(Nab[0][0]+Nab[0][1])*(double)(Nab[0][0]+Nab[1][0])
          *(double)(Nab[1][1]+Nab[0][1])*(double)(Nab[1][1]+Nab[1][0]);
 if (bunbo>0.0) CC = bunshi/sqrt(bunbo); else CC = 0.0;
 return(CC);

} /* end of Point_Correlation_Coefficient() */



void Make_Difference_CHAR3DMAP_MultiScale(Amap,Bmap,Nab)
 struct CHAR3DMAP *Amap; /* target 3D map */
 struct CHAR3DMAP *Bmap; /* reference 3D map */
 double Nab[2][2]; 
{
 int x,y,z,vA,vB,a,b;
 double CC;

 printf("#Make_Difference_CHAR3DMAP(Xmap,Pmap)\n");
 printf("#Amap %d %d %d\n",Amap->N[0],Amap->N[1],Amap->N[2]);
 printf("#Bmap %d %d %d\n",Bmap->N[0],Bmap->N[1],Bmap->N[2]);
 for (a=0;a<2;++a)
  for (b=0;b<2;++b) Nab[a][b] = 0;

 if ((Amap->N[0]!=Bmap->N[0])||(Amap->N[1]!=Bmap->N[1])||(Amap->N[2]!=Bmap->N[2])){
  printf("#ERROR:Grid sizes for Xmap and Ymap are different.\n");
  exit(1);
 }
 
 for (x=0;x<Amap->N[0];++x){
  /* if ((x%10)==0) printf("#x %d/%d\n",x,Amap->N[0]); */
  for (y=0;y<Amap->N[1];++y){
   for (z=0;z<Amap->N[2];++z){
    vA = Amap->map[x][y][z];
    vB = Bmap->map[x][y][z];
    /* if (vA!=0) vA = 1;  
       if (vB!=0) vB = 1;  */ 
         if ((vA==0)&&(vB==0)) Amap->map[x][y][z] = 0; 
    else if ((vA==1)&&(vB==1)) Amap->map[x][y][z] = 1; 
    else if ((vA==0)&&(vB==1)) Amap->map[x][y][z] = 2; 
    else if ((vA==1)&&(vB==0)) Amap->map[x][y][z] = 3; 
    else if (vA == vB)         Amap->map[x][y][z] = 4; 
    else if (vA >  vB)         Amap->map[x][y][z] = 5; 
    else if (vA <  vB)         Amap->map[x][y][z] = 6; 
    /* printf("vA %d vB %d Amap %d\n",vA,vB,Amap->map[x][y][z]);  */
   }
  }
 } 
 
 CC = Point_Correlation_Coefficient(Nab);
 printf("#Nab 00 %8.0lf 01 %8.0lf\n",Nab[0][0],Nab[0][1]);
 printf("#Nab 10 %8.0lf 11 %8.0lf\n",Nab[1][0],Nab[1][1]);
 printf("#CC  %lf\n",CC);

} /* end of Make_Difference_CHAR3DMAP_MultiScale() */


void Resize_CHAR3DMAP_by_MinMax(M,minL,maxL,newOrigPos,newNgrid,Nshift)
 struct CHAR3DMAP *M; 
 float  minL[3],maxL[3]; /* given */
 float  newOrigPos[3];  /* to be calculated */
 int    newNgrid[3];    /* to be calculated */
 int    Nshift[3];      /* to be calculated (Nshift[i] <=0 )*/
{
 int k;
 float minM[3],maxM[3]; /* for old grid */
 float minN[3],maxN[3]; /* for new grid */

 printf("#Resize_CHAR3DMAP_by_MinMax(minL %f %f %f maxL %f %f %f)\n",minL[0],minL[1],minL[2],maxL[0],maxL[1],maxL[2]);
 for (k=0;k<3;++k) {
  minM[k] = M->OrigPos[k];
  maxM[k] = M->OrigPos[k] + M->grid_width * M->N[k];
 }

 for (k=0;k<3;++k) {
   if (minM[k]<minL[k]) minN[k] = minM[k]; else minN[k] = minL[k];
   if (maxM[k]>maxL[k]) maxN[k] = maxM[k]; else maxN[k] = maxL[k];
  }

 printf("#new min %f %f %f max %f %f %f\n",minN[0],minN[1],minN[2],maxN[0],maxN[1],maxN[2]);

 for (k=0;k<3;++k){
   Nshift[k] =  (int)floor((minN[k] - M->OrigPos[k])/M->grid_width);
   newOrigPos[k] = M->OrigPos[k] + Nshift[k] * M->grid_width;
   newNgrid[k] = (int)ceil((maxN[k]-newOrigPos[k])/M->grid_width);
 }
 printf("#new newNgrid %d %d %d Nshift %d %d %d orig %f %f %f max %f %f %f\n",
 newNgrid[0], newNgrid[1], newNgrid[2], 
 Nshift[0], Nshift[1], Nshift[2], newOrigPos[0],newOrigPos[1],newOrigPos[2],
 newOrigPos[0] + M->grid_width * newNgrid[0],
 newOrigPos[1] + M->grid_width * newNgrid[1],
 newOrigPos[2] + M->grid_width * newNgrid[2]);

} /* end of Resize_CHAR3DMAP_by_MinMax() */



void Copy_CHAR3DMAP(A,B)
 struct CHAR3DMAP *A,*B;  /* A := B */
{
 int x,y,z;

 A->N[0] = B->N[0];
 A->N[1] = B->N[1];
 A->N[2] = B->N[2];
 A->OrigPos[0] = B->OrigPos[0];
 A->OrigPos[1] = B->OrigPos[1];
 A->OrigPos[2] = B->OrigPos[2];
 A->grid_width = B->grid_width;

 for (x=0;x<A->N[0];++x){
   for (y=0;y<A->N[1];++y){
     for (z=0;z<A->N[2];++z){
       A->map[x][y][z] = B->map[x][y][z];
     }
    }
  }

} /* end of Copy_CHAR3DMAP(A,B) */


void Copy_CHAR3DMAP_with_Nshift(A,B,Nshift)
 struct CHAR3DMAP *A,*B;  /* A := B */
 int Nshift[3];  /* Nshift = (origA - origB)/grid_width */
{
 int x,y,z,xx,yy,zz;

 printf("#Copy_CHAR3DMAP_with_Nshift(A,B,Nshift %d %d %d)\n",Nshift[0],Nshift[1],Nshift[2]);
 for (x=0;x<B->N[0];++x){
   for (y=0;y<B->N[1];++y){
     for (z=0;z<B->N[2];++z){
       xx = x - Nshift[0];
       yy = y - Nshift[1];
       zz = z - Nshift[2];
       /*
       xx = x + Nshift[0];
       yy = y + Nshift[1];
       zz = z + Nshift[2];
       */

       if ((xx>=0)&&(xx<A->N[0])&&(yy>=0)&&(yy<A->N[1])&&(zz>=0)&&(zz<A->N[2])){
         A->map[xx][yy][zz] = B->map[x][y][z];
         /* printf("%d %d %d = %d %d %d val %d\n",xx,yy,zz,x,y,z,A->map[xx][yy][zz]); */
       }
     }
    }
  }

} /* end of Copy_CHAR3DMAP_with_Nshift(A,B) */



void Compare_Two_CHAR3DMAP(P, minP, maxP, L, minL, maxL, NP,NL,NPL,Rec,Pre,Fme,Tani)
  struct CHAR3DMAP *P;         /* pocket (prediction) */
  unsigned char minP,maxP;
  struct CHAR3DMAP *L;         /* ligand (correct)   */
  unsigned char minL,maxL;
  int *NP,*NL,*NPL;
  float *Rec,*Pre,*Fme,*Tani;
{
 int x,y,z;
 unsigned char p,l;
 int Np,Nl,Npl;
 float rec,pre,fme,tani;

 printf("#Compare_Two_CHAR3DMAP(A,minP %d maxP %d B minL %d maxL %d)\n",minP,maxP,minL,maxL);
 Np = Nl = Npl = 0;

 for (x=0;x<P->N[0];++x){
   for (y=0;y<P->N[1];++y){
     for (z=0;z<P->N[2];++z){
       p = l = 0;
       if ((minP <= P->map[x][y][z]) && (P->map[x][y][z]<=maxP)){ p = 1;}
       if ((minL <= L->map[x][y][z]) && (L->map[x][y][z]<=maxL)){ l = 1;}
       if (p==1){ Np += 1;}
       if (l==1){ Nl += 1;}
       if ((p==1)&&(l==1)){ Npl += 1;}
     }
   }
 }

 rec = pre = fme = tani = 0.0;
 if (Nl>0){ rec = (float)Npl/(float)Nl;}
 if (Np>0){ pre = (float)Npl/(float)Np;}
 if ((Nl>0)&&(Np>0)){
   fme = 2.0/(1.0/rec + 1.0/pre);
   tani = (float)Npl/((float)Nl + (float)Np - (float)Npl);
 }

/*
 printf("#Na %d Nb %d Nab %d rec %f pre %f fme %f tani %f\n",Np,Nl,Npl,rec,pre,fme,tani);
*/

 *NP  = Np;  *NL  = Nl; *NPL = Npl; 
 *Rec = rec; *Pre = pre; *Fme = fme; *Tani = tani;

} /* end of Compare_Two_CHAR3DMAP() */


void Compare_Two_Bit_in_CHAR3DMAP(M, P_bitmask, L_bitmask, NP,NL,NPL,Rec,Pre,Fme,Tani)
  struct CHAR3DMAP *M;         /* Map */
  unsigned char P_bitmask;   /* bit for pocket prediction */
  unsigned char L_bitmask;   /* bit for correct ligand */
  int *NP,*NL,*NPL;
  float *Rec,*Pre,*Fme,*Tani;
{
 int x,y,z;
 unsigned char p,l;
 int Np,Nl,Npl;
 float rec,pre,fme,tani;

 printf("#Compare_Two_Bit_in_CHAR3DMAP(P_bitmask %d L_bitmask %d)\n",P_bitmask,L_bitmask);
 Np = Nl = Npl = 0;

 for (x=0;x<M->N[0];++x){
   for (y=0;y<M->N[1];++y){
     for (z=0;z<M->N[2];++z){
       p = l = 0; 
       if ((M->map[x][y][z]&P_bitmask) == P_bitmask){ p = 1;}
       if ((M->map[x][y][z]&L_bitmask) == L_bitmask){ l = 1;}
       if (p==1){ Np += 1;}
       if (l==1){ Nl += 1;}
       if ((p==1)&&(l==1)){ Npl += 1;}
     }
   }
 }

 /* printf("#Np %d Nl %d Npl %d\n",Np,Nl,Npl); */
 rec = pre = fme = tani = 0.0;
 if (Nl>0){ rec = (float)Npl/(float)Nl;}
 if (Np>0){ pre = (float)Npl/(float)Np;}
 if ((Nl>0)&&(Np>0)){
   fme = 2.0/(1.0/rec + 1.0/pre);
   tani = (float)Npl/((float)Nl + (float)Np - (float)Npl);
 }

 *NP  = Np;  *NL  = Nl; *NPL = Npl; 
 *Rec = rec; *Pre = pre; *Fme = fme; *Tani = tani;

/*
 printf("#Na %d Nb %d Nab %d rec %f pre %f fme %f tani %f\n",Np,Nl,Npl,rec,pre,fme,tani);
*/
} /* end of Compare_Two_Bit_in_CHAR3DMAP() */




void assign_CHAR3DMAP_from_VOXEL(cmap,vox,CutoffDensity,res_bitmask)
  struct CHAR3DMAP *cmap;
  struct VOXEL *vox;
  float  CutoffDensity;
  unsigned char res_bitmask; /* result bitmask (1 or 2 or 4 or 8 or ... or 128)*/
{
  int x,y,z,i;

  if ((cmap->N[0] != vox->N[0]) || (cmap->N[1] != vox->N[1]) || (cmap->N[2] != vox->N[2])){
    printf("#assign_CHAR3DMAP_from_VOXEL():size_of_cmap(%d %d %d) are not equal to size_of_vox(%d %d %d)\n",
      cmap->N[0], cmap->N[1], cmap->N[2], vox->N[0], vox->N[1], vox->N[2]);
    exit(1);
  }
 
  cmap->grid_width = vox->grid_width;
  for (i=0;i<3;++i){
    cmap->OrigPos[i] = vox->OrigPos[i];
  }
  
  for (x=0;x<vox->N[0];++x){
    for (y=0;y<vox->N[1];++y){
      for (z=0;z<vox->N[2];++z){
        if (vox->dat[x][y][z] >= CutoffDensity){
          cmap->map[x][y][z] = cmap->map[x][y][z]|res_bitmask;
        }
      }
    }
  }


} /* end of assign_CHAR3DMAP_from_VOXEL() */



void Foreground_MinMax_CHAR3DMAPL(cmap,bitmask,minN,maxN)
  struct CHAR3DMAP *cmap; /* (input) */
  unsigned char bitmask; /* result bitmask (1 or 2 or 4 or 8 or ... or 128) (input) */
  int minN[3],maxN[3];  /*  min / max of number of grids (to be calculated) */
{
  int x,y,z,i;
  unsigned char init; 

  printf("#Foreground_MinMax_CHAR3DMAPL(cmap,vox,bitmask %d)\n",bitmask);
  for (i=0;i<3;++i){
    minN[i] = maxN[i] = 0;
  }
 
  init = 1; 
  for (x=0;x<cmap->N[0];++x){
    for (y=0;y<cmap->N[1];++y){
      for (z=0;z<cmap->N[2];++z){
        if ((cmap->map[x][y][z]&bitmask) == bitmask){
          if (init==1){
            minN[0] = maxN[0] = x;
            minN[1] = maxN[1] = y;
            minN[2] = maxN[2] = z;
            init = 0;
          }
          if (x<minN[0]){minN[0] = x;}
          if (x>maxN[0]){maxN[0] = x;}
          if (y<minN[1]){minN[1] = y;}
          if (y>maxN[1]){maxN[1] = y;}
          if (z<minN[2]){minN[2] = z;}
          if (z>maxN[2]){maxN[2] = z;}
        }
      }
    }
  }

  for (i=0;i<3;++i){
    printf("[%d] min %d max %d minw %f maxw %f\n",i,minN[i],maxN[i],cmap->grid_width*minN[i], cmap->grid_width*maxN[i]);
  }

} /* end of Foreground_MinMax_CHAR3DMAPL() */


void resize_CHAR3DMAP_from_MinMax_and_MarginWidth(cmap,minN,maxN,margin_width,newN,newOrigPos)
  struct CHAR3DMAP *cmap; /* (input) */
  int    minN[3],maxN[3];   /* (input) */
  float  margin_width;    /* (input) */
  int    newN[3];         /* (to be calculated) */
  float  newOrigPos[3];   /* (to be calculated) */
{
  int i,marginN,shiftNmin[3],shiftNmax[3];
  
  marginN = (int)ceil(margin_width/cmap->grid_width);
  printf("#resize_CHAR3DMAP_from_MinMax_and_MarginWidth(minN %d %d %d maxN %d %d %d margin_width %f marginN %d)\n",
    minN[0], minN[1], minN[2],
    maxN[0], maxN[1], maxN[2], margin_width, marginN);
  for (i=0;i<3;++i){
    shiftNmin[i] = shiftNmax[i] = 0;
  } 

  for (i=0;i<3;++i){
    if (minN[i]<marginN){
       shiftNmin[i] = marginN - minN[i];
    }
    if ((cmap->N[i]-maxN[i])<marginN){
       shiftNmax[i] = marginN - (cmap->N[i]-maxN[i]);
    }
  } 

  printf("#shiftNmin %d %d %d shiftNmax %d %d %d\n",
  shiftNmin[0], shiftNmin[1], shiftNmin[2], shiftNmax[0], shiftNmax[1], shiftNmax[2]);

  for (i=0;i<3;++i){
    newN[i] = shiftNmin[i] + cmap->N[i] + shiftNmax[i];
    newOrigPos[i]  = cmap->OrigPos[i] - shiftNmin[i] * cmap->grid_width;
  }
  printf("#newN %d %d %d\n",newN[0],newN[1],newN[2]);

} /* end of resize_CHAR3DMAP_from_MinMax_and_MarginWidth() */


void copy_CHAR3DMAP_Same_GRID_and_Diff_OrigPos(New,Org)
  struct CHAR3DMAP *New;
  struct CHAR3DMAP *Org;
{
  int i,shiftN[3],x,y,z,xx,yy,zz;

  for (i=0;i<3;++i){
    shiftN[i] = (int)round((Org->OrigPos[i] - New->OrigPos[i])/Org->grid_width);
  } 

  printf("#copy_CHAR3DMAP_Same_GRID_and_Diff_OrigPos(shiftN %d %d %d)\n",shiftN[0],shiftN[1],shiftN[2]);
  printf("#OrigPos New %f %f %f Org %f %f %f gw %f\n",
    New->OrigPos[0], New->OrigPos[1], New->OrigPos[2], Org->OrigPos[0], Org->OrigPos[1], Org->OrigPos[2],Org->grid_width);
  for (x=0;x<Org->N[0];++x){
    xx = x + shiftN[0];
    for (y=0;y<Org->N[1];++y){
      yy = y + shiftN[1];
      for (z=0;z<Org->N[2];++z){
        zz = z + shiftN[2];
        if ((xx>=0) && (xx<New->N[0]) && (yy>=0) && (yy<New->N[1]) && (zz>=0) && (zz<New->N[2])){
          New->map[xx][yy][zz] = Org->map[x][y][z];
        }
      }
    }
 }

} /* end of copy_CHAR3DMAP_Same_GRID_and_Diff_OrigPos() */
