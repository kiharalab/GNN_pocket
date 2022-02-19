/*
 
 <ConSpheres.c>

=============================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.
This software is released under
the GNU Lesser General Public License (LGPL), see LICENSE.txt.
=============================================================

 functions for making contacting spheres

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include "globalvar.h"
#include "pdbstruct.h"
#include "PdbIO.h"
#include "Radius.h"
#include "Grid3D.h"
#include "VecFunc.h"
#include "Jacobi3.h"



/*** FUNCTIONS (GLOBAL) ***/
void Make_Sphere_Centers_from_PCaxis();
void Make_One_Sphere_Center_on_Gravity_Center();
void Setup_Radius_of_Spheres_Tangent_to_ATOMs();
void Remove_ATOM_less_than_Rmin_from_ATOMs();
void Setup_PCvar_and_PCaxis();
void Cal_Mean_and_CovM_from_ATOMs();
void write_ATOMs_as_Spheres_in_VRML();
void assign_radius_by_tFactor();
void assign_radius_by_Occup();
void split_string_to_float_array();
void Exclude_Spheres_Escape_to_Outside();


/*** FUNCTIONS (LOCAL) ***/
static double Dot_Prod3D();
static void Cross_Prod3D();
static float distance_float3D();
static void renumbering_ATOMs();
static void Multiply_Matrix3D_Vec3D_float();


void Make_Sphere_Centers_from_PCaxis(ConSpheHead,AtomHead,Nsphere,M,PCaxis)
  struct ATOM *ConSpheHead;  /* Contact Sheres ( to be generated) */
  struct ATOM *AtomHead;     /* Protein Atoms (input) */
  int    Nsphere;            /* (input) */
  double M[3];               /* (input) */
  double PCaxis[3];          /* (input) */
{
  int n,i,init;
  struct ATOM *an;
  double PfrmM[3],projPC,min_projPC, max_projPC;
  /** [1] decide min and max along PC axis */ 
  min_projPC = max_projPC = 0.0;
  an = AtomHead;
  init = 1;
  while (an->next != NULL){
    an = an->next;
    for (i=0;i<3;++i){
      PfrmM[i] = (double)an->Pos[i] - M[i];
    }
    projPC = Dot_Prod3D(PfrmM,PCaxis);
    if (init==1){
      min_projPC = projPC;
      max_projPC = projPC;
      init = 0;
    }
    else{
      if (projPC < min_projPC){min_projPC = projPC;}
      if (projPC > max_projPC){max_projPC = projPC;}
    }
  }

  printf("#min_projPC %lf max_projPC %lf\n",min_projPC, max_projPC);
  /** [2] Generate ConSpheHead  */ 
  an = ConSpheHead;
  for (n=0;n<Nsphere;++n){
    an->next = (struct ATOM*)malloc(sizeof(struct ATOM));   
    an->next->prev = an;
    an->next->next = NULL;
    an = an->next;

    for (i=0;i<3;++i){
      /*
      an->Pos[i] = M[i] + (-RPCaxis + (2.0*RPCaxis*n)/Nsphere) * PCaxis[i];
      */
      an->Pos[i] = M[i] + (min_projPC + n*(max_projPC - min_projPC)/Nsphere) * PCaxis[i];
      /* printf("#%d projPC %lf\n",n,min_projPC + n*(max_projPC - min_projPC)/Nsphere); */
    }
    an->AHtype    = 'H';
    an->Chain     = 'A';
    an->R         = 1.4;
    an->model_num = 1;
    sprintf(an->Anum,"%d",n+1);
    sprintf(an->Atom,"CA");
    sprintf(an->Resi,"UNK");
    sprintf(an->Rnum,"%3d",n+1);
  }  


} /* end of Make_Sphere_Centers_from_PCaxis() */




void Make_One_Sphere_Center_on_Gravity_Center(ConSpheHead,AtomHead,M)
  struct ATOM *ConSpheHead;  /* Contact Sheres ( to be generated) */
  struct ATOM *AtomHead;     /* Protein Atoms (input) */
  double M[3];               /* (input) */
{
  int i;
  struct ATOM *an;

  printf("#Make_One_Sphere_Center_on_Gravity_Center(ConSpheHead,AtomHead,M:%f %f %f)\n",M[0],M[1],M[2]);
  an = ConSpheHead;
  an->next = (struct ATOM*)malloc(sizeof(struct ATOM));   
  an->next->prev = an;
  an->next->next = NULL;
  an = an->next;

  for (i=0;i<3;++i){
      an->Pos[i] = M[i];
  }
  an->AHtype    = 'H';
  an->Chain     = 'A';
  an->R         = 1.4;
  an->model_num = 1;
  sprintf(an->Anum,"%d",1);
  sprintf(an->Atom,"CA");
  sprintf(an->Resi,"UNK");
  sprintf(an->Rnum,"%3d",1);

} /* end of Make_One_Sphere_Center_on_Gravity_Center() */








void Setup_Radius_of_Spheres_Tangent_to_ATOMs(ConSpheHead,AtomHead)
  struct ATOM *ConSpheHead;  /* Contact Sheres ( *->R will be assigned) */
  struct ATOM *AtomHead;     /* Protein Atoms (input) */
{
  struct ATOM *sph,*atm;
  float D, Dvdw, Dmin;
  int init;
 
  sph = ConSpheHead;
  while (sph->next != NULL){
    sph = sph->next;
    Dmin = 0.0;
    atm = AtomHead;
    init = 1;
    while (atm->next != NULL){
      atm = atm->next;
      D =  distance_float3D(atm->Pos,sph->Pos);
      Dvdw = D - atm->R;
      if ((init==1)||(Dvdw<Dmin)){
         Dmin = Dvdw;
      }
      init = 0;
    }
    sph->R = Dmin;
  }  

} /* end of Setup_Radius_of_Spheres_Tangent_to_ATOMs() */


void Remove_ATOM_less_than_Rmin_from_ATOMs(ConSpheHead,Rmin)
  struct ATOM *ConSpheHead;
  float  Rmin; 
{
  struct ATOM *sph;

  sph = ConSpheHead;
  while (sph->next != NULL){
    sph = sph->next;
    if (sph->R < Rmin){
      sph->prev->next = sph->next;
    }
  }

} /* end of Remove_ATOM_less_than_Rmin_from_ATOMs() */


void Setup_PCvar_and_PCaxis(CovM,PCvar,PCaxis)
  double CovM[3][3];   /* covariance matrix (input) */
  double PCvar[3];     /* variance for each PC axis (to be calculated) */
  double PCaxis[3][3]; /* PC axis [axis_num][x,y,z] (to be calculated) */
{
 int i,j,index[3];
 double A[3][3],U[3][3],PC0xPC1[3];

 /* printf("#Setup_PCvar_and_PCaxis()\n"); */
 for (i=0;i<3;++i){
   for (j=0;j<3;++j){ A[i][j] = CovM[i][j];}
 }
 Jacobi_Wilkinson3(A,U);

 Sort_Eigen_Value3(A,index,'D'); /* sort by decreasing order of eigen values */
 for (i=0;i<3;++i){
   PCvar[i]  = A[index[i]][index[i]];
   /* printf("#PCvar[%d] %lf\n",i,PCvar[i]); */
   for (j=0;j<3;++j){
     PCaxis[i][j] = U[j][index[i]];
   }
 }

  /* Check handedness, force the right-handed system */
  Cross_Prod3D(PC0xPC1,PCaxis[0],PCaxis[1]);
  if (Dot_Prod3D(PC0xPC1,PCaxis[2])<0.0){
    for (i=0;i<3;++i){ PCaxis[2][i] = -PCaxis[2][i];}
  }

} /* end of Setup_PCvar_and_PCaxis() */



void Cal_Mean_and_CovM_from_ATOMs(Ahead,M,CovM)
 struct ATOM *Ahead;
 double  M[3];
 double CovM[3][3];
{
 struct ATOM *an;
 int Natom;
 int i,j;
 double Dvec[3];
 double sumDensity;

 /* Caliculate Mean,Min,Max */
 sumDensity = 0.0;
 for (i=0;i<3;++i){ M[i] = 0.0;}
 Natom = 0;

 an = Ahead;
 while (an->next != NULL){
  an = an->next;
  sumDensity += an->Weight;
  for (i=0;i<3;++i){ M[i] += an->Weight * an->Pos[i];}
   ++Natom;
 }

 for (i=0;i<3;++i){
   M[i] /= sumDensity;
 }

 /* Caliculate Covariance Matrix */

 an = Ahead;
 for (i=0;i<3;++i){
   for (j=0;j<3;++j){CovM[i][j] = 0.0;}
 }

 while (an->next != NULL){
   an = an->next;
   for (i=0;i<3;++i){Dvec[i] = an->Pos[i] - M[i];}

   for (i=0;i<3;++i){
     for (j=i;j<3;++j){
       CovM[i][j] += an->Weight * Dvec[i]*Dvec[j];
     }
   }
 }

 for (i=0;i<3;++i){
   for (j=i;j<3;++j){
     CovM[i][j] /= sumDensity;
     CovM[j][i] = CovM[i][j];
   }
 }
} /* end of Cal_Mean_and_CovM_from_ATOMs() */
                                                    


void write_ATOMs_as_Spheres_in_VRML(ofname, AtomHead,RGBTstr)
 char *ofname;
 struct ATOM *AtomHead;
 char   *RGBTstr;  /* "[R]:[G]:[B]:[T]" */
{
  FILE *fpo;
  struct ATOM *atm;
  float RGBT[4];

  split_string_to_float_array(RGBTstr,':',RGBT,4);
  fpo = fopen(ofname,"w");
  printf("#write_ATOMs_as_Spheres_in_VRML()-->'%s'\n",ofname);
  fprintf(fpo,"#VRML V2.0 utf8\n");
  fprintf(fpo,"\n");
  fprintf(fpo,"#COMMAND %s\n",PAR.COMMAND);
  fprintf(fpo,"\n");

  atm = AtomHead;
  while (atm->next != NULL){
    atm = atm->next;
    fprintf(fpo,"#%s %s %s %s %c\n",atm->Anum,atm->Atom,atm->Resi,atm->Rnum,atm->Chain);
    fprintf(fpo,"Transform{\n");
    fprintf(fpo," translation %f %f %f\n",atm->Pos[0], atm->Pos[1], atm->Pos[2]);
    fprintf(fpo," scale       1.0 1.0 1.0\n");
    fprintf(fpo," children[\n");
    fprintf(fpo," Shape{\n");
    fprintf(fpo,"  appearance Appearance{\n");
    fprintf(fpo,"      material Material{\n");
    fprintf(fpo,"        diffuseColor %f %f %f\n",RGBT[0],RGBT[1],RGBT[2]);
    fprintf(fpo,"        transparency %f\n",RGBT[3]);
    fprintf(fpo,"      }\n");
    fprintf(fpo,"  }\n");
    fprintf(fpo,"  geometry Sphere{ radius %f} }\n",atm->R);
    fprintf(fpo," ]\n");
    fprintf(fpo,"}\n");
  }

} /* end of write_one_sphere_in_VRML() */



void assign_radius_by_tFactor(AtomHead)
  struct ATOM *AtomHead;
{
  struct ATOM *an;

  printf("#assign_radius_by_tFactor(AtomHead)\n");
  an = AtomHead;
  while (an->next != NULL){
    an = an->next;
    /* printf("#%s %s %s tFactor %f\n",an->Atom,an->Resi,an->Rnum,an->tFactor); */
    an->R  = an->tFactor;
    an->RR = an->R * an->R;
  }
} /* end of assign_radius_by_tFactor() */



void assign_radius_by_Occup(AtomHead)
  struct ATOM *AtomHead;
{
  struct ATOM *an;

  printf("#assign_radius_by_Occup(AtomHead)\n");
  an = AtomHead;
  while (an->next != NULL){
    an = an->next;
    /* printf("#%s %s %s Occup %f\n",an->Atom,an->Resi,an->Rnum,an->Occup); */
    an->R  = an->Occup;
    an->RR = an->R * an->R;
  }
} /* end of assign_radius_by_Occup() */




double Dot_Prod3D(a,b)
 double a[3],b[3];
{
   return(a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
} /* end of Dot_Prod3D() */



void Cross_Prod3D(c,a,b) /* c = a x b */
 double c[3],a[3],b[3];
{
 c[0] = a[1]*b[2] - b[1]*a[2];
 c[1] = a[2]*b[0] - b[2]*a[0];
 c[2] = a[0]*b[1] - b[0]*a[1];
} /* end of Cross_Prod3D() */
 

float distance_float3D(a,b)
  float a[3],b[3];
{
   float D;
   int i;
   D = 0.0;
   for (i=0;i<3;++i){
     D += (a[i]-b[i])*(a[i]-b[i]);
   } 
   return(sqrt(D));
}


void split_string_to_float_array(str,splsym,X,lengthX)
  char *str;   /* (input) */
  char splsym; /* (input) */
  float *X;     /* (ouput) array of float. size of array should be lengthX */
  int   lengthX;  /* (input) */
{
  int Nword,*Wsta,*Wend,i;
  char *buff;

  Wsta = (int *)malloc(sizeof(int)*lengthX);
  Wend = (int *)malloc(sizeof(int)*lengthX);
  buff = (char *)malloc(sizeof(char)*strlen(str));
  Split_to_Words(str,splsym,&Nword,Wsta,Wend,lengthX);

  if (Nword==lengthX){
    for (i=0;i<Nword;++i){
      Get_Part_Of_Line(buff,str,Wsta[i],Wend[i]);
      X[i] = atof(buff);
    }
  }
  else{
    printf("#ERROR(split_string_to_float_array): Nword(%d) for ('%s') is not %d.\n",Nword,str,lengthX); exit(1);
  }
  free(Wend);
  free(Wsta);
} /* end of split_string_to_float_array() */



 
void Exclude_Spheres_Escape_to_Outside(ConSpheHead,AtomHead,PCaxis,PCaxis_num)
  struct ATOM *ConSpheHead;  /* Contact Sheres ( to be generated) */
  struct ATOM *AtomHead;     /* Protein Atoms (input) */
  double PCaxis[3][3];        /* (input) */
  int    PCaxis_num;          /* PCaxis_num for contact sphere axis */ 
{
  int i;
  struct ATOM *sph,*an;
  float XX[3],YY[3],ZZ[3],Rmat[3][3];
  float x[3],xprj[3],dprj;
  int   Nconatm_plusX, Nconatm_minusX;
  float tolerance;

  printf("#Check_Sphere_can_Escape_or_Not(ConSpheHead,AtomHead,PCaxis,PCaxis_num %d)\n",PCaxis_num);
  tolerance  = 0.00001;

  for (i=0;i<3;++i){
    XX[i] = (float)PCaxis[PCaxis_num][i];
    YY[i] = (float)PCaxis[(PCaxis_num+1)%3][i];
    ZZ[i] = (float)PCaxis[(PCaxis_num+2)%3][i];
  }  

  Rmat[0][0] = XX[0]; Rmat[0][1] = YY[0]; Rmat[0][2] = ZZ[0]; 
  Rmat[1][0] = XX[1]; Rmat[1][1] = YY[1]; Rmat[1][2] = ZZ[1]; 
  Rmat[2][0] = XX[2]; Rmat[2][1] = YY[2]; Rmat[2][2] = ZZ[2]; 

  sph = ConSpheHead;

  while (sph->next != NULL){
    sph = sph->next;
    Nconatm_plusX  = 0; 
    Nconatm_minusX = 0;
 
    an = AtomHead;
    while (an->next != NULL){
      an = an->next;
      for (i=0;i<3;++i){
        x[i] = an->Pos[i] - sph->Pos[i];
      }

      Multiply_Matrix3D_Vec3D_float(xprj,Rmat,x);
      /* printf("x %f %f %f xprj %f %f %f\n",x[0],x[1],x[2],xprj[0],xprj[1],xprj[2]); */
      dprj = sqrt(xprj[1]*xprj[1] + xprj[2]*xprj[2]) + an->R;
      if (dprj <= sph->R + tolerance){
        if (xprj[0]>=0.0){
          Nconatm_plusX += 1;
        }
        else{
          Nconatm_minusX += 1;
        }
      }
    }

    printf("#CONATM %s %s %s R %f Ncon_plus %d Ncon_minus %d",sph->Anum,sph->Rnum,sph->Resi,sph->R,Nconatm_plusX, Nconatm_minusX);
 
    if ( ((Nconatm_minusX==0)&&(Nconatm_plusX>0)) || ((Nconatm_minusX>0)&&(Nconatm_plusX==0)) ){
      printf(" REMOVE!");
      sph->R = 0.0;
      if (sph->prev != NULL){
        sph->prev->next = sph->next;
      }
      if (sph->next != NULL){
        sph->next->prev = sph->prev;
      }
    }
    printf("\n");
 }

  renumbering_ATOMs(ConSpheHead);

} /* end of Exclude_Spheres_Escape_to_Outside() */


void renumbering_ATOMs(AtomHead)
  struct ATOM *AtomHead;
{
  struct ATOM *atm;
  int n;

  atm = AtomHead;
  n = 0;
  while (atm->next != NULL){
    atm = atm->next;
    sprintf(atm->Anum,"%d",n+1);
    sprintf(atm->Rnum,"%3d",n+1);
    n += 1;
  }
} /* end of renumbering_ATOMs() */

void Multiply_Matrix3D_Vec3D_float(y,A,x)    /* y = A*x */
 float y[3],A[3][3], x[3];
{
 int i,k;

 for (i=0;i<3;++i){
   y[i] = 0.0;
   for (k=0;k<3;++k){y[i] += A[i][k]*x[k];}
 }

} /* end of Multiply_Matrix3D_Vec3D_float() */


