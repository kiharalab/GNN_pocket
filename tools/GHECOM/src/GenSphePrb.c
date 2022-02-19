/*

 <GenSphePrb.c>

 
  for generating spherical probes by 3-contacting geometry 

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
#include "GenSphePrb.h" 
#include "Grid3D.h" 


/*** FUNCTION (GLOBAL) ***/
int  Cal_Sphere_On_Three_Atoms();
void Make_Probe_Spheres_Three_Contact();
void Make_Probe_Spheres_Three_Contact_NEIGHBOR_list();
int  Check_Crash_Pos_and_Atoms();
int  Check_Crash_Pos_and_NEIGHBOR_Atoms();
void Malloc_FLOAT2DMAP();
void Free_FLOAT2DMAP();
void Cal_Distance_FLOAT2DMAP();
void Remove_Probe_Atoms_Larger_Specified_Rinaccess();
void Remove_Probe_Atoms_Smaller_Specified_Pocketness();
void Remove_Probe_Atoms_Smaller_Specified_tFactor();
void Sort_Probe_By_Pocket_Cluster_Number();
void Write_Spherical_Probes_in_UCSF_DOCK_format();
void Label_Probe_Atoms_By_CHAR3DMAP();


/*** FUNCTION (LOCAL) ***/
static void Add_Atom();
static void Cal_Tij_on_Base_Plane();
static int Normalize();
static int Normalize_Double();
static float Dot_Prod();  
static double Dot_Prod_Double();  
static void Cross_Prod();
static void Cross_Prod_Double();
static void Make_NEIGHBOR_list();
static void Add_NEIGHBOR_to_Atom();







void Make_Probe_Spheres_Three_Contact(Phead,Ahead,Rprobe,Dmat)
 struct ATOM *Phead;  /* Probe   atoms to be generated */
 struct ATOM *Ahead;  /* Protein atoms */
 float  Rprobe;
 struct FLOAT2DMAP *Dmat;
{
 /*
  This functions is the basic one.
  The answer is correct, but it needs much computational time 
  for many atoms and for a large probe sphere. 
 */
 
 struct ATOM *an,*bn,*cn;
 float Rprobe2;
 float Ppos[3],Qpos[3];
 int Nprobe,Ncheck;
 char spok;

 printf("#Make_Probe_Spheres_Three_Contact(Phead,Ahead,Rprobe,Dmat)\n");
 Rprobe2 = Rprobe *2;
 Nprobe = Ncheck = 0;
 Phead->next = NULL;
 an = Ahead;

 while (an->next != NULL){
  an = an->next;
  bn = an;
  while (bn->next != NULL){
   bn = bn->next;
   /* printf("%d %d\n",an->num,bn->num); fflush(stdout); */   
   if (Dmat->m[an->num][bn->num] <= an->R+bn->R+Rprobe2){
    cn = bn;
    while (cn->next != NULL){
     cn = cn->next;
     if( (Dmat->m[an->num][cn->num] <= an->R+cn->R+Rprobe2) &&
         (Dmat->m[bn->num][cn->num] <= bn->R+cn->R+Rprobe2) ){
       spok = Cal_Sphere_On_Three_Atoms(Ppos,Qpos,Rprobe,an,bn,cn);
       if (spok==1){
        ++Ncheck;
         if (Check_Crash_Pos_and_Atoms(Ppos,Rprobe,Ahead,an,bn,cn)==0) 
        { Add_Atom(Phead,Ppos,Rprobe,' ',Nprobe,an,bn,cn); ++Nprobe;}
         if (Check_Crash_Pos_and_Atoms(Qpos,Rprobe,Ahead,an,bn,cn)==0) 
         {Add_Atom(Phead,Qpos,Rprobe,' ',Nprobe,an,bn,cn); ++Nprobe;}
       }
      } 

     } /* cn */ 

    }

   } /* bn */ 

  } /* an */

  printf("#Ncheck %d Nprobe %d\n",Ncheck,Nprobe);

} /* end of Make_Probe_Spheres_Three_Contact() */



void  Make_Probe_Spheres_Three_Contact_NEIGHBOR_list(Phead,Ahead,Rprobe,Dmat)
 struct ATOM *Phead;  /* Probe   atoms to be generated */
 struct ATOM *Ahead;  /* Protein atoms */
 float  Rprobe;
 struct FLOAT2DMAP *Dmat;
{
 /*
  This functions is the basic one.
  The answer is correct, but it needs much computational time 
  for many atoms and for a large probe sphere. 
 */
 
 struct ATOM *an,*bn,*cn;
 struct NEIGHBOR *nei_b, *nei_c;
 float Rprobe2;
 float Ppos[3],Qpos[3];
 int Nprobe,Ncheck;
 char spok;

 printf("#Make_Probe_Spheres_Three_Contact_NEIGHBOR_list(Phead,Ahead,Rprobe,Dmat)\n");
 Make_NEIGHBOR_list(Ahead,Dmat,Rprobe);

 Rprobe2 = Rprobe *2;
 Nprobe = 0;
 Phead->next = NULL;
 an = Ahead;
 Ncheck = 0;

 while (an->next != NULL){
  an = an->next;
  nei_b = &(an->NeiHead); 
  while (nei_b->next != NULL){
   nei_b = nei_b->next;
   bn  = nei_b->atom;
    if (an->num < bn->num){
    nei_c = nei_b;
    while (nei_c->next != NULL){
     nei_c = nei_c->next;
     cn = nei_c->atom;
     if ((bn->num < cn->num) && (Dmat->m[bn->num][cn->num] <= bn->R+cn->R+Rprobe2) ){
       spok = Cal_Sphere_On_Three_Atoms(Ppos,Qpos,Rprobe,an,bn,cn);
       if (spok==1){
         ++Ncheck;
         if (Check_Crash_Pos_and_NEIGHBOR_Atoms(Ppos,Rprobe,&(an->NeiHead),an,bn,cn)==0)   
        { Add_Atom(Phead,Ppos,Rprobe,' ',Nprobe,an,bn,cn); ++Nprobe;}
         if (Check_Crash_Pos_and_NEIGHBOR_Atoms(Qpos,Rprobe,&(an->NeiHead),an,bn,cn)==0)   
         {Add_Atom(Phead,Qpos,Rprobe,' ',Nprobe,an,bn,cn); ++Nprobe;}
       }
      } 

     } /* nei_c */ 

    }

   } /* nei_b */ 

  } /* an */

  printf("Ncheck %d Nprobe %d\n",Ncheck,Nprobe);

} /* end of Make_Probe_Spheres_Three_Contact_NEIGHBOR_list() */














int Cal_Sphere_On_Three_Atoms(Ppos,Qpos,Rprobe,I,J,K)
 float  Ppos[3],Qpos[3]; /* Probe Atom */
 float  Rprobe;
 struct ATOM *I,*J,*K;   /* Three atoms tangent to the probe atom */
{
 double Tij[3],Tjk[3],Tik[3];
 double X[3],Y[3],Z[3];  /** Local Coordinates **/
 double IK[3],dotIK,B[3],v1[3],v2[3];
 int m;
 double u,height;
 char ok;

 /************
  <Float or Double ?>
 
  For precise calculation all the real numbers are calculated using "double"
  in this function, although all atom Position and distance maps are 
  stored in "float". 

  <Conditions for having 3-contact sphere with R for three atoms I, J, K>

 (1) Simple distance condition

    |I - J| <= RI + RJ + 2*Rprobe   
    |J - K| <= RJ + RK + 2*Rprobe   
    |K - I| <= RK + RI + 2*Rprobe   

  This function "Cal_Sphere_On_Three_Atoms()" does not check these conditions.
  If you want efficient calculation, check these conditions before apply the function.

 (2) Even if the conditions (1) are satisfied, another condition is required. 

   We have to define some terminology.  
  
   base_plane : triangle of the three atom center I, J, K.
   base_point : point on the base_plane. The line between the base point 
                and center of probe, is perpendicular to the base_plane.  

    |Ri + Rp|  >= |base_point - I|
    |Rj + Rp|  >= |base_point - J|
    |Rk + Rp|  >= |base_point - K|

  In the program, these conditions are directly related to 
  calculation of 'height'. 

 ************/

 /*
 printf(">> Cal_Sphere_On_Three_Atoms(Ppos,Qpos,Rprobe,I,J,K) <<\n");
 */ 
 /*** Set X,Y,Z ***/
 for (m=0;m<3;++m)  X[m] = J->Pos[m] - I->Pos[m]; 
 ok = Normalize_Double(X); 
  if (ok==0) {printf("#WARNING NORM JI IS ZERO\n"); 
  /*
  printf("I %s %s %s %s %f %f %f\n",I->Anum,I->Atom,I->Resi,I->Rnum,I->Pos[0],I->Pos[1],I->Pos[2]);
  printf("J %s %s %s %s %f %f %f\n",J->Anum,J->Atom,J->Resi,J->Rnum,J->Pos[0],J->Pos[1],J->Pos[2]);
  printf("K %s %s %s %s %f %f %f\n",K->Anum,K->Atom,K->Resi,K->Rnum,K->Pos[0],K->Pos[1],K->Pos[2]);
 */
  return(0);}
 
 for (m=0;m<3;++m)  IK[m] = K->Pos[m] - I->Pos[m]; 
 dotIK = Dot_Prod_Double(X,IK); 
 for (m=0;m<3;++m) Y[m] = IK[m] - dotIK * X[m]; 
 ok = Normalize_Double(Y);  
 if (ok==0) {printf("#WARNING NORM KI IS ZERO\n"); 
 /*
 printf("I %s %s %s %s %f %f %f\n",I->Anum,I->Atom,I->Resi,I->Rnum,I->Pos[0],I->Pos[1],I->Pos[2]);
 printf("J %s %s %s %s %f %f %f\n",J->Anum,J->Atom,J->Resi,J->Rnum,J->Pos[0],J->Pos[1],J->Pos[2]);
 printf("K %s %s %s %s %f %f %f\n",K->Anum,K->Atom,K->Resi,K->Rnum,K->Pos[0],K->Pos[1],K->Pos[2]);
 */
 return(0);}
 Cross_Prod_Double(Z,X,Y);

 /*
 printf("I %s %s %s %s %f %f %f\n",I->Anum,I->Atom,I->Resi,I->Rnum,I->Pos[0],I->Pos[1],I->Pos[2]);
 printf("J %s %s %s %s %f %f %f\n",J->Anum,J->Atom,J->Resi,J->Rnum,J->Pos[0],J->Pos[1],J->Pos[2]);
 printf("K %s %s %s %s %f %f %f\n",K->Anum,K->Atom,K->Resi,K->Rnum,K->Pos[0],K->Pos[1],K->Pos[2]);
 printf("X %f %f %f\n",X[0],X[1],X[2]);
 printf("Y %f %f %f\n",Y[0],Y[1],Y[2]);
 printf("Z %f %f %f\n",Z[0],Z[1],Z[2]);
 */
 
 /*** Calculate Tij,Tjk,Tik ***/
 Cal_Tij_on_Base_Plane(Tij,I,J,Rprobe);
 Cal_Tij_on_Base_Plane(Tjk,J,K,Rprobe);
 Cal_Tij_on_Base_Plane(Tik,I,K,Rprobe);

 /* 
 printf("Tij %f %f %f\n",Tij[0],Tij[1],Tij[2]);
 printf("Tjk %f %f %f\n",Tjk[0],Tjk[1],Tjk[2]);
 printf("Tik %f %f %f\n",Tik[0],Tik[1],Tik[2]);
 */
 
 /*** Calculate point b (center of base-plane IJK) **/
 for (m=0;m<3;++m) 
  { v1[m] = Tik[m] - I->Pos[m]; 
    v2[m] = Tik[m] - Tij[m]; }

 u = Dot_Prod_Double(v1,v2)/Dot_Prod_Double(v1,Y);
 for (m=0;m<3;++m)  B[m] = Tij[m] + u * Y[m];

 /* printf("B %f %f %f\n",B[0],B[1],B[2]); */


 /*** Cal Probe center 'P' **/
 for (m=0;m<3;++m)  v1[m] = B[m] - I->Pos[m];
 height = (I->R + Rprobe)*(I->R + Rprobe) - (v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);

 /* Case for improper obtuse(donkaku) triangle */ 
 if (height <=0.0)  
 { 
   /*
   printf("#WARNING:height is zero. Improper 'obtuse' triangle.\n"); 
   */
   return(0); }

 height = sqrt(height);

 /* printf("height %f\n",height); */
 for (m=0;m<3;++m)  
  { Ppos[m] = B[m] + height * Z[m];
    Qpos[m] = B[m] - height * Z[m]; }

 return(1);

} /* end of Cal_Sphere_On_Three_Atoms() */







int Cal_Sphere_On_Three_Atoms_Float(Ppos,Qpos,Rprobe,I,J,K)
 float  Ppos[3],Qpos[3]; /* Probe Atom */
 float  Rprobe;
 struct ATOM *I,*J,*K;   /* Three atoms tangent to the probe atom */
{
 float Tij[3],Tjk[3],Tik[3];
 float X[3],Y[3],Z[3];  /** Local Coordinates **/
 float IK[3],dotIK,B[3],v1[3],v2[3];
 int m;
 float u,height;
 char ok;

 /************

  <Conditions for having 3-contact sphere with R for three atoms I, J, K>

 (1) Simple distance condition

    |I - J| <= RI + RJ + 2*Rprobe   
    |J - K| <= RJ + RK + 2*Rprobe   
    |K - I| <= RK + RI + 2*Rprobe   

  This function "Cal_Sphere_On_Three_Atoms()" does not check these conditions.
  If you want efficient calculation, check these conditions before apply the function.

 (2) Even if the conditions (1) are satisfied, another condition is required. 

   We have to define some terminology.  
  
   base_plane : triangle of the three atom center I, J, K.
   base_point : point on the base_plane. The line between the base point 
                and center of probe, is perpendicular to the base_plane.  

    |Ri + Rp|  >= |base_point - I|
    |Rj + Rp|  >= |base_point - J|
    |Rk + Rp|  >= |base_point - K|

  In the program, these conditions are directly related to 
  calculation of 'height'. 

 ************/

 /*
 printf(">> Cal_Sphere_On_Three_Atoms(Ppos,Qpos,Rprobe,I,J,K) <<\n");
 */ 
 /*** Set X,Y,Z ***/
 for (m=0;m<3;++m)  X[m] = J->Pos[m] - I->Pos[m]; 
 ok = Normalize(X); 
 if (ok==0) {printf("#WARNING NORM JI IS ZERO\n"); 
  printf("#I %f %f %f\n",I->Pos[0],I->Pos[1],I->Pos[2]);
  printf("#J %f %f %f\n",I->Pos[0],I->Pos[1],I->Pos[2]);
  printf("#K %f %f %f\n",I->Pos[0],I->Pos[1],I->Pos[2]);
  return(0);}
 
 for (m=0;m<3;++m)  IK[m] = K->Pos[m] - I->Pos[m]; 
 dotIK = Dot_Prod(X,IK); 
 for (m=0;m<3;++m) Y[m] = IK[m] - dotIK * X[m]; 
 ok = Normalize(Y);  if (ok==0) {printf("#WARNING NORM KI IS ZERO\n"); return(0);}
 Cross_Prod(Z,X,Y);

 /*
 printf("I %s %s %s %s %f %f %f\n",I->Anum,I->Atom,I->Resi,I->Rnum,I->Pos[0],I->Pos[1],I->Pos[2]);
 printf("J %s %s %s %s %f %f %f\n",J->Anum,J->Atom,J->Resi,J->Rnum,J->Pos[0],J->Pos[1],J->Pos[2]);
 printf("K %s %s %s %s %f %f %f\n",K->Anum,K->Atom,K->Resi,K->Rnum,K->Pos[0],K->Pos[1],K->Pos[2]);
 printf("X %f %f %f\n",X[0],X[1],X[2]);
 printf("Y %f %f %f\n",Y[0],Y[1],Y[2]);
 printf("Z %f %f %f\n",Z[0],Z[1],Z[2]);
 */
 
 /*** Calculate Tij,Tjk,Tik ***/
 Cal_Tij_on_Base_Plane(Tij,I,J,Rprobe);
 Cal_Tij_on_Base_Plane(Tjk,J,K,Rprobe);
 Cal_Tij_on_Base_Plane(Tik,I,K,Rprobe);

 /* 
 printf("Tij %f %f %f\n",Tij[0],Tij[1],Tij[2]);
 printf("Tjk %f %f %f\n",Tjk[0],Tjk[1],Tjk[2]);
 printf("Tik %f %f %f\n",Tik[0],Tik[1],Tik[2]);
 */
 
 /*** Calculate point b (center of base-plane IJK) **/
 for (m=0;m<3;++m) 
  { v1[m] = Tik[m] - I->Pos[m]; 
    v2[m] = Tik[m] - Tij[m]; }

 u = Dot_Prod(v1,v2)/Dot_Prod(v1,Y);
 for (m=0;m<3;++m)  B[m] = Tij[m] + u * Y[m];

 /* printf("B %f %f %f\n",B[0],B[1],B[2]); */


 /*** Cal Probe center 'P' **/
 for (m=0;m<3;++m)  v1[m] = B[m] - I->Pos[m];
 height = (I->R + Rprobe)*(I->R + Rprobe) - (v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);

 /* Case for improper obtuse(donkaku) triangle */ 
 if (height <=0.0)  
 { 
   /*
   printf("#WARNING:height is zero. Improper 'obtuse' triangle.\n"); 
   */
   return(0); }

 height = sqrt(height);

 /* printf("height %f\n",height); */
 for (m=0;m<3;++m)  
  { Ppos[m] = B[m] + height * Z[m];
    Qpos[m] = B[m] - height * Z[m]; }

 return(1);

} /* end of Cal_Sphere_On_Three_Atoms_Float() */



void Cal_Tij_on_Base_Plane(Tij,I,J,Rp)
 double Tij[3];
 struct ATOM *I,*J;
 float Rp;
{
 double RRij,D[3];
 double t,A,B;
 int m;

 for (m=0;m<3;++m) D[m] = I->Pos[m] - J->Pos[m];
 RRij = D[0]*D[0] + D[1]*D[1] + D[2]*D[2];
 if (RRij > 0.0)
 {
   A = I->R + Rp;
   B = J->R + Rp;
   t = 0.5*(A*A - B*B)/RRij;
   for (m=0;m<3;++m) 
    Tij[m] = 0.5*(I->Pos[m] + J->Pos[m])+ t * (J->Pos[m] - I->Pos[m]);
 
  }
 else for (m=0;m<3;++m) Tij[m] = I->Pos[m];

} /* end of Cal_Tij_on_Base_Plane() */





void Cal_Tij_on_Base_Plane_Float(Tij,I,J,Rp)
 float Tij[3];
 struct ATOM *I,*J;
 float Rp;
{
 float RRij,D[3];
 double t,A,B;
 int m;

 for (m=0;m<3;++m) D[m] = I->Pos[m] - J->Pos[m];
 RRij = D[0]*D[0] + D[1]*D[1] + D[2]*D[2];
 if (RRij > 0.0)
 {
   A = I->R + Rp;
   B = J->R + Rp;
   t = 0.5*(A*A - B*B)/RRij;
   for (m=0;m<3;++m) 
    Tij[m] = 0.5*(I->Pos[m] + J->Pos[m])+ t * (J->Pos[m] - I->Pos[m]);
 
  }
 else for (m=0;m<3;++m) Tij[m] = I->Pos[m];

} /* end of Cal_Tij_on_Base_Plane_Float() */





int Normalize(X)
 float X[3];
{
 float D;
 D = X[0]*X[0] + X[1]*X[1] + X[2]*X[2];
 if (D>0.0) 
 { D = sqrt(D);
   X[0] /= D; X[1] /= D; X[2] /= D; 
   return(1); }
 else return(0);

} /* end of Normalize() */


int Normalize_Double(X)
 double X[3];
{
 double D;
 D = X[0]*X[0] + X[1]*X[1] + X[2]*X[2];
 if (D>0.0) 
 { D = sqrt(D);
   X[0] /= D; X[1] /= D; X[2] /= D; 
   return(1); }
 else return(0);

} /* end of Normalize_Double() */





float Dot_Prod(A,B)  
 float A[3],B[3];
{ return(A[0]*B[0] + A[1]*B[1] + A[2]*B[2]); }

double Dot_Prod_Double(A,B)  
 double A[3],B[3];
{ return(A[0]*B[0] + A[1]*B[1] + A[2]*B[2]); }

void Cross_Prod(C,A,B)  
 float C[3],A[3],B[3];
{ C[0] = A[1]*B[2] - A[2]*B[1];
  C[1] = A[2]*B[0] - A[0]*B[2];
  C[2] = A[0]*B[1] - A[1]*B[0]; }

void Cross_Prod_Double(C,A,B)  
 double C[3],A[3],B[3];
{ C[0] = A[1]*B[2] - A[2]*B[1];
  C[1] = A[2]*B[0] - A[0]*B[2];
  C[2] = A[0]*B[1] - A[1]*B[0]; }




void Add_Atom(Ahead,pos,R,chain,Natom,c1,c2,c3)
 struct ATOM *Ahead;
 float pos[3],R;
 char chain; 
 int Natom;
 struct ATOM *c1,*c2,*c3; /* Contacted Protein Atoms */
{
 struct ATOM *an,*nn;
 int i;

  nn = Ahead;
  while (nn->next!=NULL) nn = nn->next;
  nn->next = (struct ATOM*)malloc(sizeof(struct ATOM));
  an = nn->next;
  an->next = NULL;
  an->prev = nn;

  for (i=0;i<3;++i) an->Pos[i] = pos[i];
  an->R = R;
  an->num = an->rnum = Natom;
  sprintf(an->Rnum,"%4d ",an->rnum%10000); 
  sprintf(an->Atom,"CA ");
  sprintf(an->Resi,"PRB");
  an->AHtype = 'H';
  if ((an->num %500)==0) {printf("#Add_Atom %d\n",an->num);}
  sprintf(an->Anum,"%5d",an->num);
  an->Chain = chain;
  an->tFactor = 0.0;
  an->mark = 0;
  an->cluster_num = 0; an->contact = 0;
  an->conatm1 = c1; an->conatm2 = c2; an->conatm3 = c3; 
  an->rnext = an->rprev = NULL;
  an->res = NULL;

} /* end of Add_Atom() */









void Sort_Probe_By_Pocket_Cluster_Number(Phead)
  struct ATOM *Phead;
{
 struct ATOM *an,*bn,*aprev, TmpHead;
 int c,Natom;
 char chain_str[128];

 sprintf(chain_str,"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
 printf("#Sort_Probe_By_Pocket_Cluster_Number(Phead)\n");

 /** (1) Move ATOM with cluster c into Phead **/
 TmpHead.next = Phead->next;
 if (Phead->next != NULL)  TmpHead.next->prev = &TmpHead;
 Phead->next = NULL;
 for (c=1;c<=254;++c){
   an = &TmpHead;
   while (an->next != NULL){ 
     an = an->next;
     if (an->cluster_num==c){
       bn = Phead;
       while (bn->next!=NULL) bn = bn->next;
       aprev = an->prev;
       if (an->prev!=NULL) an->prev->next = an->next;
       if (an->next!=NULL) an->next->prev = an->prev;
       bn->next = an;
       an->prev = bn;
       an->next = NULL;
       an = aprev;
     }
  }
 }

 Natom = 0;
 /** (2) Renumber atom number **/
 an = Phead;
 while (an->next !=NULL){
   an = an->next;
   ++Natom;
   an->num = Natom;
   sprintf(an->Anum,"%5d",an->num);
   if (an->cluster_num < strlen(chain_str)) an->Chain = chain_str[an->cluster_num]; else an->Chain = ' ';
    sprintf(an->Rnum,"%4d ",an->cluster_num%10000);  
 }
} /* end of Sort_Probe_By_Cluster_Number() */




int Check_Crash_Pos_and_Atoms(pos,R,Ahead,an,bn,cn)
 float pos[3],R;
 struct ATOM *Ahead;
 struct ATOM *an,*bn,*cn;
{
 /*
   Naive crash check between pos[3] and protein atoms contained in Ahead.
   if crash pos and atoms, then return 1, otherwise return 0.
 */
 struct ATOM *atom;
 int m; 
 float d[3],D,Dcrash_permit;

 Dcrash_permit = 0.0001;

 atom = Ahead;
 while (atom->next != NULL){
   atom = atom->next;
    if ((atom->num != an->num) && (atom->num != bn->num) && (atom->num != cn->num)){ 
     for (m=0;m<3;++m) d[m] = pos[m] - atom->Pos[m]; 
     D = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
     if (D>0.0) D = sqrt(D);
     if ((R + atom->R - D)> Dcrash_permit){ return(1); } 
   }  
 }

  return(0);

} /* end of Check_Crash_Pos_and_Atoms() */








int Check_Crash_Pos_and_NEIGHBOR_Atoms(pos,R,NeiHead,an,bn,cn)
 float pos[3],R;
 struct NEIGHBOR *NeiHead;
 struct ATOM *an,*bn,*cn;
{
 /*
   Naive crash check between pos[3] and protein atoms contained in NeiHead.
   if crash pos and atoms, then return 1, otherwise return 0.
 */
 struct NEIGHBOR *nei; 
 struct ATOM *atom;
 int m; 
 float d[3],D,Dcrash_permit;

 /* printf("#Check_Crash_Pos_and_NEIGHBOR_Atoms()\n"); */

 Dcrash_permit = 0.0001;

 nei = NeiHead;
 while (nei->next != NULL){
   nei  = nei->next;
   atom = nei->atom;
   if ((atom->num != an->num) && (atom->num != bn->num) && (atom->num != cn->num)){ 
     for (m=0;m<3;++m) d[m] = pos[m] - atom->Pos[m]; 
     D = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
     if (D>0.0) D = sqrt(D);
     if ((R + atom->R - D)> Dcrash_permit) { return(1); }  
   }

  }

  return(0);

} /* end of Check_Crash_Pos_and_NEIGHBOR_Atoms() */







void Malloc_FLOAT2DMAP(M,Nmalloc)
 struct FLOAT2DMAP *M;
 int Nmalloc;
{ int i,NN;
  double Mbyte;

  M->N = Nmalloc;
  NN = Nmalloc + 1;
  Mbyte = (double)sizeof(float)*M->N * M->N/1000.0/1000.0;
  printf("#Malloc_FLOAT2DMAP(N %d) --> %.2f Mbyte\n",M->N,Mbyte);

  if (Mbyte > PAR.MAX_MEMORY_MEGA){
   printf("#ERROR:Malloc_FLOAT2DMAP():%.2f Mbyte is larger than MAX (%f Mbyte).\n",
      Mbyte,PAR.MAX_MEMORY_MEGA);
   exit(1);
  }

  M->m = (float **)malloc(sizeof(float *)*NN);

  for (i=0;i<NN;++i){
     M->m[i] = (float *)malloc(sizeof(float)*NN);
  }

} /* end of Malloc_FLOAT2DMAP() */


void Free_FLOAT2DMAP(M)
 struct FLOAT2DMAP *M;
{
 int i;
 for (i=0;i<M->N;++i) free(M->m[i]);
 free(M->m);
 M->N = 0;
}



void Cal_Distance_FLOAT2DMAP(Ahead,Dmat)
 struct ATOM   *Ahead;
 struct FLOAT2DMAP *Dmat;
{
 struct ATOM *an,*bn;
 float d[3],DD,D;
 int i,a,b;

 printf("#Cal_Distance_MATRIX()\n"); 
 /*
 an = Ahead;
 while (an->next != NULL){
  an = an->next;
  printf("%d\n",an->num); fflush(stdout);
 } 
 */
 an = Ahead;
 while (an->next != NULL){
  an = an->next;
  bn = an;
  while (bn->next != NULL)
  {
   bn = bn->next;
   for (i=0;i<3;++i) d[i] = an->Pos[i] - bn->Pos[i];
   DD = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
   if (DD>0.0) D = sqrt(DD); else D = 0.0;
   a = an->num;
   b = bn->num;
   if ((a<0)|| (a>Dmat->N))
   { printf("#ERROR:Atom %d is out of matrix range [%d..%d]\n",a,0,Dmat->N-1);
     exit(1);  }
   if ((b<0)|| (b>Dmat->N))
   { printf("#ERROR:Atom %d is out of matrix range [%d..%d]\n",b,0,Dmat->N-1);
     exit(1);  }

   Dmat->m[a][b] = Dmat->m[b][a] = D;
   /* printf("Dmat %d %d %f\n",a,b,Dmat->m[a][b]);  */
  } /* bn */

 } /* an */

} /* end of Cal_Distance_FLOAT2DMAP() */



void Make_NEIGHBOR_list(Ahead,Dmap,Rprobe)
 struct ATOM *Ahead;
 struct FLOAT2DMAP *Dmap;
 float  Rprobe;
{
 struct ATOM *an,*bn;
 float Dthre,Rprobe2;
 /*
   Criteria for neighbor 
    Distance is less than (Ra + 2*Rprobe + Rb)
  */

 /* (1) Initialize */
 Rprobe2 = Rprobe * 2.0;
 an = Ahead;
 while (an->next != NULL){
  an = an->next;
  an->NeiHead.prev = an->NeiHead.next = NULL;
 }


 /* (2) make NEIGHBOR list */
 an = Ahead;
 while (an->next != NULL){
  an = an->next;
  Dthre = Rprobe2 + an->R;
  bn = an;
  while (bn->next != NULL){
    bn = bn->next;
    if ((an->num != bn->num)&&(Dmap->m[an->num][bn->num] <= (an->R + bn->R + Rprobe2))){
      Add_NEIGHBOR_to_Atom(an,bn);
      Add_NEIGHBOR_to_Atom(bn,an);   
    }
  } /* bn */

 } /* an */

} /* end of Make_NEIGHBOR_list() */



void Add_NEIGHBOR_to_Atom(an,neighbor)
 struct ATOM *an;
 struct ATOM *neighbor;
{
 struct NEIGHBOR *nn;

 nn = &(an->NeiHead);
 while (nn->next != NULL) nn = nn->next;

 nn->next = (struct NEIGHBOR*)malloc(sizeof(struct NEIGHBOR));
 nn->next->next = NULL;
 nn->next->prev = nn;
 nn->next->atom = neighbor;

} /* end of Add_NEIGHBOR_to_Atom() */



void Remove_Probe_Atoms_Larger_Specified_Rinaccess(Phead,threRinacc)
 struct ATOM *Phead;  /* Probe atoms */
 float  threRinacc;
{
 struct ATOM *pn;
 printf("#Remove_Probe_Atoms_Larger_Specified_Rinaccess(threRinacc %f)\n",threRinacc);
 pn = Phead; 
 while (pn->next != NULL){
   pn = pn->next;
   if (pn->Rinacc > threRinacc){
     pn->prev->next = pn->next;
     if (pn->next != NULL) pn->next->prev = pn->prev; 
   }
 }
} /* end of Remove_Probe_Atoms_Larger_Specified_Rinaccess() */



void Remove_Probe_Atoms_Smaller_Specified_Pocketness(Phead,threPocketness)
 struct ATOM *Phead;  /* Probe atoms */
 float  threPocketness;
{
 struct ATOM *pn;
 printf("#Remove_Probe_Atoms_Smaller_Specified_Pocketness(threRinacc %f)\n",threPocketness);
 pn = Phead; 
 while (pn->next != NULL){
   pn = pn->next;
   if (pn->pocketness <= threPocketness){
     pn->prev->next = pn->next;
     if (pn->next != NULL) pn->next->prev = pn->prev; 
   }
 }
} /* end of Remove_Probe_Atoms_Larger_Specified_threPocketness() */



void Remove_Probe_Atoms_Smaller_Specified_tFactor(Phead,thre_tFactor)
 struct ATOM *Phead;  /* Probe atoms */
 float  thre_tFactor;
{
 struct ATOM *pn;
 printf("#Remove_Probe_Atoms_Smaller_Specified_tFactor(thre_tFactor %f)\n",thre_tFactor);
 pn = Phead; 
 while (pn->next != NULL){
   pn = pn->next;
   if (pn->tFactor < thre_tFactor){
     pn->prev->next = pn->next;
     if (pn->next != NULL) pn->next->prev = pn->prev; 
   }
 }
} /* end of Remove_Probe_Atoms_Smaller_Specified_tFactor() */




void Write_Spherical_Probes_in_UCSF_DOCK_format(ofname,Phead)
 char *ofname;
 struct ATOM *Phead;
{
  FILE *fp;
  struct ATOM *an; 
  int  c,init,Natom_cluster[256];
  printf("#Write_Spherical_Probes_in_UCSF_DOCK_format()-->'%s'\n",ofname);
  
/*
>> FILE FORMAT EXAMPLE <<
DOCK 3.5 receptor_spheres
cluster     1   number of spheres in cluster   160
   23  12.41121  29.98135  57.55564   1.818  277 0  0
   24  12.84669  28.55700  60.03390   1.902  431 0  0
   67   5.58404  50.91005  59.97029   1.411   96 0  0
*/

  fp = fopen(ofname,"w");
  fprintf(fp,"DOCK 3.5 receptor_spheres  (generated by ghecom program)\n");
  /** count Natom_cluster **/
  for (c=0;c<256;++c) Natom_cluster[c] = 0;
  an = Phead;
  while (an->next!=NULL){
    an = an->next;
    if (an->cluster_num<256) Natom_cluster[an->cluster_num] += 1; 
  } 
  /** output **/
  an = Phead;
  init = 1;
  while (an->next!=NULL){
    an = an->next;
    if ((init==1) || (an->prev->cluster_num != an->cluster_num)){
      fprintf(fp,"cluster   %3d   number of spheres in cluster",an->cluster_num);
      if (an->cluster_num<256) fprintf(fp,"%6d",Natom_cluster[an->cluster_num]);
      fprintf(fp,"\n");
    }
    fprintf(fp,"%5s%10.5f%10.5f%10.5f%8.3f%5s 0  0\n",
      an->conatm1->Anum,an->Pos[0],an->Pos[1],an->Pos[2],an->R,an->conatm2->Anum);
    init = 0;
  }
  fclose(fp);
} /* end of Write_Spherical_Probes_in_UCSF_DOCK_format() */




void Label_Probe_Atoms_By_CHAR3DMAP(AtmHead,X,tar_bitmask,ValType,AssignType,hit_val,nonhit_val)
 struct ATOM *AtmHead;
 struct CHAR3DMAP *X;
 unsigned char tar_bitmask; /* if tar_bitmask==0, then true if map is non-zero. */
 char  ValType;   /* 'T':tFactor, R'esidue number,'C'luster_num */
 char AssignType;  /* 'M':from map, 'A':assigned value by 'hit_val' */
 float hit_val;    /* assigned value for hit probes */ 
 float nonhit_val; /* assigned value for non hit probes */
{
 struct ATOM *an;
 int i,m,maxNmap,maxm,x,y,z,minN[3],maxN[3];
 char hit;
 float D[3],DD[3],ddis,RR;
 int Nmap_count[256];

 printf("#Label_Probe_Atoms_By_CHAR3DMAP(AtmHead,X,tar_bitmask %d AssignType %c ValType %c)\n",tar_bitmask,AssignType,ValType);
 an = AtmHead;
 while (an->next != NULL){ 
   an = an->next;
   an->tFactor = 0.0;
   RR = an->R * an->R;

   for (i=0;i<3;++i) {
      minN[i] = (int)floor((an->Pos[i] - an->R  - X->OrigPos[i])/X->grid_width) - 1;
      maxN[i] = (int)ceil((an->Pos[i]  + an->R  - X->OrigPos[i])/X->grid_width) + 1;
      if (minN[i]<0) minN[i] = 0;
      if (maxN[i]>=X->N[i]) maxN[i] = X->N[i]-1;
   }

   for (m=0;m<256;++m){ Nmap_count[m] = 0; }

   hit = 0;
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
           /* printf("%d %d %d map %d\n",x,y,z,X->map[x][y][z]); */
           if( ((tar_bitmask==0)&&(X->map[x][y][z]>0))||((tar_bitmask>0)&&(X->map[x][y][z]&tar_bitmask) == tar_bitmask ) ) { 
             hit = 1; 
             Nmap_count[X->map[x][y][z]] += 1;
             /* x = X->N[0]; y = X->N[1]; z = X->N[2];  */
           }
        }
       } /* z */
      } /* y */
    } /* x */

    if (AssignType=='M'){
      if (hit==1){
        maxNmap = 0;
        maxm = 0;
        for (m=0;m<256;++m){ 
          if (Nmap_count[m]>maxNmap){maxm = m; maxNmap = Nmap_count[m];} 
        }
        if (ValType=='T'){ an->tFactor    = (float)maxm;}
        if (ValType=='R'){ sprintf(an->Rnum,"%3d ",(int)maxm);}
        if (ValType=='C'){ 
          an->cluster_num = maxm;
          sprintf(an->Rnum,"%3d ",(int)maxm);
        }
     }
     else{
       if (ValType=='T'){ an->tFactor    = nonhit_val;}
       if (ValType=='R'){ sprintf(an->Rnum,"%3.0f ",nonhit_val);}
       if (ValType=='C'){ 
         an->cluster_num = (int)nonhit_val;
         sprintf(an->Rnum,"%3.0f ",nonhit_val);
       }
     }
    }

    if (AssignType=='A'){
      if (ValType=='T'){
        if (hit==1) an->tFactor = hit_val; else  an->tFactor = nonhit_val; 
      }

      if (ValType=='R'){
        if (hit==1) an->Rinacc = hit_val; else an->Rinacc = nonhit_val; 
      }
    }


  } /* an */



} /* end of Label_Probe_Atoms_By_CHAR3DMAP() */
