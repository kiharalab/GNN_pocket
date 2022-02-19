/*
 <Radius.c>
 
 Functions for assigning radius to atoms


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
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include "globalvar.h" 
#include "pdbstruct.h" 
#include "defradius_Chothia1976.h" 

/*** FUNCTIONS (GLOBAL) ***/
void Assign_Radius_to_ATOMs();
void Assign_Weight_to_ATOMs();


/*** FUNCTIONS (LOCAL) ***/
static int Match();
static void Get_Part_Of_Line();
static float vdW_Radius_Bondi1964();
static float vdW_Radius_Chothia1976();
static float Atomic_Weight();




void Assign_Radius_to_ATOMs(Ahead,RadiusType,Runified)
  struct ATOM   *Ahead;
  char  RadiusType; /* 'C'hothia, 'B'ond, 'U'nified, 'O'ccupancy */
  float Runified;
{
 struct ATOM *an;

 an = Ahead;
 while (an->next != NULL){
  an = an->next;
  an->R = 0.0;
  if (RadiusType=='B'){
    an->R =  vdW_Radius_Bondi1964(an->element);
  }
  else if (RadiusType=='C'){
    an->R =  vdW_Radius_Chothia1976(an->Atom, an->Resi);
  }
  else if (RadiusType=='U'){
    an->R  = Runified;
  }
  else if (RadiusType=='O'){
    an->R  = an->Occup;
  }
  an->RR = an->R * an->R;
 }

} /* end of Assign_Radius_to_ATOMs() */




int Match(Pat,plen,Str,offset)
 char *Pat;
 int  plen;
 char *Str;
 int offset;
{
  /*** Strict Pattern Matching ***/
  int i,j,ok,slen;

  slen = strlen(Str);
  if ((plen==0)||(slen==0)) return(0);

  if (slen>=plen){
   for (i= offset;i<(slen - plen+1);++i){ ok = 1; j = 0;
      while ((ok==1) && (j<plen)){
       if ((Pat[j]!='x')&&(Pat[j] != Str[i+j])) { ok = 0; }
       ++j; }

     if (ok==1) return(1);
   }

  }
 return(0);

} /* end of Match */




void Get_Part_Of_Line(part,line,s,e)
  char *part;
  char *line;
  int  s,e;
{
 int i,E,L;
 L = strlen(line)-1;
 if (line[L] == '\n') L -= 1;
 if (e>L) E = L; else E = e;
 for (i=s;i<=E;++i) part[i-s] = line[i];
 part[E-s+1] = '\0';

} /* end of Get_Part_of_Line() */





float vdW_Radius_Bondi1964(ele)
 char *ele;
{
  int L;
  float radius;
/*
 * >> Van der Waals radii taken from Bondi's compilation (1964). <<
 * Bondi, A. (1964). "Van der Waals Volumes and Radii".  J. Phys. Chem., 1964, 68 (3), pp 441–451 DOI: 10.1021/j100785a001
 *
 * Hydrogen(H) 1.20 
 * Carbon(C) 1.70 
 * Nitrogen(N) 1.55 
 * Oxygen(O) 1.52 
 * Fluorine(F) 1.47 
 * Phosphorus(P) 1.80 
 * Sulfur(S) 1.80 
 * Chlorine(Cl) 1.75 
 * Copper(Cu) 1.4 
 * */
  L = strlen(ele);
       if ((L==1)&&(strncmp(ele,"H",1)==0)){  radius = 1.20;}
  else if ((L==1)&&(strncmp(ele,"C",1)==0)){  radius = 1.70;}
  else if ((L==1)&&(strncmp(ele,"N",1)==0)){  radius = 1.55;}
  else if ((L==1)&&(strncmp(ele,"O",1)==0)){  radius = 1.52;}
  else if ((L==1)&&(strncmp(ele,"F",1)==0)){  radius = 1.47;}
  else if ((L==1)&&(strncmp(ele,"S",1)==0)){  radius = 1.80;}
  else if ((L==1)&&(strncmp(ele,"P",1)==0)){  radius = 1.80;}
  else if ((L==1)&&(strncmp(ele,"I",1)==0)){  radius = 1.98;}
  else if ((L==1)&&(strncmp(ele,"K",1)==0)){  radius = 2.75;}
  else if ((L==2)&&(strncmp(ele,"Cl",2)==0)){ radius = 1.75;}
  else if ((L==2)&&(strncmp(ele,"CL",2)==0)){ radius = 1.75;}
  else if ((L==2)&&(strncmp(ele,"Cu",2)==0)){ radius = 1.40;}
  else if ((L==2)&&(strncmp(ele,"CU",2)==0)){ radius = 1.40;}
  else if ((L==2)&&(strncmp(ele,"Br",2)==0)){ radius = 1.85;}
  else if ((L==2)&&(strncmp(ele,"BR",2)==0)){ radius = 1.85;}
  else if ((L==2)&&(strncmp(ele,"Zn",2)==0)){ radius = 1.39;}
  else if ((L==2)&&(strncmp(ele,"ZN",2)==0)){ radius = 1.39;}
  else if ((L==2)&&(strncmp(ele,"Mg",2)==0)){ radius = 1.73;}
  else if ((L==2)&&(strncmp(ele,"MG",2)==0)){ radius = 1.73;}
  else if ((L==2)&&(strncmp(ele,"Na",2)==0)){ radius = 2.27;}
  else if ((L==2)&&(strncmp(ele,"NA",2)==0)){ radius = 2.27;}
  else{ radius= 1.70; /* radius for Carbon */ }
  return(radius);

} /* end of Set_radius() */



float vdW_Radius_Chothia1976(Atom, Resi)
  char *Atom;  /* atom name of PDB     */
  char *Resi;  /* residue name of PDB */
{
 char line[1024],pat_atom[5],pat_resi[4],pat_val[32];
 int i,L,hit;
 float radius;
/*
           1 
 01234567890123 
" N   ALA 1.65",
" CA  ALA 1.87",
" C   ALA 1.76",
" O   ALA 1.40",
" CB  ALA 1.87",
" N   ARG 1.65",
" CB  ARG 1.87",
:
" O1D HEM 1.35",
" O2D HEM 1.35",
" O   HOH 1.40",
"xCxx xxx 1.87",
"xNxx xxx 1.65",
"xOxx xxx 1.40",
"xPxx xxx 1.90",
"xSxx xxx 1.85",
"xxxx xxx 1.50"
*/
 i = 0;
 hit = 0;
 radius  = 0.0;
 while ((i<NRADLINE_CHOTHIA1976) && (hit==0)){ 
   line[0] = '\0';
   sprintf(line,"%s",RADLINE_CHOTHIA1976[i]);
   L = strlen(line);
   if ((L>5) && (line[0] != '#')){
     /* printf("#'%s'\n",line); */
     Get_Part_Of_Line(pat_atom,line,0,3);
     Get_Part_Of_Line(pat_resi,line,5,7);
     Get_Part_Of_Line(pat_val, line,8,L-1);
     if ((Match(pat_resi,3,Resi,0)==1) && (Match(pat_atom,3,Atom,0)==1) ){
        hit  = 1;
        radius = atof(pat_val);
     }
   } 
   i += 1;  
 } 
 return(radius);

} /* end of vdW_Radius_Chothia1976() */


void Assign_Weight_to_ATOMs(Ahead)
 struct ATOM   *Ahead;
{
 struct ATOM   *an;
 char ele[4];

 an = Ahead;
 while (an->next != NULL){
   an = an->next;
   ele[0] = an->Atom[1];
   ele[1] = '\0';
   an->Weight = Atomic_Weight(ele);
   /* printf("#'%s' '%s' ele '%s' W %f\n",an->Resi,an->Atom,ele,an->Weight);  */
 }
} /* end of Assign_Weight() */



float Atomic_Weight(ele)
  char *ele; /* element name */
{
 /*
 *   These values are taken from 
      Urabyoushi I. Genshiryou-hyou (1985)
       in "Rikagakujiten dai-4-ban".
  */
       if (strncmp(ele,"C",1)==0){ return(12.011);}
  else if (strncmp(ele,"N",1)==0){ return(14.00674);}
  else if (strncmp(ele,"O",1)==0){ return(15.9994);}
  else if (strncmp(ele,"S",1)==0){ return(32.066);}
  else if (strncmp(ele,"P",1)==0){ return(30.973762);}
  /*
      else if (strncmp(ele," H",2)==0){ return(1.0);}
  */
 /**  
       NOTE: Weight of the Hydrogen atom is assigned to ZERO.
            Because, GMM does not take account for hydrogens.
                                     
 **/
  else if (strncmp(ele," H",2)==0){ return(0.0);}
  else { return(12.011);}
} /* end of Atomic_Weight() */


