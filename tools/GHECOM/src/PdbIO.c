/*

<PdbIO.c>

 functions for input/output of PDB files
 using list structure (ATOM and RESIDUE node). 

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
#include <math.h>
#include <time.h>
#include "globalvar.h" 
#include "pdbstruct.h" 
#include "ReadOpt.h" 
#include "PdbIO.h" 

/*** Functions (GLOBAL) ***/
void Read_PDB_File();
void Make_Residue_List();
void Find_Filename_Core();
void Find_Filename_Core_with_Directory();
void Get_Part_Of_Line();
char *Get_Date_String_PDB();
int  Number_Of_ATOMs();
int  Number_Of_RESIDUEs();
int  Number_Of_CHAINs();
int  Number_Of_MODELs();
int  Free_ATOMs();
void Set_Region_Of_Atoms();
void Set_Constant_tFactor();
void Renumber_Atom_Number();
void Renumber_Residue_Number();
void Write_PDB_File();
void Write_PDB_File_Residue();
void Add_Atom_To_Residue_Head();
void Add_Atom_To_Residue_Tail();
void Split_to_Words();
void Remove_Space();
void Add_String_to_STRING_NODE();
void Show_STRING_NODE();
void Write_One_Point_In_PDB_Format();
int  Find_XYZ_By_Pattern();
int  Find_XYZ_By_Pattern_In_PDB_File();
void Make_Copied_Atoms();
double Set_Molecule_Gcenter_to_Origin();
void   Rotate_Molecule();
void   InvRotate_Molecule();
void   Translate_Molecule();
void Delete_Doubling_altLoc_Atoms();
void Change_N_CA_C_O_Hetatom_to_Atom();
void Delete_HETATMs_or_ATOMs();


/*** Functions (LOCAL) ***/
static int Find_Same_Res_Rnum_Atom_Before();
static void Renumber_Atom_num();
static int ChainIDStr_Match();
/* static struct ATOM *Find_Atom_with_Anum(); */
static void Element_from_PDB_Atomname();
static int Check_Element_String_Or_Not();


#define MAX_NELEMENT         87 

static char element_list[MAX_NELEMENT][3] =
{"H",
 "HE","LI","BE","B","C","N","O","F","NE",
 "NA","MG","AL","SI","P","S","CL","AR","K",
 "CA","SC","TI","V","CR","MN","FE","CO","NI","CU","ZN","GA","GE","AS","SE","BR",
 "KR","RB","SR","Y","ZR","NB","MO","TC","RU","RH","PD","AG","CD","IN","SN","SB","TE","I","XE",
 "CS","BA",
 "LA","CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER","TM","YB","LU",
 "HF","TA","W", "RE","OS","IR","PT","AU","HG","TL","PB","BI","PO","AT","RN",
 "?"
};


static char element_list_lower[MAX_NELEMENT][3] =
{"H",
 "He","Li","Be","B","C","N","O","F","Ne",
 "Na","Mg","Al","Si","P","S","Cl","Ar","K",
 "Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br",
 "Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
 "Cs","Ba",
 "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
 "Hf","Ta","W", "Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
 "?"
};


void Read_PDB_File(fname, Ahead,pdbid, AtomSelect, ChainIDStr, AtomHetSelect, ReadTER, macro_exclude, MODELtype)
 char *fname;
 struct ATOM *Ahead;
 char  *pdbid;       /* pdbid (read from HEADER of fname) */
 char   AtomSelect;    /* 'A':only Calpha  'B':only Cbeta 'H':read only heavy atom */
 char  *ChainIDStr;    /* Chain Identifier String. Multiple Chain can be assigned. (ex). "-", "A", "AB", "ABCD", "1234") */
 char  AtomHetSelect; /* 'A' : Only read ATOM, 'H':Only read HETATM  'B'oth */
                      /* 'h':read HETATM or non 'ChainID' matched ATOMs */
 char  ReadTER;       /* Read and store 'TER' line  'T' or 'F' */
 char  macro_exclude; /* exclude macromolecule('pro','dna','rna')'T' or 'F' */
 char  MODELtype;     /* 'S':read single model 'M':read multiple model */
{
 FILE *fp;
 char line[256],B[32],Aname[32],Rname[32],atomtype,chain,altloc, command[MAX_FILE_NAME];
 char atminfo[32],atminfo0[32],end;
 int accept,anum, L, model_num;
 struct ATOM *an;

 printf("#Read_PDB_File('%s' AtomSelect %c ChainIDStr '%s' AtomHetSelect %c ReadTER %c macro_exclude %c)\n",
  fname,AtomSelect,ChainIDStr,AtomHetSelect,ReadTER,macro_exclude);
 anum = 0;  
 Ahead->next  = NULL;
 Ahead->prev  = NULL;
 Ahead->AHtype = 'S';
 an = Ahead; end = 0;

 L = strlen(fname);
 if ((L>2)&&(fname[L-2]=='.')&&(fname[L-1]=='Z')){
   sprintf(command,"uncompress -c %s",fname);
   fp = popen(command,"r");
 }
 else if ((L>2)&&(fname[L-3]=='.')&&(fname[L-2]=='g')&&(fname[L-1]=='z')){
   sprintf(command,"zcat %s",fname);
   fp = popen(command,"r");
 }
 else{
   fp = fopen(fname,"r");
 }

 if (fp==NULL){printf("#ERROR:Can't open pdbfile \"%s\"\n",fname); exit(1);} 

 pdbid[0] = '\0';

 model_num = 0;
 
 while ((feof(fp)==0)&&(end==0)){
  line[0] = '\0';
  fgets(line,255,fp);
/*
          1         2         3         4         5         6         7
01234567890123456789012345678901234567890123456789012345678901234567890123456789
HEADER    GLYCOSYLTRANSFERASE                     07-JUN-95   1XYZ
*/
  if (strncmp(line,"HEADER",5)==0){ Get_Part_Of_Line(pdbid,line,62,65);}

       if (strncmp(line,"ATOM",4)==0){   atomtype = 'A';}
  else if (strncmp(line,"HETATM",6)==0){ atomtype = 'H';}
  else if (strncmp(line,"TER",3)==0){    atomtype = 'T';}
  else atomtype = ' ';
  Get_Part_Of_Line(Rname,line,17,19);

  if (strncmp(line,"MODEL",5)==0){
/*
          1
01234567890123
MODEL        1
MODEL        2
:
MODEL       29
*/
    Get_Part_Of_Line(B,line,5,14);
    model_num = atoi(line); 
  }


  if ((MODELtype=='S')&&(strncmp(line,"ENDMDL",6)==0)){ end = 1; }

  if ((ReadTER == 'T')&&(atomtype=='T')){
     an->next = (struct ATOM *)malloc(sizeof(struct ATOM));
     an->next->prev = an;
     an = an->next;
     an->next = NULL;
     an->Pos[0]=an->Pos[1] = an->Pos[2] = 0.0;
     an->Anum[0] = an->Atom[0] = an->Resi[0] = an->Rnum[0] = '\0'; 
     an->Occup = an->tFactor = 0.0; 
     an->Chain  = ' '; an->altLoc = ' '; an->num   = 0;
     an->res = NULL;
     an->rnext = an->rprev = NULL;
     an->AHtype = atomtype;
     an->shellAcc = an->Rinacc = 0.0;
     an->conatm1 =  an->conatm2 =  an->conatm3 =  NULL; an->contact = 0;
     ++anum; 
     accept = 0;
  }

  accept = 1; 

  if ((atomtype != 'A')&&(atomtype!='H')){ accept = 0; } 
/*
  if ((AtomHetSelect == 'A')&&(atomtype == 'H')){ accept = 0;}
  if ((AtomHetSelect == 'H')&&(atomtype == 'A')){ accept = 0;}
 */
  if ((macro_exclude=='T') && 
      ((strcmp(Rname,"dna")==0)||(strcmp(Rname,"rna")==0)||(strcmp(Rname,"pro")==0))){ accept = 0; }
 
  if (accept == 1){ 
    Get_Part_Of_Line(Aname,line,12,15);
    Get_Part_Of_Line(atminfo,line,12,25);
    chain = line[21]; 
    atminfo[4] = ' '; 
    /* printf("\"%s\"\n",atminfo); */
    altloc = line[16]; 
    
    /*** Decide Accept or Not ***/
/*
    if ((AtomHetSelect == 'h')&&(atomtype == 'A')&&(ChainIDStr_Match(ChainIDStr,chain)==1)){ accept = 0;}
 */
    if ((Aname[1]=='H')&&(AtomSelect=='H')){ accept = 0;}

    if (strncmp(Rname,"HOH",3)==0) accept = 0;
    if (strncmp(Rname,"DOD",3)==0) accept = 0;

    if (AtomSelect == 'A')
    { if (strncmp(Aname," CA ",4)==0) accept = 1; else accept = 0; }

    if (AtomSelect == 'B'){ 
      if (strncmp(Aname," CB ",4)==0) accept = 1; else accept = 0; 
      if ((strncmp(Rname,"GLY",3)==0)&&(strncmp(Aname," CA ",4)==0)){ accept = 1; }
     }

    if (ChainIDStr_Match(ChainIDStr,chain)!=1){ accept = 0;}
/*
    if ((AtomHetSelect=='A')&&(ChainIDStr_Match(ChainIDStr,chain)!=1)){ accept = 0;}
    if ((AtomHetSelect=='H')&&(ChainIDStr_Match(ChainIDStr,chain)!=1)){ accept = 0;}
 */
    if (accept==1){ 
/*    
      012345678901234567890123456789012345678901234567890123456789012345678901234567890
               1         2         3         4         5         6         7         8  
      HETATM 2225 BR1  HMD A 400      23.120  26.871  28.136  1.00 63.13          BR  
      HETATM 2226  O1  HMD A 400      25.184  31.320  30.286  1.00 45.50           O  
      HETATM 2227  O2  HMD A 400      30.537  29.959  25.467  1.00 43.12           O  
      HETATM 2228  N1  HMD A 400      24.542  29.157  28.792  1.00 48.35           N  
      HETATM 5274  C4  GCP A1749      14.439  21.481  -2.931  1.00  7.13           C  
      HETATM 5275 MG    MG A1750      10.612  12.184  -8.148  1.00  5.01          MG  
      HETATM 5276  PG  GCP D1749      21.687  11.358  17.910  1.00 11.10           P  
      HETATM 5277  O1G GCP D1749      20.774  10.642  16.983  1.00  7.76           O  
      HETATM 5278  O2G GCP D1749      21.638  11.013  19.350  1.00 14.21           O  
      HETATM 5279  O3G GCP D1749      23.127  11.469  17.544  1.00 11.34           O  
      HETATM 5280  C3B GCP D1749      21.121  13.108  17.927  1.00  9.79           C  
      HETATM 5281  PB  GCP D1749      19.444  13.422  18.583  1.00  8.41           P  
      HETATM 2050  FE  670 A   1      20.305  11.132  13.458  1.00  0.00  B   300 FE
      HETATM 2051  NAA 670 A   1      17.394   8.763  12.431  1.00  0.00  B   300  N
      HETATM 2052  NAB 670 A   1      16.774   7.874  13.019  1.00  0.00  B   300  N

 */


      an->next = (struct ATOM *)malloc(sizeof(struct ATOM));
      an->next->prev = an;
      an = an->next;
      an->next = NULL;
      Get_Part_Of_Line(B,line,30,37); an->Pos[0] = atof(B);
      Get_Part_Of_Line(B,line,38,45); an->Pos[1] = atof(B);
      Get_Part_Of_Line(B,line,46,53); an->Pos[2] = atof(B);
      Get_Part_Of_Line(an->Anum,line,6,10);
      Get_Part_Of_Line(an->Atom,line,12,15);
      Get_Part_Of_Line(an->Resi,line,17,19);
      Get_Part_Of_Line(an->Rnum,line,22,26);
      Get_Part_Of_Line(B,line,54,59); an->Occup  = atof(B);
      an->tFactor = 0.000;
      Get_Part_Of_Line(B,line,60,65); an->tFactor = atof(B);
      an->Chain  = chain; 
      an->altLoc = altloc; 
      an->num   = anum;
      an->model_num = model_num;

      Get_Part_Of_Line(B,line,76,77);
      if ((B[0]==' ')&&(B[1]!=' ')){ B[0] = B[1]; B[1] = '\0'; }
      if (Check_Element_String_Or_Not(B)==1){
        sprintf(an->element,"%s",B);
      }
      else{
        Element_from_PDB_Atomname(an->element,an->Atom);
      }
      an->res = NULL;
      an->rnext = an->rprev = NULL;
      an->AHtype = atomtype;
      an->shellAcc =  an->Rinacc = 0.0;
      an->conatm1 = an->conatm2 = an->conatm3 = NULL;  an->contact = 0;
      ++anum; 
/*
      printf("#AHtype %c %s %s %s %s element '%s' model_num %d\n",an->AHtype,an->Anum, an->Atom,an->Resi,an->Rnum,an->element,an->model_num); 
 */
    }

   sprintf(atminfo0,"%s",atminfo); 
  }

 } /* while */

 fclose(fp);

 Delete_Doubling_altLoc_Atoms(Ahead);


 if (PAR.Change_N_CA_C_O_Hetatom_to_Atom == 'T'){ Change_N_CA_C_O_Hetatom_to_Atom(Ahead);}

 if (AtomHetSelect == 'A'){
   Delete_HETATMs_or_ATOMs(Ahead,'H');
 }

 
} /* end of Read_PDB_File() */








void Write_PDB_File(fname,mode,Ahead,comment,CommentHead)
 char *fname;
 char mode; /* 'a'ppend, 'w'rite */
 struct ATOM *Ahead;
 char *comment;
 struct STRING_NODE *CommentHead;
{
 FILE *fp;
 struct ATOM *an;
 char core[MAX_FILE_NAME],buff[MAX_FILE_NAME],buff2[MAX_FILE_NAME],init;
 struct STRING_NODE *sn; 
 printf("#Write_PDB_File(mode %c) -> \"%s\"\n",mode,fname);
 if (fname[0]=='-') fp = stdout;
 else{
    if (mode=='a') fp = fopen(fname,"a");
             else  fp = fopen(fname,"w");
  if (fp==NULL) { printf("#ERROR:Can't write to \"%s\"\n",fname); exit(1);}
 }

 Find_Filename_Core(core,fname);
 Get_Part_Of_Line(buff,fname,0,40);
 Get_Part_Of_Line(buff2,core,0,9);
 fprintf(fp,"HEADER    %-40s%s   %-10s0XXX   1\n",buff,Get_Date_String_PDB(),buff2);

 fprintf(fp,"REMARK    DATE %s\n",Get_Date_String());
 if (PAR.COMMAND[0] != '\0') fprintf(fp,"REMARK  COMMAND %s\n",PAR.COMMAND);
 if (comment[0]!='\0') fprintf(fp,"REMARK  COMMENT %s\n",comment);
 
 sn = CommentHead;
 while (sn->next != NULL){ 
   sn = sn->next;
   fprintf(fp,"REMARK  COMMENT %s\n",sn->line); 
 }



 fprintf(fp,"REMARK    Occupancy [55-60] : Radius of the atom\n");

 /*** Write ATOM/HETATM ***/

/*
            1
  01234567890123
  MODEL        1
  MODEL        2
    :
  MODEL       29
*/

 init = 1;
 an = Ahead;
 while (an->next != NULL){
   an = an->next;
   if ((init==1) || (an->model_num != an->prev->model_num )){
     if (init==0){
       fprintf(fp,"ENDMDL\n");
     }
     fprintf(fp,"MODEL%9d\n",an->model_num);
   }
   init = 0;

        if (an->AHtype=='A'){fprintf(fp,"ATOM  ");}
   else if (an->AHtype=='H'){fprintf(fp,"HETATM");}
   else { fprintf(fp,"ATOM? ");}
   fprintf(fp,"%5s %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
        an->Anum,an->Atom,an->Resi,an->Chain,an->Rnum,
        an->Pos[0],an->Pos[1],an->Pos[2],an->R,an->tFactor);
 }
 fprintf(fp,"TER\n");
 fprintf(fp,"ENDMDL\n");

 if (fp!=stdout){
   fclose(fp);
 }
} /* end of Write_PDB_File() */













void Write_PDB_File_Residue(fname,Rhead,mode,comment,command)
 char *fname;
 struct RESIDUE *Rhead;
 char mode; /* 'a'ppend, 'w'rite */
 char *comment;
 char *command;
{
 FILE *fp;
 struct RESIDUE *rn;
 struct ATOM *an;
 int i,j,L;
 char line[512],subline[100];

 printf("#Write_PDB_File_Residue -> \"%s\"\n",fname);
 if (fname[0]=='-') fp = stdout;
 else
 {
    if (mode=='a') fp = fopen(fname,"a");
             else  fp = fopen(fname,"w");
  if (fp==NULL) { printf("#ERROR:Can't write to \"%s\"\n",fname); exit(1);}
 }

 sprintf(line,"FILENAME \"%s\"",fname);
 Get_Part_Of_Line(subline,line,0,40);
 fprintf(fp,"HEADER    %-40s%s   0XXX      0XXX   1\n",subline,Get_Date_String_PDB());

 if (comment[0] != '\0'){
  L = strlen(comment);
  i = j = 0;
  fprintf(fp,"REMARK #   ");
  while (i<L){ 
    if ((j>60)||(comment[i]=='\n')) { fprintf(fp,"\nREMARK #   "); j = 0;}
    if (comment[i]!='\n') fprintf(fp,"%c",comment[i]);
     ++i; ++j;
  } 
  fprintf(fp,"\n");
 }

 if (command[0] != '\0'){fprintf(fp,"REMARK COMMAND %s\n",command);}

 rn = Rhead;
 while (rn->next != NULL){
   rn = rn->next;
   an = &(rn->Ahead);
   while (an->rnext != NULL){
     an = an->rnext;
          if (an->AHtype=='A')  fprintf(fp,"ATOM  ");
     else if (an->AHtype=='H')  fprintf(fp,"HETATM");
     else fprintf(fp,"ATOM? ");
     fprintf(fp,"%5s %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
          an->Anum,an->Atom,an->Resi,an->Chain,an->Rnum,
          an->Pos[0],an->Pos[1],an->Pos[2],an->R,an->tFactor);
    } /* an */

 } /* rn */

 fprintf(fp,"TER\n");
 if (fp!=stdout) fclose(fp);

} /* end of Write_PDB_File_Residue() */




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



/*

char *Get_Date_String()
{
 time_t      now_t;
 struct tm  *loc_t;
 static char Mon[][4] = {"Jan","Feb","Mar","Apr","May","Jun",
                         "Jul","Aug","Sep","Oct","Nov","Dec"};
 static char string[64];
 now_t = time(NULL);
 loc_t = localtime(&now_t);

 sprintf(string,"%s %d,%d %d:%d:%d",
  Mon[loc_t->tm_mon],loc_t->tm_mday,loc_t->tm_year+1900,
  loc_t->tm_hour,loc_t->tm_min,loc_t->tm_sec);
 return(string);

} 
*/

char *Get_Date_String_PDB()
{
 time_t      now_t;
 struct tm  *loc_t;
 static char Mon[][4] = {"JAN","FEB","MAR","APR","MAY","JUN",
                         "JUL","AUG","SEP","OCT","NOV","DEC"};
 static char string[64];
 char day[16],year[16];
 int  Nyear;  
 now_t = time(NULL);
 loc_t = localtime(&now_t);

 if (loc_t->tm_mday <10) sprintf(day,"0%d",loc_t->tm_mday);
 else sprintf(day,"%2d",loc_t->tm_mday);

 if (loc_t->tm_year<100) Nyear = loc_t->tm_year;
   else                  Nyear = loc_t->tm_year - 100;

 if (Nyear < 10) sprintf(year,"0%d",Nyear);
            else sprintf(year,"%2d",Nyear);

 sprintf(string,"%s-%3s-%s",day, Mon[loc_t->tm_mon],year);
 return(string);

} /* end of Get_Date_String_PDB() */


void Find_Filename_Core(core,fname)
 char *core,*fname;
{
 int Wsta[100],Wend[100],Nword;
 char lastfile[128],hit;
 int i,dot_pos;

 Split_to_Words(fname,'/',&Nword,Wsta,Wend,100);
 Get_Part_Of_Line(lastfile,fname,Wsta[Nword-1],Wend[Nword-1]);
 i = strlen(lastfile)-1; 
 hit = 0;
 while ((i>0)&&(hit==0)){
  if ((lastfile[i] == '.') && (lastfile[i-1] != '.')){
    hit = 1; 
    dot_pos = i-1;
   }
  --i; 
 }

 if (hit==1){ 
   Get_Part_Of_Line(core,lastfile,0,dot_pos);
 }
 else{
   sprintf(core,"%s",lastfile);
 }
} /* end of Find_Filename_Core() */
 

void Find_Filename_Core_with_Directory(core,fname)
 char *core,*fname;
{
 char hit;
 int i,dot_pos;

 i = strlen(fname)-1; 
 hit = 0;
 while ((i>0)&&(hit==0)){
  if ((fname[i] == '.') && (fname[i-1] != '.')){
    hit = 1; 
    dot_pos = i-1;
  }
  --i; 
 }

 if (hit==1){ 
   Get_Part_Of_Line(core,fname,0,dot_pos);
 }
 else{
   sprintf(core,"%s",fname);
 }
} /* end of Find_Filename_Core_with_Directory() */

 
void Split_to_Words(str,splsym,Nword,Wsta,Wend,Nwordmax)
 char *str;          /* Input String */
 char splsym;        /* Symbol for split */
 int *Nword;         /* Number of words  */
 int Wsta[];         /* Start point of str for a wowd */
 int Wend[];         /* End point of str for a wowd */
 int Nwordmax;       /* Maxium number of Nword  */
{
 /* [Example]
  str = "abc/d/ef//ghi//"
  splsym = '/'
   -->
  Nword 4
  (Wsta,Wend) = {(0,2), (4,4), (6,7), (10,12)}
 */
 int i,L;
 
 L = strlen(str);
 *Nword = 0; i = 0;

 while ((i<L)&&(*Nword < Nwordmax)){
  if (str[i]!=splsym)
   { Wsta[*Nword] = i;
     while ((str[i]!=splsym)&&(i<=(L-1))) { ++i; }
     Wend[*Nword] = i-1;
     ++(*Nword); }
  ++i;
 }

} /* end of Split_to_Words() */

                                                                                                    
void Remove_Space(rsp,orig)
 char *rsp,*orig;
{
 int i,j;
 j =0;
 for (i=0;i<strlen(orig);++i)
  {
   if (orig[i]!=' ') { rsp[j] = orig[i];  ++j;}
  }
 rsp[j] = '\0';
                                                                                                    
} /* end of Remove_Space() */

int Number_Of_ATOMs(Ahead)
 struct ATOM *Ahead;
{ struct ATOM *an;
  int Natom;
  Natom = 0; an = Ahead;
  while (an->next != NULL){ 
    an = an->next; ++Natom; 
   }
  return(Natom);
} /* end of  Number_Of_ATOMs() */


int Number_Of_RESIDUEs(Rhead)
 struct RESIDUE *Rhead;
{ struct RESIDUE *rn;
  int Nres;

  Nres = 0; rn = Rhead;
  while (rn->next != NULL){
    rn = rn->next; 
    Nres += 1; 
  }

  return(Nres);

} /* end of  Number_Of_RESIDUEs() */




int Number_Of_CHAINs(Ahead)
 struct ATOM *Ahead;
{
 struct ATOM *an;
 int Nchain;

 Nchain = 0;
 an = Ahead;

 while (an->next != NULL){
   an = an->next;
   if (Nchain!=0){ 
     if ((an->Chain != an->prev->Chain)||(an->model_num != an->prev->model_num)){
       Nchain += 1; 
     }
   }
   else{
     Nchain += 1;
   } 
 }

 return(Nchain);

} /* end of  Number_Of_CHAINs() */


int Number_Of_MODELs(Ahead)
 struct ATOM *Ahead;
{
 struct ATOM *an;
 int Nmodel;

 Nmodel = 0;
 an = Ahead;

 while (an->next != NULL){
   an = an->next;
   if (Nmodel!=0){ 
     if (an->model_num != an->prev->model_num){
        Nmodel += 1; 
     }
   }
   else{
     Nmodel += 1;
   } 
 }

 return(Nmodel);

} /* end of  Number_Of_MODELs() */






int Free_ATOMs(Ahead)
 struct ATOM *Ahead;
{
 /* return "Number of free atom" */
 struct ATOM *an,*pn;
 int Nfree;

 /* printf("#Free_ATOMs()\n"); */

 if (Ahead->next == NULL) return(0);
                                                                                                        
 an = Ahead; pn = NULL; Nfree = 0;
 while (an->next != NULL){
   pn = an;
   an = an->next;
   if ((pn != Ahead)&&(pn != NULL)) { free(pn); pn = NULL; ++Nfree;}
  }
                                                                                                        
 if (an!=NULL) { free(an); an = NULL; ++Nfree;}
 Ahead->next = NULL;
 return(Nfree);

} /* end of Free_ATOMs() */





void Set_Region_Of_Atoms(Ahead,chainA,chainB,resA,resB)
 struct ATOM *Ahead;
 char chainA,chainB,*resA,*resB;
{
 struct ATOM *an;
 int NatomA,NatomB;

 NatomA = NatomB = 0;

 an = Ahead;
 while (an->next != NULL){
  an = an->next;
       if ( ((chainA=='-')||(an->Chain==chainA))  &&
            ((strncmp(resA,"xxx",3)==0)||(strncmp(resA,an->Resi,3)==0)))
        { an->region = 'A'; ++NatomA;}
  else if ( ((chainB=='-')||(an->Chain==chainB))  &&
            ((strncmp(resB,"xxx",3)==0)||(strncmp(resB,an->Resi,3)==0)))
 { an->region = 'B'; ++NatomB;}
  else  an->region = ' ';
 }

 if ((NatomA==0)||(NatomB==0)){
  printf("#ERROR : regions are not proper(NatomA %d NatomB %d).\n",NatomA, NatomB);
  exit(1);
 }
 else 
  { 
    printf("#REGION_ASSIGNED NatomA %d NatomB %d\n",NatomA,NatomB);
   }
} /* end of Set_Region_Of_Atoms() */



void Make_Residue_List(Ahead,Rhead)
 struct ATOM    *Ahead;
 struct RESIDUE *Rhead;
{
 struct ATOM *an;
 struct RESIDUE *rn;
 int  Nresidue,init;

 printf("#Make_Residue_List(Ahead,Rhead)\n");
 Rhead->next = NULL;

 an = Ahead;
 rn = Rhead;
 Nresidue = 0; init = 1;
 while (an->next != NULL){
   an = an->next;
 
   if ((init==1)  ||
        ((an->prev!=NULL)&&((strncmp(an->Rnum,an->prev->Rnum,5)!=0)||(an->Chain != an->prev->Chain)))){

/*
    printf("MAKE_NEW_RESIDUE an ATOM %s %s RES %s %s %c\n",an->Anum,an->Atom,an->Rnum,an->Resi,an->Chain);
 */
     init = 0; 
     rn->next = (struct RESIDUE *)malloc(sizeof(struct RESIDUE));
     rn->next->prev = rn;
     rn->next->next = NULL;
     rn = rn->next;
     rn->Chain = an->Chain;
     sprintf(rn->Resi,"%s",an->Resi);
     sprintf(rn->Rnum,"%s",an->Rnum);
     rn->Natom = 0; 
     rn->Natom_contact = 0; 
     ++Nresidue;  
     rn->num = Nresidue;
     rn->rep_atom = an; 
     rn->Ahead.rnext = NULL;  
   }

   if (strncmp(an->Atom," CA ",4)==0) rn->rep_atom = an;

   Add_Atom_To_Residue_Tail(rn,an); 
 }

 printf("#Nresidue %d\n",Nresidue);
} /* end of Make_Residue_List() */




void Add_Atom_To_Residue_Tail(rn,an)
 struct RESIDUE *rn;
 struct ATOM    *an;
{
 struct ATOM *bn;

 /*
 printf("#Add_Atom_To_Residue_Tail(rnum %d anum %d)\n",rn->num,an->num); fflush(stdout);
 */
 
 bn = &(rn->Ahead);
 while (bn->rnext != NULL) bn = bn->rnext;
 bn->rnext = an;
 an->rnext = NULL; 
 ++(rn->Natom);
 an->res = rn;
 an->rnum = rn->num; 
 sprintf(an->Rnum,"%s",rn->Rnum);
 sprintf(an->Resi,"%s",rn->Resi);

} /* end of Add_Atom_To_Residue_Tail() */




void Add_Atom_To_Residue_Head(rn,an)
 struct RESIDUE *rn;
 struct ATOM    *an;
{
 struct ATOM *bn;

 bn = rn->Ahead.rnext;
 rn->Ahead.rnext = an; 
 if (bn != NULL) { an->rnext = bn; } 
 ++(rn->Natom);
 an->res = rn;
 an->rnum = rn->num; 
 sprintf(an->Rnum,"%s",rn->Rnum);
 sprintf(an->Resi,"%s",rn->Resi);

} /* end of Add_Atom_To_Residue_Head() */






void Set_Constant_tFactor(Phead,val)
 struct ATOM *Phead;
 float val;
{  
 struct ATOM *an;
   
 an = Phead;
 while (an->next != NULL){
   an = an->next;
   an->tFactor = val;
 } 

} /* end of Set_Constant_tFactor() */



void Renumber_Atom_Number(Phead)
 struct ATOM *Phead;
{  
 struct ATOM *an;
 int Natom;
 char buff[128];  
 an = Phead;
 Natom = 0;
 while (an->next != NULL){
   an = an->next;
   ++Natom;
   an->num = Natom;
   /* sprintf(an->Anum,"%5d",an->num); */
   sprintf(buff,"%5d",an->num);
   Get_Part_Of_Line(an->Anum, buff, 0,6);
   if (an->res != NULL) { 
      sprintf(an->Rnum,"%s",an->res->Rnum);
      an->rnum = an->res->num;}
   else{ 
     /*sprintf(an->Rnum,"%4d ",an->num%10000); */ 
     sprintf(buff,"%4d ",an->num%10000);  
     Get_Part_Of_Line(an->Rnum, buff, 0,6);
     an->rnum = an->num;
   }
 } 

} /* end of Renumber_Atom_Number() */



void Renumber_Residue_Number(Rhead)
 struct RESIDUE *Rhead;
{  
 struct RESIDUE *rn;
 int Nres;
   
 rn = Rhead;
 Nres = 0;
 while (rn->next != NULL){
   rn = rn->next;
   ++Nres;
   rn->num = Nres;
   sprintf(rn->Rnum,"%4d ",rn->num%10000);
 } 

} /* end of Renumber_Residue_Number() */






void Delete_Doubling_altLoc_Atoms(HeadAtom)
 struct ATOM *HeadAtom;
/*
  NOTE: This functions does not work for molecules generated by mmCIF assembly !!
*/
{
 struct ATOM     *an;
 char hit;
 int Nex;

 Nex = 0;
 an = HeadAtom;
 while (an->next != NULL){
   an = an->next;
   if (an->altLoc != ' '){
     hit = Find_Same_Res_Rnum_Atom_Before(HeadAtom,an);
     if (hit==1){ 
       ++Nex;
       if (an->prev != NULL){ an->prev->next   = an->next;}
       if (an->next != NULL){ an->next->prev   = an->prev;}
     }
   }

 } 

 Renumber_Atom_num(HeadAtom);

} /* end of Delete_Doubling_altLoc_Atoms() */




void Renumber_Atom_num(HeadAtom)
 struct ATOM *HeadAtom;
{
 int Natom;
 struct ATOM     *an;

 Natom = 0;
 an = HeadAtom;
 while (an->next != NULL)
 { an = an->next; 
   an->num = Natom;
   ++Natom; }

} /* end of Renumber_Atom_num() */




int Find_Same_Res_Rnum_Atom_Before(HeadAtom,focus)
 struct ATOM  *HeadAtom;
 struct ATOM  *focus;
{
 struct ATOM *an;
 char hit;

 hit = 0;
 an = HeadAtom;
 while ((an->next != NULL)&&(hit==0)&&(an->next->num!=focus->num)){
  an = an->next;
  if ((an->altLoc != ' ')&&(an->num != focus->num)){
    if ( (focus->Chain == an->Chain) &&
         (strncmp(focus->Resi, an->Resi,3)==0) &&
         (strncmp(focus->Rnum, an->Rnum,5)==0) &&
         (strncmp(focus->Atom, an->Atom,4)==0)) { hit = 1; return(hit);}
   }

 } /* an */

  return(hit);
} /* end of Find_Same_Res_Rnum_Atom_Before() */







void Change_N_CA_C_O_Hetatom_to_Atom(HeadAtom)
 struct ATOM *HeadAtom;
{
 struct ATOM *atom,*satom;
 char Nok,Cok,CAok,Ook,init;
 int Nhetatm;

 printf("#Change_N_CA_C_O_Hetatom_to_Atom()\n");

 atom = HeadAtom;
 Nok = Cok = CAok = Ook = 0; 
 Nhetatm = 0; 
 satom = NULL; 
 init = 1;

 while (atom->next != NULL){
  
   if ((init==0)&&((atom->Chain != atom->next->Chain)||(strcmp(atom->Rnum,atom->next->Rnum)!=0))){
     /* printf("#Chain %c Resi %s Rnum %s Nhetatm %d N C CA O %d %d %d %d \n",
        atom->Chain,atom->Resi,atom->Rnum,Nhetatm,Nok,Cok,CAok,Ook); */
     if (((Nok*Cok*CAok*Ook)==1)&&(Nhetatm>0)){
       /* printf("#Change_to_ATOM\n"); */
       while ((satom!=NULL)&&(satom->num <= atom->num)){
         if (satom->AHtype != 'A') { satom->AHtype = 'A'; }
          satom = satom->next; 
       }
     }
      satom = NULL;
      Nok = Cok = CAok = Ook = 0; 
      Nhetatm = 0;
   }

   atom = atom->next;
 /*
   printf("#AHtype %c Chain %c Resi %s Rnum %s\n",atom->AHtype,atom->Chain,atom->Resi,atom->Rnum);
  */
   if (satom==NULL){ satom = atom;}
   if (atom->AHtype != 'A'){ ++Nhetatm;}
   if (strncmp(atom->Atom," CA ",4)==0){ CAok = 1;}
   if (strncmp(atom->Atom," N  ",4)==0){ Nok  = 1;}
   if (strncmp(atom->Atom," C  ",4)==0){ Cok  = 1;}
   if (strncmp(atom->Atom," O  ",4)==0){ Ook  = 1;}
   init = 0;
 }

   if (((Nok*Cok*CAok*Ook)==1)&&(Nhetatm>0)){
     /* printf("#Change_to_ATOM\n"); */
     while ((satom!=NULL)&&(satom->num <= atom->num)){
       if (satom->AHtype != 'A') { satom->AHtype = 'A'; }
        satom = satom->next; 
     }
   }


} /* end of  Change_N_CA_C_O_Hetatom_to_Atom() */



void Delete_HETATM_ATOMs(HeadAtom)
 struct ATOM *HeadAtom;
{
  struct ATOM *atm;
  int Nhetatm,Ndelatm;

  Nhetatm = Ndelatm = 0;
  atm = HeadAtom;
  while (atm->next != NULL){
    atm = atm->next;
    if (atm->AHtype == 'H'){
      Nhetatm += 1;
    }
  }

  atm = HeadAtom;
  while (atm->next != NULL){
    atm = atm->next;
    if (atm->AHtype == 'H'){
      if (atm->prev != NULL){
        atm->prev->next = atm->next;
      }
      if (atm->next != NULL){
        atm->next->prev = atm->prev;
      }
      Ndelatm += 1;
    }
  }

 /*
  printf("#Nhetatm %d\n",Nhetatm);
  printf("#Nhetatm %d Ndelatm %d\n",Nhetatm, Ndelatm);
 */


} /* end of Delete_HETATM_ATOMs() */



void Delete_HETATMs_or_ATOMs(HeadAtom,DeleteAtmHet)
 struct ATOM *HeadAtom;
 char   DeleteAtmHet;  /* 'A':delete ATOM, 'H':delete HETATM */
{
  struct ATOM *atm;
  int Ndelatm_suppose,Ndelatm;

  Ndelatm = Ndelatm_suppose = 0;
  atm = HeadAtom;
  while (atm->next != NULL){
    atm = atm->next;
    if (atm->AHtype == DeleteAtmHet){
      Ndelatm_suppose += 1;
    }
  }

  atm = HeadAtom;
  while (atm->next != NULL){
    atm = atm->next;
    if (atm->AHtype == DeleteAtmHet){
      if (atm->prev != NULL){
        atm->prev->next = atm->next;
      }
      if (atm->next != NULL){
        atm->next->prev = atm->prev;
      }
      Ndelatm += 1;
    }
  }

 printf("#Ndelatm_suppose %d Ndelatm %d\n",Ndelatm_suppose, Ndelatm);

} /* end of Delete_HETATMs_or_ATOMs() */







void Add_String_to_STRING_NODE(Shead,string)
  struct STRING_NODE *Shead; 
  char *string;
{
 struct STRING_NODE *sn;
 /* printf("#Add_String_to_STRING_NODE('%s')\n",string); */

 sn = Shead;
 while (sn->next != NULL){sn = sn->next;}
 sn->next = (struct STRING_NODE*)malloc(sizeof(struct STRING_NODE));
 sn->next->prev = sn;
 sn->next->next = NULL; 
 Get_Part_Of_Line(sn->next->line,string,0,MAX_STRING_NODE_LINE-1);
} /* end of Add_String_to_STRING_NODE() */


void Show_STRING_NODE(Shead)
  struct STRING_NODE *Shead;
{
 struct STRING_NODE *sn;
 
 sn = Shead;
 while (sn->next != NULL){
   sn = sn->next;
   printf("%s\n",sn->line);
 }
} /* end of Show_STRING_NODE() */




void Write_One_Point_In_PDB_Format(fname,Pos,mode,comment)
 char *fname;
 float Pos[3];
 char  mode; /* 'a'ppend, 'w'rite */
 char *comment;
{
 FILE *fp;
 char core[MAX_FILE_NAME],buff[MAX_FILE_NAME],buff2[MAX_FILE_NAME];

 printf("#Write_One_Point_In_PDB_Format -> \"%s\"\n",fname);
 if (fname[0]=='-') fp = stdout;
 else
 { if (mode=='a') fp = fopen(fname,"a");
             else  fp = fopen(fname,"w");
  if (fp==NULL) { printf("#ERROR:Can't write to \"%s\"\n",fname); exit(1);}
 }

 Find_Filename_Core(core,fname);
 Get_Part_Of_Line(buff,fname,0,40);
 Get_Part_Of_Line(buff2,core,0,9);
 fprintf(fp,"HEADER    %-40s%s   %-10s0XXX   1\n",buff,Get_Date_String_PDB(),buff2);
 fprintf(fp,"REMARK    DATE %s\n",Get_Date_String());
 if (PAR.COMMAND[0] != '\0') fprintf(fp,"REMARK  COMMAND %s\n",PAR.COMMAND);
 if (comment[0]!='\0') fprintf(fp,"REMARK  COMMENT %s\n",comment);
 
 fprintf(fp,"HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
       1," C  ","CEN",1, Pos[0],Pos[1],Pos[2],0.0,0.0);
 fprintf(fp,"TER\n");

 if (fp!=stdout) fclose(fp);

} /* end of Write_One_Point_In_PDB_Format() */






int Find_XYZ_By_Pattern(Pos,Ahead,pattern)
 float  Pos[3];
 struct ATOM *Ahead;
 char *pattern;  /* '[ATOM]:[RES]:[CHAIN]:[RNUM]' */
{
 struct ATOM *an;
 int Nword, Ws[10],We[10];
 char Pres[10],Prnum[10],Pchain[10],Patm[10]; /* Pattern str for various properties */
 char buff[10],rnumR[10];
                                                                                                    
 Split_to_Words(pattern,':',&Nword,Ws,We,10);
 if (Nword>0) Get_Part_Of_Line(Patm,  pattern,Ws[0],We[0]); else sprintf(Patm,"xxx");
 if (Nword>1) Get_Part_Of_Line(Pres,  pattern,Ws[1],We[1]); else sprintf(Pres,"xxx");
 if (Nword>2) Get_Part_Of_Line(Pchain,pattern,Ws[2],We[2]); else sprintf(Pchain,"x");
 if (Nword>3) {Get_Part_Of_Line(buff, pattern,Ws[3],We[3]); Remove_Space(Prnum,buff); }
         else sprintf(Prnum,"xxx");

 printf("#Find_Atom_By_Pattern(\"%s\" -> Patm '%s' Pres '%s' Pchain '%s' Prnum '%s')\n",
  pattern,Patm, Pres, Pchain,Prnum);
 
 an = Ahead;
 while (an->next != NULL){
  an = an->next;
                                                                                                    
   Remove_Space(rnumR,an->Rnum);
   if ( 
        ((Patm[0] == 'x')||(strncmp(Patm,an->Atom,4)==0))&&
        ((Pres[0]=='x')||(strncmp(Pres,an->Resi,3)==0)) &&
        ((Pchain[0]=='-')||(Pchain[0]=='x')||(Pchain[0]==an->Chain)) &&
        ((Prnum[0]=='x')||(strncmp(Prnum,rnumR,3)==0))
        )
   { 
    Pos[0] = an->Pos[0]; Pos[1] = an->Pos[1]; Pos[2] = an->Pos[2];
    printf("%s %s %c %s %f %f %f\n",
    an->Atom, an->Resi, an->Chain, an->Rnum,an->Pos[0], an->Pos[1], an->Pos[2]);

    return(1);    
    }
 } /* an */

 printf("#ERROR:Can't find atom for the pattern '%s'\n",pattern);
 exit(1); 
 return(0);

} /* end of Find_Atom_By_Pattern() */







int Find_XYZ_By_Pattern_In_PDB_File(Pos,ifname,pattern)
 float  Pos[3];
 char *ifname; 
 char *pattern;  /* '[ATOM]:[RES]:[CHAIN]:[RNUM]' */
{
 FILE *fp;
 int Nword, Ws[10],We[10];
 char Pres[10],Prnum[10],Pchain[10],Patm[10]; /* Pattern str for various properties */
 char line[255],buff[16],rnumR[8];
 char Anum[8],Atom[8],Resi[8],Rnum[8],Chain;

 Split_to_Words(pattern,':',&Nword,Ws,We,10);
 if (Nword>0) Get_Part_Of_Line(Patm,  pattern,Ws[0],We[0]); else sprintf(Patm,"xxx");
 if (Nword>1) Get_Part_Of_Line(Pres,  pattern,Ws[1],We[1]); else sprintf(Pres,"xxx");
 if (Nword>2) Get_Part_Of_Line(Pchain,pattern,Ws[2],We[2]); else sprintf(Pchain,"x");
 if (Nword>3) {Get_Part_Of_Line(buff, pattern,Ws[3],We[3]); Remove_Space(Prnum,buff); }
         else sprintf(Prnum,"xxx");

 printf("#Find_Atom_By_Pattern_In_PDB_File('%s' \"%s\" -> Patm '%s' Pres '%s' Pchain '%s' Prnum '%s')\n",
  ifname,pattern,Patm, Pres, Pchain,Prnum);
 
 fp = fopen(ifname,"r");
 if (fp==NULL){printf("#ERROR:Can't open pdbfile \"%s\"\n",ifname); exit(1);}
 
 while (feof(fp)==0){
  line[0] = '\0';
  fgets(line,255,fp);
  if ((strncmp(line,"ATOM",4)==0)||(strncmp(line,"HETATM",6)==0)){ 
   Get_Part_Of_Line(Anum,line,6,10);
   Get_Part_Of_Line(Atom,line,12,15);
   Get_Part_Of_Line(Resi,line,17,19);
   Get_Part_Of_Line(Rnum,line,22,26);
   Remove_Space(rnumR,Rnum);
   Chain = line[21];
   if ( 
        ((Patm[0] == 'x')||(strncmp(Patm,Atom,4)==0))&&
        ((Pres[0]=='x')  ||(strncmp(Pres,Resi,3)==0)) &&
        ((Pchain[0]=='-')||(Pchain[0]=='x')||(Pchain[0]==Chain)) &&
        ((Prnum[0]=='x') ||(strncmp(Prnum,rnumR,3)==0))
        ){ 
    Get_Part_Of_Line(buff,line,30,37); Pos[0] = atof(buff);
    Get_Part_Of_Line(buff,line,38,45); Pos[1] = atof(buff);
    Get_Part_Of_Line(buff,line,46,53); Pos[2] = atof(buff);
    printf("#%s",line);
    return(1);    
    }
  } 
 
 } /* while */

 printf("#ERROR:Can't find atom for the pattern '%s'\n",pattern);
 exit(1); 
 return(0);

} /* end of Find_Atom_By_Pattern_In_PDB_File() */
 


void Make_Copied_Atoms(CopyHead, OrigHead)
 struct ATOM *CopyHead, *OrigHead;
{
 struct ATOM *cn, *on;
 int k;

 cn = CopyHead;
 on = OrigHead;
 while (on->next!=NULL){
  on = on->next;
  cn->next = (struct ATOM *)malloc(sizeof(struct ATOM));
  cn->next->prev = cn;
  for (k=0;k<3;++k) cn->next->Pos[k] = on->Pos[k];
  cn->next->R  = on->R;
  cn->next->RR = on->RR;
  cn->next->num = on->num;
  cn->next->rnum = on->rnum;
  sprintf(cn->next->Anum,"%s",on->Anum);
  sprintf(cn->next->Atom,"%s",on->Atom);
  sprintf(cn->next->Resi,"%s",on->Resi);
  sprintf(cn->next->Rnum,"%s",on->Rnum);
  cn->next->Chain = on->Chain;
  cn->next->altLoc = on->altLoc;
  cn->next->Occup = on->Occup;
  cn->next->tFactor = on->tFactor;
  cn->next->AHtype = on->AHtype;
  cn->next->region = on->region;
  cn->next->mark = on->mark;
  cn->next->next = NULL;
  cn = cn->next;
 }

} /* end of Make_Copied_Atoms() */


double Set_Molecule_Gcenter_to_Origin(AtomHead)
 struct ATOM *AtomHead;
 /* return maximum distance from the origin */
{
 struct ATOM *an;
 float Gpos[3],DD,D,maxDori;
 int i,Natom;
 Natom = 0;
 for (i=0;i<3;++i) Gpos[i] = 0.0;
 /* (1) Cal Gpos */
 an = AtomHead;
 while (an->next != NULL){
  an = an->next;
  for (i=0;i<3;++i) Gpos[i] += an->Pos[i];
  ++Natom;
 }
 for (i=0;i<3;++i) Gpos[i] /= Natom;
 printf("#Set_Molecule_Gcenter_to_Origin() --> Gpos %f %f %f\n",Gpos[0],Gpos[1],Gpos[2]);
 /* (2) Substract Gpos and Decide maximum Distace from origin (maxDori)*/
 maxDori = 0.0;
 an = AtomHead;
 while (an->next != NULL){
  an = an->next;
  DD = 0.0;
  for (i=0;i<3;++i) {
   an->Pos[i] -= Gpos[i];
   DD += an->Pos[i] * an->Pos[i];
  }
                                                                                                                               
  if (DD>0.0) D = sqrt(DD); else D = 0.0;
  D += an->R;
  if (D>maxDori) maxDori = D;
 }
 return(maxDori);
}  /* end of Set_Molecule_Gcenter_to_Origin() */



void Rotate_Molecule(AtomHead,Rmat)
 struct ATOM *AtomHead;
 float Rmat[3][3];
{
 struct ATOM *an;
 float Rpos[3];
 int i,j;
 an = AtomHead;
 while (an->next != NULL){
  an = an->next;
  for (i=0;i<3;++i){
   Rpos[i] = 0.0;
   for (j=0;j<3;++j) Rpos[i] += Rmat[i][j] * an->Pos[j];
  }
  for (i=0;i<3;++i) an->Pos[i] = Rpos[i];
 }
}  /* end of Rotate_Molecule() */
     
void InvRotate_Molecule(AtomHead,Rmat)
 struct ATOM *AtomHead;
 float Rmat[3][3];
{
 struct ATOM *an;
 float Rpos[3];
 int i,j;
 an = AtomHead;
 while (an->next != NULL){
  an = an->next;
  for (i=0;i<3;++i){
   Rpos[i] = 0.0;
   for (j=0;j<3;++j) Rpos[i] += Rmat[j][i] * an->Pos[j];
  }
  for (i=0;i<3;++i) an->Pos[i] = Rpos[i];
 }
}  /* end of InvRotate_Molecule() */
                                                                                                                   

                                                                                                                          
void Translate_Molecule(AtomHead,Tvec)
 struct ATOM *AtomHead;
 float Tvec[3];
{
 struct ATOM *an;
 int i;
 an = AtomHead;
 while (an->next != NULL)
 {
  an = an->next;
  for (i=0;i<3;++i) an->Pos[i] += Tvec[i];
 }
}  /* end of Translate_Molecule() */


int ChainIDStr_Match(ChainIDStr,chain)
   char *ChainIDStr;
   char chain;
{
 int i,L;
 L = strlen(ChainIDStr);
 for (i=0;i<L;++i){
   if ((ChainIDStr[i]=='-')||(ChainIDStr[i]==chain)) return(1);
 }
 return(0);
} /* end of ChainIDStr_Match() */

/*
>> This function is only for testing <<

struct ATOM *Find_Atom_with_Anum(anum,Ahead)
 char   *anum;
 struct ATOM *Ahead;
{
 struct ATOM *an;
 an = Ahead;
 while (an->next != NULL){
  an = an->next;
  if (strncmp(an->Anum,anum,6)==0) return(an);
  }

 printf("#ERROR:Can't find anum \"%s\"\n",anum); exit(1);

}
*/




void Element_from_PDB_Atomname(element,atomname)
  char *element;
  char *atomname; 
{
/*    
>>> EXAMPLES OF PDB FILES <<<
>>pdb1aa9.ent<<
HETATM 2727 MG    MG A 173     129.369  -1.143  10.024  1.00  3.45          MG
HETATM 2728  PB  GDP A 180     132.580  -0.884  10.946  1.00  1.28           P
HETATM 2729  O1B GDP A 180     131.602  -1.788  11.591  1.00  2.05           O
HETATM 2730  O2B GDP A 180     131.979   0.121  10.041  1.00  1.89           O
HETATM 2734  O1A GDP A 180     132.281  -0.115  13.953  1.00  2.00           O
HETATM 2736  O5' GDP A 180     134.673   0.349  13.942  1.00  1.40           O
HETATM 2757 H5'' GDP A 180     134.197  -0.301  15.838  1.00  1.76           H
HETATM 2758  H4' GDP A 180     134.640   1.663  17.022  1.00  1.33           H
HETATM 2760 HO3' GDP A 180     132.224   1.890  16.425  1.00  1.95           H
HETATM 2761  H2' GDP A 180     134.224   3.920  13.786  1.00  1.22           H
HETATM 2762 HO2' GDP A 180     134.163   5.712  15.429  1.00  1.71           H
HETATM 2763  H1' GDP A 180     135.817   4.485  16.154  1.00  1.00           H
HETATM 2765  HN1 GDP A 180     140.538   7.338  14.020  1.00  1.32           H
HETATM 2766 HN21 GDP A 180     140.460   7.807  16.187  1.00  1.44           H
HETATM 2767 HN22 GDP A 180     139.262   7.113  17.256  1.00  1.31           H
>> pdb4aai.ent <<
ATOM      9  H   MET A   1       9.517  18.716 -26.058  1.00  0.00           H
ATOM     10  HA  MET A   1       9.836  17.425 -24.009  1.00  0.00           H
ATOM     57  HA  SER A   4       1.442  13.939 -19.192  1.00  0.00           H
ATOM     58  HB2 SER A   4       0.861  14.036 -22.152  1.00  0.00           H
ATOM     59  HB3 SER A   4       0.707  15.457 -21.119  1.00  0.00           H
ATOM     60  HG  SER A   4      -0.721  13.020 -20.965  1.00  0.00           H
ATOM   2266  HA  PRO B  60      -6.046  15.537  10.707  1.00  0.00           H
ATOM   2267  HB2 PRO B  60      -6.055  17.031  12.445  1.00  0.00           H
ATOM   2268  HB3 PRO B  60      -4.429  17.409  11.865  1.00  0.00           H
ATOM   2269  HG2 PRO B  60      -4.832  16.612  14.376  1.00  0.00           H
ATOM   2270  HG3 PRO B  60      -3.409  16.072  13.464  1.00  0.00           H
ATOM   2271  HD2 PRO B  60      -4.254  14.006  14.051  1.00  0.00           H
ATOM   2272  HD1 PRO B  60      -5.925  14.594  14.109  1.00  0.00           H
ATOM   2301  CG2 ILE B  62      -5.021  10.552   6.917  1.00  0.00           C
ATOM   2302  CD1 ILE B  62      -2.865  10.208  10.115  1.00  0.00           C
ATOM   2303  HA  ILE B  62      -2.225  11.381   7.606  1.00  0.00           H
ATOM   2304  HB  ILE B  62      -4.964  11.559   8.787  1.00  0.00           H
ATOM   2305 HG12 ILE B  62      -3.147   9.202   8.279  1.00  0.00           H
ATOM   2306 HG13 ILE B  62      -4.567   9.244   9.319  1.00  0.00           H
ATOM   2307 HG21 ILE B  62      -4.334  10.482   6.087  1.00  0.00           H
ATOM   2308 HG22 ILE B  62      -5.803  11.258   6.680  1.00  0.00           H
ATOM   2309 HG23 ILE B  62      -5.458   9.581   7.104  1.00  0.00           H
ATOM   2310 HD11 ILE B  62      -3.462  10.736  10.845  1.00  0.00           H
ATOM   2311 HD12 ILE B  62      -2.064  10.850   9.772  1.00  0.00           H
>> pdb3arc.ent <<
HETATM41688  O1  OEX A 601     -26.510 -36.601 204.050  1.00 22.48           O
HETATM41689 CA1  OEX A 601     -27.969 -36.452 202.243  1.00 24.57          CA
HETATM41690 MN1  OEX A 601     -24.977 -35.548 203.849  1.00 22.62          MN
HETATM41691  O2  OEX A 601     -28.560 -34.524 203.696  1.00 25.22           O
HETATM41692 MN2  OEX A 601     -27.388 -35.240 205.326  1.00 22.88          MN
HETATM41693  O3  OEX A 601     -25.808 -34.096 204.540  1.00 23.27           O
HETATM41694 MN3  OEX A 601     -27.257 -33.264 203.226  1.00 23.49          MN
HETATM41695  O4  OEX A 601     -28.740 -32.584 201.922  1.00 26.62           O
HETATM41696 MN4  OEX A 601     -27.599 -33.237 200.276  1.00 26.26          MN
>> pdb1aa5.ent <<
HETATM  133  HN1 3FG A   7       1.904   2.986   5.854  1.00  5.72           H
HETATM  134  HA  3FG A   7       0.691   0.779   4.801  1.00  6.88           H
HETATM  135  HD1 3FG A   7       4.075   1.650   1.428  1.00  5.99           H
HETATM  136  HZ  3FG A   7       6.148   0.117   3.505  1.00  5.61           H
HETATM  137  HD2 3FG A   7       6.605  -0.752   5.609  1.00  7.09           H
HETATM  138  HG2 3FG A   7       3.497  -0.129   6.511  1.00  6.76           H
HETATM  139  HOT 3FG A   7       1.007   0.445   8.470  1.00  9.61           H
>> pdb3a0w.ent <<
HETATM 2502 HG  AEMC A 801      49.554  58.731  15.205  0.30 29.81          HG
HETATM 2503 HG  BEMC A 801      49.086  63.531  15.237  0.30 44.78          HG
HETATM 2504  C1 AEMC A 801      48.824  56.764  16.257  0.30 27.68           C
HETATM 2505  C1 BEMC A 801      49.487  64.170  17.491  0.30 43.48           C
HETATM 2506  C2 AEMC A 801      49.010  56.925  17.750  0.30 28.92           C
HETATM 2507  C2 BEMC A 801      51.053  63.827  17.815  0.30 43.57           C
HETATM 2508 HG   EMC A 802      50.971  48.709  -4.694  0.80 25.19          HG
>> pdb1cc8.ent <<
HETATM  569 HG    HG A  74       5.897   0.005   3.783  1.00 10.60          HG
HETATM  570  C1  BEN A 186       6.778  10.510  20.665  1.00  4.96           C
HETATM  571  C2  BEN A 186       5.994   9.710  19.842  1.00  5.48           C
>> pdb1o7t.ent <<
HETATM21403  OBD HF5 A1310      12.209  52.345 -17.389  1.00 55.42           O
HETATM21404 HFA  HF5 A1310      10.810  51.834 -13.339  1.00 33.27          HF
HETATM21405  O00 HF5 A1310      12.522  51.019 -14.388  1.00 66.11           O
HETATM21406  OAB HF5 A1310      10.204  50.460 -15.311  1.00 42.43           O
HETATM21407  OAC HF5 A1310      12.821  52.802 -12.248  1.00 34.32           O
HETATM21408  OAD HF5 A1310      10.622  53.390 -14.955  1.00 86.30           O
HETATM21409  OAE HF5 A1310      10.540  53.809 -11.877  1.00 56.64           O
HETATM21410  OA1 HF5 A1310       8.751  50.422 -13.681  1.00 37.07           O
HETATM21411 HFB  HF5 A1310      12.391  50.336 -16.504  1.00 43.82          HF
HETATM21412  OBC HF5 A1310      14.747  50.922 -16.256  1.00 64.16           O
HETATM21413  OB1 HF5 A1310      12.516  47.750 -16.386  1.00 39.29           O
HETATM21414  OB2 HF5 A1310       9.974  50.462 -17.864  1.00 25.83           O
HETATM21415  OB3 HF5 A1310      14.284  50.097 -18.387  1.00 46.74           O
HETATM21416 HFC  HF5 A1310      14.480  51.988 -14.021  1.00 44.44          HF
HETATM21417  OCD HF5 A1310      14.371  54.160 -14.611  1.00 53.36           O
HETATM21418  OC1 HF5 A1310      16.837  52.625 -14.715  1.00 63.62           O
HETATM21419  OC2 HF5 A1310      15.120  54.166 -11.785  1.00 40.71           O
HETATM21420 HFD  HF5 A1310      12.141  54.214 -15.784  1.00 65.07          HF
HETATM21421  ODE HF5 A1310      11.786  56.117 -14.349  1.00 99.50           O
HETATM21422  OD1 HF5 A1310      13.605  55.217 -17.842  1.00 89.38           O
HETATM21423  OD2 HF5 A1310      11.175  56.187 -17.319  1.00 68.46           O
HETATM21424 HFE  HF5 A1310      10.581  56.275 -12.625  1.00 99.50          HF
HETATM21425  OE1 HF5 A1310      12.578  56.675 -11.378  1.00 85.88           O
HETATM21426  OE2 HF5 A1310       8.798  55.140 -13.958  1.00 99.50           O
HETATM21427  OE3 HF5 A1310       9.748  58.496 -13.264  1.00 84.21           O
>>pdb1tmn.ent<<
HETATM 2454  N1  0ZN E 317      52.952  15.765  -3.983  1.00 18.44           N
HETATM 2455  CA2 0ZN E 317      53.530  14.448  -4.298  1.00 23.01           C
HETATM 2456  C1  0ZN E 317      52.643  13.531  -5.133  1.00 26.89           C
HETATM 2457  O1  0ZN E 317      51.572  14.027  -5.605  1.00 22.68           O
HETATM 2458  CB3 0ZN E 317      54.101  13.735  -3.037  1.00 23.46           C
HETATM 2459  CG2 0ZN E 317      52.943  13.258  -2.164  1.00 19.83           C
HETATM 2460 CD11 0ZN E 317      51.832  13.908  -1.823  1.00 24.81           C
HETATM 2461 CD21 0ZN E 317      52.798  11.977  -1.642  1.00 24.95           C
HETATM 2462  NE1 0ZN E 317      50.892  12.989  -1.314  1.00 25.46           N
HETATM 2463 CE21 0ZN E 317      51.442  11.840  -1.286  1.00 25.54           C
HETATM 2464  CE3 0ZN E 317      53.723  10.946  -1.480  1.00 30.14           C
HETATM 2465 CZ21 0ZN E 317      50.904  10.638  -0.855  1.00 28.63           C
HETATM 2466  CZ3 0ZN E 317      53.213   9.741  -0.964  1.00 30.85           C
HETATM 2467  CH2 0ZN E 317      51.831   9.584  -0.690  1.00 30.40           C
HETATM 2468  OXT 0ZN E 317      53.142  12.469  -5.652  1.00 24.16           O
HETATM 2469 CA    CA E 318      61.303  22.116   4.786  1.00 13.41          CA
HETATM 2470 CA    CA E 319      61.569  23.253   8.472  1.00 14.54          CA
HETATM 2471 CA    CA E 320      35.745  32.533  15.823  1.00 12.78          CA
HETATM 2472 CA    CA E 321      57.297  12.904  10.398  1.00 16.76          CA
HETATM 2473 ZN    ZN E 322      54.519  20.316  -7.288  1.00 14.25          ZN
HETATM 2474  O   HOH E 331      68.539  16.913 -21.003  1.00 18.54           O
>>pdb1cll.ent<<
HETATM 1135 CA    CA A 149       1.708  47.776  24.197  1.00 12.07          CA
HETATM 1136 CA    CA A 150      -2.112  44.438  13.360  1.00 13.42          CA
HETATM 1137 CA    CA A 151      30.463  11.962  10.473  1.00 25.88          CA
HETATM 1138 CA    CA A 152      29.181  17.963   0.804  1.00 29.77          CA
HETATM 1139  C1  EOH A 153       7.868  37.805  18.104  1.00 63.73           C
HETATM 1140  C2  EOH A 153       6.465  38.277  17.751  1.00 64.50           C
HETATM 1141  H21 EOH A 153       6.499  39.538  17.106  1.00 64.47           H
>>pdb4fik.ent<<
ATOM   1677  CB  CYS A 206      -4.229   8.811  22.820  1.00  7.14           C
ATOM   1678  SG  CYS A 206      -4.184   7.846  21.278  1.00  7.76           S
HETATM 2175  S1 AD2S A 304      -8.717   9.120   6.669  0.71 13.76           S
HETATM 2176  C2 AD2S A 304     -10.591  10.818   7.506  0.71 12.18           C
HETATM 2177  S2 AD2S A 304      -9.076   7.758   8.098  0.71 14.67           S
HETATM 2178  C3 AD2S A 304      -9.260  10.594   7.339  0.71 13.01           C
HETATM 2179  O3 AD2S A 304     -10.490  14.182   8.981  0.71 15.27           O

>> pdb4e7r.ent <<
HETATM 4678 CL2  0NW H 302      23.679   6.056  26.437  1.00 29.41          CL
HETATM 4679  C36 0NW H 302      23.660   3.431  25.705  1.00 30.99           C
HETATM 4680 CL1  0NW H 302      23.546   0.734  25.512  1.00 48.22          CL

>> pdb1pxi.ent <<
HETATM 2367 CL7AACK1 A 500      14.171  44.275  29.729  0.73 13.63          CL
HETATM 2370  S2AACK1 A 500      14.041  42.285  27.535  0.73 22.98           S
HETATM 2372 CL6AACK1 A 500      13.557  41.890  24.844  0.73 13.46          CL
HETATM 2373  C5AACK1 A 500      13.388  44.250  25.907  0.73 16.88           C

>> pdb2gde.ent <<
HETATM 2360  C30 SN3 H 401      18.158 -16.049  26.995  1.00 44.63           C
HETATM 2361 CL31 SN3 H 401      19.334 -14.269  24.819  1.00 47.30          CL
HETATM 2362  N32 SN3 H 401      16.043 -15.948  23.241  1.00 50.01           N

>> pdb4ht0.ent <<
HETATM 2155  F13 V50 A 304      -2.858   1.546  16.530  1.00 33.38           F
HETATM 2163  F11 V50 A 304      -6.992   3.162  15.322  1.00 35.74           F
HETATM 2164  C1  V50 A 304      -5.067   4.374  14.600  1.00 31.51           C

>> pdb3rz7.ent <<
HETATM 2098 F2'1 RZ7 A 262      20.534   7.261  14.299  1.00 25.62           F
HETATM 2099 F3'1 RZ7 A 262      18.846   9.517  12.136  1.00 27.91           F
HETATM 2100 F4'1 RZ7 A 262      20.745   9.802  15.258  1.00 29.96           F
HETATM 2101 F5'1 RZ7 A 262      21.189  11.328  12.368  1.00 32.43           F
HETATM 2102 F6'1 RZ7 A 262      22.957  11.252  15.526  1.00 43.46           F
HETATM 2103 F2'2 RZ7 A 262      18.702   8.275  14.563  1.00 21.80           F

>> pdb1mbd.ent <<
HETATM 1237 FE   HEM A 155      15.271  27.962   0.622  1.00  7.86          FE
HETATM 1238  CHA HEM A 155      16.869  31.030   0.421  1.00  4.49           C
HETATM 1242  NA  HEM A 155      16.206  29.133   1.939  1.00  7.65           N


*/

/*
>>> SUMMARY <<<
ATOM      9  H   MET A   1      -> 'H'
ATOM     10  HA  MET A   1      -> 'H'
HETATM 2678  HN5 X37 A1480      -> 'H'
HETATM  139  HOT 3FG A   7      -> 'H'
HETATM 2679 HN5A X37 A1480      -> 'H'
ATOM   2305 HG12 ILE B  62      -> 'H'
ATOM     60  HG  SER A   4      -> 'H'
HETATM 2766 HN21 GDP A 180      -> 'H'
HETATM 2767 HN22 GDP A 180      -> 'H'
HETATM 2373 HG    HG A 601      -> 'HG'
HETATM21404 HFA  HF5 A1310      -> 'HF'
HETATM21424 HFE  HF5 A1310      -> 'HF'
HETATM 2374  C1  BNZ A 602      -> 'C'
HETATM 2375  C2  BNZ A 602      -> 'C'
HETATM 5274  C4  GCP A1749      -> 'C'
HETATM 5280  C3B GCP D1749      -> 'C'
HETATM 2460 CD11 0ZN E 317      -> 'C'
HETATM 2463 CE21 0ZN E 317      -> 'C'
HETATM 2465 CZ21 0ZN E 317      -> 'C'
HETATM 4680 CL1  0NW H 302      -> 'CL'
HETATM 3647 CL    CL A 800      -> 'CL'
HETATM 2372 CL6AACK1 A 500      -> 'CL'
HETATM 2367 CL7AACK1 A 500      -> 'CL'
HETATM 2361 CL31 SN3 H 401      -> 'CL'
HETATM 2228  N1  HMD A 400      -> 'N'
HETATM 2677  N5  X37 A1480      -> 'N'
HETATM 2226  O1  HMD A 400      -> 'O'
HETATM 3646  O5  PG4 A 500      -> 'O'
HETATM 5277  O1G GCP D1749      -> 'O'
HETATM 3648  O   HOH A 501      -> 'O'
HETATM    1  BR  MOL A   1      -> 'BR'
HETATM 2225 BR1  HMD A 400      -> 'BR'
HETATM 5275 MG    MG A1750      -> 'MG'
HETATM 5276  PG  GCP D1749      -> 'P'
HETATM 5281  PB  GCP D1749      -> 'P'
HETATM 7746  K     K A1525      -> 'K'
HETATM 7747  P   JLN A1526      -> 'P'
HETATM41689 CA1  OEX A 601      -> 'CA'
HETATM 1135 CA    CA A 149      -> 'CA'
HETATM41690 MN1  OEX A 601      -> 'MN'
ATOM   1678  SG  CYS A 206      -> 'S'
HETATM 2175  S1 AD2S A 304      -> 'S'
HETATM 2177  S2 AD2S A 304      -> 'S'
HETATM 2155  F13 V50 A 304      -> 'F'
HETATM 2098 F2'1 RZ7 A 262      -> 'F'
HETATM 1237 FE   HEM A 155      -> 'FE'



>> Rule to distinguish hydrogens(H) or carbons(C) and others (HG and HF, CA and CL) << .
 (1) If the first symbol is 'space' and the second is 'H' or 'C' (' H  ', ' HA  ',  ' HN5', ' C  ',' CA ') 
      -> element := 'H' or 'C' 
 (2) If the first symbol is 'H' or 'C'

    (2a) if it starts with 'CL', ('CL  ', 'CL7A','CL7A') 
         element := 'CL'
    (2b) else if the length of atomname == 4 ('HN5A', 'HG12','CD11','CZ21') 
      -> element := 'H' or 'C'
    (3c) else if the length of atomname < 4  ('HG  ', 'HFA ',  'HFE ', 'CL  ' 'CA  ','CA1 ') 
      -> element := the first and second symobol. (such as 'HG', 'HF')
*/
  int i,len;

  len = 0;
  for (i=0;i<strlen(atomname);++i){
    if ((atomname[i]!=' ')&&(atomname[i] != '\0')){ len += 1;}
  }

  if (atomname[0]!=' '){
    if ((atomname[0]=='H')||(atomname[0]=='C')){
      if ((atomname[0]=='C')||(atomname[0]=='L')){
        element[0] = 'C'; element[1] = 'L'; element[2] = '\0';
      }
      else if (len==4){
       /* (2) the first symbol is 'H' or 'C' and the length of atomname == 4 ('HN5A', 'HG12','CD11','CE12')  */
        element[0] = atomname[0];
        element[1] = '\0';
      }
      else if (len<4){
       /* (3) the first symbol is 'H' or 'C' and the length of atomname < 4  ('HG  ', 'HFA ',  'HFE ', 'CA  ','CA1 ')  */
        element[0] = atomname[0]; 
        element[1] = atomname[1];
        element[2] = '\0';
      }
    }
    else{
      if (isalpha(atomname[1])){
        /* 'MN1 ','BR1 ' */
        element[0] = atomname[0];
        element[1] = atomname[1];
        element[2] = '\0';
      }
      else{
        /* 'F2'1' */
        element[0] = atomname[0];
        element[1] = '\0';
      }
    }
  }
  else if (atomname[0]==' '){  
        if ((atomname[1]=='B') && (atomname[2]=='R'))  sprintf(element,"BR"); /* ' BR ' */
   else if ((atomname[1]=='C') && (atomname[2]=='L'))  sprintf(element,"CL"); /* ' CL ' */
   else{ /* ' H  ', ' HA  ', ' C1 ', ' O  ', ' O1G ', ' PG ', ' PB ' */
     element[0] = atomname[1];
     element[1] = '\0';
   }
 }

 /* printf("#atomname '%s' --> element '%s'\n",atomname,element); */

} /* end of Element_from_PDB_Atomname() */


int Check_Element_String_Or_Not(ele)
  char *ele;
{
 int i,Nele;
 char hit;
 Nele = 54;
 i = 0;
 hit = 0;
 while ((i<Nele) && (hit==0)){
  if (strcmp(ele,element_list[i])==0)       {hit = 1; return(1);}
  if (strcmp(ele,element_list_lower[i])==0) {hit = 1; return(1);}
  i += 1;
 }
 return(0);
} /* end of Check_Element_String_Or_Not() */

