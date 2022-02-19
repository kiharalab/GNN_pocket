/*

<Amino.c>

coded by Takeshi Kawabata.

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


char Amino_Number();
int Number_Amino();
char TriAA_to_OneAA();
void Get_TriAmino();
void Set_AAnum();

/******** Amino Acid Table (globally defined) **********/
char  AAnum[21];     /* Number -> Amino  table */
int   numAA[256];    /* Amino  -> Number table */


char Amino_Number(num)
 int num;
{
 char i;

 i = 'A';
 switch(num){
   case 0:i = 'A'; break;
   case 1:i = 'I'; break;
   case 2:i = 'L'; break;
   case 3:i = 'M'; break;
   case 4:i = 'F'; break;
   case 5:i = 'P'; break;
   case 6:i = 'V'; break;
   case 7:i = 'R'; break;
   case 8:i = 'D'; break;
   case 9:i = 'E'; break;
   case 10:i ='K'; break;
   case 11:i = 'N'; break;
   case 12:i = 'C'; break;
   case 13:i = 'Q'; break;
   case 14:i = 'H'; break;
   case 15:i = 'S'; break;
   case 16:i = 'T'; break;
   case 17:i = 'W'; break;
   case 18:i = 'Y'; break;
   case 19:i = 'G'; break;
  }
 return(i);

} /* end of Amino_Number() */







int Number_Amino(rsym)
 char rsym;
{
 int i;
 switch(rsym)
  {
   case 'A':i = 0; break;
   case 'I':i = 1; break;
   case 'L':i = 2; break;
   case 'M':i = 3; break;
   case 'm':i = 3; break;
   case 'F':i = 4; break;
   case 'P':i = 5; break;
   case 'p':i = 5; break;
   case 'V':i = 6; break;
   case 'R':i = 7; break;
   case 'D':i = 8; break;
   case 'E':i = 9; break;
   case 'K':i = 10; break;
   case 'N':i = 11; break;
   case 'C':i = 12; break;
   case 'Q':i = 13; break;
   case 'H':i = 14; break;
   case 'S':i = 15; break;
   case 'T':i = 16; break;
   case 'W':i = 17; break;
   case 'Y':i = 18; break;
   case 'G':i = 19; break;
   case 'B':i =  8; break;
   case 'Z':i =  9; break;
   case 'X':i = 0; break;
   default:i = 0; break;
  }
 return(i);

} /* end of Number_Amino() */




char TriAA_to_OneAA(t)
 char *t;
{ char r;
      if (strncmp(t,"ALA",3)==0)  r = 'A';
 else if (strncmp(t,"GLU",3)==0)  r = 'E';
 else if (strncmp(t,"GLN",3)==0)  r = 'Q';
 else if (strncmp(t,"ASP",3)==0)  r = 'D';
 else if (strncmp(t,"ASN",3)==0)  r = 'N';
 else if (strncmp(t,"LEU",3)==0)  r = 'L';
 else if (strncmp(t,"GLY",3)==0)  r = 'G';
 else if (strncmp(t,"LYS",3)==0)  r = 'K';
 else if (strncmp(t,"SER",3)==0)  r = 'S';
 else if (strncmp(t,"VAL",3)==0)  r = 'V';
 else if (strncmp(t,"ARG",3)==0)  r = 'R';
 else if (strncmp(t,"THR",3)==0)  r = 'T';
 else if (strncmp(t,"PRO",3)==0)  r = 'P';
 else if (strncmp(t,"ILE",3)==0)  r = 'I';
 else if (strncmp(t,"MET",3)==0)  r = 'M';
 else if (strncmp(t,"PHE",3)==0)  r = 'F';
 else if (strncmp(t,"TYR",3)==0)  r = 'Y';
 else if (strncmp(t,"CYS",3)==0)  r = 'C';
 else if (strncmp(t,"TRP",3)==0)  r = 'W';
 else if (strncmp(t,"HIS",3)==0)  r = 'H';
 else if (strncmp(t,"UNK",3)==0)  r = 'u';
 else if (strncmp(t,"ACE",3)==0)  r = 'a';
 else if (strncmp(t,"FOR",3)==0)  r = 'f';
 else r = '?';
 return(r);

} /* end of TriAA_to_OneAA() */





void Get_TriAmino(r,t)
 char r,*t;
{     if (r == 'A') sprintf(t,"ALA");
 else if (r == 'E') sprintf(t,"GLU");
 else if (r == 'Q') sprintf(t,"GLN");
 else if (r == 'D') sprintf(t,"ASP");
 else if (r == 'N') sprintf(t,"ASN");
 else if (r == 'L') sprintf(t,"LEU");
 else if (r == 'G') sprintf(t,"GLY");
 else if (r == 'K') sprintf(t,"LYS");
 else if (r == 'S') sprintf(t,"SER");
 else if (r == 'V') sprintf(t,"VAL");
 else if (r == 'R') sprintf(t,"ARG");
 else if (r == 'T') sprintf(t,"THR");
 else if (r == 'P') sprintf(t,"PRO");
 else if (r == 'I') sprintf(t,"ILE");
 else if (r == 'M') sprintf(t,"MET");
 else if (r == 'F') sprintf(t,"PHE");
 else if (r == 'Y') sprintf(t,"TYR");
 else if (r == 'C') sprintf(t,"CYS");
 else if (r == 'W') sprintf(t,"TRP");
 else if (r == 'H') sprintf(t,"HIS");
 else if (r == 'u') sprintf(t,"UNK");
 else if (r == 'a') sprintf(t,"ACE");
 else if (r == 'f') sprintf(t,"FOR");
 else sprintf(t,"XXX");
} /* end of Get_TriAmino() */





void Set_AAnum()
{ int i;
  for (i=0;i<256;++i) numAA[i] = 20;
  numAA['A'] = 0;  AAnum[0]  = 'A';
  numAA['I'] = 1;  AAnum[1]  = 'I';
  numAA['L'] = 2;  AAnum[2]  = 'L';
  numAA['M'] = 3;  AAnum[3]  = 'M';
  numAA['m'] = 3;  
  numAA['F'] = 4;  AAnum[4]  = 'F';
  numAA['P'] = 5;  AAnum[5]  = 'P';
  numAA['p'] = 5;  
  numAA['V'] = 6;  AAnum[6]  = 'V';
  numAA['R'] = 7;  AAnum[7]  = 'R';
  numAA['D'] = 8;  AAnum[8]  = 'D';
  numAA['E'] = 9;  AAnum[9]  = 'E';
  numAA['K'] = 10; AAnum[10] = 'K';
  numAA['N'] = 11; AAnum[11] = 'N';
  numAA['C'] = 12; AAnum[12] = 'C';
  numAA['Q'] = 13; AAnum[13] = 'Q';
  numAA['H'] = 14; AAnum[14] = 'H';
  numAA['S'] = 15; AAnum[15] = 'S';
  numAA['T'] = 16; AAnum[16] = 'T';
  numAA['W'] = 17; AAnum[17] = 'W';
  numAA['Y'] = 18; AAnum[18] = 'Y';
  numAA['G'] = 19; AAnum[19] = 'G';
  numAA['X'] = 20; AAnum[20] = 'X';
}

