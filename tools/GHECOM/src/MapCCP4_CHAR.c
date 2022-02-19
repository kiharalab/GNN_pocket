/*

 <MapCCP4_CHAR.c>

 Functions for changing CCP4 map file to a simple text file

=========== "ghecom" program ===========

Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under
the GNU Lesser General Public License (LGPL) version 3, see LICENSE.txt.

=========== Installation and Usage ===========

See html file under the "doc/" directory.

=========== Citing "ghecom" program ===========

Please cite:

1) Kawabata T. (2010) Detection of multi-scale pockets on protein surfaces using mathematical morphology. Proteins,78, 1195-1121.



>> FILE FORMAT CCP4 MAP DATA <<
 (taken from http://xtal.pharm.nwu.edu/info/ccp4/maplib.html)
 The term "word" may be 4 characters in following text.


2) DETAILED DESCRIPTION OF THE MAP FORMAT
The overall layout of the file is as follows:

a) File header (256 longwords) 
b) Symmetry information 
c) Map, stored as a 3-dimensional array 

The files are read & written using the diskio package, which allows the file to be treated as a direct-access byte stream, essentially as in C fread & fwrite calls.

The header is organised as 56 words followed by space for ten 80 character text labels as follows:

 
 1      NC              # of Columns    (fastest changing in map)
 2      NR              # of Rows
 3      NS              # of Sections   (slowest changing in map)
 4      MODE            Data type
                          0 = envelope stored as signed bytes (from
                              -128 lowest to 127 highest)
                          1 = Image     stored as Integer*2
                          2 = Image     stored as Reals
                          3 = Transform stored as Complex Integer*2
                          4 = Transform stored as Complex Reals
                          5 == 0	
 
                          Note: Mode 2 is the normal mode used in
                                the CCP4 programs. Other modes than 2 and 0
                                may NOT WORK
 
 5      NCSTART         Number of first COLUMN  in map
 6      NRSTART         Number of first ROW     in map
 7      NSSTART         Number of first SECTION in map
 8      NX              Number of intervals along X
 9      NY              Number of intervals along Y
10      NZ              Number of intervals along Z
11      X length        Cell Dimensions (Angstroms)
12      Y length                     "
13      Z length                     "
14      Alpha           Cell Angles     (Degrees)
15      Beta                         "
16      Gamma                        "
17      MAPC            Which axis corresponds to Cols.  (1,2,3 for X,Y,Z)
18      MAPR            Which axis corresponds to Rows   (1,2,3 for X,Y,Z)
19      MAPS            Which axis corresponds to Sects. (1,2,3 for X,Y,Z)
20      AMIN            Minimum density value
21      AMAX            Maximum density value
22      AMEAN           Mean    density value    (Average)
23      ISPG            Space group number
24      NSYMBT          Number of bytes used for storing symmetry operators
25      LSKFLG          Flag for skew transformation, =0 none, =1 if foll
26-34   SKWMAT          Skew matrix S (in order S11, S12, S13, S21 etc) if
                        LSKFLG .ne. 0.
35-37   SKWTRN          Skew translation t if LSKFLG .ne. 0.
                        Skew transformation is from standard orthogonal
                        coordinate frame (as used for atoms) to orthogonal
                        map frame, as
 
                                Xo(map) = S * (Xo(atoms) - t)
 
38      future use       (some of these are used by the MSUBSX routines
 .          "              in MAPBRICK, MAPCONT and FRODO)
 .          "   (all set to zero by default)
 .          "
52          "

53	MAP	        Character string 'MAP ' to identify file type
54	MACHST		Machine stamp indicating the machine type
			which wrote file
55      ARMS            Rms deviation of map from mean density
56      NLABL           Number of labels being used
57-256  LABEL(20,10)    10  80 character text labels (ie. A4 format)


Symmetry records follow - if any - stored as text as in International Tables, operators separated by * and grouped into 'lines' of 80 characters (ie. symmetry operators do not cross the ends of the 80-character 'lines' and the 'lines' do not terminate in a *).

Map data array follows. 

---------------------------------------------------------------------------------------
                                                                                                            
>> FILE FORMAT "MRC" MAP DATA <<
                                                                                                            
The following information is taken from
  http://splorg.org:8080/~tobin/projects/downing/mrc/mrc_format.html
  http://ami.scripps.edu/prtl_data/mrc_specification.htm
                                                                                                            
* The term "word" may be 4 characters in following text.
* Basically, it is the same format as the CCP4 format.


 format for MRC image files
                                                                                                            
*                Map/Image Header Format for imsubs2000                 *
*        Length = 1024 bytes, organized as 56 LONG words followed       *
*       by space for 10 80 byte text labels.                            *
*                                                                       *
*       1       NX      number of columns (fastest changing in map)     *
*       2       NY      number of rows                                  *
*       3       NZ      number of sections (slowest changing in map)    *
*       4       MODE    data type :                                     *
*                       0       image : signed 8-bit bytes range -128   *
*                                       to 127                          *
*                       1       image : 16-bit halfwords                *
*                       2       image : 32-bit reals                    *
*                       3       transform : complex 16-bit integers     *
*                       4       transform : complex 32-bit reals        *
*       5       NXSTART number of first column in map (Default = 0)     *
*       6       NYSTART number of first row in map       "              *
*       7       NZSTART number of first section in map   "              *
*       8       MX      number of intervals along X                     *
*       9       MY      number of intervals along Y                     *
*       10      MZ      number of intervals along Z                     *
*       11-13   CELLA   cell dimensions in angstroms                    *
*       14-16   CELLB   cell angles in degrees                          *
*       17      MAPC    axis corresp to cols (1,2,3 for X,Y,Z)          *
*       18      MAPR    axis corresp to rows (1,2,3 for X,Y,Z)          *
*       19      MAPS    axis corresp to sections (1,2,3 for X,Y,Z)      *
*       20      DMIN    minimum density value                           *
*       21      DMAX    maximum density value                           *
*       22      DMEAN   mean density value                              *
*       23      ISPG    space group number 0 or 1 (default=0)           *
*       24      NSYMBT  number of bytes used for symmetry data (0 or 80)*
*       25-49   EXTRA   extra space used for anything  - 0 by default   *
*       50-52   ORIGIN  origin in X,Y,Z used for transforms             *
*       53      MAP     character string 'MAP ' to identify file type   *
*       54      MACHST  machine stamp                                   *
*       55      RMS     rms deviation of map from mean density          *
*       56      NLABL   number of labels being used                     *
*       57-256  LABEL(20,10) 10 80-character text labels                *
*                                                                       *
*       Symmetry records follow - if any - stored as text as in         *
*       International Tables, operators separated by * and grouped into *
*       'lines' of 80 characters (ie. symmetry operators do not cross   *
*       the ends of the 80-character 'lines' and the 'lines' do not     *
*       terminate in a *).                                              *
*                                                                       *
*       Data records follow.                                            *
---------------------------------------------------------------------------------------


*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "Grid3D.h" 

/*** FUNCTIONS (GLOBAL) ***/
void Read_MapCCP4_File_CHAR3DMAP();
void Write_MapCCP4_File_CHAR3DMAP();

/*** FUNCTIONS (LOCAL) ***/
static void Read_String();
static int   Read_Int();
static float Read_Float();
static void Inverse_4bytes();


void Read_MapCCP4_File_CHAR3DMAP(fname,cmap)
  char   *fname; 
  struct CHAR3DMAP *cmap;
{
 FILE *fp;
 char Order,buff[128];
 int i,j,x,y,z;
 int NC,NR,NS,MODE,NCSTART,NRSTART,NSSTART,NX,NY,NZ,MAPC,MAPR,MAPS,ISPG,NSYMBT,LSKFLG,NLABL;
 float Xlength, Ylength, Zlength,Alpha,Beta,Gamma,AMIN,AMAX,AMEAN,ARMS;
 float SKWMAT[3][3],SKWTRN[3];
 char MAP[5],MACHST[5],LABEL[10][81];

 printf("#Read_MapCCP4_File(\"%s\")\n",fname);
 /* Checking byte number for int and float */

 if (sizeof(int)!=4) {printf("#ERROR:Size of int(%d) is not 4\n",(int)sizeof(int));}
 if (sizeof(float)!=4) {printf("#ERROR:Size of float(%d) is not 4\n",(int)sizeof(float));}

 fp = fopen(fname,"r");
 if (fp==NULL) {printf("#ERROR:Can't open mapfile \"%s\"\n",fname); exit(1);}

 /** Deciding "Byte Order" from the value of NC,**/ 
 Order = 'N'; /* First try Normal byte order */
 NC = Read_Int(fp,Order,"NC"); 
 if ((NC > 1024)||(NC<0))  /* If NC is too large or negative, try Inverse byte order */
 { Order = 'I';
   Inverse_4bytes(&NC); 
   printf("NC %d\n",NC); }
 /** Read Headers **/
 NR = Read_Int(fp,Order,"NR"); 
 NS = Read_Int(fp,Order,"NS"); 
 MODE = Read_Int(fp,Order,"MODE"); 
 NCSTART = Read_Int(fp,Order,"NCSTART"); 
 NRSTART = Read_Int(fp,Order,"NRSTART"); 
 NSSTART = Read_Int(fp,Order,"NSSTART"); 
 NX = Read_Int(fp,Order,"NX"); 
 NY = Read_Int(fp,Order,"NY"); 
 NZ = Read_Int(fp,Order,"NZ"); 
 Xlength = Read_Float(fp,Order,"Xlength"); 
 Ylength = Read_Float(fp,Order,"Ylength"); 
 Zlength = Read_Float(fp,Order,"Zlength"); 
 Alpha = Read_Float(fp,Order,"Alpha"); 
 Beta = Read_Float(fp,Order,"Beta"); 
 Gamma = Read_Float(fp,Order,"Gamma"); 
 MAPC = Read_Int(fp,Order,"MAPC"); 
 MAPR = Read_Int(fp,Order,"MAPR"); 
 MAPS = Read_Int(fp,Order,"MAPS"); 
 AMIN = Read_Float(fp,Order,"AMIN"); 
 AMAX = Read_Float(fp,Order,"AMAX"); 
 AMEAN = Read_Float(fp,Order,"AMEAN"); 
 ISPG = Read_Int(fp,Order,"ISPG"); 
 NSYMBT = Read_Int(fp,Order,"NSYMBT"); 
 LSKFLG = Read_Int(fp,Order,"LSKFLG"); 
  
 for (i=0;i<3;++i) 
  for (j=0;j<3;++j) 
   { sprintf(buff,"SKWMAT%d%d",i+1,j+1); 
     SKWMAT[i][j] = Read_Float(fp,Order,buff); }

  for (i=0;i<3;++i) 
  { sprintf(buff,"SKWTRN%d",i+1); 
    SKWTRN[i] = Read_Float(fp,Order,buff); }
  
  for (i=0;i<15;++i) Read_Int(fp,Order,"-"); 

  fread(MAP,1,4,fp); MAP[4] = '\0';
  printf("MAP \"%s\"\n",MAP);
 
  fread(MACHST,1,4,fp); MACHST[4] = '\0';
  printf("MACHST \"%s\"\n",MACHST);
  
  ARMS = Read_Float(fp,Order,"ARMS"); 
  NLABL = Read_Int(fp,Order,"NLABL"); 

  for (i=0;i<10;++i) 
  { sprintf(buff,"LABEL%d",i+1); 
    Read_String(fp,LABEL[i],80,buff); }

  /** Malloc Voxel **/
  cmap->N[0] = NC;   cmap->N[1] = NR;   cmap->N[2] = NS;
 
  Malloc_CHAR3DMAP(cmap,cmap->N[0],cmap->N[1],cmap->N[2]);   
  cmap->grid_width = (float)Xlength/(float)NX;

  cmap->OrigPos[0] = NCSTART * cmap->grid_width;
  cmap->OrigPos[1] = NRSTART * cmap->grid_width;
  cmap->OrigPos[2] = NSSTART * cmap->grid_width;
 
  /** Read Voxel **/
  for (z=0;z<cmap->N[2];++z){ 
    for (y=0;y<cmap->N[1];++y){
      for (x=0;x<cmap->N[0];++x){
        cmap->map[x][y][z] =  Read_Float(fp,Order,"-"); 
      }
    }
  }
 
  fclose(fp);

} /* end of Read_MapCCP4_File_CHAR3DMAP() */




void Write_MapCCP4_File_CHAR3DMAP(fname, cmap, Char2Float, Cutoff,Nradius, Rarray)
  char   *fname;
  struct CHAR3DMAP *cmap;
  char   Char2Float;  /* 'D'irect, 'C'utoff 'R'array*/
  float  Cutoff;   /* (only for Char2Float == 'C') */ 
  int    Nradius;  /* (only for Char2Float == 'R') */
  float  *Rarray;  /* (only for Char2Float == 'R') */
{
 FILE *fp;
 int i,x,y,z,init;
 float v;
 int NCSTART,NRSTART,NSSTART;
 float AMIN,AMAX,AMEAN,ARMS;
 char MAP[5],MACHST[5],LABEL[10][81];
 double SUM;
 

 printf("#Write_MapCCP4_File_CHAR3DMAP(%d %d %d Char2Float %c)-->\"%s\"\n",cmap->N[0],cmap->N[1],cmap->N[2],Char2Float,fname);

 fp = fopen(fname,"w");
 if (fp==NULL) {printf("#ERROR:Can't write to CCP4file \"%s\"\n",fname); exit(1);}

 /** Calculate AMIN, AMAX, AMEAN,ARMS **/
 SUM = 0.0;
 init = 1;
 for (x=0;x<cmap->N[0];++x){
   for (y=0;y<cmap->N[1];++y){
     for (z=0;z<cmap->N[2];++z){
       v = 0.0;
       if (Char2Float=='D'){
         v = (float)cmap->map[x][y][z]; 
       }
       else if (Char2Float=='C'){
         if (v>=Cutoff){ v = 1.0;}
         else { v = 0.0;}
       }
       else if (Char2Float=='R'){
              if (cmap->map[x][y][z]==255){v = 0.0;}
         else if (cmap->map[x][y][z]==0) {v = 0.0;}
         else {v = (float)(1.0/Rarray[cmap->map[x][y][z]]);}
       }
/*
         printf("v %f AMIN %f AMAX %f init %d\n",v,AMIN,AMAX,init); 
*/
       if ((v<AMIN)||(init==1)){ AMIN =v;}
       if ((v>AMAX)||(init==1)){ AMAX =v;}
       init = 0;
       SUM += v;  
      }
   }
 }
 
 AMEAN = SUM/cmap->N[0]/cmap->N[1]/cmap->N[2]; 

 SUM = 0.0;
 for (x=0;x<cmap->N[0];++x){
   for (y=0;y<cmap->N[1];++y){
     for (z=0;z<cmap->N[2];++z){
       v = cmap->map[x][y][z];
       SUM += (v-AMEAN)*(v-AMEAN);   
     }
   }
 }
 
 ARMS = SUM/cmap->N[0]/cmap->N[1]/cmap->N[2]; 
 if (ARMS>0.0){ ARMS = sqrt(ARMS);}

 /* NC, NR, NS */
 fwrite(&(cmap->N[0]),1,4,fp); 
 fwrite(&(cmap->N[1]),1,4,fp); 
 fwrite(&(cmap->N[2]),1,4,fp); 
 /* Mode */
 i = 2; fwrite(&i,1,4,fp); 
 /* NCSTART, NRSTART, NSSTART */
 /*  
 cmap->OrigPos[0] = NCSTART * cmap->grid_width;
 cmap->OrigPos[1] = NRSTART * cmap->grid_width;
 cmap->OrigPos[2] = NSSTART * cmap->grid_width;
 */
 NCSTART = (int)(cmap->OrigPos[0]/cmap->grid_width);
 fwrite(&NCSTART,1,4,fp); 
 NRSTART = (int)(cmap->OrigPos[1]/cmap->grid_width);
 fwrite(&NRSTART,1,4,fp); 
 NSSTART = (int)(cmap->OrigPos[2]/cmap->grid_width);
 fwrite(&NSSTART,1,4,fp); 

 
 /* NX, NY, NZ */
 fwrite(&(cmap->N[0]),1,4,fp); 
 fwrite(&(cmap->N[1]),1,4,fp); 
 fwrite(&(cmap->N[2]),1,4,fp); 

 /* X length, Y length, Z length */
 v = cmap->N[0]*cmap->grid_width; fwrite(&v,1,4,fp); 
 v = cmap->N[1]*cmap->grid_width; fwrite(&v,1,4,fp); 
 v = cmap->N[2]*cmap->grid_width; fwrite(&v,1,4,fp); 
 
 /* Alpha, Beta, Gamma */
 v = 90.0;
 fwrite(&v,1,4,fp); 
 fwrite(&v,1,4,fp); 
 fwrite(&v,1,4,fp); 

 /* MAPC, MAPR, MAPS */
 i = 1; fwrite(&i,1,4,fp); 
 i = 2; fwrite(&i,1,4,fp); 
 i = 3; fwrite(&i,1,4,fp); 

 /* AMIN, AMAX, AMEAN */
 printf("#AMIN %f AMAX %f AMEAN %f\n",AMIN,AMAX,AMEAN);
 fwrite(&AMIN,1,4,fp); 
 fwrite(&AMAX,1,4,fp); 
 fwrite(&AMEAN,1,4,fp); 

 /* ISPG, NSYMBT, LSKFLG */
 i = 1; fwrite(&i,1,4,fp); 
 i = 0; fwrite(&i,1,4,fp); 
 i = 0; fwrite(&i,1,4,fp); 
 
 /* SKWMAT11, SKWMAT12, ..., SKWMAT33 */
 v = 0.0;
 for (x=1;x<=3;++x){
   for (y=1;y<=3;++y){ 
     fwrite(&v,1,4,fp);
   }
 }
 
 /* SKWTRN1, SKWTRN2, SKWTRN2 */
 /*
 v = 0.0;
 fwrite(&v,1,4,fp); 
 fwrite(&v,1,4,fp); 
 fwrite(&v,1,4,fp); 
 */
 fwrite(&(cmap->OrigPos[0]),1,4,fp); 
 fwrite(&(cmap->OrigPos[1]),1,4,fp); 
 fwrite(&(cmap->OrigPos[2]),1,4,fp); 
 /*
 printf("cmap->min %f %f %f\n",cmap->min[0],cmap->min[1],cmap->min[2]);
 */
 
 /* future use (from 38 to 52 words) */
 i = 0;
 for (x=38;x<=52;++x) fwrite(&i,1,4,fp); 

 /* MAP, MACHST */
 sprintf(MAP,"MAP "); fwrite(MAP,1,4,fp);
 sprintf(MACHST,"DA"); fwrite(MACHST,1,4,fp);
 
 /* ARMS */
 fwrite(&ARMS,1,4,fp); 

 /* NLABL */
 i = 10; fwrite(&i,1,4,fp); 

 /* LABEL */
 for (x=0;x<10;++x)
 for (y=0;y<80;++y) LABEL[x][y] = ' ';

 sprintf(LABEL[0],"THIS IS MADE BY ghecom program");
 fwrite(LABEL[0],1,80,fp);
 for (x=1;x<10;++x) fwrite(LABEL[x],1,80,fp);
 
 /* write Voxel values */

 for (z=0;z<cmap->N[2];++z){
   for (y=0;y<cmap->N[1];++y){
     for (x=0;x<cmap->N[0];++x){
       v = 0.0;
       if (Char2Float=='D'){
         v = (float)cmap->map[x][y][z]; 
       }
       else if (Char2Float=='C'){
         if (v>=Cutoff){ v = 1.0;}
         else { v = 0.0;}
       }
       else if (Char2Float=='R'){
             if (cmap->map[x][y][z]==255){v = 0.0;}
        else if (cmap->map[x][y][z]==0) {v = 0.0;}
        else {v = (float)(1.0/Rarray[cmap->map[x][y][z]]);}
      }
      fwrite(&v,1,4,fp);
     }
   }
 }

 fclose(fp);                   

} /* end of Write_MapCCP4_File_CHAR3DMAP() */



void Read_String(fp,instr,length,str)
 FILE *fp;
 char *instr;
 int  length;
 char *str;
{ int i;
  int ret;

  ret = fread(instr,1,length,fp); 
  if (ret<length) {printf("#ERROR:Can't read %d chars.\n",length); exit(1);}
  for (i=0;i<length;++i)
    if (iscntrl(instr[i])!=0) instr[i] = ' ';
  instr[length] = '\0';
  if (str[0]!='-') printf("#%s \"%s\"\n",str,instr);

} /* end of Read_String() */



int Read_Int(fp,Order,str)
 FILE *fp;
 char Order; /* 'I'nverse */
 char *str;
{ int i;
  int ret; 
 ret = fread(&i,1,4,fp); 
 if (ret<4) {printf("#ERROR:Can't read 4 chars.\n"); exit(1);}
 if (Order=='I') Inverse_4bytes(&i); 
 if (str[0]!='-') printf("#%s %d\n",str,i);
 return(i); 
} /* end of Read_Int() */



float Read_Float(fp,Order,str)
 FILE *fp;
 char Order; /* 'I'nverse */
 char *str;
{ float f;
  int ret; 

 ret = fread(&f,1,4,fp); 
 if (ret<4) {printf("#ERROR:Can't read 4 chars.\n"); exit(1);}
 if (Order=='I') Inverse_4bytes(&f); 
 if (str[0] !='-') printf("#%s %f\n",str,f);
 return(f); 
} /* end of Read_Int() */




void Inverse_4bytes(X)
 char X[4];
{
 /* 
 X0 X1 X2 X3 
    ->
 X3 X2 X1 X0 
 */
 char buff;
 /*
 printf("%x %x %x %x\n",X[0],X[1],X[2],X[3]);
 */ 
 buff = X[3]; X[3] = X[0]; X[0] = buff;
 buff = X[2]; X[2] = X[1]; X[1] = buff;

}  /* end of Inverse_4bytes() */
