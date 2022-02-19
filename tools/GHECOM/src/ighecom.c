/*

 <ighecom.c> 
 "Grid-based HECOMi finder"

 Program for pocket detection based on the operations
 of Mathematical Morphology


 The source code 'ighecom' is for the minimum functions for detecting pockets.

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
#include "OpeMMorph.h"
#include "OpeStd3DPr.h"
#include "OpeClus3DPr.h"

static char LastModDate[] = "Oct 25, 2007";
struct PARAMETERS PAR;

/*** FUNCTIONS (LOCAL) ***/


main(argc,argv)
 int argc;
 char **argv;
{
 char MODE; 
 int k,L, NatomPdb;
 char corename[MAX_FILE_NAME];
 char oVdWpdbfile[MAX_FILE_NAME],ores_gdbfile[MAX_FILE_NAME];
 char iprbSpdbfile[MAX_FILE_NAME]; 
 char AtomSelect;
 char comment[512];
 struct STRING_NODE ComStrHead;
 struct RADIUS   HeadRadius; 
 float Rall; 
 int SeedRand;
 /* for probe spheres */
 float Rprobe,RprobeS,RprobeL;
 int    Ngrid_probe, Ngrid_probeS,Ngrid_probeL; 
 struct CHAR3DMAP MapP;
 struct CHAR3DMAP MapPS;
 struct CHAR3DMAP MapPL;
 float  margin_width;
 struct ATOM PrbSAtomHead;
 float  ScaleRadiusPrbS;
 
 /* for protein molecules  */
 struct CHAR3DMAP MapX;  /* 3D map for protein shape */
 struct CHAR3DMAP MapY;  /* Additional 3D map for multi-scale probes */
 char ipdbfile[MAX_FILE_NAME];
 struct ATOM ProAtomHead;
 struct RESIDUE ProResHead;
 float MinP[3],MaxP[3];  /* Min and Maximum Values for XYZ cordinates of protein */
 char  ChainID;
 int   Ngrid_VdW,Ngrid;
 char  ReadAtmHetType; 

 /** For detected pockets **/
 int  Ngrid_pocket;
 int  Ncluster, Nbit_cluster[256];
 char  buffstr[128];

 /* for ligand molecule to be docked */
 int    Nrotate_docklig;
 struct ROT_MATRIX *RmatArray;  /* [Nrotate_docklig] (malloc later)*/
 float  Dorig_max;


 /*** Set Default Value ***/
 PAR.SubType[0] = '\0';
 Rprobe = 4.0; 
 RprobeL = 6.0; 
 RprobeS = 1.87;
 MODE = 'P';
 SeedRand = 0;
 ChainID = '-';
 AtomSelect = '-';
 ReadAtmHetType = 'A'; 
 oVdWpdbfile[0] = ores_gdbfile[0] = '\0';
 PAR.grid_width = 0.8;
 iprbSpdbfile[0] = '\0'; 
 ComStrHead.next = NULL;
 PAR.radfile[0] = '\0';
 PAR.NeighborNum = 26;
 Nrotate_docklig = 10; 
 ScaleRadiusPrbS = 1.0;
 
 if (argc<2){
  printf("ighecom [pdbfile] <options>\n");
  printf(" for surface/volume/pocket/cavity in 3D grid represention.\n");
  printf(" coded by T.Kawabata. LastModified:%s\n",LastModDate);
  printf(" <options>\n");
  printf(" -M   : MODE 'D'ilation 'E'rosion molceular 'C'losing 'O'pening\n");
  printf("             'P'ocket(masuya_doi) 'p'ocket(kawabata_go)\n");
  printf("             'A'symmetric probe pocket [%c]\n",MODE);
  printf(" -ip  : Input pdbfile [%s]\n",ipdbfile); 
  printf(" -ch  : ChainID for target pdbfile [%c]\n",ChainID); 
  printf(" -ah  : Read only ATOM 'a', only HETATM 'h', both 'b' [%c]\n",ReadAtmHetType); 
  printf(" -gw  : Grid width [%f]\n",PAR.grid_width); 
  printf(" -ovw : Output VdW protein 3D grid PDB files[%s]\n",oVdWpdbfile); 
  printf(" -org : Output final result 3D grid PDB files[%s]\n",ores_gdbfile); 
  printf(" -rp  : Radius for probe spheres [%f]\n",Rprobe); 
  printf(" -rs  : Radius for small probe spheres [%f]\n",RprobeS); 
  printf(" -rl  : Radius for large probe spheres [%f]\n",RprobeL); 
  printf(" -isp : Input small Probe PDB file [%s]\n",iprbSpdbfile); 
  printf(" -ssp : Scale for radius of small Probe PDB file [%f]\n",ScaleRadiusPrbS); 
  printf(" -sub : SubType for various purpose [%s]\n",PAR.SubType); 
  exit(1);
 }


 /****** READ ARGUMENTS ***********/
 PAR.COMMAND[0] = '\0';
 for (k=0;k<argc;++k)
 { if (k>0) strcat(PAR.COMMAND," ");
    strcat(PAR.COMMAND,argv[k]);   }

 k = 0;
 if (argv[1][0]!='-')
 { sprintf(ipdbfile,"%s",argv[1]);
   printf("#ipdbfile '%s'\n",ipdbfile);
  Find_Filename_Core(corename,ipdbfile); k = 1;}

 while (k<argc){
    if (argv[k][0]=='-'){
      L = strlen(argv[k]);
         if ((L==3)&&(argv[k][1]=='c')&&(argv[k][2]=='h')) {++k; ChainID = argv[k][0];}
    else if ((argv[k][1]=='g')&&(argv[k][2]=='w')) {++k; PAR.grid_width = atof(argv[k]); }
    else if ((L==3)&&(argv[k][1]=='i')&&(argv[k][2]=='p')) {++k; sprintf(ipdbfile,"%s",argv[k]);}
    else if ((L==4)&&(argv[k][1]=='o')&&(argv[k][2]=='r')&&(argv[k][3]=='g')) {++k; sprintf(ores_gdbfile,"%s",argv[k]);}
    else if ((L==4)&&(argv[k][1]=='o')&&(argv[k][2]=='v')&&(argv[k][3]=='w')) {++k; sprintf(oVdWpdbfile,"%s",argv[k]);}
    else if ((L==3)&&(argv[k][1]=='r')&&(argv[k][2]=='p')) {++k; Rprobe  = atof(argv[k]);}
    else if ((L==3)&&(argv[k][1]=='r')&&(argv[k][2]=='s')) {++k; RprobeS = atof(argv[k]);}
    else if ((L==3)&&(argv[k][1]=='r')&&(argv[k][2]=='l')) {++k; RprobeL = atof(argv[k]);}
    else if ((L==2)&&(argv[k][1]=='M')) {++k; MODE = argv[k][0]; }
    else if ((L==2)&&(argv[k][1]=='A')) {++k; AtomSelect = argv[k][0]; }
    else if ((L==3)&&(argv[k][1]=='a')&&(argv[k][2]=='h')) {++k; ReadAtmHetType = argv[k][0];}
    else if ((L==3)&&(argv[k][1]=='n')&&(argv[k][2]=='n')) {++k; PAR.NeighborNum = atoi(argv[k]);}
    else if ((L==4)&&(argv[k][1]=='i')&&(argv[k][2]=='s')&&(argv[k][3]=='p')) 
       {++k; sprintf(iprbSpdbfile,"%s",argv[k]);}
    else if ((L==4)&&(argv[k][1]=='n')&&(argv[k][2]=='r')&&(argv[k][3]=='d')) 
         {++k; Nrotate_docklig  = atoi(argv[k]);}
    else if ((L==4)&&(argv[k][1]=='s')&&(argv[k][2]=='r')&&(argv[k][3]=='a')) 
         {++k; SeedRand  = atoi(argv[k]); }
    else if ((L==4)&&(argv[k][1]=='s')&&(argv[k][2]=='u')&&(argv[k][3]=='b')) 
         {++k; sprintf(PAR.SubType,"%s",argv[k]); }
    else if ((L==4)&&(argv[k][1]=='s')&&(argv[k][2]=='s')&&(argv[k][3]=='p')) 
         {++k; ScaleRadiusPrbS = atof(argv[k]); }
    else { printf("#ERROR:Can't understand option %s\n",argv[k]); exit(1);}
 }
   
  ++k;

 } /* while k */

 srand(SeedRand);
 MapX.grid_width = MapY.grid_width = MapP.grid_width = MapPS.grid_width 
 = MapPL.grid_width = PAR.grid_width;
 
 printf("#COMMAND %s\n",PAR.COMMAND);

  
 /*** [0] Set Radius  ***/

     if ((AtomSelect=='A')||(AtomSelect=='B')) Set_Unified_Radius(&HeadRadius,Rall);
 else if (PAR.radfile[0] != '\0') Read_Radius_File(PAR.radfile,&HeadRadius);
 else Set_Default_Radius(&HeadRadius);

 /*** [1] Read PDB file  ***/
 if (ipdbfile[0]!='\0'){
  Read_PDB_File(ipdbfile,&ProAtomHead,'H',ChainID,ReadAtmHetType,'F');
  NatomPdb = Number_Of_Atom(&ProAtomHead); 
  Assign_Radius(&ProAtomHead,&HeadRadius);
  MinMaxXYZ_Of_Atoms(&ProAtomHead,MinP,MaxP); 
  printf("#MinP %f %f %f MaxP %f %f %f\n", MinP[0], MinP[1], MinP[2], MaxP[0], MaxP[1], MaxP[2]);
  Make_Residue_List(&ProAtomHead,&ProResHead);

 /*** [2] Malloc Huge 3D grids for the PDB molecule  ***/
   if ((MODE=='D')||(MODE=='E')||(MODE=='C')||(MODE=='O'))     margin_width = 2.0*Rprobe; 
  else if ((MODE=='P')||(MODE=='p')||(MODE=='A')||(MODE=='R')) margin_width = 2.0*RprobeL; 
  
 margin_width += MapX.grid_width*2;
 
 Decide_Num_Grid_3D(MapX.N,MapX.OrigPos,MapX.grid_width, margin_width, MinP,MaxP, '-');

 Malloc_CHAR3DMAP(&MapX, MapX.N[0], MapX.N[1], MapX.N[2]);
 
 Ngrid_VdW = Set_Bit_3DMAP_Inside_Atoms(&MapX,&ProAtomHead,1,1.0);
 Ngrid_VdW = Count_Specific_3DMAP(&MapX,1);
 printf("Ngrid_VdW %d\n",Ngrid_VdW);
 sprintf(comment,"pdbfile %s Natom %d", ipdbfile,NatomPdb);
 Add_String_to_STRING_NODE(&ComStrHead,comment); 
 sprintf(comment,"grid_width %6.3f A VdW volume: %d grids, %.2f AAA",
          MapX.grid_width,Ngrid_VdW,Ngrid_VdW*MapX.grid_width*MapX.grid_width*MapX.grid_width);
 Add_String_to_STRING_NODE(&ComStrHead,comment); 

  if (oVdWpdbfile[0]!='\0') 
  { Output_CHAR3DMAP_PDB_format(oVdWpdbfile,&MapX,1,"Protein VdW Volume",&ComStrHead); }

 }


 /*************************/ 
 /** MODE 'D' : Dilation **/
 /*************************/ 

   if (MODE=='D'){ 
    Ngrid_probe = Make_Sphere_Probe_Map(&MapP,Rprobe,MapP.grid_width);
    Add_String_to_STRING_NODE(&ComStrHead,"woops1"); 
    Add_String_to_STRING_NODE(&ComStrHead,"woops2"); 
    Output_CHAR3DMAP_PDB_format("probe.gdb",&MapP,1,"Probe Sphere",&ComStrHead);
    Dilation(&MapX,&MapP,1,4);
    if (ores_gdbfile[0]=='\0') sprintf(ores_gdbfile,"%s_dila.gdb",corename);
    Output_CHAR3DMAP_PDB_format(ores_gdbfile,&MapX,1,"original(1) and Dilated (4) protein",&ComStrHead); 
  }

 /************************/ 
 /** MODE 'E' : Erosion **/
 /************************/ 
 
 else if (MODE=='E')
  { Ngrid_probe = Make_Sphere_Probe_Map(&MapP,Rprobe,MapP.grid_width);
    Output_CHAR3DMAP_PDB_format("probe.gdb",&MapP,1,"Probe Sphere",&ComStrHead);
    Erosion(&MapX,&MapP,1,4);
    if (ores_gdbfile[0]=='\0') sprintf(ores_gdbfile,"%s_eros.gdb",corename);
    Output_CHAR3DMAP_PDB_format(ores_gdbfile,&MapX,1,"original(1) and Erosion(4)  protein",&ComStrHead); 
 }
 
 /**********************************************/
 /** MODE 'C': Closing (or moleculer Surface )**/
 /**********************************************/

 else if (MODE=='C'){ 
    Ngrid_probe = Make_Sphere_Probe_Map(&MapP,Rprobe,MapP.grid_width);
    Output_CHAR3DMAP_PDB_format("probe.gdb",&MapP,1,"Probe Sphere",&ComStrHead);
    
    Dilation(&MapX,&MapP,1,4);
    
    Ngrid = Erosion(&MapX,&MapP,4,8);

    Output_CHAR3DMAP_PDB_format("clos.gdb",&MapX,1,"Closing(molecular surace)",&ComStrHead); 
    
    sprintf(comment,"Molecular Volume %d grids, %.2f AAA",
       Ngrid,Ngrid*MapX.grid_width*MapX.grid_width*MapX.grid_width);
    Add_String_to_STRING_NODE(&ComStrHead,comment);

    Reset_Specific_3DMAP(&MapX,4);
   /*
    Output_CHAR3DMAP_PDB_format("eros.gdb",&MapX,1,"erosion",&ComStrHead); 
    */ 
    Extract_Specific_3DMAP(&MapX,8);
    Extract_Surface_3DMAP(&MapX,8,16);
    Extract_Specific_3DMAP(&MapX,16);
    if (ores_gdbfile[0]=='\0') sprintf(ores_gdbfile,"%s_clos.gdb",corename);
    Output_CHAR3DMAP_PDB_format(ores_gdbfile,&MapX,1,"Molecular Surface with probe",&ComStrHead); 
  }


 /*************************/
 /** MODE 'O' : Opening  **/
 /*************************/

 else if (MODE=='O'){ 
    Ngrid_probe = Make_Sphere_Probe_Map(&MapP,Rprobe,MapP.grid_width);
    Output_CHAR3DMAP_PDB_format("probe.gdb",&MapP,1,"Probe Sphere",&ComStrHead);
  
    Ngrid = Erosion(&MapX,&MapP,1,2);
    sprintf(comment,"Eroded Volume %d grids, %.2f AAA",
       Ngrid,Ngrid*MapX.grid_width*MapX.grid_width*MapX.grid_width);
    Add_String_to_STRING_NODE(&ComStrHead,comment);
    Dilation(&MapX,&MapP,2,4);
    Ngrid = Extract_Specific_3DMAP(&MapX,4);
    if (ores_gdbfile[0]=='\0') sprintf(ores_gdbfile,"%s_open.gdb",corename);
    Output_CHAR3DMAP_PDB_format(ores_gdbfile,&MapX,1,"Opening",&ComStrHead); 
  }



 /******************************************************/
 /** MODE 'P': Pocket detection  (non-contact pocket) **/
 /**  (definition of Masuya and Doi (1995))           **/
 /******************************************************/
 
  else if (MODE=='P'){ 
    Ngrid_probeS = Make_Sphere_Probe_Map(&MapPS,RprobeS,MapPS.grid_width);

    Output_CHAR3DMAP_PDB_format("probeS.gdb",&MapPS,1,"Small probe sphere",&ComStrHead);
        
    Ngrid_probeL = Make_Sphere_Probe_Map(&MapPL,RprobeL,MapPL.grid_width);
    Output_CHAR3DMAP_PDB_format("probeL.gdb",&MapPL,1,"Large probe sphere",&ComStrHead);
   
    printf("#Ngrid_probeS %4d probeL %5d\n",Ngrid_probeS,Ngrid_probeL);

    Dilation(&MapX,&MapPL,1,2);
    Erosion(&MapX,&MapPL,2,4);

    Extract_A_and_notB_3DMAP(&MapX,4,1,8);
    
    Erosion(&MapX,&MapPS,8,16);
    Dilation(&MapX,&MapPS,16,64);
    
    if (ores_gdbfile[0]=='\0') sprintf(ores_gdbfile,"%s_pock.gdb",corename);
     Output_CHAR3DMAP_PDB_format(ores_gdbfile,&MapX,64,"Pocket Volume",&ComStrHead);  
 }
 

 /******************************************************/
 /** MODE 'p': Pocket detection  (contact pocket) **/
 /**  (definition of Kawbata and Go (2006))           **/
 /******************************************************/
 else if (MODE=='p'){

    Ngrid_probeS = Make_Sphere_Probe_Map(&MapPS,RprobeS,MapPS.grid_width);
    Output_CHAR3DMAP_PDB_format("probeS.gdb",&MapPS,1,"Small probe sphere",&ComStrHead);

    Ngrid_probeL = Make_Sphere_Probe_Map(&MapPL,RprobeL,MapPL.grid_width);
    Output_CHAR3DMAP_PDB_format("probeL.gdb",&MapPL,1,"Large probe sphere",&ComStrHead);
   
    printf("#Ngrid_probeS %4d probeL %5d\n",Ngrid_probeS,Ngrid_probeL);
 
    Dilation(&MapX,&MapPL,1,2);
    Erosion(&MapX,&MapPL,2,4);

    Dilation(&MapX,&MapPS,1,8);
    Dilation(&MapX,&MapPS,8,16);
   
    Extract_A_and_B_and_notC_3DMAP(&MapX,4,16,1,32);

    Erosion(&MapX,&MapPS,32,64);
    Dilation(&MapX,&MapPS,64,128);
  
    sprintf(comment,"Pocket detected with Rsmall: %.2f A Rlarge: %.2f A. Volume: %d grids, %.2f AAA",
     RprobeS, RprobeL, Ngrid_pocket, Ngrid_pocket*MapX.grid_width*MapX.grid_width*MapX.grid_width);  
    Add_String_to_STRING_NODE(&ComStrHead,comment);  

    if (ores_gdbfile[0]=='\0') sprintf(ores_gdbfile,"%s_pock.gdb",corename);
     Output_CHAR3DMAP_PDB_format(ores_gdbfile,&MapX,128,"Pocket Volume",&ComStrHead);  

 }


 /***********************************************/
 /** MODE 'A': Asymmetric (small) probe pocket **/
 /***********************************************/
 else if (MODE=='A'){
   if (iprbSpdbfile[0] != '\0'){
    /* Read Docking Ligand file */ 
    Read_PDB_File(iprbSpdbfile,&PrbSAtomHead,    'H','-','B','F');
    Assign_Radius(&PrbSAtomHead,&HeadRadius);
    MinMaxXYZ_Of_Atoms(&PrbSAtomHead,MinP,MaxP); 
    printf("#SmallProbeAtoms MinP %f %f %f MaxP %f %f %f\n", 
         MinP[0], MinP[1], MinP[2], MaxP[0], MaxP[1], MaxP[2]);
    Decide_Num_Grid_3D(MapPS.N,MapPS.OrigPos,MapPS.grid_width,0.0,MinP,MaxP,'O');
    Malloc_CHAR3DMAP(&MapPS, MapPS.N[0], MapPS.N[1], MapPS.N[2]);
    Ngrid_VdW = Set_Bit_3DMAP_Inside_Atoms(&MapPS,&PrbSAtomHead,1,ScaleRadiusPrbS);
   }   
   else{ 
    Ngrid_probeS = Make_Sphere_Probe_Map(&MapPS,RprobeS,MapPS.grid_width);
  } 

   Output_CHAR3DMAP_PDB_format("probeS.gdb",&MapPS,1,"Small probe sphere",&ComStrHead);

   Ngrid_probeL = Make_Sphere_Probe_Map(&MapPL,RprobeL,MapPL.grid_width);
   Output_CHAR3DMAP_PDB_format("probeL.gdb",&MapPL,1,"Large probe sphere",&ComStrHead);

   Dilation(&MapX,&MapPL,1,2);
   Erosion(&MapX,&MapPL,2,4);

   Extract_A_and_notB_3DMAP(&MapX,4,1,8);

   Opening(&MapX,&MapPS,8,64);
   
   if (ores_gdbfile[0]=='\0') sprintf(ores_gdbfile,"%s_pock.gdb",corename);
   Output_CHAR3DMAP_PDB_format(ores_gdbfile,&MapX,64,"Pocket using asymmetric probe",&ComStrHead); 
 }

 else {
  printf("#ERROR: MODE '%c' is improper.\n",MODE);
 } 

}/* end of main() */
