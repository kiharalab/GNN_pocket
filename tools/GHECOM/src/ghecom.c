/*

 <ghecom.c> 

 "Grid-based HECOMi finder"


 Program for pocket detection based on the operations
 of Mathematical Morphology

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
#include "MscProbe.h"
#include "MscShPrb.h"
#include "MscLabelAtm.h"
#include "MscLabelRes.h"
#include "MscPockClus.h"
#include "MscSkelClus.h"
#include "HeapSort.h"
#include "OpnRotPrb.h"
#include "TauOpnRotPrb.h"
#include "MapCCP4_CHAR.h"
#include "MapCCP4.h"
#include "LigsitePSP.h"
#include "GenSphePrb.h"
#include "ClusSphePrb.h"
#include "TreeSphePrb.h"
#include "FldFillQue.h"
#include "ReadOpt.h"
#include "ClusBinPock.h"
#include "C3mapMergeSort.h"
#include "ConSpheres.h"
#include "ATOMs_ghecom_from_PDBmmCIF.h"



struct PARAMETERS PAR;

/*** FUNCTIONS (LOCAL) ***/

int main(argc,argv)
 int argc;
 char **argv;
{
 char MODE[8]; 
 int i,k, NatomPdb;
 char oVdWpdbfile[MAX_FILE_NAME],oVdWmapfile[MAX_FILE_NAME];
 char  oMoVpdbfile[MAX_FILE_NAME], oMoVmapfile[MAX_FILE_NAME];
 char opockpdbfile[MAX_FILE_NAME], opockmapfile[MAX_FILE_NAME];
 char omprobfile[MAX_FILE_NAME],omapfile[MAX_FILE_NAME]; 
 char oprbSfile[MAX_FILE_NAME],oprbLfile[MAX_FILE_NAME]; 
 char comment[512];
 struct STRING_NODE ComStrHead;
 char   SurfacePocket;
 /* for probe spheres */
 float RprobeS,RprobeL,RprobeLmin,RprobeLmax;
 int Ngrid_probe, Ngrid_probeS,Ngrid_probeL; 
 struct CHAR3DMAP MapP;
 struct CHAR3DMAP MapPS;
 struct CHAR3DMAP MapPL;
 float margin_width, voxel_volume;
 char iprbSpdbfile[MAX_FILE_NAME]; 

 /* for multiscale probe sphere */
 float Bin_of_RprobeL;
 struct MULSC_SHELL_PRB  MSprb;
 float  Wshell;   /* Width of shell. (default 2.0*RprobeS) */
 
 /* for multiscale pocket clustering  */
 struct GRID_CLUSTER Gclus;
 int    threCluster_char_min;
 int    threCluster_char_max;
 int    RankKeepCluster;
 char   tFactorType; 
 char   ClusPocketType;

 /* for protein molecules  */
 struct CHAR3DMAP MapX;  /* 3D map for protein shape */
 struct CHAR3DMAP MapY;  /* Additional 3D map for multi-scale probes    */
 struct CHAR3DMAP MapZ;  /* Additional 3D map for clustering of pockets */
 struct CHAR3DMAP ClusterMapHead;  /* Head of linked list of Cluster (small) maps  */
 char ipdbfile[MAX_FILE_NAME];
 char  iciffile[MAX_FILE_NAME], assembly_id[MAX_FILE_NAME];
 int   MaxAtomForAll;
 struct ATOM ProAtomHead;
 struct RESIDUE ProResHead;
 float MinP[3],MaxP[3];  /* Min and Maximum XYZ Values for ordinates of protein */
 int   MinN[3],MaxN[3];  /* Min and Maximum Grid Int   for foreground of the map */
 float MinL[3],MaxL[3];  /* Min and Maximum XYZ Values for ordinates of ligand  */
 char  ChainID_Str[32], asym_id_Str[32];
 int   Ngrid_VdW,Ngrid;
 char  opdbfile[MAX_FILE_NAME],oresfile[MAX_FILE_NAME],orpdbfile[MAX_FILE_NAME],oallpdbfile[MAX_FILE_NAME];
 char  ReadAtmHetType; 
 struct FLOAT2DMAP DisMapPro;
 char   pdbidPro[8],pdbidProA[8],pdbidProB[8];
 char  MODELtype;
 char   select_model_num_str[64];
 double M[3],CovM[3][3],PCaxis[3][3],PCvar[3];

 /* for ligand molecule  */
 char  iligpdbfile[MAX_FILE_NAME],oligpdbfile[MAX_FILE_NAME], iligciffile[MAX_FILE_NAME];
 struct ATOM LigAtomHead;
 char   vdwRadiusType, iligRadiusType;
 float  Runified;
 char   pdbidLig[8];

 /* for 3D density map  */
 struct VOXEL Voxel;
 char imapfile[MAX_FILE_NAME];
 float CutoffDensity;

 /* for spherical probes  */
 struct ATOM ProbeAtomHead;
 struct RESIDUE ProbeClusHead;
 char   GenSphePrb,ospheprbfile[MAX_FILE_NAME],ospheprbfile_dock[MAX_FILE_NAME];
 float  Dclus_SpheProbe,Dtree_SpheProbe;

 /** For detected pockets **/
 int  Ngrid_pocket;
 int  Ncluster, Nbit_cluster[MAX_UNSIGNED_CHAR];
 float  Min_Pocket_ProbeS; 
 int    Min_Pocket_GridNum;
 char  cluster_success; 
 char  OrigPosString[128], NgridString[128];
 int   Nword, Wsta[10],Wend[10];
 char  buffstr[128];

 /** For Interface Between Two Chains  **/
 char ipdbfileA[MAX_FILE_NAME],ipdbfileB[MAX_FILE_NAME];
 struct ATOM ProAtomHeadA;
 struct ATOM ProAtomHeadB;
 int NatomPdbA, NatomPdbB;
 float MinA[3],MaxA[3];  /* Min and Maximum Values for XYZ cordinates of protein A */ 
 float MinB[3],MaxB[3];  /* Min and Maximum Values for XYZ cordinates of protein B */
 int   Ngrid_VdW_A,Ngrid_VdW_B;
 int   Ngrid_gapAB;
 char ChainID_StrA[32],ChainID_StrB[32];
 int Natom_interA, Natom_interB;

 /*** for Ray-based LIGSITE-like inaccessibiliy */
 int  Nraydir;      /* Number of ray direction */
 char TermRayType; /* Ray Teminal Type. Count 'B'oth-end(PSP),Count 'O'ne-end */

 /*** MODE=ConSphe (Contacting Spheres) ***/
 int Nconspheres;
 struct ATOM ConSpheAtomHead;
 char  ConSphe_RGBTstr[MAX_FILE_NAME];
 char icen_pdbfile[MAX_FILE_NAME], isphe_pdbfile[MAX_FILE_NAME],osphe_pdbfile[MAX_FILE_NAME],osphe_wrlfile[MAX_FILE_NAME];
 int PCnum_ConSphe;
 float Rmin_ConSphe; 
 char ExcludeEscapeConSphe;
 Nconspheres = 10;
 sprintf(ConSphe_RGBTstr,"0:1:0:0.5");
 icen_pdbfile[0] = '\0';
 isphe_pdbfile[0] = '\0';
 sprintf(osphe_pdbfile,"consphe.pdb");
 sprintf(osphe_wrlfile,"consphe.wrl");
 PCnum_ConSphe = 0;
 Rmin_ConSphe = -1.0;
 ExcludeEscapeConSphe = 'F';

 /*** Other variables ***/
 char igridfile[MAX_FILE_NAME], igridfileA[MAX_FILE_NAME], igridfileB[MAX_FILE_NAME],ogridfile[MAX_FILE_NAME];
 int  thre_min_gridA, thre_max_gridA;
 char chain_select_grid,chain_select_gridA;
 char DetailOptHelp; 
 float Vvoxel; 
 int Nshift[3]; 
 int  NP,NL,NPL;
 float Rec,Pre,Fme,Tani;

 struct OPTION_VALUE *ov;

 /*** Set Default Value ***/
 sprintf(PAR.START_DATE,"%s",Get_Date_String());
 sprintf(PAR.HOSTNAME,"%s",string_getenv("HOSTNAME"));
 PAR.START_TIME_SEC = Get_Time_in_Second_by_Double();
 PAR.COMP_TIME_SEC          = 0.0;
 PAR.COMP_TIME_SEC_PROGRESS = 0.0;
 PAR.RecursiveQueue = 'Q';

 PAR.MAX_MEMORY_MEGA = 2000.0;

 PAR.N_ATOM    = 0;
 PAR.N_RESIDUE = 0;
 PAR.N_ATOM    = 0;
 PAR.N_CHAIN   = 0;

 PAR.Change_N_CA_C_O_Hetatom_to_Atom = 'T';
 GenSphePrb = 'F';
 ospheprbfile[0] = ospheprbfile_dock[0] = '\0';

 Nraydir = 7;
 TermRayType = 'B';

 thre_min_gridA = 1; thre_max_gridA = 254;
 chain_select_grid = chain_select_gridA = '-';
 PAR.SubType[0] = '\0';
 RprobeL = 6.0; RprobeS = 1.87;
 RprobeLmin = 2.0; RprobeLmax = 10.0;
 Bin_of_RprobeL  = 0.5;
 Wshell = RprobeS * 2.0;

 /* MODE = 'P'; */
 sprintf(MODE,"P");
 sprintf(ChainID_Str, "-");
 sprintf(ChainID_StrA,"-");
 sprintf(ChainID_StrB,"-");
 asym_id_Str[0] = '\0';
 ReadAtmHetType = 'A'; 
 oVdWpdbfile[0] = oVdWmapfile[0] = '\0';
 PAR.grid_width = 0.8;
 Min_Pocket_ProbeS = 0.0; Min_Pocket_GridNum = 0; 
 SurfacePocket = 'F';
 ClusPocketType = 'T';
 opockpdbfile[0] = opockmapfile[0] = oMoVpdbfile[0] = oMoVmapfile[0] = '\0';
 OrigPosString[0] = NgridString[0] = '\0';
 iprbSpdbfile[0] = '\0'; 
 MSprb.Nradius = 0;
 ipdbfile[0] = iciffile[0] = ipdbfileA[0] = ipdbfileB[0] = imapfile[0] = '\0';
 assembly_id[0] = '\0';
 MaxAtomForAll = -1;
 oprbSfile[0] = oprbLfile[0] = '\0'; 
 opdbfile[0] = oresfile[0] = orpdbfile[0] =  omapfile[0] = oallpdbfile[0] = '\0';
 ComStrHead.next = NULL; ComStrHead.prev = NULL; ComStrHead.line[0] = '\0';
 iligpdbfile[0] = iligciffile[0] = oligpdbfile[0] = omprobfile[0] = '\0';
 vdwRadiusType  = 'C';
 iligRadiusType = 'C';
 Runified = 2.5;
 igridfile[0] = igridfileA[0] = igridfileB[0] = ogridfile[0] = '\0';
 PAR.radfile[0] = '\0';
 /* PAR.NeighborNum = 6; */
 PAR.NeighborNum = 26; 
 tFactorType = 'P'; 
 threCluster_char_min  = 0; 
 threCluster_char_max  = 254; 
 RankKeepCluster = -1;
 PAR.Nmax_cluster = 5;
 Dclus_SpheProbe = 0.5;
 Dtree_SpheProbe = -1.0;
 MODELtype = 'M';
 CutoffDensity = 0.0;
 sprintf(select_model_num_str,"-1");

 if((argc>1) && ((strncmp(argv[1],"-h",2)==0) || (strncmp(argv[1],"-H",2)==0)) ) {
   DetailOptHelp = 1; 
 }
 else{
   DetailOptHelp = 0;
 }
 
 if ((argc<2)||(DetailOptHelp==1)){
  printf("ghecom [input_pdbfile] <options>\n");
  printf(" for surface/volume/pocket/cavity in 3D grid represention.\n");
  printf(" coded by T.Kawabata. LastModified:%s\n",LAST_MOD_DATE);
  printf(" <'P':standard usage for simple pocket detection>\n");
  printf("   ghecom P -ipdb [input_receptor_file] -opocpdb [pocket_in_pdb] -opdb [out_receptor_file_with_ClusterNum]\n");
  printf(" <'M':standard usage for multi-scale pocket detection>\n");
  printf("   ghecom M -ipdb [input_receptor_file] -opocpdb [pocket_in_pdb] -opdb [out_receptor_file_with_Rinaccess]\n");
  printf(" <options>\n");
  printf(" -h   : Show all of the options \n");
  printf(" -M  or the first argument without hyphen.:\n");
  printf("       'D'ilation 'E'rosion, 'C'losing(molecular surface), 'O'pening.\n");
  printf("       'P'ocket(masuya_doi),'p'ocket(kawabata_go),'V':ca'V'ity, 'e'roded pocket.\n");
  printf("       'V':cavity(center-connected) 'CV':Connolly-void (probe-connected cavity) \n");
  printf("       'CP':Cave-pocket 'CPKG:Cave-pocket(KG-like definition) 'CorP':CavityOrPocket\n");
  printf("       'M'ultiscale_closing/pocket,'I'nterface_pocket_bwn_two_chains.\n");
  printf("       'G'rid_comparison_binary 'g'rid_comparison_mutiscale.\n");
  printf("       'ConSphe' Contacting-spheres\n");
  printf("       'GcmpS' Grid-vs-Sphere comparison\n");
  printf("       'GcmpL' Grid-vs-Ligand comparison\n");
  printf("       'GcmpG' Grid-vs-Ligand comparison\n");
  printf("       'GtoV'  Grid-to-Volume_map conversion\n");
  printf("       'Gsurf' Surfacizaton of Grid\n");
  printf("       'L'igand-grid comparison (-ilg and -igA are required)[%s]\n",MODE);
  printf("       'X' do nothing. just read pdb/mmCIF file\n");
  printf("<options for input PDB file>\n");
  printf(" -ipdb : input PDB   file [%s]\n",ipdbfile); 
  printf(" -icif : input mmCIF file [%s]\n",iciffile); 
  printf(" -ch       : ChainID for target pdbfile [%s]\n",ChainID_Str); 
  printf(" -asym_id  : asym_id for target pdbfile (only for -icif) [%s]\n",asym_id_Str); 
  printf(" -atmhet  : Read only ATOM 'A', only HETATM 'H', both 'B' [%c]\n",ReadAtmHetType); 
  printf(" -hetpep2atm  : Change HETATM residues with N_CA_C_O to ATOM ('T' or 'F') [%c]\n", PAR.Change_N_CA_C_O_Hetatom_to_Atom);
  printf(" -assembly: assembly_id for mmCIF file (only for -icif) [%s]\n",assembly_id);
  printf(" -maxatm : maximum allowed number of atoms for '-atmsel A'. If over '-maxatm', then change '-atmsel R'.[%d]\n",MaxAtomForAll);
  printf("<options for grid pocket calculation>\n");
  printf(" -gw  : Grid width [%f]\n",PAR.grid_width); 
  printf(" -rs  : Radius for small probe spheres [%f]\n",RprobeS); 
  printf(" -rl  : Radius for large probe spheres [%f]\n",RprobeL); 
  printf("<options for input 3D density map>\n");
  printf(" -imap   : Input 3D density map (*.map|*.mrc) [%s]\n",imapfile); 
  printf(" -cutoff : cutoff density [%f]\n",CutoffDensity); 
  printf("<other options>\n");
  printf(" -opocpdb: Output file for Grid Points of Pocket in PDB[%s]\n",opockpdbfile); 
  printf(" -opocmap: Output file for Pocket in 3D density map (*.map) [%s]\n",opockmapfile); 
  printf(" -oall: Output PDB file contains '-opoc','-opdb', '-olg', and '-oprb' [%s]\n",oallpdbfile); 
  printf(" -opdb: Output receptor PDB file with Cluster Number in tFactor[%s]\n",opdbfile);
  printf(" -model :'S':read only single model (for NMR). 'M':read multiple models (for biological unit) [%c]\n",MODELtype);
  printf(" -model_num : select model number  '-1':no-select '1':model=1 'le2' model<=2[%s]\n",select_model_num_str);
  printf(" <options only for MODE 'M'>\n"); 
  printf(" -rli : Radius for min_large probe spheres [%f]\n",RprobeLmin); 
  printf(" -rlx : Radius for max_large probe spheres [%f]\n",RprobeLmax); 
  printf(" -br  : bin of large probe radius for MODE 'M' [%f]\n",Bin_of_RprobeL); 
  printf(" -opdb: Output receptor PDB file with Rinacess/Pocketness in tFactor[%s]\n",opdbfile); 
  printf(" -ores: Output Rediue-based property file with calculated Rinaceess/Pocketness[%s]\n",oresfile); 
  printf(" -clus: clustering pockets ('T' or 'F') [%c]\n",ClusPocketType);
  printf(" -iligpdb : Input ligand PDB file for Rinaccess calculation  (only for MODE 'M') [%s]\n",iligpdbfile); 
  printf(" -oligpdb : Output ligand PDB file with calcualted Rinaccess (only for MODE 'M') [%s]\n",oligpdbfile); 
  printf(" -oprb: Output PDB file for 3-contacting spherical probes [%s]\n",ospheprbfile);
  printf(" -oprd: Output PDB file for 3-contacting spherical probes in DOCK sphere format [%s]\n",ospheprbfile_dock);
  printf(" <options only for MODE 'ConSphe'>\n"); 
  printf(" -ipdb    :Input pdbfile [%s]\n",ipdbfile); 
  printf(" -Nconsphe:number of contacting spheres [%d]\n",Nconspheres); 
  printf(" -icenpdb  :Input  center positions file in PDB [%s]\n",icen_pdbfile); 
  printf(" -osphepdb :Output spheres in PDB [%s]\n",osphe_pdbfile); 
  printf(" -osphewrl :Output spheres in VRML [%s]\n",osphe_wrlfile); 
  printf(" -csRGBT   :RGBT for -osphewrl [%s]\n",ConSphe_RGBTstr); 
  printf(" -PCconsphe: PCA axis number for generating spheres 0 or 1 or 2. if -1, put one sphere on the center. [%d]\n",PCnum_ConSphe);
  printf(" -csRmin   : minimum Radius for contact spheres [%f]\n",Rmin_ConSphe);
  printf(" -Xescsphe : Exclude escaping contacting spheres ('T' or 'F') [%c]\n",ExcludeEscapeConSphe);
  printf(" <options only for MODE 'GcmpS'>\n"); 
  printf(" -igridpdb  :Input grid pdb file [%s]\n",igridfile);
  printf(" -isphepdb  :Input sphere file in PDB [%s]\n",icen_pdbfile); 
  printf(" -chaingrid : Cluster_num(chainID) select for grid. '1':1,'2':1,2,'3':1,2,3 '-':all. [%c]\n",chain_select_grid);
  printf(" <options only for MODE 'GcmpL'>\n"); 
  printf(" -igridpdb  :Input grid pdb file [%s]\n",igridfile);
  printf(" -iligpdb   :Input ligand atom file in PDB [%s]\n",iligpdbfile); 
  printf(" -iligcif   :Input ligand atom file in mmCIF [%s]\n",iligciffile); 
  printf(" -assembly  : assembly_id for mmCIF file (-iligcif) [%s]\n",assembly_id);
  printf(" -chaingrid : Cluster_num(chainID) select for grid. '1':1,'2':1,2,'3':1,2,3 '-':all. [%c]\n",chain_select_grid);
  printf(" <options only for MODE 'GcmpG'>\n"); 
  printf(" -igridpdbA  :Input grid pdb fileA [%s]\n",igridfileA);
  printf(" -igridpdbB  :Input grid pdb fileA [%s]\n",igridfileB);
  printf(" <options only for MODE 'GtoV'>\n"); 
  printf(" -igridpdb  :Input grid pdb file [%s]\n",igridfile);
  printf(" -omap      :Output map file [%s]\n",omapfile);
  printf(" -chaingrid : Cluster_num(chainID) select for grid. '1':1,'2':1,2,'3':1,2,3 '-':all. [%c]\n",chain_select_grid);
  printf(" -recque    : 'R'ecursive or 'Q'ueue for flood-fill or labeling [%c]\n",PAR.RecursiveQueue);
  printf(" <options only for MODE 'Gsurf'>\n"); 
  printf(" -ogridpdb  :Output grid pdb file [%s]\n",ogridfile);

 if (DetailOptHelp==1){ 
  printf(" <standard usage for MODE 'L'>\n");
  printf("    ghecom -M L -ilg [input_ligand_PDB_file] -igA [input_pocket_GRID_file]\n");
  printf(" <options>\n");
  printf(" -ovdwpdb : Output VdW protein 3D grid PDB file[%s]\n",oVdWpdbfile); 
  printf(" -ovdwmap : Output VdW protein 3D density map file[%s]\n",oVdWmapfile); 
  printf(" -omovpdb : Output Molecular Volume 3D grid in PDB file[%s]\n",oMoVpdbfile); 
  printf(" -omovmap : Output Molecular Volume 3D grid in 3D density mapfile[%s]\n",oMoVpdbfile); 
  printf(" -thc : max thereshold char value for clustering for MODE 'M' and 'L'. 255 for everything[%d]\n",threCluster_char_max);
  printf(" -thci: min thereshold char value for clustering for MODE 'R'[%d]\n",threCluster_char_min);
  printf(" -mp  : Min small probe number of pocket cluster [%f]\n",Min_Pocket_ProbeS); 
  printf(" -mg  : Min grid size          of pocket cluster [%d]\n",Min_Pocket_GridNum); 
  printf(" -su  : surfacing  detected pockets ('T' or 'F') [%c]\n",SurfacePocket);
  printf(" -opos: origin point (OrigPos) position (x):(y):(z) [%s]\n",OrigPosString); 
  printf(" -ngri: number of grids (x):(y):(z) [%s]\n",NgridString); 
  printf(" -rvdw : vwRadius type for input PDB    file     'C'hothia1976 'B'ondi1964, 'O'ccupancy, 'U'nified [%c]\n",vdwRadiusType); 
  printf(" -rlig : vwRadius type for input lignad PDB file 'C'hothia1976 'B'ondi1964, 'O'ccupancy, 'U'nified [%c]\n",iligRadiusType); 
  printf(" -runi: Radius_unified for '-rvdw U' or '-rlig U' [%f]\n",Runified); 
  printf(" -nnei : Number of Neighbor pixel for flooding and labeling  (6,18,or 26) [%d]\n",PAR.NeighborNum); 
  printf(" -ops : output small probe grid file [%s]\n",oprbSfile);  
  printf(" -opl : output large probe grid file [%s]\n",oprbLfile);  
  printf(" <options only for MODE 'M'>\n"); 
  printf(" -wsh : width of accessible shell (default is 2.0*RprobeS)[%f]\n",Wshell); 
  printf(" -rkc : rank keep for large cluster for MODE 'M' and 'R' [%d]\n",RankKeepCluster);
  printf(" -nr  : Number of radius for MODE 'M' [%d]\n",MSprb.Nradius); 
  printf(" -omp : Output multiscale probe PDB file [%s]\n",omprobfile); 
  printf(" -tfac: tFactor value for opdb 'R'inacc,'P'ocketness 'Z'ero value[%c]\n",tFactorType); 
  printf(" -omap: Output CCP4 style 3D density map[%s]\n",omapfile); 
  printf(" -orpd: Output PDB file with Rediue-based Rinaceess/Pocketness[%s]\n",orpdbfile); 
  printf(" <options for spherical probes (MODE 'M')>\n"); 
  printf(" -sprb: Generate 3-contacting spherical probe (T or F) [%c]\n",GenSphePrb);
  printf(" -dprb: dis thre for reducing spherical probes by single linkage clustering [%f]\n",Dclus_SpheProbe);
  printf(" -dwhc: dis thre for weighted hierarchical clustering of spherical probes[%f]\n",Dtree_SpheProbe);
  printf(" <options only for MODE 'R'>\n"); 
  printf(" ## protein: '-rs'-closed vdW, target: '-rl'-dilated '-rs'-closed vdW.\n"); 
  printf(" -nry : Number of Ray directions (7 or 13) [%d]\n",Nraydir); 
  printf(" -try : Terminal Type of ray 'B'oth(PSP) or 'O'ne.[%c]\n",TermRayType); 
  printf(" <options only for MODE 'G'>\n"); 
  printf(" -igA : Input grid pdb file A [%s]\n",igridfileA);
  printf(" -igB : Input grid pdb file B [%s]\n",igridfileB);
  printf(" <options only for MODE 'I'>\n"); 
  printf(" -ipdbA : Input pdb file A [%s]\n",ipdbfileA);
  printf(" -ipdbB : Input pdb file B [%s]\n",ipdbfileB);
  printf(" -chA  : ChainID for pdb file A [%s]\n",ChainID_StrA); 
  printf(" -chB  : ChainID for pdb file B [%s]\n",ChainID_StrB); 
  printf(" <options for MODE 'L'>\n"); 
  printf(" -tgAx: max thre value for grid A for mode 'L' [%d]\n",thre_min_gridA);
  printf(" -tgAi: min thre value for grid A for mode 'L' [%d]\n",thre_max_gridA);
  printf(" -cgA : Cluster_num(chainID) select for gridA. '1':1,'2':1,2,'3':1,2,3 '-':all. [%c]\n",chain_select_gridA);
  printf(" <other options>\n"); 
  printf(" -isp : Input small Probe PDB file[%s]\n",iprbSpdbfile); 
  printf(" -sub : SubType for various purpose [%s]\n",PAR.SubType); 
  printf(" -max_memory  : maximum memory size in mega byte [%.1f]\n",PAR.MAX_MEMORY_MEGA);

  } 
  exit(1);
 }


 /****** READ ARGUMENTS ***********/
 PAR.OptValHead.next = NULL;
 Read_Options_From_Arguments(argc,argv,PAR.COMMAND,&(PAR.OptValHead));

 if (argv[1][0]!='-'){ 
   sprintf(MODE,"%s",argv[1]);
 }


 ov = &(PAR.OptValHead);
 while (ov->next != NULL){
   ov = ov->next;
        if (strcmp(ov->opt,"ch")==0){ sprintf(ChainID_Str,"%s",ov->val);}
   else if (strcmp(ov->opt,"asym_id")==0){ sprintf(asym_id_Str,"%s",ov->val);}
   else if (strcmp(ov->opt,"opos")==0){ sprintf(OrigPosString,"%s",ov->val);}
   else if (strcmp(ov->opt,"ngri")==0){ sprintf(NgridString,"%s",ov->val);}
   else if (strcmp(ov->opt,"opdb")==0){ sprintf(opdbfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"oall")==0){ sprintf(oallpdbfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"ores")==0){ sprintf(oresfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"orpd")==0){ sprintf(orpdbfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"tfac")==0){ tFactorType = ov->val[0];}
   else if (strcmp(ov->opt,"gw")==0){ PAR.grid_width = atof(ov->val); }
   else if (strcmp(ov->opt,"gw")==0){ PAR.grid_width = atof(ov->val); }
   else if (strcmp(ov->opt,"wsh")==0){ Wshell = atof(ov->val);}
   else if (strcmp(ov->opt,"ipdb")==0){ sprintf(ipdbfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"icif")==0){ sprintf(iciffile,"%s",ov->val);}
   else if (strcmp(ov->opt,"assembly")==0){ sprintf(assembly_id,"%s",ov->val);}
   else if (strcmp(ov->opt,"ipdbA")==0){ sprintf(ipdbfileA,"%s",ov->val);}
   else if (strcmp(ov->opt,"ipdbB")==0){ sprintf(ipdbfileB,"%s",ov->val);}
   else if (strcmp(ov->opt,"icenpdb")==0){ sprintf(icen_pdbfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"osphepdb")==0){ sprintf(osphe_pdbfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"osphewrl")==0){ sprintf(osphe_wrlfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"csRGBT")==0){ sprintf(ConSphe_RGBTstr,"%s",ov->val);}
   else if (strcmp(ov->opt,"PCconsphe")==0){ PCnum_ConSphe = atoi(ov->val);}
   else if (strcmp(ov->opt,"csRmin")==0){ Rmin_ConSphe = atof(ov->val);}
   else if (strcmp(ov->opt,"imap")==0){ sprintf(imapfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"chA")==0){ sprintf(ChainID_StrA,"%s",ov->val);}
   else if (strcmp(ov->opt,"chB")==0){ sprintf(ChainID_StrB,"%s",ov->val);}
   else if (strcmp(ov->opt,"igA")==0) { sprintf(igridfileA,"%s",ov->val);}
   else if (strcmp(ov->opt,"igB")==0) { sprintf(igridfileB,"%s",ov->val);}
   else if (strcmp(ov->opt,"igridpdbA")==0) { sprintf(igridfileA,"%s",ov->val);}
   else if (strcmp(ov->opt,"igridpdbB")==0) { sprintf(igridfileB,"%s",ov->val);}
   else if (strcmp(ov->opt,"igridpdb")==0) { sprintf(igridfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"ogridpdb")==0) { sprintf(ogridfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"isphepdb")==0) { sprintf(isphe_pdbfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"opo")==0) { sprintf(opockpdbfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"opocpdb")==0) { sprintf(opockpdbfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"opocmap")==0) { sprintf(opockmapfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"ovdwpdb")==0) { sprintf(oVdWpdbfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"ovdwmap")==0) { sprintf(oVdWmapfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"omovpdb")==0) { sprintf(oMoVpdbfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"omovmap")==0) { sprintf(oMoVmapfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"iligpdb")==0) { sprintf(iligpdbfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"iligcif")==0) { sprintf(iligciffile,"%s",ov->val);}
   else if (strcmp(ov->opt,"oligpdb")==0) { sprintf(oligpdbfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"rvdw")==0) { vdwRadiusType  = ov->val[0];}
   else if (strcmp(ov->opt,"rlig")==0) { iligRadiusType = ov->val[0];}
   else if (strcmp(ov->opt,"runi")==0) { Runified  = atof(ov->val);}
   else if (strcmp(ov->opt,"tgAx")==0) { thre_max_gridA = atoi(ov->val);}
   else if (strcmp(ov->opt,"tgAi")==0) { thre_min_gridA = atoi(ov->val);}
   else if (strcmp(ov->opt,"cgA")==0) { chain_select_gridA = ov->val[0];}
   else if (strcmp(ov->opt,"chaingrid")==0) { chain_select_grid = ov->val[0];}
   else if (strcmp(ov->opt,"omap")==0) {sprintf(omapfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"omp")==0) { sprintf(omprobfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"rs")==0)  { RprobeS = atof(ov->val); }
   else if (strcmp(ov->opt,"rl")==0)  { RprobeL = atof(ov->val); }
   else if (strcmp(ov->opt,"rli")==0)  { RprobeLmin = atof(ov->val); }
   else if (strcmp(ov->opt,"rlx")==0)  { RprobeLmax = atof(ov->val); }
   else if (strcmp(ov->opt,"nry")==0)  { Nraydir = atoi(ov->val);}
   else if (strcmp(ov->opt,"try")==0)  { TermRayType = atoi(ov->val);}
   else if (strcmp(ov->opt,"M")==0)   { sprintf(MODE,"%s",ov->val); }
   else if (strcmp(ov->opt,"mp")==0) { Min_Pocket_ProbeS  = atof(ov->val);}
   else if (strcmp(ov->opt,"cl")==0)   { ClusPocketType = ov->val[0];}
   else if (strcmp(ov->opt,"clus")==0) { ClusPocketType = ov->val[0];}
   else if (strcmp(ov->opt,"su")==0) { SurfacePocket = ov->val[0];}
   else if (strcmp(ov->opt,"nr")==0) { MSprb.Nradius = atoi(ov->val);}
   else if (strcmp(ov->opt,"br")==0) { Bin_of_RprobeL    = atof(ov->val);}
   else if (strcmp(ov->opt,"atmhet")==0) { ReadAtmHetType = ov->val[0];}
   else if (strcmp(ov->opt,"hetpep2atm")==0) { PAR.Change_N_CA_C_O_Hetatom_to_Atom = ov->val[0]; }
   else if (strcmp(ov->opt,"nnei")==0) { PAR.NeighborNum = atoi(ov->val);}
   else if (strcmp(ov->opt,"isp")==0) { sprintf(iprbSpdbfile,"%s",ov->val);}
   else if (strcmp(ov->opt,"sub")==0) { sprintf(PAR.SubType,"%s",ov->val); }
   else if (strcmp(ov->opt,"thc")==0) { threCluster_char_max = atoi(ov->val); }
   else if (strcmp(ov->opt,"thci")==0) { threCluster_char_min = atoi(ov->val); }
   else if (strcmp(ov->opt,"rkc")==0) {; RankKeepCluster = atoi(ov->val); }
   else if (strcmp(ov->opt,"ops")==0) { sprintf(oprbSfile,"%s",ov->val); }
   else if (strcmp(ov->opt,"opl")==0) { sprintf(oprbLfile,"%s",ov->val); }
   else if (strcmp(ov->opt,"sprb")==0) { GenSphePrb= ov->val[0];}
   else if (strcmp(ov->opt,"oprb")==0) { sprintf(ospheprbfile,"%s",ov->val); GenSphePrb='T';}
   else if (strcmp(ov->opt,"oprd")==0) { sprintf(ospheprbfile_dock,"%s",ov->val); GenSphePrb = 'T';}
   else if (strcmp(ov->opt,"dprb")==0) {  Dclus_SpheProbe = atof(ov->val); }
   else if (strcmp(ov->opt,"dwhc")==0) { Dtree_SpheProbe = atof(ov->val); }
   else if (strcmp(ov->opt,"model")==0) {  MODELtype  = ov->val[0]; }
   else if (strcmp(ov->opt,"model_num")==0) { sprintf(select_model_num_str,"%s",ov->val); }
   else if (strcmp(ov->opt,"cutoff")==0) { CutoffDensity = atof(ov->val); }
   else if (strcmp(ov->opt,"Nconsphe")==0) { Nconspheres = atoi(ov->val); }
   else if (strcmp(ov->opt,"recque")==0) { PAR.RecursiveQueue = ov->val[0]; }
   else if (strcmp(ov->opt,"Xescsphe")==0) { ExcludeEscapeConSphe = ov->val[0]; }
   else if (strcmp(ov->opt,"max_memory")==0) { PAR.MAX_MEMORY_MEGA = atof(ov->val); }
   else { printf("#ERROR:Can't understand option %s\n",ov->opt); exit(1);}
 }

 /*** [0] General procedures commonly used for all the MODE ***/

 MapX.grid_width = MapY.grid_width = MapZ.grid_width  = PAR.grid_width;
 MapP.grid_width = MapPS.grid_width = MapPL.grid_width =  PAR.grid_width;
 MapX.malloced = MapY.malloced = MapZ.malloced = MapP.malloced = MapPS.malloced = MapPL.malloced =  0; 
 
 printf("#COMMAND %s\n",PAR.COMMAND);

 if (  (opockpdbfile[0]=='\0')
     &&( (strcmp(MODE,"P")==0)||(strcmp(MODE,"p")==0)||(strcmp(MODE,"F")==0)||
         (strcmp(MODE,"r")==0)||(strcmp(MODE,"s")==0)||(strcmp(MODE,"u")==0)||(strcmp(MODE,"V")==0)||(strncmp(MODE,"CP",2)==0)) 
   ){ sprintf(opockpdbfile,"%s.pock.pdb",ipdbfile);}
 

 ProAtomHead.next = NULL;
 ProAtomHeadA.next = NULL;
 ProAtomHeadB.next = NULL;

 /** [1] Setup margin_width **/
  margin_width = 2.0*RprobeS; 
  if ((strcmp(MODE,"D")==0)||(strcmp(MODE,"E")==0)||(strcmp(MODE,"C")==0)|| 
      (strcmp(MODE,"v")==0)||(strcmp(MODE,"O")==0)||(strcmp(MODE,"o")==0)){
       margin_width = 2.0*RprobeS; 
  }

  else if ((strcmp(MODE,"P")==0)||(strcmp(MODE,"p")==0)||(strcmp(MODE,"m")==0)||(strcmp(MODE,"V")==0)||
           (strcmp(MODE,"e")==0)|| (strncmp(MODE,"CP",2)==0)||(strncmp(MODE,"CPKG",4)==0)||(strncmp(MODE,"CorP",4)==0)||(strcmp(MODE,"X")==0)){
       margin_width = 2.0*RprobeL; 
  }
  else if (strcmp(MODE,"M")==0){ margin_width = 2.0*RprobeLmax;}
  
  margin_width += MapX.grid_width*2;


 /*** [2-1] Read PDB file  ***/
 if ((ipdbfile[0]!='\0')||(iciffile[0]!='\0')){

   if (ipdbfile[0] != '\0'){
     Read_PDB_File(ipdbfile,&ProAtomHead,pdbidPro,'H',ChainID_Str,ReadAtmHetType,'F','F',MODELtype);
   }
   else if (iciffile[0] != '\0'){
     Read_mmCIF_File(iciffile,&ProAtomHead,'A',ReadAtmHetType,ChainID_Str,asym_id_Str,assembly_id,MaxAtomForAll);
   }

   if (pdbidPro[0]!='\0') {sprintf(PAR.pdbidPro,"%s",pdbidPro);} else {sprintf(PAR.pdbidPro,"0XXX");}

   NatomPdb = Number_Of_ATOMs(&ProAtomHead); 
   Assign_Radius_to_ATOMs(&ProAtomHead,vdwRadiusType,Runified);
   MinMaxXYZ_Of_Atoms(&ProAtomHead,MinP,MaxP); 
   printf("#MinP %f %f %f MaxP %f %f %f\n", MinP[0], MinP[1], MinP[2], MaxP[0], MaxP[1], MaxP[2]);
   Make_Residue_List(&ProAtomHead,&ProResHead);
   PAR.N_ATOM    = Number_Of_ATOMs(&ProAtomHead);
   PAR.N_RESIDUE = Number_Of_RESIDUEs(&ProResHead);
   PAR.N_MODEL   = Number_Of_MODELs(&ProAtomHead);
   PAR.N_CHAIN   = Number_Of_CHAINs(&ProAtomHead);
   printf("N_ATOM    %d\n",PAR.N_ATOM);
   printf("N_RESIDUE %d\n",PAR.N_RESIDUE);
   printf("N_MODEL   %d\n",PAR.N_MODEL);
   printf("N_CHAIN   %d\n",PAR.N_CHAIN);

   sprintf(comment,"INPUT_RECEPTOR_PDB_FILE:%s",ipdbfile);
   Add_String_to_STRING_NODE(&ComStrHead,comment); 
   sprintf(comment,"NATOM_OF_RECEPTOR:%d",NatomPdb);
   Add_String_to_STRING_NODE(&ComStrHead,comment); 
  
   if (opdbfile[0]!='\0'){ 
     Write_PDB_File(opdbfile,'w',&ProAtomHead,"",&ComStrHead);
   }

    /*** [2-1-1] Generate Spherical Probes  ***/
    if (GenSphePrb=='T'){
      Malloc_FLOAT2DMAP(&DisMapPro,Number_Of_ATOMs(&ProAtomHead)); 
      Cal_Distance_FLOAT2DMAP(&ProAtomHead,&DisMapPro);
      Make_Probe_Spheres_Three_Contact_NEIGHBOR_list(&ProbeAtomHead,&ProAtomHead,RprobeS,&DisMapPro);
      if (Dclus_SpheProbe>0.0) {
        Make_Single_Linkage_Cluster_Of_Atoms(&ProbeAtomHead,&ProbeClusHead,Dclus_SpheProbe);
        Extract_Rep_Atoms_in_Cluster_Center(&ProbeAtomHead,&ProbeClusHead);
       }

      /* Write_PDB_File("prbN.pdb",'w',&ProbeAtomHead,"",&ComStrHead); */
    }

    /*** [2-1-2] Malloc Huge 3D grids for the PDB molecule  ***/
 
   if ((OrigPosString[0] !='\0') && (NgridString[0] != '\0')){
     Split_to_Words(OrigPosString,':',&Nword,Wsta,Wend,10);
     Get_Part_Of_Line(buffstr,OrigPosString,Wsta[0],Wend[0]); MapX.OrigPos[0] = atof(buffstr);
     Get_Part_Of_Line(buffstr,OrigPosString,Wsta[1],Wend[1]); MapX.OrigPos[1] = atof(buffstr);
     Get_Part_Of_Line(buffstr,OrigPosString,Wsta[2],Wend[2]); MapX.OrigPos[2] = atof(buffstr);
  
     Split_to_Words(NgridString,':',&Nword,Wsta,Wend,10);
     Get_Part_Of_Line(buffstr,NgridString,Wsta[0],Wend[0]); MapX.N[0] = atoi(buffstr);
     Get_Part_Of_Line(buffstr,NgridString,Wsta[1],Wend[1]); MapX.N[1] = atoi(buffstr);
     Get_Part_Of_Line(buffstr,NgridString,Wsta[2],Wend[2]); MapX.N[2] = atoi(buffstr);
    }
 
   else{
    Decide_Num_Grid_3D(MapX.N,MapX.OrigPos,MapX.grid_width, margin_width, MinP,MaxP, '-');
   }


   Malloc_CHAR3DMAP(&MapX, MapX.N[0], MapX.N[1], MapX.N[2]);
 
   Ngrid_VdW = Set_Bit_3DMAP_Inside_Atoms(&MapX,&ProAtomHead,1,1.0);
   Ngrid_VdW = Count_Specific_3DMAP(&MapX,1);
   printf("#Ngrid_VdW %d\n",Ngrid_VdW);
   PAR.VOL_VDW = Ngrid_VdW * MapX.grid_width * MapX.grid_width * MapX.grid_width;
   printf("VOL_VDW %f\n",PAR.VOL_VDW);
   sprintf(comment,"grid_width %6.3f A VdW volume: %d grids, %.2f AAA",
          MapX.grid_width,Ngrid_VdW,Ngrid_VdW*MapX.grid_width*MapX.grid_width*MapX.grid_width);
   Add_String_to_STRING_NODE(&ComStrHead,comment); 
 }

 /*** [2-2] Read 3D density map  ***/
 else if (imapfile[0]!='\0'){
   Read_MapCCP4_File_VOXEL(&Voxel,imapfile,'-',-1.0, 'F', 'F',-1, 0.0, 0.0);
   Malloc_CHAR3DMAP(&MapX, Voxel.N[0], Voxel.N[1], Voxel.N[2]);
   Initialize_CHAR3DMAP(&MapX);
   assign_CHAR3DMAP_from_VOXEL(&MapX,&Voxel,CutoffDensity,1);
   Foreground_MinMax_CHAR3DMAPL(&MapX,1,MinN,MaxN);
   printf("#margin_width %f margin_voxel %.1f\n",margin_width, margin_width/MapX.grid_width);
   resize_CHAR3DMAP_from_MinMax_and_MarginWidth(&MapX,MinN,MaxN,margin_width,MapY.N,MapY.OrigPos);

   if ((MapX.N[0] != MapY.N[0])|| (MapX.N[1] != MapY.N[1])|| (MapX.N[2] != MapY.N[2])){
      MapY.grid_width = MapX.grid_width;
      Malloc_CHAR3DMAP(&MapY, MapY.N[0], MapY.N[1], MapY.N[2]);
      Initialize_CHAR3DMAP(&MapY);
      copy_CHAR3DMAP_Same_GRID_and_Diff_OrigPos(&MapY,&MapX);
      Free_CHAR3DMAP(&MapX);
      Malloc_CHAR3DMAP(&MapX, MapY.N[0], MapY.N[1], MapY.N[2]);
      Copy_CHAR3DMAP(&MapX, &MapY);
      Free_CHAR3DMAP(&MapY);
   }
 }

 if (oVdWpdbfile[0]!='\0'){ 
   Output_CHAR3DMAP_PDB_format(oVdWpdbfile,'w',&MapX,1,"Protein VdW Volume",&ComStrHead); 
 }

 if (oVdWmapfile[0]!='\0'){ 
   Write_MapCCP4_File_CHAR3DMAP(oVdWmapfile,&MapX,'D');
 }

 /*************************/ 
 /** MODE 'D' : Dilation **/
 /*************************/ 

  if (strcmp(MODE,"D")==0){ 
    Ngrid_probe = Make_Sphere_Probe_Map(&MapP,RprobeS,MapP.grid_width);
    /* Output_CHAR3DMAP_PDB_format("probe.gdb",'w',&MapP,1,"Probe Sphere",&ComStrHead); */
    Dilation(&MapX,&MapP,1,4,'O');
    if (opockpdbfile[0] != '\0'){
      Output_CHAR3DMAP_PDB_format(opockpdbfile,'w',&MapX,1,"original(1) and Dilated (4) protein",&ComStrHead); 
    }
  }

 /************************/ 
 /** MODE 'E' : Erosion **/
 /************************/ 
 
 if (strcmp(MODE,"E")==0){
    Ngrid_probe = Make_Sphere_Probe_Map(&MapP,RprobeS,MapP.grid_width);
    /* Output_CHAR3DMAP_PDB_format("probe.gdb",'w',&MapP,1,"Probe Sphere",&ComStrHead); */
    Erosion(&MapX,&MapP,1,4);
    if (opockpdbfile[0] != '\0'){
      Output_CHAR3DMAP_PDB_format(opockpdbfile,'w',&MapX,1,"original(1) and Erosion(4)  protein",&ComStrHead); 
    }
 }
 
 /***************************************************/
 /** MODE 'C' : Closing or moleculer Surface **/
 /***************************************************/

 if (strcmp(MODE,"C")==0){
    Ngrid_probe = Make_Sphere_Probe_Map(&MapP,RprobeS,MapP.grid_width);
    Output_CHAR3DMAP_PDB_format("probe.gdb",'w',&MapP,1,"Probe Sphere",&ComStrHead);
    
    Dilation(&MapX,&MapP,1,4,'O');
    
    Ngrid = Erosion(&MapX,&MapP,4,8);

    Output_CHAR3DMAP_PDB_format("clos.gdb",'w',&MapX,1,"Closing(molecular surace)",&ComStrHead); 
    
    sprintf(comment,"Molecular Volume %d grids, %.2f AAA",
       Ngrid,Ngrid*MapX.grid_width*MapX.grid_width*MapX.grid_width);
    Add_String_to_STRING_NODE(&ComStrHead,comment);

    Reset_Specific_3DMAP(&MapX,4);
    Extract_Specific_3DMAP(&MapX,8);
    Extract_Surface_3DMAP(&MapX,8,16);
    Extract_Specific_3DMAP(&MapX,16);
    if (oMoVpdbfile[0]!= '\0'){
      Output_CHAR3DMAP_PDB_format(oMoVpdbfile,'w',&MapX,1,"Molecular Surface with probe",&ComStrHead); 
    }
    if (oMoVmapfile[0]!= '\0'){
      Write_MapCCP4_File_CHAR3DMAP(oMoVmapfile,&MapX,'D');
    }
  }


 /*************************/
 /** MODE 'O' : Opening  **/
 /*************************/

 if (strcmp(MODE,"O")==0){
    Ngrid_probe = Make_Sphere_Probe_Map(&MapP,RprobeS,MapP.grid_width);
    Output_CHAR3DMAP_PDB_format("probe.gdb",'w',&MapP,1,"Probe Sphere",&ComStrHead);
  
    Ngrid = Erosion(&MapX,&MapP,1,2);
    sprintf(comment,"Eroded Volume %d grids, %.2f AAA",
       Ngrid,Ngrid*MapX.grid_width*MapX.grid_width*MapX.grid_width);
    Add_String_to_STRING_NODE(&ComStrHead,comment);
    Dilation(&MapX,&MapP,2,4,'O');
    Ngrid = Extract_Specific_3DMAP(&MapX,4);
    sprintf(comment,"Opening Volume %d grids, %.2f AAA",
       Ngrid,Ngrid*MapX.grid_width*MapX.grid_width*MapX.grid_width);
    Add_String_to_STRING_NODE(&ComStrHead,comment);
  }


 /*************************/
 /** MODE 'o' : Opening (one-step opening) **/
 /*************************/

 if (strcmp(MODE,"o")==0){
    Ngrid_probe = Make_Sphere_Probe_Map(&MapP,RprobeS,MapP.grid_width);
    Output_CHAR3DMAP_PDB_format("probe.gdb",'w',&MapP,1,"Probe Sphere",&ComStrHead);
  
    Ngrid = Opening(&MapX,&MapP,1,2);
    Extract_Specific_3DMAP(&MapX,2);
 }


 /*********************************************************************/
 /** MODE 'P': Pocket detection  (non-contact pocket)                **/
 /**  (definition of Masuya and Doi (1995))                          **/
 /*    PocketMD = ((X closing P) and noX)) opening S                 **/
 /** MODE 'p': Pocket detection  (contact pocket)                    **/
 /**  (definition of Kawbata and Go (2006))                          **/
 /*    PocketKG = ((X closing P) and noX and (X + 2S))) opening S    **/
 /*********************************************************************/

  if ((strcmp(MODE,"P")==0)||(strcmp(MODE,"p")==0)){

    Ngrid_pocket = 0;

    Ngrid_probeS = Make_Sphere_Probe_Map(&MapPS,RprobeS,MapPS.grid_width);

   if ((Min_Pocket_GridNum==0)&&(Min_Pocket_ProbeS>0.0)){
    Min_Pocket_GridNum =  (int)(Min_Pocket_ProbeS * Ngrid_probeS); 
    printf("Min_Pocket_ProbeS %f Min_Pocket_GridNum %d\n",Min_Pocket_ProbeS,Min_Pocket_GridNum);
   }
    if (oprbSfile[0] != '\0'){
      Output_CHAR3DMAP_PDB_format(oprbSfile,'w',&MapPS,1,"Small probe sphere",&ComStrHead);
    }
 
    Ngrid_probeL = Make_Sphere_Probe_Map(&MapPL,RprobeL,MapPL.grid_width);
    
    if (oprbLfile[0] != '\0'){
     Output_CHAR3DMAP_PDB_format(oprbLfile,'w',&MapPL,1,"Large probe sphere",&ComStrHead);
    }
 
    printf("#Ngrid_probeS %4d probeL %5d\n",Ngrid_probeS,Ngrid_probeL);

     if (strcmp(MODE,"P")==0){
      Dilation(&MapX,&MapPL,1,2,'O');
      Erosion(&MapX,&MapPL,2,4);
/*
      if (oMoVpdbfile[0]!= '\0'){
       Output_CHAR3DMAP_PDB_format(oMoVpdbfile,'w',&MapX,4,"Molecular Volume",&ComStrHead); 
      }
*/

      Extract_A_and_notB_3DMAP(&MapX,4,1,8);
      Ngrid_pocket = Count_Specific_3DMAP(&MapX,8);
      printf("#Ngrid_pocket_before_openingS %d\n",Ngrid_pocket);
      Erosion(&MapX,&MapPS,8,16);
      Dilation(&MapX,&MapPS,16,64,'O');

      if (oMoVpdbfile[0]!= '\0'){
       Output_CHAR3DMAP_PDB_format(oMoVpdbfile,'w',&MapX,4,"Molecular Volume",&ComStrHead); 
      }

      Reset_Specific_3DMAP(&MapX,2);
      Reset_Specific_3DMAP(&MapX,4); 
      Reset_Specific_3DMAP(&MapX,8); 
      Reset_Specific_3DMAP(&MapX,16); 
      Reset_Specific_3DMAP(&MapX,32); 
      Ngrid_pocket = Count_Specific_3DMAP(&MapX,64);
      printf("#Ngrid_pocket %d\n",Ngrid_pocket);
      sprintf(comment,"Pocket detected with Rsmall: %.2f A Rlarge: %.2f A. Volume: %d grids, %.2f AAA",
       RprobeS, RprobeL, Ngrid_pocket, Ngrid_pocket*MapX.grid_width*MapX.grid_width*MapX.grid_width);  
      Add_String_to_STRING_NODE(&ComStrHead,comment);  
   }


    if (strcmp(MODE,"p")==0){
      Ngrid_probeS = Make_Sphere_Probe_Map(&MapPS,RprobeS,MapPS.grid_width);
      Ngrid_probeL = Make_Sphere_Probe_Map(&MapPL,RprobeL,MapPL.grid_width);
      Output_CHAR3DMAP_PDB_format("probeL.gdb",'w',&MapPL,1,"Large probe sphere",&ComStrHead);
   
      printf("#Ngrid_probeS %4d probeL %5d\n",Ngrid_probeS,Ngrid_probeL);
 
      Dilation(&MapX,&MapPL,1,2,'O');
      Erosion(&MapX,&MapPL,2,4);
    
      if (oMoVpdbfile[0]!= '\0'){
       Output_CHAR3DMAP_PDB_format(oMoVpdbfile,'w',&MapX,4,"Molecular Volume",&ComStrHead); 
      }
      Dilation(&MapX,&MapPS,1,8,'O');
      Dilation(&MapX,&MapPS,8,16,'O');
   
      Extract_A_and_B_and_notC_3DMAP(&MapX,4,16,1,32);
      Extract_Specific_3DMAP(&MapX,32);

      Erosion(&MapX,&MapPS,32,16);
      Dilation(&MapX,&MapPS,16,64,'O');
  
      sprintf(comment,"Pocket detected with Rsmall: %.2f A Rlarge: %.2f A. Volume: %d grids, %.2f AAA",
         RprobeS, RprobeL, Ngrid_pocket, Ngrid_pocket*MapX.grid_width*MapX.grid_width*MapX.grid_width);  
      Add_String_to_STRING_NODE(&ComStrHead,comment);  
   }

    printf("#opockmapfile '%s'\n",opockmapfile);
    printf("#ClusPocketType %c\n",ClusPocketType); 
    /** clustering pocket **/
    cluster_success = 0;
    if (ClusPocketType != 'F'){
      printf("#ClusPocketType %c\n",ClusPocketType); 
      ClusterMapHead.next = NULL; 
      ClusterMapHead.prev = NULL; 
      if (PAR.RecursiveQueue=='R'){
        Make_Clusters_for_Binary_CHAR3DMAP_Recursive(&MapX,64,128,2,&ClusterMapHead);
      }
      else{
        Make_Clusters_for_Binary_CHAR3DMAP_Queue(&MapX,64,128,2,&ClusterMapHead);
      }
      Merge_Sort_Double_Linked_List_CHAR3DMAP(&ClusterMapHead,'D');

      if (opdbfile[0]!='\0'){ Write_PDB_File(opdbfile,'w',&ProAtomHead,"",&ComStrHead);}
      
      if (opockpdbfile[0] != '\0'){
        Output_Clusters_of_CHAR3DMAP_in_PDB_format(opockpdbfile,'w',&MapX,&ClusterMapHead,64);
      }
      if (opockmapfile[0] != '\0'){
        Output_Clusters_of_CHAR3DMAP_in_CCP4_format(opockmapfile,&ClusterMapHead,128);
      }

      if (oallpdbfile[0] != '\0'){
        Write_PDB_File(oallpdbfile,'w',&ProAtomHead,"",&ComStrHead);
        Output_Clusters_of_CHAR3DMAP_in_PDB_format(oallpdbfile,'a',&MapX,&ClusterMapHead,64);
      }
    }

    /** not clustering pocket **/
    if (ClusPocketType == 'F'){
      printf("#ClusPocketType %c\n",ClusPocketType); 
 
      Assign_Map_Value_to_tFactor_of_Atoms(&ProAtomHead,&MapX,'B',64); 
      sprintf(comment,"The tFactor of ATOM indicates a value of its neighboring grid.");
      Add_String_to_STRING_NODE(&ComStrHead,comment);  

      if (opdbfile[0]!='\0'){
        Write_PDB_File(opdbfile,'w',&ProAtomHead,PAR.COMMAND,&ComStrHead);
      }
 
      if (SurfacePocket=='T'){ 
         Extract_Surface_3DMAP(&MapX,64,128);
         Extract_Specific_3DMAP(&MapX,128); 
         Output_CHAR3DMAP_PDB_format(opockpdbfile,'w',&MapX,1,"Surface of Pocket Volume",&ComStrHead);  
         if (oallpdbfile[0]!='\0') Output_CHAR3DMAP_PDB_format(oallpdbfile,'w',&MapX,1,"Surface of Pocket Volume",&ComStrHead);  
       }
       else{
         Extract_Specific_3DMAP(&MapX,64);
         if (opockpdbfile[0] != '\0'){
           Output_CHAR3DMAP_PDB_format(opockpdbfile,'w',&MapX,1,"Pocket Volume",&ComStrHead);  
         }
         if (opockmapfile[0] != '\0'){
           Write_MapCCP4_File_CHAR3DMAP(opockmapfile,&MapX,'D');
         }
       
         if (oallpdbfile[0]!='\0'){ 
           Output_CHAR3DMAP_PDB_format(oallpdbfile,'w',&MapX,1,"Pocket Volume",&ComStrHead); 
         }
       }
    }

      /*** output the spherical probes **/ 
     if (GenSphePrb=='T'){
        Label_Probe_Atoms_By_CHAR3DMAP(&ProbeAtomHead,&MapX,0,'T','A',1.0,0.0); 
        Remove_Probe_Atoms_Smaller_Specified_tFactor(&ProbeAtomHead,0.5); 
        Label_Probe_Atoms_By_CHAR3DMAP(&ProbeAtomHead,&MapX,0,'C','M',1.0,0.0); 
        Sort_Probe_By_Pocket_Cluster_Number(&ProbeAtomHead); 
        if (ospheprbfile[0]!='\0'){
         Write_PDB_File(ospheprbfile,'w',&ProbeAtomHead, "PROBES: Pocket for spherical probe atoms",&ComStrHead);
        }
        if (ospheprbfile_dock[0]!='\0'){
         Write_Spherical_Probes_in_UCSF_DOCK_format(ospheprbfile_dock,&ProbeAtomHead);
        }

        if (oallpdbfile[0]!='\0'){
         Write_PDB_File_with_Multiscale_MolVol_and_Pocket(oallpdbfile,'a',&ProbeAtomHead,'P',&MSprb,'R',
                                    "PROBES:Multiscale MolVol/Pocket for spherical probe atoms",&ComStrHead,'-',ClusPocketType);
        }
      }
 }




 /********************************************************************/
 /** MODE 'CP': Cave Pocket detection                                */
 /* Pocket.KGlike  = (Outer[(X^c erosion P)^c] dilation P)^c and  X^c) opening S */
 /*                               OR                                 */ 
 /* Pocket.Manak  = (Outer[(X dilation P)^c] dilation P and  X)^c opening S */
 /*******************************************************/

  if ((strncmp(MODE,"CP",2)==0) ||(strncmp(MODE,"CPKG",2)==0)){

    Ngrid_pocket = 0;

    Ngrid_probeS = Make_Sphere_Probe_Map(&MapPS,RprobeS,MapPS.grid_width);

   if ((Min_Pocket_GridNum==0)&&(Min_Pocket_ProbeS>0.0)){
    Min_Pocket_GridNum =  (int)(Min_Pocket_ProbeS * Ngrid_probeS); 
    printf("Min_Pocket_ProbeS %f Min_Pocket_GridNum %d\n",Min_Pocket_ProbeS,Min_Pocket_GridNum);
   }

    if (oprbSfile[0] != '\0'){ 
      Output_CHAR3DMAP_PDB_format(oprbSfile,'w',&MapPS,1,"Small probe sphere",&ComStrHead);
    }
 
    Ngrid_probeL = Make_Sphere_Probe_Map(&MapPL,RprobeL,MapPL.grid_width);
    
    if (oprbLfile[0] != '\0'){ 
     Output_CHAR3DMAP_PDB_format(oprbLfile,'w',&MapPL,1,"Large probe sphere",&ComStrHead);
    }

    printf("#Ngrid_probeS %4d probeL %5d\n",Ngrid_probeS,Ngrid_probeL);




 /* Pocket.KGlike  = (Out[ (X closing P)^c - P)^c and  X^c) opn S */
   if (strcmp(MODE,"CPKG")==0){
     Dilation(&MapX,&MapPL,1,2,'O');
     Erosion(&MapX,&MapPL,2,4);

     printf("#Ngrid_ClosingP %d\n",Count_Specific_3DMAP(&MapX,4));

     Inverse_Specific_3DMAP(&MapX,4);
     Reset_Specific_3DMAP(&MapX,2);

     Bounded_Erosion(&MapX,&MapPL,4,2,'O');
     printf("#Ngrid_ClosingP_comp_erosionP %d\n",Count_Specific_3DMAP(&MapX,2));
     Output_CHAR3DMAP_with_tar_bitmask_PDB_format("XclsP_c_erP.pdb",&MapX,2);
     if (PAR.RecursiveQueue == 'R'){
       Extract_Outerspace_Region_Recursive(&MapX,2,8);
     }
     else{
       Extract_Outerspace_Region_Queue(&MapX,2,8);
     }
     printf("#Ngrid_ClosingP_comp_erosionP_outer %d\n",Count_Specific_3DMAP(&MapX,8));
     Output_CHAR3DMAP_with_tar_bitmask_PDB_format("XclsP_c_erP_out.pdb",&MapX,8);

     Extract_A_and_notB_3DMAP(&MapX,2,8,128);

     Output_CHAR3DMAP_with_tar_bitmask_PDB_format("diff.pdb",&MapX,128);

     Dilation(&MapX,&MapPL,8,16,'O');
     Extract_notA_and_notB_3DMAP(&MapX,16,1,32);
     Ngrid_pocket = Count_Specific_3DMAP(&MapX,32);
     printf("#Ngrid_pocket_before_openingS %d\n",Ngrid_pocket);
     Erosion(&MapX,&MapPS,32,64);
     Dilation(&MapX,&MapPS,64,128,'O');
  }


 /* Pocket.Manak  = (Out[(X dilation P)^c] dilation P and  X)^c opening S */
   else {
    Dilation(&MapX,&MapPL,1,2,'O');
    Inverse_Specific_3DMAP(&MapX,2);
    if (PAR.RecursiveQueue == 'R'){
     Extract_Outerspace_Region_Recursive(&MapX,2,8);
    }
    else{
      Extract_Outerspace_Region_Queue(&MapX,2,8); 
    }

    Dilation(&MapX,&MapPL,8,16,'O');
    Extract_A_or_B_3DMAP(&MapX,1,16,32);
    Inverse_Specific_3DMAP(&MapX,32);

    Ngrid_pocket = Count_Specific_3DMAP(&MapX,32);
    printf("#Ngrid_pocket_before_openingS %d\n",Ngrid_pocket);

    Erosion(&MapX,&MapPS,32,64);
    Dilation(&MapX,&MapPS,64,128,'O');
   }

   
    if (oMoVpdbfile[0]!= '\0'){
     Output_CHAR3DMAP_PDB_format(oMoVpdbfile,'w',&MapX,8,"Molecular Volume",&ComStrHead); 
    }

    Reset_Specific_3DMAP(&MapX,2);
    Reset_Specific_3DMAP(&MapX,4); 
    Reset_Specific_3DMAP(&MapX,8); 
    Reset_Specific_3DMAP(&MapX,16); 
    Reset_Specific_3DMAP(&MapX,32); 
    Reset_Specific_3DMAP(&MapX,64); 
    Ngrid_pocket = Count_Specific_3DMAP(&MapX,128);
    printf("#Ngrid_pocket %d\n",Ngrid_pocket);
    sprintf(comment,"Pocket detected with Rsmall: %.2f A Rlarge: %.2f A. Volume: %d grids, %.2f AAA",
     RprobeS, RprobeL, Ngrid_pocket, Ngrid_pocket*MapX.grid_width*MapX.grid_width*MapX.grid_width);  
    Add_String_to_STRING_NODE(&ComStrHead,comment);  


    /** clustering pocket **/
    cluster_success = 0;
    if (ClusPocketType != 'F'){
      ClusterMapHead.next = NULL; 
      ClusterMapHead.prev = NULL; 
      if (PAR.RecursiveQueue=='R'){
        Make_Clusters_for_Binary_CHAR3DMAP_Recursive(&MapX,128,64,32,&ClusterMapHead);
      }
      else{
        Make_Clusters_for_Binary_CHAR3DMAP_Queue(&MapX,128,64,32,&ClusterMapHead);
      }
      Merge_Sort_Double_Linked_List_CHAR3DMAP(&ClusterMapHead,'D');
      printf("#opockpdbfile '%s'\n",opockpdbfile);
      if (opdbfile[0]!='\0'){ Write_PDB_File(opdbfile,'w',&ProAtomHead,"",&ComStrHead);}
      
      if (opockpdbfile[0] != '\0'){
        Output_Clusters_of_CHAR3DMAP_in_PDB_format(opockpdbfile,'w',&MapX,&ClusterMapHead,128);
      } 
      if (opockmapfile[0] != '\0'){
        Output_Clusters_of_CHAR3DMAP_in_CCP4_format(opockmapfile,&ClusterMapHead,128);
      }

      if (oallpdbfile[0] != '\0'){
        Write_PDB_File(oallpdbfile,'w',&ProAtomHead,"",&ComStrHead);
        Output_Clusters_of_CHAR3DMAP_in_PDB_format(oallpdbfile,'a',&MapX,&ClusterMapHead,128);
      }
      cluster_success = 1;
    }

    /** not clustering pocket **/
    if ((ClusPocketType == 'F')||(cluster_success==0)){
 
     Assign_Map_Value_to_tFactor_of_Atoms(&ProAtomHead,&MapX,'B',128); 
     sprintf(comment,"The tFactor of ATOM indicates a value of its neighboring grid.");
     Add_String_to_STRING_NODE(&ComStrHead,comment);  

     if (opdbfile[0]!='\0'){
        Write_PDB_File(opdbfile,'w',&ProAtomHead,PAR.COMMAND,&ComStrHead);
     }
 
     if (SurfacePocket=='T'){ 
       Extract_Surface_3DMAP(&MapX,128,64);
       Extract_Specific_3DMAP(&MapX,64); 
       Output_CHAR3DMAP_PDB_format(opockpdbfile,'w',&MapX,1,"Surface of Pocket Volume",&ComStrHead);  
       if (oallpdbfile[0]!='\0') Output_CHAR3DMAP_PDB_format(oallpdbfile,'w',&MapX,1,"Surface of Pocket Volume",&ComStrHead);  
      }
     else{
       Extract_Specific_3DMAP(&MapX,128);
       if (opockpdbfile[0] != '\0'){
         Output_CHAR3DMAP_PDB_format(opockpdbfile,'w',&MapX,1,"Pocket Volume",&ComStrHead);  
       }
       if (opockmapfile[0] != '\0'){
         Write_MapCCP4_File_CHAR3DMAP(opockmapfile,&MapX,'D');
       }
       if (oallpdbfile[0]!='\0'){ Output_CHAR3DMAP_PDB_format(oallpdbfile,'w',&MapX,1,"Pocket Volume",&ComStrHead); }
      }
    }

      /*** output the spherical probes **/ 
     if (GenSphePrb=='T'){
       Label_Probe_Atoms_By_CHAR3DMAP(&ProbeAtomHead,&MapX,0,'T','A',1.0,0.0); 
       Remove_Probe_Atoms_Smaller_Specified_tFactor(&ProbeAtomHead,0.5); 
       Label_Probe_Atoms_By_CHAR3DMAP(&ProbeAtomHead,&MapX,0,'C','M',1.0,0.0); 
       Sort_Probe_By_Pocket_Cluster_Number(&ProbeAtomHead); 
       if (ospheprbfile[0]!='\0'){
        Write_PDB_File(ospheprbfile,'w',&ProbeAtomHead, "PROBES: Pocket for spherical probe atoms",&ComStrHead);
       }
       if (ospheprbfile_dock[0]!='\0')
        Write_Spherical_Probes_in_UCSF_DOCK_format(ospheprbfile_dock,&ProbeAtomHead);

       if (oallpdbfile[0]!='\0')
        Write_PDB_File_with_Multiscale_MolVol_and_Pocket(oallpdbfile,'a',&ProbeAtomHead,'P',&MSprb,'R',
                                   "PROBES:Multiscale MolVol/Pocket for spherical probe atoms",&ComStrHead,'-',ClusPocketType);
      }
 }



 /*****************************************************/
 /** MODE 'CorP': Cavity-or-Pocket                   **/ 
 /**  (((X closing P) and X^c) or X.cavity ) opening S **/    
 /**= (((X closing P) and X^c) or (fi[X^c-P]+P))   opening S **/    
 /**= (((X closing P) and X^c) or (fi[(X+P)^c]+P)) opening S **/    
 /*****************************************************/

 if (strcmp(MODE,"CorP")==0){
    Ngrid_probeS = Make_Sphere_Probe_Map(&MapPS,RprobeS,MapPS.grid_width);
    Ngrid_probeL = Make_Sphere_Probe_Map(&MapPL,RprobeL,MapPS.grid_width);
    /* Output_CHAR3DMAP_PDB_format("probe.gdb",'w',&MapP,1,"Probe Sphere",&ComStrHead); */
   
    /* 2: X+P */
    Dilation(&MapX,&MapPL,1,2,'O');
    /* 2: (X+P)^c */
    Inverse_Specific_3DMAP(&MapX,2);
    /* 4: fo[(X+P)^c] */
    Extract_Outerspace_Region_Queue(&MapX,2,4);
    /* 8: fi[(X+P)^c] */
    Extract_A_and_notB_3DMAP(&MapX,2,4,8);
    /* 16: fi[(X+P)^c]+P */
    Dilation(&MapX,&MapPL,8,16);
    Output_CHAR3DMAP_PDB_format("cavity.pdb",'w',&MapX,16,"",&ComStrHead); 

    /* 2: X+P */
    Inverse_Specific_3DMAP(&MapX,2);
    /* 32: (X+P)-P  = X clsng P */
    Erosion(&MapX,&MapPL,2,32,'O');

    /* 64: (((X clsing P) and X^c) ) **/    
    Extract_A_and_notB_3DMAP(&MapX,32,1,64);
    Output_CHAR3DMAP_PDB_format("XclsngP_notX.pdb",'w',&MapX,64,"",&ComStrHead); 
 
    /* 128: (((X clsing P) and X^c) or (fi[(X+P)^c]+P)) **/    
    Extract_A_or_B_3DMAP(&MapX,16,64,128);
    Output_CHAR3DMAP_PDB_format("CorP.pdb",'w',&MapX,128,"",&ComStrHead); 
    /* 64: (((X closing P) and X^c) or (fi[(X+P)^c]+P)) opening S **/    
    Reset_Specific_3DMAP(&MapX,2);
    Reset_Specific_3DMAP(&MapX,4);
    Erosion(&MapX,&MapPS,128,2,'O');
    Dilation(&MapX,&MapPS,2,4,'O');

    Ngrid_pocket = Extract_Specific_3DMAP(&MapX,4);

    sprintf(comment,"Cavity detected with Radius: %.2f A. Volume: %d grids, %.2f AAA",
     RprobeS,Ngrid_pocket,Ngrid_pocket*MapX.grid_width*MapX.grid_width*MapX.grid_width);  

    Add_String_to_STRING_NODE(&ComStrHead,comment);

    if (opockpdbfile[0] != '\0'){
      Output_CHAR3DMAP_PDB_format(opockpdbfile,'w',&MapX,4,"cavity volume",&ComStrHead); 
    }

    if (opockmapfile[0] != '\0'){
      Write_MapCCP4_File_CHAR3DMAP(opockmapfile,&MapX,'D');
    }


  }








 /**********************************************************/
 /** MODE 'e': Eroded-Pocket detection                    **/
 /*  StandardPocket = ((X closing P) and noX)) opening S  **/
 /*                 = (((X closing P) and noX)) - S) + S) **/
 /*  ErodedPocket   = (((X closing P) and noX)) - S)      **/
 /*                 = StandardPocket - S                  **/
 /**********************************************************/

   if (strcmp(MODE,"e")==0){

    Ngrid_pocket = 0;

    Ngrid_probeS = Make_Sphere_Probe_Map(&MapPS,RprobeS,MapPS.grid_width);

   if ((Min_Pocket_GridNum==0)&&(Min_Pocket_ProbeS>0.0)){
    Min_Pocket_GridNum =  (int)(Min_Pocket_ProbeS * Ngrid_probeS); 
    printf("Min_Pocket_ProbeS %f Min_Pocket_GridNum %d\n",Min_Pocket_ProbeS,Min_Pocket_GridNum);
   }
    if (oprbSfile[0] != '\0') 
      Output_CHAR3DMAP_PDB_format(oprbSfile,'w',&MapPS,1,"Small probe sphere",&ComStrHead);
    
    Ngrid_probeL = Make_Sphere_Probe_Map(&MapPL,RprobeL,MapPL.grid_width);
    
    if (oprbLfile[0] != '\0') 
     Output_CHAR3DMAP_PDB_format(oprbLfile,'w',&MapPL,1,"Large probe sphere",&ComStrHead);
   
    printf("#Ngrid_probeS %4d probeL %5d\n",Ngrid_probeS,Ngrid_probeL);

    Dilation(&MapX,&MapPL,1,2,'O');
    Erosion(&MapX,&MapPL,2,4);
/*
      if (oMoVpdbfile[0]!= '\0'){
       Output_CHAR3DMAP_PDB_format(oMoVpdbfile,'w',&MapX,4,"Molecular Volume",&ComStrHead); 
      }
*/

    Extract_A_and_notB_3DMAP(&MapX,4,1,8);
  
    /*
    Erosion(&MapX,&MapPS,8,16);
    Dilation(&MapX,&MapPS,16,64,'O');
    */

    Erosion(&MapX,&MapPS,8,64);

    if (oMoVpdbfile[0]!= '\0'){
     Output_CHAR3DMAP_PDB_format(oMoVpdbfile,'w',&MapX,4,"Molecular Volume",&ComStrHead); 
    }

     Reset_Specific_3DMAP(&MapX,2);
     Reset_Specific_3DMAP(&MapX,4); 
     Reset_Specific_3DMAP(&MapX,8); 
     Reset_Specific_3DMAP(&MapX,16); 
     Reset_Specific_3DMAP(&MapX,32); 
     Ngrid_pocket = Count_Specific_3DMAP(&MapX,64);
     sprintf(comment,"Pocket detected with Rsmall: %.2f A Rlarge: %.2f A. Volume: %d grids, %.2f AAA",
     RprobeS, RprobeL, Ngrid_pocket, Ngrid_pocket*MapX.grid_width*MapX.grid_width*MapX.grid_width);  
     Add_String_to_STRING_NODE(&ComStrHead,comment);  


    /** clustering pocket **/
    cluster_success = 0;
    if (ClusPocketType != 'F'){
      cluster_success = Label_Connected_Bit_Clusters(&MapX,64,Min_Pocket_GridNum,&Ncluster,Nbit_cluster);
      Assign_Map_Value_to_tFactor_of_Atoms(&ProAtomHead,&MapX,'C',0);
      sprintf(comment,"The tFactor of ATOM indicates Cluster numbers of its neighboring grid.");
      Add_String_to_STRING_NODE(&ComStrHead,comment);  
 
      if (opdbfile[0]!='\0'){ Write_PDB_File(opdbfile,'w',&ProAtomHead,"",&ComStrHead);}
      
      if (opockpdbfile[0] != '\0'){
        Output_CHAR3DMAP_PDB_format_Grouped_By_Char_Value(opockpdbfile,'w',&MapX,Ncluster,Nbit_cluster,
                 "Clustered Pockets",&ComStrHead);
      }


      if (oallpdbfile[0] != '\0'){
        Output_CHAR3DMAP_PDB_format_Grouped_By_Char_Value(oallpdbfile,'w',&MapX,Ncluster,Nbit_cluster,
                 "Clustered Pockets",&ComStrHead);
      }
    }

    /** not clustering pocket **/
    if ((ClusPocketType == 'F')||(cluster_success==0)){
 
     Assign_Map_Value_to_tFactor_of_Atoms(&ProAtomHead,&MapX,'B',64); 
     sprintf(comment,"The tFactor of ATOM indicates a value of its neighboring grid.");
     Add_String_to_STRING_NODE(&ComStrHead,comment);  

     if (opdbfile[0]!='\0'){
        Write_PDB_File(opdbfile,'w',&ProAtomHead,PAR.COMMAND,&ComStrHead);
     }
 
     if (SurfacePocket=='T'){ 
       Extract_Surface_3DMAP(&MapX,64,128);
       Extract_Specific_3DMAP(&MapX,128); 
       Output_CHAR3DMAP_PDB_format(opockpdbfile,'w',&MapX,1,"Surface of Pocket Volume",&ComStrHead);  
       if (oallpdbfile[0]!='\0'){ Output_CHAR3DMAP_PDB_format(oallpdbfile,'w',&MapX,1,"Surface of Pocket Volume",&ComStrHead);}
      }
     else{
       Extract_Specific_3DMAP(&MapX,64);
       Output_CHAR3DMAP_PDB_format(opockpdbfile,'w',&MapX,1,"Pocket Volume",&ComStrHead);  
       if (oallpdbfile[0]!='\0'){ Output_CHAR3DMAP_PDB_format(oallpdbfile,'w',&MapX,1,"Pocket Volume",&ComStrHead); } 
      }
    }
 }















 /************************************/
 /** MODE 'm': real masuya and doi  **/
 /************************************/
  
  if (strcmp(MODE,"m")==0){

    Ngrid_probeS = Make_Sphere_Probe_Map(&MapPS,RprobeS,MapPS.grid_width);

   if ((Min_Pocket_GridNum==0)&&(Min_Pocket_ProbeS>0.0)){
    Min_Pocket_GridNum =  (int)(Min_Pocket_ProbeS * Ngrid_probeS); 
    printf("Min_Pocket_ProbeS %f Min_Pocket_GridNum %d\n",Min_Pocket_ProbeS,Min_Pocket_GridNum);
   }
    Output_CHAR3DMAP_PDB_format("probeS.gdb",'w',&MapPS,1,"Small probe sphere",&ComStrHead);
     
    Ngrid_probeL = Make_Sphere_Probe_Map(&MapPL,RprobeL,MapPL.grid_width);
    Output_CHAR3DMAP_PDB_format("probeL.gdb",'w',&MapPL,1,"Large probe sphere",&ComStrHead);
   
    printf("#Ngrid_probeS %4d probeL %5d\n",Ngrid_probeS,Ngrid_probeL);

    Dilation(&MapX,&MapPS,1,2,'O');
    Erosion(&MapX,&MapPS,2,4,'O');
    
    Dilation(&MapX,&MapPL,1,8,'O');
    Erosion(&MapX,&MapPL,8,16,'O');

    Extract_A_and_notB_3DMAP(&MapX,16,4,32);

    Erosion(&MapX,&MapPS,32,64,'O');
    Dilation(&MapX,&MapPS,64,128,'O');
    if (opockpdbfile[0] != '\0'){
     Output_CHAR3DMAP_PDB_format(opockpdbfile,'w',&MapX,128,"real masuya and doi pocket",&ComStrHead); 
    }

  }





 /************************************************/
 /** MODE 'V': Cavity detection (center-conneced)**/
 /** U_i{ fi[X^c-P]+P } = U_i{fi[(X+P)^c]+P } **/
 /************************************************/

 if (strcmp(MODE,"V")==0){
    if (opdbfile[0]!='\0'){ Write_PDB_File(opdbfile,'w',&ProAtomHead,"",&ComStrHead);}

    Ngrid_probe = Make_Sphere_Probe_Map(&MapP,RprobeL,MapP.grid_width);
    /* Output_CHAR3DMAP_PDB_format("probe.gdb",'w',&MapP,1,"Probe Sphere",&ComStrHead); */
   
    Inverse_Specific_and_Set_3DMAP(&MapX,1,2);
    Bounded_Erosion(&MapX,&MapP,2,4,'O');

    if (PAR.RecursiveQueue == 'R'){
       Extract_Outerspace_Region_Recursive(&MapX,4,8);
    }
    else{
      Extract_Outerspace_Region_Queue(&MapX,4,8);
    }

    Extract_A_and_notB_3DMAP(&MapX,4,8,16);

    Dilation(&MapX,&MapP,16,32);

    Extract_Specific_3DMAP(&MapX,32);

    Ngrid_pocket = Extract_Specific_3DMAP(&MapX,32);

    sprintf(comment,"Cavity detected with Radius: %.2f A. Volume: %d grids, %.2f AAA",
     RprobeS,Ngrid_pocket,Ngrid_pocket*MapX.grid_width*MapX.grid_width*MapX.grid_width);  

    Add_String_to_STRING_NODE(&ComStrHead,comment);
    cluster_success = 0;

    if (ClusPocketType != 'F'){
      ClusterMapHead.next = NULL; 
      ClusterMapHead.prev = NULL; 
      if (PAR.RecursiveQueue=='R'){
        Make_Clusters_for_Binary_CHAR3DMAP_Recursive(&MapX,32,64,128,&ClusterMapHead);
      }
      else{
        Make_Clusters_for_Binary_CHAR3DMAP_Queue(&MapX,32,64,128,&ClusterMapHead);
      }
      Merge_Sort_Double_Linked_List_CHAR3DMAP(&ClusterMapHead,'D');
      if (opockpdbfile[0] != '\0'){
        Output_Clusters_of_CHAR3DMAP_in_PDB_format(opockpdbfile,'w',&MapX,&ClusterMapHead,32);
      } 
      if (opockmapfile[0] != '\0'){
        Output_Clusters_of_CHAR3DMAP_in_CCP4_format(opockmapfile,&ClusterMapHead,32);
      }
   }
   else{
     if (opockpdbfile[0] != '\0'){
       Output_CHAR3DMAP_PDB_format(opockpdbfile,'w',&MapX,32,"cavity volume",&ComStrHead); 
     }

     if (opockmapfile[0] != '\0'){
       Write_MapCCP4_File_CHAR3DMAP(opockmapfile,&MapX,'D');
     }
   }
  }


 /************************************************/
 /** MODE 'CV': Connoly-void Cavity (probe-conneced)**/
 /** Union_i fi[(X closing P)^c] = Union_fi[X^c opening P]    **/
 /************************************************/
 if (strcmp(MODE,"v")==0){
    Ngrid_probe = Make_Sphere_Probe_Map(&MapP,RprobeL,MapP.grid_width);
    if (oprbLfile[0] != '\0'){
     Output_CHAR3DMAP_PDB_format(oprbLfile,'w',&MapP,1,"Large probe sphere",&ComStrHead);
    }
    /*
    Ngrid_probe = Make_Sphere_Probe_Map(&MapP,RprobeL,MapP.grid_width);
    */   
 /* Output_CHAR3DMAP_PDB_format("probe.gdb",'w',&MapP,1,"Probe Sphere",&ComStrHead); */
   
    /* X closing P --> '8' */  
    Dilation(&MapX,&MapP,1,4,'O');
    Erosion(&MapX,&MapP,4,8);
    Extract_Specific_3DMAP(&MapX,8);

    if (oMoVpdbfile[0]!='\0'){
     Output_CHAR3DMAP_PDB_format(oMoVpdbfile,'w',&MapX,1,"Molecular Volume",&ComStrHead);
    }
    if (oMoVmapfile[0]!='\0'){
      Write_MapCCP4_File_CHAR3DMAP(oMoVmapfile,&MapX,'D');
    }

    /* Cavity for target '8' -->  blank:16, cavity:32 */  
    Extract_Cavity_Region(&MapX,8,16,32);
    Ngrid_pocket = Extract_Specific_3DMAP(&MapX,32);

    sprintf(comment,"Cavity detected with Radius: %.2f A. Volume: %d grids, %.2f AAA",
     RprobeL,Ngrid_pocket,Ngrid_pocket*MapX.grid_width*MapX.grid_width*MapX.grid_width);  

    Add_String_to_STRING_NODE(&ComStrHead,comment);
   
    if (ClusPocketType != 'F'){
      ClusterMapHead.next = NULL; 
      ClusterMapHead.prev = NULL; 
      if (PAR.RecursiveQueue=='R'){
        Make_Clusters_for_Binary_CHAR3DMAP_Recursive(&MapX,32,64,128,&ClusterMapHead);
      }
      else{
        Make_Clusters_for_Binary_CHAR3DMAP_Queue(&MapX,32,64,128,&ClusterMapHead);
      }
      Merge_Sort_Double_Linked_List_CHAR3DMAP(&ClusterMapHead,'D');
      if (opockpdbfile[0] != '\0'){
        Output_Clusters_of_CHAR3DMAP_in_PDB_format(opockpdbfile,'w',&MapX,&ClusterMapHead,32);
      }
    }
    else{
      if (opockpdbfile[0] != '\0'){
        Output_CHAR3DMAP_PDB_format(opockpdbfile,'w',&MapX,32,"cavity volume",&ComStrHead); 
      }
      if (opockmapfile[0] != '\0'){
        Write_MapCCP4_File_CHAR3DMAP(opockmapfile,&MapX,'D');
      }
    }

   /*** [OPTION] output the spherical probes **/ 
   if (GenSphePrb=='T'){
      Label_Probe_Atoms_By_CHAR3DMAP(&ProbeAtomHead,&MapX,32,'T','A',1.0,0.0);
      Remove_Probe_Atoms_Smaller_Specified_tFactor(&ProbeAtomHead,0.5);
    
      if (ospheprbfile[0]!='\0'){
       Write_PDB_File(ospheprbfile,'w',&ProbeAtomHead,"CAVITY PROBES",&ComStrHead);
      }
   }

  }







 /*******************************************/
 /** MODE 'M': Multi Scale Closing /Pocket **/
 /*******************************************/
 
 if (strcmp(MODE,"M")==0){
   
   /** [1] Set Radius Array **/
   if (Bin_of_RprobeL>0.0) 
    Set_Rarray_By_Min_Max_Bin(MSprb.Rarray,&(MSprb.Nradius),RprobeLmin,RprobeLmax,Bin_of_RprobeL);
   else if (MSprb.Nradius==1) MSprb.Rarray[0] = RprobeLmax;
   else if (MSprb.Nradius > 0)
    Set_Rarray_By_Nradius_Min_Max(MSprb.Rarray,MSprb.Nradius,RprobeLmin,RprobeLmax,Bin_of_RprobeL);
   else
   { printf("#ERROR:MSprb.Nradius or Bin_of_RprobeL were not properly assigned.\n");
     exit(1);}

   sprintf(comment,"MIN_RLARGE:%f",RprobeLmin);
   Add_String_to_STRING_NODE(&ComStrHead,comment); 
   sprintf(comment,"MAX_RLARGE:%f",RprobeLmax);
   Add_String_to_STRING_NODE(&ComStrHead,comment); 
   sprintf(comment,"BIN_RLARGE:%f",Bin_of_RprobeL);
   Add_String_to_STRING_NODE(&ComStrHead,comment); 

   /*
   sprintf(comment,"MULSC_PROBE: %3d th for Radius %6.3f angstrom (VdW Volume)",0.0,0,0);
   Add_String_to_STRING_NODE(&ComStrHead,comment); 
   for (k=1;k<=MSprb.Nradius;++k)
    { sprintf(comment,"MULSC_PROBE: %3d th for Radius %6.3f angstrom",k,MSprb.Rarray[k]); 
      Add_String_to_STRING_NODE(&ComStrHead,comment); }
   */
 
  for (k=1;k<=MSprb.Nradius;++k){ 
    printf("#k %d R %f\n",k,MSprb.Rarray[k]);
  }
 
   /** [2] Make MultiScale Shell Sphere Probe **/
   MSprb.grid_width = MapX.grid_width;
   Make_MultiScale_Sphere_Shell_Probe_Map(&MSprb);
   if (omprobfile[0] != '\0')
   Write_PDB_File_MultiSc_Shell_Spheres("mprobe.gdb",&MSprb,"MultiScale Shell Sphere Probe",&ComStrHead);
   
   /** [3] Malloc MapY(addtional CHAR3DMAP for multiscale probe calculations) */
   for (k=0;k<3;++k){ 
     MapY.N[k] = MapX.N[k]; 
     MapY.OrigPos[k] = MapX.OrigPos[k]; 
   }
   MapY.grid_width = MapX.grid_width; 
   Malloc_CHAR3DMAP(&MapY, MapY.N[0], MapY.N[1], MapY.N[2]);

 
   /** [4] MultiScale Closing (Dilation and Erosion)   */
   Dilation_MultiScaleShell(&MapX,&MapY,&MSprb);
   Erosion_MultiScaleShell(&MapY,&MapX,&MSprb);
  /* >> Final states  <<
     MapX : multiscale closing (molecular volume) 
     MapY:  multiscale dilation */ 
   
   /** [5] Closing by Small Probe for the "POCKET" detection  */
    Ngrid_probeS = Make_Sphere_Probe_Map(&MapPS,RprobeS,MapP.grid_width);
    sprintf(comment,"RSMALL:%6.3f",RprobeS);
    Add_String_to_STRING_NODE(&ComStrHead,comment); 
    Opening_MultiScaleShell(&MapX,&MapY,&MapPS);
  
   /* >> Final states  <<
     MapX : multiscale closing (molecular volume) 
     MapY:  multiscale pocket 
   */ 
 
   if (oMoVpdbfile[0]!='\0'){
     Output_CHAR3DMAP_PDB_format_with_MULSC_info(oMoVpdbfile,'w',&MapX,1,255,
                      "Multi-scale Molecular Volume",&MSprb,&ComStrHead);
   }
 
  /*** [5B (optional) ] Chose Spherical Probe in Pocket Regions***/
  if ((ClusPocketType!='F')||(Dtree_SpheProbe>0.0)){
     for (k=0;k<3;++k) { MapZ.N[k] = MapX.N[k]; MapZ.OrigPos[k] = MapX.OrigPos[k]; }
     MapZ.grid_width = MapX.grid_width; 
     Malloc_CHAR3DMAP(&MapZ, MapZ.N[0], MapZ.N[1], MapZ.N[2]);
  }

  if (GenSphePrb=='T'){
    Label_Multiscale_MolVol_or_Pocket_to_Ligand_Atoms(&ProbeAtomHead,&MapX,&MSprb,Wshell,'M');
    Label_Multiscale_MolVol_or_Pocket_to_Ligand_Atoms(&ProbeAtomHead,&MapY,&MSprb,Wshell,'P');
    Remove_Probe_Atoms_Smaller_Specified_Pocketness(&ProbeAtomHead,0.00);

    if (Dtree_SpheProbe>0.0){
       Make_Weighted_Hierarichical_Cluster_Of_Atoms(&Gclus,&ProbeAtomHead,Dtree_SpheProbe,'C');
       Project_Sphe_Probe_Cluster_into_Grid(&Gclus,&ProbeAtomHead,&MapY,&MapZ);
       Sort_Probe_By_Pocket_Cluster_Number(&ProbeAtomHead); 
       if (RankKeepCluster>0){
        Set_Map_to_Zero_for_Low_Ranked_Grid_Cluster(&MapY,&Gclus,RankKeepCluster);
       }
     }

  }
 
   /** [6] Clustering Multiscale pocket and filtering **/
    if ((ClusPocketType!='F')&&(Dtree_SpheProbe<=0.0)){ 
     Cluster_MultiScale_Pocket_Single_Linkage(&MapY,&MapZ,&Gclus,threCluster_char_max,'X',&MSprb,'I');
     if (GenSphePrb=='T'){ 
       Label_Multiscale_Each_Pocket_Cluster_to_Ligand_Atoms(&ProbeAtomHead,&MapY,&MapZ,&MSprb);
       Sort_Probe_By_Pocket_Cluster_Number(&ProbeAtomHead); 
     }
     if (RankKeepCluster>0){
       Set_Map_to_Zero_for_Low_Ranked_Grid_Cluster(&MapY,&Gclus,RankKeepCluster);
     } 
   }

   /** [7]  Labeling Protein Atoms by Multiscale MolVol and pocket  */
   
    Label_Multiscale_MolVol_or_Pocket_to_Protein_Atom_Shells(&ProAtomHead,&MapX,&MSprb,Wshell,'M');
    Label_Multiscale_MolVol_or_Pocket_to_Protein_Atom_Shells(&ProAtomHead,&MapY,&MSprb,Wshell,'P');
    if (ClusPocketType!='F'){
      Label_Multiscale_Each_Pocket_Cluster_to_Protein_Atom_Shells(&ProAtomHead,&MapY,&MapZ,&MSprb,Wshell);
    } 

   /** [8]  Labeling Ligand Atoms by Multiscale MolVol and pocket  */
 
   if (iligpdbfile[0]!='\0'){
     Read_PDB_File(iligpdbfile,&LigAtomHead,pdbidLig,'H',"-",'h','F','F',MODELtype);
     NatomPdb = Number_Of_ATOMs(&LigAtomHead); 
     printf("#NatomLig %d\n",NatomPdb);
     printf("#iligRadiusType %c\n",iligRadiusType); 
     Assign_Radius_to_ATOMs(&LigAtomHead,iligRadiusType,Runified);
 /*
    if (iligRadiusType == 'O') {
     printf("#Assign_Radius_By_Occup(&LigAtomHead);\n");
     Assign_Radius_By_Occup(&LigAtomHead);
    }
     else Assign_Radius(&LigAtomHead,&HeadRadius);
    */
  
    sprintf(comment,"INPUT_LIGAND_PDB_FILE:%s",iligpdbfile);
    Add_String_to_STRING_NODE(&ComStrHead,comment); 
    sprintf(comment,"NATOM_OF_LIGAND:%d",NatomPdb);
    Add_String_to_STRING_NODE(&ComStrHead,comment); 
    Cal_Number_of_Contacting_Atoms_with_Ligand(&ProAtomHead,&LigAtomHead);

    /* Set_Bit_3DMAP_Inside_Atoms(&MapX,&LigAtomHead,128,1.0); */ 
    Label_Multiscale_MolVol_or_Pocket_to_Ligand_Atoms(&LigAtomHead,&MapX,&MSprb,Wshell,'M');
    Label_Multiscale_MolVol_or_Pocket_to_Ligand_Atoms(&LigAtomHead,&MapY,&MSprb,Wshell,'P');
    if (ClusPocketType!='F'){
      Label_Multiscale_Each_Pocket_Cluster_to_Ligand_Atoms(&LigAtomHead,&MapY,&MapZ,&MSprb);
    }
   }
   
    /** [9] Output result files  **/

   /***[9-1] output pocket in PDB **/ 
   if (opockpdbfile[0]!='\0'){
     Output_CHAR3DMAP_PDB_format_with_MULSC_info(opockpdbfile,'w',&MapY,1,254,
                      "Pockets by MultiScale Probe Sphere",&MSprb,&ComStrHead);
    }
   if (oallpdbfile[0]!='\0'){
     Output_CHAR3DMAP_PDB_format_with_MULSC_info(oallpdbfile,'w',&MapY,1,254,
                      "Pockets by MultiScale Probe Sphere",&MSprb,&ComStrHead);
    }

   if (ClusPocketType!='F'){   
     if (opockpdbfile[0]!='\0')
       Output_CHAR3DMAP_PDB_format_By_Grid_Cluster(
         opockpdbfile,'w',&MapY,&Gclus,&MSprb,0,RankKeepCluster,"Clustered_MultiScale_Pocket",&ComStrHead);
     if (oallpdbfile[0]!='\0')
       Output_CHAR3DMAP_PDB_format_By_Grid_Cluster(
         oallpdbfile,'w',&MapY,&Gclus,&MSprb,0,RankKeepCluster,"Clustered_MultiScale_Pocket",&ComStrHead);
   }

   /***[9-2] output pocket in 3D density map **/ 
   if (omapfile[0]!='\0') {
      Write_MapCCP4_File_CHAR3DMAP(omapfile,&MapY, 'R', MSprb.Nradius, MSprb.Rarray);
   }
 
   if (opockmapfile[0]!='\0') {
      Write_MapCCP4_File_CHAR3DMAP(opockmapfile,&MapY, 'R', MSprb.Nradius, MSprb.Rarray);
   }

   /***[9-3] output the protein **/ 
    if (opdbfile[0]!='\0'){
       Write_PDB_File_with_Multiscale_MolVol_and_Pocket(opdbfile,'w',&ProAtomHead,'A',&MSprb,tFactorType,
                                  "RECEPTOR:Multiscale MolVol/Pocket for protein atoms",&ComStrHead,'-',ClusPocketType);
    }

    if (oallpdbfile[0]!='\0'){
       Write_PDB_File_with_Multiscale_MolVol_and_Pocket(oallpdbfile,'a',&ProAtomHead,'A',&MSprb,tFactorType,
                                  "RECEPTOR:Multiscale MolVol/Pocket for protein atoms",&ComStrHead,'-',ClusPocketType);
    }

   /***[9-3] output the ligand **/ 
    
   if (iligpdbfile[0] !='\0'){
    if (oligpdbfile[0] != '\0'){ 
      Write_PDB_File_with_Multiscale_MolVol_and_Pocket(oligpdbfile,'w',&LigAtomHead,'H',&MSprb,'R',
                                  "LIGAND:Multiscale MolVol/Pocket for ligand atoms",&ComStrHead,'-',ClusPocketType);
    }

     if (oallpdbfile[0] != '\0'){ 
       Write_PDB_File_with_Multiscale_MolVol_and_Pocket(oallpdbfile,'a',&LigAtomHead,'H',&MSprb,'R',
                                   "LIGAND:Multiscale MolVol/Pocket for ligand atoms",&ComStrHead,'-',ClusPocketType);
     }
    }

   /***[9-4] output the spherical probes **/ 
   if (GenSphePrb=='T'){

     if (ospheprbfile[0]!='\0'){
      Write_PDB_File_with_Multiscale_MolVol_and_Pocket(ospheprbfile,'w',&ProbeAtomHead,'P',&MSprb,'R',
                                 "PROBES:Multiscale MolVol/Pocket for spherical probe atoms",&ComStrHead,'-',ClusPocketType);
     }

     if (ospheprbfile_dock[0]!='\0'){
      Write_Spherical_Probes_in_UCSF_DOCK_format(ospheprbfile_dock,&ProbeAtomHead);
     }

     if (oallpdbfile[0]!='\0'){
       Write_PDB_File_with_Multiscale_MolVol_and_Pocket(oallpdbfile,'a',&ProbeAtomHead,'P',&MSprb,'R',
                                 "PROBES:Multiscale MolVol/Pocket for spherical probe atoms",&ComStrHead,'-',ClusPocketType);
     }
   }
   /** [10]  Labeling Protein Residues by Multiscale MolVol and pocket  */
    if ((oresfile[0]!='\0')||(orpdbfile[0]!='\0')){
     Label_Multiscale_MolVol_or_Pocket_to_Protein_Residue_Shells(&ProResHead,&MapX,&MSprb,Wshell,'M');
     Label_Multiscale_MolVol_or_Pocket_to_Protein_Residue_Shells(&ProResHead,&MapY,&MSprb,Wshell,'P');
     if (ClusPocketType!='F') 
       Label_Multiscale_Each_Pocket_Cluster_to_Protein_Residue_Shells(&ProResHead,&MapY,&MapZ,&MapX,&MSprb,Wshell);
       /*** CAUTION !! after this line, the content of MapX is destroyed !! ***/
    } 
    if (orpdbfile[0]!='\0')
       Write_PDB_File_with_Multiscale_MolVol_and_Pocket(orpdbfile,'w',&ProAtomHead,'A',&MSprb,tFactorType,
                                  "Multiscale MolVol/Pocket for protein atoms",&ComStrHead,'R',ClusPocketType);

   /** [11] Output Residue Based file with pocket 1/Rinaccess values **/
    
    if (oresfile[0]!='\0'){
     Write_Residue_File_Assigned_Multiscale_MolVol_and_Pocket(oresfile,&ProResHead,&MSprb,ClusPocketType,
         "Pocket invRinacc value for residues",&ComStrHead); }

     if (orpdbfile[0]!='\0'){
      Write_PDB_File_with_Multiscale_MolVol_and_Pocket(orpdbfile,'w',&ProAtomHead,'A',&MSprb,tFactorType,
                                  "Pocket inRinacc value for protein molecule",&ComStrHead,'R',ClusPocketType);}


 } /* MODE=='M'  */


 /***************************************************/
 /** MODE 'G' or 'g': Comparison between two grids **/
 /***************************************************/

  if ((strcmp(MODE,"G")==0)||(strcmp(MODE,"g")==0)){

   if (igridfileA[0]!='\0'){ 
     Read_CHAR3DMAP_PDB_format(igridfileA,&MapX,chain_select_gridA,select_model_num_str);
    Output_CHAR3DMAP_PDB_format("outA.pdb",'w',&MapX,1,"igridfileA",&ComStrHead); 
   }
  
   if (igridfileB[0]!='\0'){ 
    Read_CHAR3DMAP_PDB_format(igridfileB,&MapY,'-',select_model_num_str); 
    Output_CHAR3DMAP_PDB_format("outB.pdb",'w',&MapY,1,"igridfileB",&ComStrHead); 
   }
 

   printf("#CorrCoeff %f\n",Corr_Coeff_of_Two_CHAR3DMAPs(&MapX,&MapY));
   if ((igridfileA[0]!='\0')&&(igridfileB[0]!='\0')){
    Output_Histogram_of_Two_CHAR3DMAPs("hist.out",&MapX,&MapY);
    Output_Two_CHAR3DMAPs_in_PDB_format("diffmap.pdb",&MapX,&MapY,1,255,"diffence map");

   }

 } 


 /***************************************************/
 /** MODE 'GcmpS' : Comparison of grid and spheres **/
 /***************************************************/

 if (strcmp(MODE,"GcmpS")==0){
   if (igridfile[0]!='\0'){ 
     Read_CHAR3DMAP_PDB_format(igridfile,&MapX,chain_select_grid,select_model_num_str);
     Output_CHAR3DMAP_PDB_format("pock_grid.pdb",'w',&MapX,1,"igridfile",&ComStrHead); 
   }
   else{
     printf("#ERROR(MODE 'GcmpS'): -igridpdb is required.\n");
     exit(1);
   }
   if (isphe_pdbfile[0]!='\0'){ 
     Read_PDB_File(isphe_pdbfile,&ConSpheAtomHead,pdbidPro,'-',"-",'B','F','F','M');
     assign_radius_by_Occup(&ConSpheAtomHead);
   }
   else{
     printf("#ERROR(MODE 'GcmpS'): -isphepdb is required.\n");
     exit(1);
   }
   Malloc_CHAR3DMAP(&MapY, MapX.N[0], MapX.N[1], MapX.N[2]);
   for (i=0;i<3;++i){
     MapY.N[i]       = MapX.N[i];
     MapY.OrigPos[i] = MapX.OrigPos[i];
   }
   MapY.grid_width = MapX.grid_width;
   Set_Bit_3DMAP_Inside_Atoms(&MapY,&ConSpheAtomHead,1,1.0);
   Output_CHAR3DMAP_PDB_format("sphe_grid.pdb",'w',&MapY,1,"igridfile",&ComStrHead); 
   Compare_Two_CHAR3DMAP(&MapX,1,255,&MapY,1,255, &NP,&NL,&NPL,&Rec,&Pre,&Fme,&Tani);
   voxel_volume = MapX.grid_width *MapX.grid_width *MapX.grid_width;
   printf("NP %d\n",NP);
   printf("VP %f\n",NP*voxel_volume);
   printf("NL %d\n",NL);
   printf("VL %f\n",NL*voxel_volume);
   printf("NPL %d\n",NPL);
   printf("VPL %f\n",NPL*voxel_volume);
   printf("RECALL %f\n",Rec);
   printf("PRECISION %f\n",Pre);
   printf("FMEASURE %f\n",Fme);
   printf("TANIMOTO %f\n",Tani);
 }



 /********************************************************/
 /** MODE 'GcmpL' : Comparison of grid and ligand atoms **/
 /********************************************************/

 if (strcmp(MODE,"GcmpL")==0){
   if (igridfile[0]!='\0'){ 
     Read_CHAR3DMAP_PDB_format(igridfile,&MapX,chain_select_grid,select_model_num_str);
     /* Output_CHAR3DMAP_PDB_format("pock_grid.pdb",'w',&MapX,1,"igridfile",&ComStrHead);  */
   }
   else{
     printf("#ERROR(MODE 'GcmpS'): -igridpdb is required.\n");
     exit(1);
   }

   if (iligpdbfile[0]!='\0'){ 
     Read_PDB_File(iligpdbfile,&LigAtomHead,pdbidLig,'H',"-",'B','F','T',MODELtype);
   }
   else if (iligciffile[0]!='\0'){ 
     Read_mmCIF_File(iligciffile,&LigAtomHead,'A','H',ChainID_Str, asym_id_Str, assembly_id,MaxAtomForAll);
   }
   else{
     printf("#ERROR(MODE 'GcmpL'): -isphepdb is required.\n");
     exit(1);
   }

   NatomPdb = Number_Of_ATOMs(&LigAtomHead); 
   printf("#iligRadiusType %c\n",iligRadiusType); 
   Assign_Radius_to_ATOMs(&LigAtomHead,iligRadiusType,Runified);
   MinMaxXYZ_Of_Atoms(&LigAtomHead,MinL,MaxL);
   margin_width = 1.87;  
   for (i=0;i<3;++i){
     MinP[i] = MapX.OrigPos[i];
     MaxP[i] = MapX.OrigPos[i] + MapX.grid_width * MapX.N[i];
     MinL[i] = MinL[i] - margin_width;
     MaxL[i] = MaxL[i] + margin_width;
   }
   printf("#MinP %f %f %f MaxP %f %f %f\n",MinP[0],MinP[1],MinP[2],MaxP[0],MaxP[1],MaxP[2]);
   printf("#MinL %f %f %f MaxL %f %f %f\n",MinL[0],MinL[1],MinL[2],MaxL[0],MaxL[1],MaxL[2]);

   Min_Vec(MinA,MinL,MinP);
   Max_Vec(MaxA,MaxL,MaxP);
   printf("#MinA %f %f %f MaxA %f %f %f\n",MinA[0],MinA[1],MinA[2],MaxA[0],MaxA[1],MaxA[2]);

   MapY.grid_width = MapX.grid_width;
 
   for (i=0;i<3;++i){
     MapY.OrigPos[i] = MapX.grid_width * (int)floor((MinA[i]-MinP[i])/MapX.grid_width) + MinP[i]; 
     printf("MapY %d OrigPos %f MinA %f MinP %f\n",i,MapY.OrigPos[i],MinA[i],MinP[i]);
     MapY.N[i] = (int)ceil((MaxA[i]-MapY.OrigPos[i])/MapX.grid_width); 
   }

   printf("#MapY.OrigPos %f %f %f N %d %d %d\n", MapY.OrigPos[0], MapY.OrigPos[1], MapY.OrigPos[2], MapY.N[0], MapY.N[1], MapY.N[2]);

   Malloc_CHAR3DMAP(&MapY, MapY.N[0], MapY.N[1], MapY.N[2]);
 
/* 
   Write_PDB_File("lig.pdb",'w',&LigAtomHead,"",&ComStrHead);
*/
/* 
   Output_CHAR3DMAP_PDB_format("lig_grid.pdb",'w',&MapY,1,"igridfile",&ComStrHead); 
 */
   MapY.grid_width = MapX.grid_width;

   Set_Bit_By_Another_Map_in_CHAR3DMAP(&MapY,&MapX,1,0);
   Set_Bit_3DMAP_Inside_Atoms(&MapY,&LigAtomHead,2,1.0);

   Compare_Two_Bit_in_CHAR3DMAP(&MapY,1,2, &NP,&NL,&NPL,&Rec,&Pre,&Fme,&Tani);

/*
   Output_CHAR3DMAP_with_tar_bitmask_PDB_format("pocgrid.pdb",&MapY,1);
   Output_CHAR3DMAP_with_tar_bitmask_PDB_format("liggrid.pdb",&MapY,2);
*/

/*
 Compare_Two_CHAR3DMAP(&MapX,1,255,&MapY,1,255, &NP,&NL,&NPL,&Rec,&Pre,&Fme,&Tani);
 */
   voxel_volume = MapX.grid_width *MapX.grid_width *MapX.grid_width;
   printf("NATOM_LIG %d\n",NatomPdb);
   printf("NP %d\n",NP);
   printf("VP %f\n",NP*voxel_volume);
   printf("NL %d\n",NL);
   printf("VL %f\n",NL*voxel_volume);
   printf("NPL %d\n",NPL);
   printf("VPL %f\n",NPL*voxel_volume);
   printf("RECALL %f\n",Rec);
   printf("PRECISION %f\n",Pre);
   printf("FMEASURE %f\n",Fme);
   printf("TANIMOTO %f\n",Tani);
 }



 /********************************************************/
 /** MODE 'GcmpG' : Comparison of grid and grid **/
 /********************************************************/
 if (strcmp(MODE,"GcmpG")==0){
   if (igridfileA[0]!='\0'){ 
     Read_CHAR3DMAP_PDB_format(igridfileA,&MapX,chain_select_grid,select_model_num_str);
     Output_CHAR3DMAP_PDB_format("gridA.pdb",'w',&MapX,1,"igridfileA",&ComStrHead); 
   }
   else{
     printf("#ERROR(MODE 'GcmpG'): -igridpdbA is required.\n");
     exit(1);
   }

   if (igridfileB[0]!='\0'){ 
     Read_CHAR3DMAP_PDB_format(igridfileB,&MapY,chain_select_grid,select_model_num_str);
     Output_CHAR3DMAP_PDB_format("gridB.pdb",'w',&MapY,1,"igridfileB",&ComStrHead); 
   }
   else{
     printf("#ERROR(MODE 'GcmpG'): -igridpdbB is required.\n");
     exit(1);
   }
   if ( MapX.grid_width != MapY.grid_width){
     printf("#ERROR:grid_widthA %f is different from grid_widthB %f\n",MapX.grid_width,MapY.grid_width);
     exit(1);
   }
   if ((MapX.OrigPos[0]!= MapY.OrigPos[0])||(MapX.OrigPos[1]!= MapY.OrigPos[1])||(MapX.OrigPos[2]!= MapY.OrigPos[2])){
     printf("#ERROR:OrigPosA %f %f %f is different from OrigPosB %f %f %f \n",
                MapX.OrigPos[0], MapX.OrigPos[1], MapX.OrigPos[2], MapY.OrigPos[0], MapY.OrigPos[1], MapY.OrigPos[2]);
     exit(1);
   }

   Compare_Two_CHAR3DMAP(&MapX,1,255,&MapY,1,255, &NP,&NL,&NPL,&Rec,&Pre,&Fme,&Tani);
   voxel_volume = MapX.grid_width *MapX.grid_width *MapX.grid_width;
   printf("NP %d\n",NP);
   printf("VP %f\n",NP*voxel_volume);
   printf("NL %d\n",NL);
   printf("VL %f\n",NL*voxel_volume);
   printf("NPL %d\n",NPL);
   printf("VPL %f\n",NPL*voxel_volume);
   printf("RECALL %f\n",Rec);
   printf("PRECISION %f\n",Pre);
   printf("FMEASURE %f\n",Fme);
   printf("TANIMOTO %f\n",Tani);

 }

 /************************************************/
 /** MODE 'GtoV' : convert a grid to volume_map **/
 /************************************************/

 if (strcmp(MODE,"GtoV")==0){
   if (igridfile[0]!='\0'){ 
     Read_CHAR3DMAP_PDB_format(igridfile,&MapX,chain_select_grid,select_model_num_str);
     Output_CHAR3DMAP_PDB_format("pock_grid.pdb",'w',&MapX,1,"igridfile",&ComStrHead); 
   }
   else{
     printf("#ERROR(MODE 'GcmpS'): -igridpdb is required.\n");
     exit(1);
   }
   if (omapfile[0] != '\0'){
     Write_MapCCP4_File_CHAR3DMAP(omapfile,&MapX,'D');
   }
 }

 /*******************************************/
 /** MODE 'Gsurf' : surfacization of grid  **/
 /**   X - (X erosion P)                   **/
 /*******************************************/

 if (strcmp(MODE,"Gsurf")==0){
    printf("# MODE 'Gsurf' : surfacization of grid\n");
    printf("#    X - (X erosion P)                \n");

   if (igridfile[0]!='\0'){ 
     Read_CHAR3DMAP_PDB_format(igridfile,&MapX,chain_select_grid,select_model_num_str);
     Extract_NonZero_3DMAP(&MapX,1);
     Output_CHAR3DMAP_PDB_format("pock_grid.pdb",'w',&MapX,1,"igridfile",&ComStrHead); 
   }
   else{
     printf("#ERROR(MODE 'Gsurf'): -igridpdb and -ogridpdb are required.\n");
     exit(1);
   }
    MapP.grid_width = MapX.grid_width;
    printf("#Probe_radius %f\n",MapP.grid_width);
    Ngrid_probe = Make_Sphere_Probe_Map(&MapP,MapP.grid_width,MapP.grid_width);
    /* Output_CHAR3DMAP_PDB_format("probe.pdb",'w',&MapP,1,"Probe Sphere",&ComStrHead); */
    Erosion(&MapX,&MapP,1,2);
    Extract_A_and_notB_3DMAP(&MapX,1,2,4);
    Extract_Specific_3DMAP(&MapX,4); 

    if (ogridfile[0] != '\0'){
      Output_CHAR3DMAP_PDB_format(ogridfile,'w',&MapX,4,"original(1) and Erosion(4)  protein",&ComStrHead); 
    }
 }





 /******************************************************************************/
 /* MODE 'I'  Interface Between Two Chains                                     */
 /*  InterfacePocket(A,B,S)  = [(A or B) cls P - (A cls P) - (B cls P)] opn S  */
 /******************************************************************************/

 if (strcmp(MODE,"I")==0){
  if (ipdbfileA[0]=='\0') {printf("#ERROR(MODE 'I'): no option '-ipA'.\n"); exit(1);} 
  if (ipdbfileB[0]=='\0') {printf("#ERROR(MODE 'I'): no option '-ipB'.\n"); exit(1);} 

  /** (1) Read ipdbfile A **/ 
  Read_PDB_File(ipdbfileA,&ProAtomHeadA,pdbidProA,'H',ChainID_StrA,ReadAtmHetType,'F','F',MODELtype);
  NatomPdbA = Number_Of_ATOMs(&ProAtomHeadA); 
  Assign_Radius_to_ATOMs(&ProAtomHeadA,vdwRadiusType,Runified);
  MinMaxXYZ_Of_Atoms(&ProAtomHeadA,MinA,MaxA); 
  sprintf(comment,"#ipdbfileA '%s' Natom %d\n",ipdbfileA,NatomPdbA);  
  Add_String_to_STRING_NODE(&ComStrHead,comment); 
  /** (2) Read ipdbfile B **/ 
  Read_PDB_File(ipdbfileB,&ProAtomHeadB,pdbidProB,'H',ChainID_StrB,ReadAtmHetType,'F','F',MODELtype);
  NatomPdbB = Number_Of_ATOMs(&ProAtomHeadB); 
  Assign_Radius_to_ATOMs(&ProAtomHeadB,vdwRadiusType,Runified);
  MinMaxXYZ_Of_Atoms(&ProAtomHeadB,MinB,MaxB); 
  sprintf(comment,"#ipdbfileB '%s' Natom %d\n",ipdbfileB,NatomPdbB);  
  Add_String_to_STRING_NODE(&ComStrHead,comment); 

  for (k=0;k<3;++k){
    if (MinA[k]<MinB[k]) MinP[k] = MinA[k]; else MinP[k] = MinB[k];
    if (MaxA[k]>MaxB[k]) MaxP[k] = MaxA[k]; else MaxP[k] = MaxB[k];
   }   

  /** (3) Malloc Grids **/ 
 
  margin_width = 2.0*RprobeL; 
  margin_width += MapX.grid_width*2;
  Decide_Num_Grid_3D(MapX.N,MapX.OrigPos,MapX.grid_width, margin_width, MinP,MaxP,'-');
  Malloc_CHAR3DMAP(&MapX, MapX.N[0], MapX.N[1], MapX.N[2]);
  Vvoxel = MapX.grid_width * MapX.grid_width * MapX.grid_width;

  /** (4) Set grids 1 for ipdbA, 2 for ipdbB, 4 for A or B **/ 
  Ngrid_VdW_A = Set_Bit_3DMAP_Inside_Atoms(&MapX,&ProAtomHeadA,1,1.0);

  sprintf(comment,"VdW_A    Natom %d Ngrid %d Vol %.2f AAA",
             NatomPdbA,Ngrid_VdW_A,Ngrid_VdW_A*Vvoxel);
  Add_String_to_STRING_NODE(&ComStrHead,comment); 

  Ngrid_VdW_B = Set_Bit_3DMAP_Inside_Atoms(&MapX,&ProAtomHeadB,2,1.0);
  
  sprintf(comment,"VdW_B    Natom %d Ngrid %d Vol %.2f AAA",
        NatomPdbB,Ngrid_VdW_B,Ngrid_VdW_B*Vvoxel);
  Add_String_to_STRING_NODE(&ComStrHead,comment); 
  
  Ngrid = Extract_A_or_B_3DMAP(&MapX,1,2,4);
  sprintf(comment,"VdW_AorB Natom %d Ngrid %d Vol %.2f AAA",
       NatomPdbA+NatomPdbB,Ngrid,Ngrid*Vvoxel);
  Add_String_to_STRING_NODE(&ComStrHead,comment); 
 

  /** (5) Make Probe Grid for sphere with RprobeL **/
  Ngrid_probe = Make_Sphere_Probe_Map(&MapPL,RprobeL,MapP.grid_width);
  /*
  Output_CHAR3DMAP_PDB_format("probeL.gdb",'w',&MapPL,1,"ProbeL Sphere",&ComStrHead);
  */
  sprintf(comment,"RprobeL %f Ngrid_probe %d",RprobeL,Ngrid_probe);
  Add_String_to_STRING_NODE(&ComStrHead,comment); 
  if (RprobeS>0.0){
    Ngrid_probeS = Make_Sphere_Probe_Map(&MapPS,RprobeS,MapPS.grid_width);
    sprintf(comment,"RprobeS %f Ngrid_probe %d",RprobeS,Ngrid_probeS);
    Add_String_to_STRING_NODE(&ComStrHead,comment); 
  } 
 /** (6) Molecular Volume for A or B **/

  printf("#### Calcualte Molecular Volume for A or B ####\n"); 
  Dilation(&MapX,&MapPL,4,8,'O');
  Ngrid = Erosion(&MapX,&MapPL,8,16);
  sprintf(comment,"MolVol_AorB        : Ngrid %d Vol %.2f AAA",Ngrid,Ngrid*Vvoxel);
  Add_String_to_STRING_NODE(&ComStrHead,comment); 

  /** (7) Molecular Volume for A **/
  printf("#### Calcualte Molecular Volume for A  ####\n"); 
  Reset_Specific_3DMAP(&MapX,8);
  Dilation(&MapX,&MapPL,1,8,'O');
  Ngrid = Erosion(&MapX,&MapPL,8,32);
  sprintf(comment,"MolVol_A           : Ngrid %d Vol %.2f AAA",Ngrid,Ngrid*Vvoxel);
  Add_String_to_STRING_NODE(&ComStrHead,comment); 
  
  /** (8) Molecular Volume for B **/
  printf("#### Calcualte Molecular Volume for B  ####\n"); 
  Reset_Specific_3DMAP(&MapX,8);
  Dilation(&MapX,&MapPL,2,8,'O');
  Ngrid = Erosion(&MapX,&MapPL,8,64);
  sprintf(comment,"MolVol_B           : Ngrid %d Vol %.2f AAA",Ngrid,Ngrid*Vvoxel);
  Add_String_to_STRING_NODE(&ComStrHead,comment); 
  
  Reset_Specific_3DMAP(&MapX,8);

  sprintf(comment,"BIT  1  : ipdbfileA"); Add_String_to_STRING_NODE(&ComStrHead,comment); 
  sprintf(comment,"BIT  2  : ipdbfileB"); Add_String_to_STRING_NODE(&ComStrHead,comment); 
  sprintf(comment,"BIT  4  : ipdbfileA or ipdbfileB"); Add_String_to_STRING_NODE(&ComStrHead,comment); 
  sprintf(comment,"BIT 16  : MolVol for A or B "); Add_String_to_STRING_NODE(&ComStrHead,comment); 
  sprintf(comment,"BIT 32  : MolVol for A"); Add_String_to_STRING_NODE(&ComStrHead,comment); 
  sprintf(comment,"BIT 64  : MolVol for B"); Add_String_to_STRING_NODE(&ComStrHead,comment); 
 
  /** (9) Set MolVol(A or B) - MolVol(A) - MolVol(B)  **/
  Extract_A_and_notB_and_notC_3DMAP(&MapX,16,32,64,128);
  
  /** (9) Opening (MolVol(A or B) - MolVol(A) - MolVol(B)) opn probeS  **/
  if (RprobeS > 0.0){
    printf("#### Calcualte Opening by probeS   ####\n"); 
    Reset_Specific_3DMAP(&MapX,8);
    Erosion(&MapX,&MapPS,128,8);
    Reset_Specific_3DMAP(&MapX,128);
    Dilation(&MapX,&MapPS,8,128,'-');
    sprintf(comment,"BIT 128 : [MolVol(AorB)-MolVol(A)-MolVol(B)] opn probeS"); Add_String_to_STRING_NODE(&ComStrHead,comment); 
  } 
  else{
    sprintf(comment,"BIT 128 : MolVol(AorB)-MolVol(A)-MolVol(B)"); Add_String_to_STRING_NODE(&ComStrHead,comment); 
  }

  Ngrid_gapAB = Count_Specific_3DMAP(&MapX,128);

  sprintf(comment,"Interface pocket Ngrid %d Vol %.2f AAA", Ngrid_gapAB,Ngrid_gapAB*Vvoxel);
  Add_String_to_STRING_NODE(&ComStrHead,comment); 

  Natom_interA =  Assign_Map_Value_to_tFactor_of_Atoms(&ProAtomHeadA,&MapX,'B',128);
  Natom_interB =  Assign_Map_Value_to_tFactor_of_Atoms(&ProAtomHeadB,&MapX,'B',128);
  sprintf(comment,"Natom_inter A %d B %d interface_pocket_per_Natom_inter A %.3f B %.3f",
    Natom_interA, Natom_interB, (float)Ngrid_gapAB*Vvoxel/(float)Natom_interA,(float)Ngrid_gapAB*Vvoxel/(float)Natom_interB);
  Add_String_to_STRING_NODE(&ComStrHead,comment); 

  Show_STRING_NODE(&ComStrHead); 

  if (opockpdbfile[0]=='\0') sprintf(opockpdbfile,"inpock.gpdb");

  if (opockpdbfile[0]!='\0')
   Output_CHAR3DMAP_PDB_format(opockpdbfile,'w',&MapX,128,"Interface pocket between two chains",&ComStrHead);  
   
  if (opdbfile[0]!='\0'){
    sprintf(comment,"The ATOM with more than zero tFactor contacts with the interface volume.");
    Add_String_to_STRING_NODE(&ComStrHead,comment); 
    Write_PDB_File(opdbfile,'w',&ProAtomHeadA,"",&ComStrHead);
    Write_PDB_File(opdbfile,'a',&ProAtomHeadB,"",&ComStrHead);
  }
 } /* MODE 'I' */
 



 /************************************************/ 
 /** MODE 'L' : Ligand vs grid comparison      **/
 /***********************************************/ 
 
 if (strcmp(MODE,"L")==0){
   if (igridfileA[0]=='\0')  {printf("#ERROR:MODE 'L' requires -igA option"); exit(1);}
   if (iligpdbfile[0]=='\0') {printf("#ERROR:MODE 'L' requires -ilg option"); exit(1);}

   /** (1) Read grid file  **/
   Read_CHAR3DMAP_PDB_format(igridfileA,&MapX,chain_select_gridA,select_model_num_str); 
   Output_CHAR3DMAP_PDB_format_Min_Max("init_grid.pdb",'w',&MapX,1,255,"pocket grid",&ComStrHead); 

   for (k=0;k<3;++k){
    MinA[k] = MapX.OrigPos[k];
    MaxA[k] = MapX.OrigPos[k] + MapX.grid_width * MapX.N[k];
   }

    if ((threCluster_char_max>0) ||(threCluster_char_min>0)){
      for (k=0;k<3;++k) { MapZ.N[k] = MapX.N[k]; MapZ.OrigPos[k] = MapZ.OrigPos[k]; }
      MapZ.grid_width = MapX.grid_width; 
      Malloc_CHAR3DMAP(&MapZ, MapZ.N[0], MapZ.N[1], MapZ.N[2]);
    }

    if (threCluster_char_max>0) {
      Cluster_MultiScale_Pocket_Single_Linkage(&MapX,&MapZ,&Gclus,threCluster_char_max,'X',&MSprb,'N');
    }
    if (threCluster_char_min>0){
       Cluster_MultiScale_Pocket_Single_Linkage(&MapX,&MapZ,&Gclus,threCluster_char_min,'I',&MSprb,'N',PAR.SubType[0]);
    }

   if (RankKeepCluster>0){
     Set_Map_to_Zero_for_Low_Ranked_Grid_Cluster(&MapX,&Gclus,RankKeepCluster);
   }

   if (opockpdbfile[0]!='\0'){
      Output_CHAR3DMAP_PDB_format_Min_Max(opockpdbfile,'w',&MapX,1,255,"X.gdb",&ComStrHead); 
   }

   /** (2) Read Ligand file **/
   Read_PDB_File(iligpdbfile,&LigAtomHead,pdbidLig,'H',"-",'B','F','T',MODELtype);
   NatomPdb = Number_Of_ATOMs(&LigAtomHead); 
   printf("#NatomLig %d\n",NatomPdb);
   printf("#iligRadiusType %c\n",iligRadiusType); 
   Assign_Radius_to_ATOMs(&LigAtomHead,iligRadiusType,Runified);
   MinMaxXYZ_Of_Atoms(&LigAtomHead,MinP,MaxP); 
   printf("#Ligand Min  %8.3f %8.3f %8.3f Max %8.3f %8.3f %8.3f\n", MinP[0], MinP[1], MinP[2], MaxP[0], MaxP[1], MaxP[2]);
   printf("#Grid   Min  %8.3f %8.3f %8.3f Max %8.3f %8.3f %8.3f\n", MinA[0], MinA[1], MinA[2], MaxA[0], MaxA[1], MaxA[2]);

  if ((MinP[0]<MinA[0])|| (MaxP[0]>MaxA[0]) || (MinP[1]<MinA[1])|| (MaxP[1]>MaxA[1]) || 
      (MinP[2]<MinA[2])|| (MaxP[2]>MaxA[2]) ){
     printf("#ERROR:Ligand is over the grid igA\n");
     Resize_CHAR3DMAP_by_MinMax(&MapX,MinP,MaxP,MapY.OrigPos,MapY.N,Nshift);
     MapY.grid_width = MapX.grid_width; 
     Malloc_CHAR3DMAP(&MapY, MapY.N[0], MapY.N[1], MapY.N[2]);
     Copy_CHAR3DMAP_with_Nshift(&MapY,&MapX,Nshift);
    /*
     Output_CHAR3DMAP_PDB_format_Min_Max("X.gdb",'w',&MapX,1,255,"X.gdb",&ComStrHead); 
     Output_CHAR3DMAP_PDB_format_Min_Max("Y.gdb",'w',&MapY,1,255,"Y.gdb",&ComStrHead); 
    */
     Free_CHAR3DMAP(&MapX);
     Malloc_CHAR3DMAP(&MapX, MapY.N[0], MapY.N[1], MapY.N[2]);
     for (k=0;k<3;++k){ 
       MapX.N[k] = MapY.N[k]; 
       MapX.OrigPos[k] = MapY.OrigPos[k]; 
     }
     MapX.grid_width = MapY.grid_width; 
     Copy_CHAR3DMAP(&MapX,&MapY);
     Initialize_CHAR3DMAP(&MapY);

    } 
  else{
    for (k=0;k<3;++k){ 
      MapY.N[k] = MapX.N[k]; 
      MapY.OrigPos[k] = MapX.OrigPos[k]; 
     }
     MapY.grid_width = MapX.grid_width; 
     Malloc_CHAR3DMAP(&MapY, MapY.N[0], MapY.N[1], MapY.N[2]);
  }

   /** (3) Malloc MapY for ligand **/
   Ngrid_VdW = Set_Bit_3DMAP_Inside_Atoms(&MapY,&LigAtomHead,1,1.0);
   Output_CHAR3DMAP_PDB_format_Min_Max("pock_grid.pdb",'w',&MapX,1,255,"pocket grid",&ComStrHead); 
   Output_CHAR3DMAP_PDB_format_Min_Max("lig_grid.pdb",'w', &MapY,1,255,"ligand grid",&ComStrHead); 
   Compare_Two_CHAR3DMAP(&MapX,thre_min_gridA,thre_max_gridA,&MapY,1,1, &NP,&NL,&NPL,&Rec,&Pre,&Fme,&Tani);
   printf("GRID_LIG_COMP %s %s threA %d %d NP %4d NL %4d NPL %4d Rec %f Pre %f Fme %f Tani %f\n",
   igridfileA,iligpdbfile,thre_min_gridA, thre_max_gridA,NP,NL,NPL,Rec,Pre,Fme,Tani);
   if (opockpdbfile[0]!='\0')
     Output_CHAR3DMAP_PDB_format_Min_Max(opockpdbfile,'w',&MapX,1,255,"-igA map",&ComStrHead); 
   if (oligpdbfile[0]!='\0')
     Output_CHAR3DMAP_PDB_format_Min_Max(oligpdbfile,'w',&MapY,1,255,"-ilg map",&ComStrHead); 
 }


 /*****************************************/ 
 /** MODE 'ConSphe' : Contacting Spheres **/
 /*****************************************/ 
 if (strcmp(MODE,"ConSphe")==0){
   if (ipdbfile[0]!='\0'){
     Read_PDB_File(ipdbfile,&ProAtomHead,pdbidPro,'H',ChainID_Str,ReadAtmHetType,'F','F',MODELtype);
     NatomPdb = Number_Of_ATOMs(&ProAtomHead); 
     Assign_Radius_to_ATOMs(&ProAtomHead,vdwRadiusType,Runified);
     Assign_Weight_to_ATOMs(&ProAtomHead);
     Make_Residue_List(&ProAtomHead,&ProResHead);
     Cal_Mean_and_CovM_from_ATOMs(&ProAtomHead,M,CovM);
     Setup_PCvar_and_PCaxis(CovM,PCvar,PCaxis);

     if (icen_pdbfile[0] != '\0'){
       Read_PDB_File(icen_pdbfile,&ConSpheAtomHead,pdbidPro,'-',"-",'B','F','F','M');
     }
     else if ((PCnum_ConSphe==0) ||(PCnum_ConSphe==1) ||(PCnum_ConSphe==2)){
       Make_Sphere_Centers_from_PCaxis(&ConSpheAtomHead,&ProAtomHead,Nconspheres,M,PCaxis[PCnum_ConSphe],sqrt(PCvar[PCnum_ConSphe]));
     }
     else if (PCnum_ConSphe <0){ 
       Make_One_Sphere_Center_on_Gravity_Center(&ConSpheAtomHead,&ProAtomHead,M);
     }


     Setup_Radius_of_Spheres_Tangent_to_ATOMs(&ConSpheAtomHead,&ProAtomHead);

     if (ExcludeEscapeConSphe == 'T'){
       Exclude_Spheres_Escape_to_Outside(&ConSpheAtomHead, &ProAtomHead, PCaxis,PCnum_ConSphe);
     }
     if (Rmin_ConSphe>0.0){
       Remove_ATOM_less_than_Rmin_from_ATOMs(&ConSpheAtomHead,Rmin_ConSphe);
     }

     if (osphe_pdbfile[0] != '\0'){
       Write_PDB_File(osphe_pdbfile,'w',&ConSpheAtomHead,"",&ComStrHead);
     }
     if (osphe_wrlfile[0] != '\0'){
       write_ATOMs_as_Spheres_in_VRML(osphe_wrlfile,&ConSpheAtomHead,ConSphe_RGBTstr);
     }
   }
   else{
     printf("#ERROR: pdbfile is missing\n");
     exit(1);
   }

 }
 return(1);
}/* end of main() */
