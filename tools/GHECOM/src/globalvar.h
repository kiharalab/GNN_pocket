/*

<globalvar.h>

Definition of basic structures

*/

#define MAX_UNSIGNED_CHAR     256
#define MAX_FILE_NAME         1024
#define MAX_STRING_NODE_LINE  100

#define LAST_MOD_DATE     "2021/12/01"
#define CODED_BY          "coded by Takeshi Kawabata."

struct OPTION_VALUE{
  char *opt;   /* key string   of option */
  char *val;   /* value string of option */
  int  Lopt;   /* Length of opt[] string */
  int  num;
  char accept; /* accept = 1, this is already translated */
  struct OPTION_VALUE *next,*prev;
};



/** PARAMETERS: for globally defined important parameters **/
struct PARAMETERS{
 char   COMMAND[MAX_FILE_NAME];
 char   HOSTNAME[MAX_FILE_NAME];
 struct OPTION_VALUE OptValHead;  
 char   pdbidPro[8];
 char   radfile[MAX_FILE_NAME];  /* radius file */

 /** Properties for input molecules **/
 int   N_RESIDUE;    /* number of RESIDUEs */
 int   N_ATOM;       /* number of ATOMs */
 int   N_MODEL;      /* number of MODELs */
 int   N_CHAIN;      /* number of CHAINs */
 float VOL_VDW;      /* VdW volume (A^3) */

 /* Others */
 float  MAX_MEMORY_MEGA;  /* Maximum memory size (MEGA BYTE) for one voxel */
 char   SubType[64];    /* Type for Various Purpuse */
 char   Vtype;          /* Visible for progress of calculation */ 
 int    NeighborNum;    /* 6 or 18 or 26 */
 int    NeighX[26],NeighY[26],NeighZ[26];   /* for neighbor voxels */
 float  grid_width;    
 char   OutCalProcess;  
 int    Nmax_cluster;   /* Number of maximum cluster number */ 

 char   RecursiveQueue;  /* 'R'ecursive call or 'Q'ueur for 'Flood-fill' or 'labeling'  */

 char   Change_N_CA_C_O_Hetatom_to_Atom;  /* 'T' or 'F' */


 /** Date and Computation time **/
 char   START_DATE[MAX_FILE_NAME];  /* Start Date */
 char   END_DATE[MAX_FILE_NAME];    /* End   Date */
 double START_TIME_SEC;   /* Date on starting time ( in second) */
 double COMP_TIME_SEC_PROGRESS;  /* Computation time in progress (in second); */
 double COMP_TIME_SEC;           /* Final Computation time (in second); */


};


/** GLOBAL VARIABLES **/

extern struct PARAMETERS PAR;

