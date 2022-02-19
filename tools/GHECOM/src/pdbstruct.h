/*

<pdbstructs.h>

Definition of basic structures

*/


#define MAX_CLUSTER_NUM  10


struct NEIGHBOR{
 struct ATOM *atom;
 struct NEIGHBOR *next,*prev;
};




/** ATOM : for both protein atoms and probe spheres */
struct ATOM{
 float Pos[3];        /* X,Y,Z coordinates */
 float R;             /* Radius            */ 
 float RR;            /* Radius^2          */ 
 float Weight;        /* Weight            */ 
 int   num;           /* Atom Number */
 int   rnum;          /* Residue Number */
 char  Anum[7];       /* Atom Number String */
 char  Atom[5];       /* Name of Atom */
 char  Resi[5];       /* Name of Residue */
 char  Rnum[6];       /* Residue Number String */
 char  Chain;         /* Chain ID */
 char  altLoc;        /* Alternate location indicator */
 float Occup;         /* Occupancy */ 
 float tFactor;       /* Temparature Factor */ 
 char  AHtype;        /* 'A'tom, 'H'etatm */
 char  element[3];    /* chemical element, such as 'C','O','Zn','Fe' */
 int   model_num;      /* model number */
 struct ATOM *next,*prev;    /* Pointer for atom list (bi-direction) */
 struct RESIDUE *res;        /* Pointer to corresponding residues */
 struct ATOM *rnext,*rprev;  /* Pointer for the list from RESIDUE (bi-direction)*/
 char   region;
 char   mark;           /* Mark for various purpose */

 /* For Multiscape MolVol and Pocket */
 unsigned char  charRinacc;  /* index of Rinaccess value assigned MULSC_SHELL_PROBE->Rarray[] */
   /* 
    for MolVol 
      Rinacc is  maximum value among the MolVol grids around the atom.
    for Pocket 
      Rinacc is  minimum value among the Pocket grids around the atom 
   */
 float           shellAcc;     /* Accessibility of the shell (ratio[0:1]) */
 float           Rinacc;       /* Rinaccess value among the MolVol grids around the atom */
 float           rw_shellAcc;  /* Rinacc-weighted Accessibility of the shell (ratio[0:1])    */
 float           pocketness;   /* pocketness = average_invRinacc * RprobeS */
 float  pocketness_clus[MAX_CLUSTER_NUM];  /* pocketness for each pocket cluster */
 int    cluster_num;                       /* cluster number of pocket */ 
 char   contact;                /* 1:contact with ligand, 0:not contact*/
 struct NEIGHBOR  NeiHead;      /* Head of NEIGHBOR atoms */
 struct ATOM *conatm1,*conatm2,*conatm3;  /* contacting three atoms (only for spherical probes) */
};



/** RESIDUE : for both protein residues and probe clusters */
struct RESIDUE{
 int    num;           /* Residue Number */
 char   Resi[5];       /* Name of Residue */
 char   Rnum[6];       /* Residue Number String */
 char   Chain;         /* Chain ID */
 int    Natom;         /* Number of Atoms */ 
 float  Pos[3];        /* Position of Center of the residue */ 
 char   *con_clus;     /* [Ncluster] : String of Contacted Cluster Number (malloc later) */       
 struct RESIDUE *next,*prev; /* Pointer for RESIDUE list (bi-direction) */
 struct ATOM    Ahead; /* Head of ATOM list, belonged to the residue (using ATOM->rnext, rprev) */
 struct ATOM    *rep_atom;   /* Pointer to representative atoms, such as 'Calpha' */
 float  value;               /* Value for various purpose such as Sorting */  
 float  Rinacc;              /* averaged Rinaccess value among the MolVol grids around the residue */
 float  shellAcc;            /* Accessibility of the shell */
 float  rw_shellAcc;        /* Rinacc-weighted Accessibility of the shell (ratio[0:1])    */
 float  pocketness;         /* pocketness = average_invRinacc * RprobeS */
 float  pocketness_clus[MAX_CLUSTER_NUM];  /* pocketness for each pocket cluster */
 int    Natom_contact;     /* Number of contacted atoms with ligand */
};




/** STRING_NODE : for comment strings of output files **/
struct STRING_NODE{
  char   line[MAX_STRING_NODE_LINE]; 
  struct STRING_NODE *next,*prev;
};


