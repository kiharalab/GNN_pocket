/*

 <MscPockClus.h>

*/


/** for member grids of GRID_CLUSTER **/
struct GRID_MEMBER{
/*
 int    x,y,z;
 */
 short  x,y,z;
 struct GRID_MEMBER *next;
 /*
 struct GRID_MEMBER *next,*prev;
 char   value;  char value for map[x][y][z]  
 char   type;  s'K'eleton, 'D'ilated voxel 
 */
};

/** for Cluster of grids **/
struct GRID_CLUSTER{
 int    num;         /* ranges [1..Ncluster] */
 int    Nmember;     /* Number of Member Grids */ 
 float  value;       /* value for sorting */
 struct GRID_MEMBER  HeadGridMember; /* Head Pointer for Grid Member */
 struct GRID_CLUSTER *next,*prev;
 /* Various Properties for GRID_CLUSTER */
 float  Pos[3];          /* Position of center of GRID_CLUSTER */ 
 int    NRinacc[256];    /* Histogram of NRinacc */  
 float  aveRinacc;       /* Harmonic Average value for aveRinacc */
 float  minRinacc;       /* minimum value for aveRinacc */
 float  maxRinacc;       /* maximum value for aveRinacc */
 float  Volume;          /* Volume(AAA) */
 float  invRvolume;      /* 1/Rinacc-weighted Volume (1/A * AAA = AA)*/
 char   flag;            /* Flag for various purpose */
};



/*** FUNCTIONS (GLOBAL) ***/
extern void Cluster_MultiScale_Pocket_Single_Linkage();
extern void Sort_GridCluster();
extern void Output_CHAR3DMAP_PDB_format_By_Grid_Cluster();
extern void Set_Map_to_Zero_for_Low_Ranked_Grid_Cluster();
extern void Label_Pocket_ClusterNum_to_Protein_Atom_Shells();
extern void Label_Pocket_ClusterNum_to_Ligand_Atoms();
extern void Add_Grid_Member();
extern void Free_Grid_Member();
extern int Check_Overlap_of_Grid_Member();
extern int Number_Of_Grid_Member();
