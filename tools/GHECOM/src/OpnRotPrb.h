/*

 <OpnRotPrb.h>

*/

struct CHAR4DMAP{
 int N[4];
 unsigned char ****map; /* [N[0]][N[1]][N[2]][N[3]] */
};


struct GRIDPNT3D{
 int x,y,z;
 int  Nprobe;      
 struct GRIDPNT3D *next,*prev;
 char mark;
 int  value; 
};



/*** FUNCTIONS (GLOBAL) ***/
extern void Make_CHAR4DMAP_from_ProbeAtoms_and_RmatArray();
extern void Make_GRIDPNT3D_List_from_CHAR4DMAP_Multiple_Probes();
extern void Dilation_CHAR4DMAP();
extern int Opening_Using_CHAR4DMAP_and_GRIDPNT3D();
extern int Opening_Using_CHAR4DMAP_and_GRIDPNT3D_with_Fitness();
extern void Mark_Minimum_Filtering_GRIDPNT3Ds();
extern void Make_Another_GRID3DPNT_List_from_Marked_List();
extern int NonAllZero_Region_GRID3DPNT_Listed_Probes();
extern int Opening_for_Focused_Region();
extern int Malloc_CHAR4DMAP();
extern void Make_Rod_Probe();
extern void Make_Disc_Probe();
extern void Make_CHAR4DMAP_for_Rotated_RodDisc_by_RmatArray();
