/*

 <SphePnt3.h>

*/

/*** FUNCTION (GLOBAL) ***/
extern int  Cal_Sphere_On_Three_Atoms();
extern void Make_Probe_Spheres_Three_Contact();
extern void Make_Probe_Spheres_Three_Contact_NEIGHBOR_list();
extern int  Check_Crash_Pos_and_Atoms();
extern void Malloc_FLOAT2DMAP();
extern void Free_FLOAT2DMAP();
extern void Cal_Distance_FLOAT2DMAP();
extern void Remove_Probe_Atoms_Larger_Specified_Rinaccess();
extern void Remove_Probe_Atoms_Smaller_Specified_Pocketness();
extern void Remove_Probe_Atoms_Smaller_Specified_tFactor();
extern void Sort_Probe_By_Pocket_Cluster_Number();
extern void Write_Spherical_Probes_in_UCSF_DOCK_format();
extern void Label_Probe_Atoms_By_CHAR3DMAP();

/** FLOAT2DMAP : simple 2D matrix for general purpose **/
struct FLOAT2DMAP{
 int   N;     /* Number of rows and columns        */
 float **m;  /* 2D valiable: M[1..N][1..N]   */
};

