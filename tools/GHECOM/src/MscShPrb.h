/*

 <MscShPrb.h>
 
 for dealing multi scaling "shell" sphere probes  

*/


struct POS_INT{
  int    iPos[3];
  int    num;
  struct POS_INT *next,*prev;
};


#define MAX_RADIUS_TYPE 256

struct MULSC_SHELL_PRB{
  int    Nradius;                         /* Number of radius */
  float  grid_width;
  int    Ngrid[3];                        /* Number of Grids */
  int    Ngrid_half[3];                   /* Ngrid = 2 * Ngrid_halft + 1 */
  float  Rarray[MAX_RADIUS_TYPE];         /* Rarray[1..Nradius] : Array for radius */ 
  struct POS_INT  HeadPosInt[MAX_RADIUS_TYPE];
                      /* HeadPosInt[0..Nradius-1] : Head of linked list of POS_INTs */
  int    Npoint[MAX_RADIUS_TYPE];          /* Npoint[1..Nradius] : Number of points for the shell  */ 
};


/** FUNCTIONS(global) **/
extern void Set_Rarray_By_Min_Max_Bin();
extern void Set_Rarray_By_Nradius_Min_Max();
extern void Make_MultiScale_Sphere_Shell_Probe_Map();
extern void Write_PDB_File_MultiSc_Shell_Spheres();
extern void Dilation_MultiScaleShell();
extern void Erosion_MultiScaleShell();
extern void Opening_MultiScaleShell();
extern int  Label_Target_Shape_By_Multiscale_Closing();
extern int  Label_Target_Shape_By_Multiscale_Pocket();
