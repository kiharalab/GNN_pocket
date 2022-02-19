/*

 <Grid3D.h>

*/

struct CHAR3DMAP{
  int  N[3];            /* Number of grids */
  float OrigPos[3];     /* Position of Origin */
  float grid_width;     /* Grid width */
  unsigned char ***map; /* map[N[0]][N[1]][N[2]] (malloc later) */
  int  minN[3],maxN[3]; /* min and max range for the non zelo value region  */
  char  malloced;       /* 1:malloced 0:not-malloced */
  int   Nforeground;    /* number of foreground pixels */
  int   rank;            /* 1,2,3,... for a linked list for CLUSTER */
  struct CHAR3DMAP *prev, *next;  /* for a linked list for CLUSTER */
  int   Noffset[3];     /* offset grid num[3] from the original map */
  float value;          /* for sorting of CLUSTERs */
};


struct INT3DMAP{
 int N[3];
 int ***map;
 char  malloced;       /* 1:malloced 0:not-malloced */
};


struct VOXEL{
 int   N[3];           /* Number of XYZ column 0:X, 1:Y, 2:Z */
 float grid_width;     /* Grid Length  */
 float OrigPos[3];     /* Real XYZ coordinates on the grid (0,0,0) */
 float min,max,ave,sd; /* 'min'imum, 'max'imum, 'ave'rage, 'sd':standard deviation of voxel values */
 char mal;             /* 1:already malloced 0: not malloced */
 float ***dat;         /* [0..X-1][0..Y-1][0..Z-1] */
 unsigned char MODE;
/*
     MODE    0 8-bit signed integer (range -128 to 127)
             1 16-bit signed integer
             2 32-bit signed real
             3 transform : complex 16-bit integers
             4 transform : complex 32-bit reals
             6 16-bit unsigned integer       2
*/

};



extern int  Malloc_CHAR3DMAP();
extern int  Malloc_INT3DMAP();
extern int  Malloc_VOXEL();
extern void Initialize_CHAR3DMAP();
extern void Initialize_INT3DMAP();
extern void Initialize_VOXEL();
extern void Free_CHAR3DMAP();
extern void Free_INT3DMAP();
extern void Free_VOXEL();
extern void Decide_Num_Grid_3D();
extern int  Set_Bit_3DMAP_Inside_Atoms();
extern void Set_Bit_By_Another_Map_in_CHAR3DMAP();
extern void Output_CHAR3DMAP_PDB_format();
extern void Output_CHAR3DMAP_with_tar_bitmask_PDB_format();
extern void Output_CHAR3DMAP_PDB_format_Min_Max();
extern void Output_CHAR3DMAP_PDB_format_with_MULSC_info();
extern void Output_CHAR3DMAP_PDB_format_Grouped_By_Char_Value();
extern void Output_Two_CHAR3DMAPs_in_PDB_format();
extern void MinMaxXYZ_Of_Atoms();
extern void Output_Two_CHAR3DMAPs_in_PDB_format();
extern void Count_CHAR3DMAP_Insides_Overlap();
extern int  Assign_Map_Value_to_tFactor_of_Atoms();
extern void Read_CHAR3DMAP_PDB_format();
extern float Corr_Coeff_of_Two_CHAR3DMAPs();
extern void Output_Histogram_of_Two_CHAR3DMAPs();
extern void Make_Difference_CHAR3DMAP();
extern double Point_Correlation_Coefficient();
extern void Make_Difference_CHAR3DMAP_MultiScale();
extern void Resize_CHAR3DMAP_by_MinMax();
extern void Copy_CHAR3DMAP();
extern void Copy_CHAR3DMAP_with_Nshift();
extern void Compare_Two_CHAR3DMAP();
extern void Compare_Two_Bit_in_CHAR3DMAP();
extern void assign_CHAR3DMAP_from_VOXEL();
extern void Foreground_MinMax_CHAR3DMAPL();
extern void resize_CHAR3DMAP_from_MinMax_and_MarginWidth();
extern void copy_CHAR3DMAP_Same_GRID_and_Diff_OrigPos();

