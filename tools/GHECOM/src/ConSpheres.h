/*
 <ConSpheres.h>
*/

/*** FUNCTIONS (GLOBAL) ***/
extern void Make_Sphere_Centers_from_PCaxis();
extern void Make_One_Sphere_Center_on_Gravity_Center();
extern void Setup_Radius_of_Spheres_Tangent_to_ATOMs();
extern void Remove_ATOM_less_than_Rmin_from_ATOMs();
extern void Setup_PCvar_and_PCaxis();
extern void Cal_Mean_and_CovM_from_ATOMs();
extern void write_ATOMs_as_Spheres_in_VRML();
extern void assign_radius_by_tFactor();
extern void assign_radius_by_Occup();
extern void split_string_to_float_array();
extern void Exclude_Spheres_Escape_to_Outside();
