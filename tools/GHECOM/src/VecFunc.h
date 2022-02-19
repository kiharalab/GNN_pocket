/*
 <VecFunc.h>
*/

                                                                                                                   
/** ROT_MATRIX : for Array of Rotational Matrix **/
struct ROT_MATRIX{
 float R[3][3];
};


/*** FUNCTIONS (GLOBAL) */
extern void Cross_Vec();
extern float Dot_Prod();
extern void Sub_Vec();
extern void Equal_Vec();
extern void Min_Vec();
extern void Max_Vec();
extern float Normalize();
extern float Vec_Length();
extern float Distance();
extern void Print_Vec();
extern void Make_Random_Rotation_Matrix();
extern void Set_Unitary_Matrix();
extern void Make_RotMatrix_By_Quaternion();
