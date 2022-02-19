 /*
 
  <PdbIO.h>
 
 */
extern void Read_PDB_File();
extern void Make_Residue_List();
extern void Find_Filename_Core();
extern void Find_Filename_Core_with_Directory();
extern void Get_Part_Of_Line();
extern char *Get_Date_String();
extern char *Get_Date_String_PDB();
extern int  Number_Of_ATOMs();
extern int  Number_Of_RESIDUEs();
extern int  Number_Of_CHAINs();
extern int  Number_Of_MODELs();
extern int  Free_ATOMs();
extern void Set_Region_Of_Atoms();
extern void Set_Constant_tFactor();
extern void Renumber_Atom_Number();
extern void Renumber_Residue_Number();
extern void Write_PDB_File();
extern void Write_PDB_File_Residue();
extern void Add_Atom_To_Residue_Head();
extern void Add_Atom_To_Residue_Tail();
extern void Split_to_Words();
extern void Remove_Space();
extern void Add_String_to_STRING_NODE();
extern void Show_STRING_NODE();
extern void Write_One_Point_In_PDB_Format();
extern int  Find_XYZ_By_Pattern();
extern int  Find_XYZ_By_Pattern_In_PDB_File();
extern void Make_Copied_Atoms();
extern double Set_Molecule_Gcenter_to_Origin();
extern void Rotate_Molecule();
extern void InvRotate_Molecule();
extern void Translate_Molecule();
extern void Delete_Doubling_altLoc_Atoms();
extern void Change_N_CA_C_O_Hetatom_to_Atom();
extern void Delete_HETATMs_or_ATOMs();
