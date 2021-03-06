/*

 <ATOMs_ghecom_from_PDBmmCIF.c>

=============================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.
This software is released under
the GNU Lesser General Public License (LGPL), see LICENSE.txt.
=============================================================


 This file is for getting a linked list of "struct ATOM" 
 from the "struct PDB_DATA" (generated from "struct BLOCK").

 This file should be modified considering the definition 
 of "struct ATOM".

 This file is modified for the program "ghecom".

 LastModDate: 2019/06/28.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "globalvar.h" 
#include "io_mmCIF.h" 
#include "PDB_from_mmCIF.h" 
#include "ATOMs_ghecom_from_PDBmmCIF.h" 
#include "pdbstruct.h"
#include "PdbIO.h"

/** FUNCTIONS (GLOBAL) **/
void Read_mmCIF_File();
void make_ATOM_list_from_ASSEMBLY();
void make_ATOM_list_from_UNITMOLs();

/** FUNCTIONS (LOCAL) **/
static void insert_ATOM_data();


void Read_mmCIF_File(iciffile,HeadAtom,AtomSelect,AtomHetSelect,ChainIDStr,asym_id_Str,assembly_id,MaxAtomForAll)
  char *iciffile;
  struct ATOM *HeadAtom;
  char   AtomSelect; /* Atom selection. 'A'll atom except hydrogen, 'R'esidue-based (only ' CA ' and ' P ')  */
  char   AtomHetSelect; /* 'A' : Only read ATOM, 'H':Only read HETATM  'B'oth */
  char   *ChainIDStr;  /* auth_asym_id */
  char   *asym_id_Str;
  char   *assembly_id;
  int    MaxAtomForAll;
{
  struct BLOCK BLK;
  struct PDB_DATA PDB;
  int    Natom;

  printf("#Read_mmCIF_File('%s' AtomSelect:%c AtomHetSelect %c ChainID:'%s' assembly_id:'%s' MaxTomForAll %d)\n",
     iciffile,AtomSelect,AtomHetSelect,ChainIDStr,assembly_id,MaxAtomForAll); fflush(stdout);
 
  HeadAtom->next = NULL;

  read_mmCIF_file_into_BLOCK(iciffile,&BLK);

  printf("#Ntable %d\n",BLK.Ntable);

  /*
  tab = &(BLK.head_table);
  while (tab->next != NULL){
    tab = tab->next;
    printf("#TABLE '%s' Ncolumn %d Nrow %d\n",tab->table_name, tab->Ncolumn, tab->Nrow); 
    word = &(tab->head_column_name); 
    while (word->next != NULL){
      word = word->next;
      printf("'%s' c %d\n",word->word,Integer_from_key(tab->colnum_hash_tab,word->word));
    }
    rl = &(tab->head_rowline);
    while (rl->next != NULL){
      rl = rl->next;
      printf("#oligomeric_details '%s'\n",data_from_rowline(tab,rl,"oligomeric_details"));
      for (c=0;c<tab->Ncolumn;++c){
        printf(" c %d '%s'\n",c,rl->data[c]);
      }
    }
  }
  */

  PDB.head_assembly.next  = NULL;
  PDB.head_operation.next = NULL;
  PDB.head_unitmol.next   = NULL;
  PDB.head_asmblmol.next  = NULL;
 

  make_ASSEMBLY_from_BLOCK(&PDB,&BLK);
  make_OPERATION_from_BLOCK(&PDB,&BLK);
  make_UNITMOL_from_BLOCK(&PDB,&BLK);

  if (assembly_id[0] != '\0'){
    make_ASMBLMOL_from_BLOCK(&PDB,&BLK,assembly_id);
    make_ATOM_list_from_ASSEMBLY(HeadAtom,&PDB,&BLK,assembly_id,AtomSelect,ChainIDStr,asym_id_Str);
  }
  else{
    make_ATOM_list_from_UNITMOLs(HeadAtom,&PDB,&BLK,AtomSelect,ChainIDStr,asym_id_Str);
  }

  Natom  = Number_Of_ATOMs(HeadAtom);
  printf("#Number_Of_Atom:%d\n",Natom);

  /** If Natom>MaxAtomForAll -> change AtomSelect='R' */ 
  if ((AtomSelect=='A') && (MaxAtomForAll > 0) && (Natom > MaxAtomForAll)){
    Free_ATOMs(HeadAtom);
    if (assembly_id[0] != '\0'){
      make_ATOM_list_from_ASSEMBLY(HeadAtom,&PDB,&BLK,assembly_id,'R',ChainIDStr, asym_id_Str);
    }
    else{
      make_ATOM_list_from_UNITMOLs(HeadAtom,&PDB,&BLK,'R',ChainIDStr,asym_id_Str);
    }
    Natom  = Number_Of_ATOMs(HeadAtom);
    printf("#Number_Of_Atom for AtomSelect 'R' :%d\n",Natom);
  } 

  free_BLOCK(&BLK);
  free_PDB_DATA(&PDB);

 /*
 Delete_Doubling_altLoc_Atoms(HeadAtom);
 */
 if (PAR.Change_N_CA_C_O_Hetatom_to_Atom=='T'){ Change_N_CA_C_O_Hetatom_to_Atom(HeadAtom);}

  if (AtomHetSelect == 'A'){
   Delete_HETATMs_or_ATOMs(HeadAtom,'H');
  }
  if (AtomHetSelect == 'H'){
   Delete_HETATMs_or_ATOMs(HeadAtom,'A');
  }

  Natom  = Number_Of_ATOMs(HeadAtom);

  printf("#Natom %d\n",Natom);
  printf("#Number_Of_Atom:%d\n",Number_Of_ATOMs(HeadAtom));

} /* end of Read_mmCIF_File() */





void make_ATOM_list_from_ASSEMBLY(head_atom,P,blk,assembly_id,AtomSelect,ChainIDStr,asym_id_Str)
  struct ATOM *head_atom;
  struct PDB_DATA *P;
  struct BLOCK *blk;
  char *assembly_id;
  char   AtomSelect; /* Atom selection. 'A'll atom except hydrogen, 'R'esidue-based (only ' CA ' and ' P ')  */
  char   *ChainIDStr;      /* (auth_asym_id) */
  char   *asym_id_Str;  
{
  struct TABLE *tb;
  struct ROWLINE *rn;
  int group_PDB,id,auth_atom_id,label_alt_id,auth_comp_id,auth_asym_id,auth_seq_id,pdbx_PDB_ins_code;
  int Cartn_x,Cartn_y,Cartn_z,occupancy,B_iso_or_equiv,type_symbol,label_asym_id, pdbx_PDB_model_num;
  struct ASSEMBLY *A;
  struct UNITMOL  *um;
  struct ASMBLMOL *am;
  struct ATOM *atom;
  int i,k,Natom,model_num;
  char accept,representative,pre_asym_id[16];
  int curr_pdbx_PDB_model_num, prev_pdbx_PDB_model_num, Npdbx_PDB_model_num;

  printf("#make_ATOM_list_from_ASSEMBLY()\n");
  tb = TABLE_from_key(blk->table_hash_tab,"atom_site");
  if (tb==NULL){
     printf("#ERROR:Can't find table '%s'.\n","atom_site");
     exit(1);
  }

   setup_column_numbers_for_ATOM_HETATM_line(
    tb, &group_PDB,&id, &auth_atom_id,&label_alt_id, &auth_comp_id,&auth_asym_id,&auth_seq_id,&pdbx_PDB_ins_code,
    &Cartn_x, &Cartn_y, &Cartn_z,&occupancy, &B_iso_or_equiv,&type_symbol,&label_asym_id,&pdbx_PDB_model_num);

  head_atom->next = NULL;
  head_atom->num  = -1;
  atom = head_atom;
  A = &(P->head_assembly);
  Natom = 0; 
  model_num = 0;
  pre_asym_id[0] = '\0';
  Npdbx_PDB_model_num = 0;
  curr_pdbx_PDB_model_num = 0;
  prev_pdbx_PDB_model_num    = -1;

  while (A->next != NULL){
    A = A->next;
    if ((assembly_id[0]=='\0')||(strcmp(A->id,assembly_id)==0)){
      printf("#**assembly_id '%s' Nasmblmol %d\n",A->id,A->Nasmblmol); fflush(stdout);
      for (i=0;i<A->Nasmblmol;++i){
        um = find_unitmol(P,A->ASYM_ID_LIST[i]);
        am = find_asmblmol(P,A->ASYM_ID_LIST[i],A->OPER_EXPRESSION_LIST[i]);
        if ((um != NULL)&&(am != NULL)){ 
          rn = um->head_rowline_pointer;
          k = 0;
          /* while ((rn->next != NULL)&&(rn->next != um->tail_rowline_pointer)){ */
           while ((rn != NULL)&&(rn->next != NULL)&&(rn->next != um->tail_rowline_pointer)){ 
            rn = rn->next;

            if (strcmp(rn->data[label_asym_id],pre_asym_id)!=0){model_num += 1;} /* Different "asym_id" is represented by different "model_num".*/

            accept = 1;
            if (strncmp(rn->data[auth_comp_id],"HOH",3)==0){accept = 0;}
            if (strncmp(rn->data[auth_comp_id],"DDD",3)==0){accept = 0;}

            if ( (ChainIDStr[0]!='\0')&&(ChainIDStr[0] != '-') && (strcmp(ChainIDStr,rn->data[auth_asym_id])!=0)){ 
              accept = 0;
            }
            if ( (asym_id_Str[0]!='\0')&&(asym_id_Str[0] != '-') && (strcmp(asym_id_Str,rn->data[label_asym_id])!=0)){ 
              accept = 0;
            }

            /** only take the first model  **/ 
            curr_pdbx_PDB_model_num = atoi(rn->data[pdbx_PDB_model_num]);
            if (curr_pdbx_PDB_model_num != prev_pdbx_PDB_model_num){
             Npdbx_PDB_model_num += 1;
            }
            prev_pdbx_PDB_model_num = curr_pdbx_PDB_model_num;
            if (Npdbx_PDB_model_num > 1) { accept = 0;}

            if (AtomSelect == 'R'){
              representative = 0;
              if ((strcmp(rn->data[type_symbol],"C")==0)&& (strcmp(rn->data[auth_atom_id],"CA")==0)){representative = 1;}
              if ((strcmp(rn->data[type_symbol],"P")==0)&& (strcmp(rn->data[auth_atom_id],"P")==0)) {representative = 1;}
              if (representative==0) { accept = 0;}
            }

            if (accept==1){ 
              atom->next = (struct ATOM*)malloc(sizeof(struct ATOM)); 
              Natom += 1;
              atom->next->prev = atom;
              atom = atom->next;
              atom->next = NULL;
              atom->rnext = NULL;
              atom->rprev = NULL;
              atom->num = atom->prev->num + 1;

              insert_ATOM_data(
                atom,
                rn->data[group_PDB], rn->data[id], rn->data[auth_atom_id], rn->data[label_alt_id], rn->data[auth_comp_id], 
                rn->data[auth_asym_id], rn->data[auth_seq_id], rn->data[pdbx_PDB_ins_code],
                am->Cartn_x[k], am->Cartn_y[k], am->Cartn_z[k], 
                rn->data[occupancy], rn->data[B_iso_or_equiv],rn->data[type_symbol]
              );
              atom->model_num = model_num;
              atom->model_num = i + 1; 
            }
            k += 1;
            sprintf(pre_asym_id,"%s",rn->data[label_asym_id]);
           }
         }
       }
     }
   }

  printf("#Natom %d\n",Natom);

} /* end of make_ATOM_list_from_ASSEMBLY() */




void make_ATOM_list_from_UNITMOLs(head_atom,P,blk,AtomSelect,ChainIDStr, asym_id_Str)
  struct ATOM *head_atom;
  struct PDB_DATA *P;
  struct BLOCK *blk;
  char   AtomSelect; /* Atom selection. 'A'll atom except hydrogen, 'R'esidue-based (only ' CA ' and ' P ')  */
  char   *ChainIDStr;      /* (auth_asym_id) */
  char   *asym_id_Str;  
{
  struct TABLE *tb;
  struct ROWLINE *rn;
  int group_PDB,id,auth_atom_id,label_alt_id,auth_comp_id,auth_asym_id,auth_seq_id,pdbx_PDB_ins_code;
  int Cartn_x,Cartn_y,Cartn_z,occupancy,B_iso_or_equiv,type_symbol,label_asym_id, pdbx_PDB_model_num;
  struct UNITMOL  *um;
  struct ATOM *atom;
  int Natom,model_num;
  char accept,representative,pre_asym_id[16];
  int curr_pdbx_PDB_model_num, prev_pdbx_PDB_model_num, Npdbx_PDB_model_num;

  printf("#make_ATOM_list_from_UNITMOLs()\n");
  tb = TABLE_from_key(blk->table_hash_tab,"atom_site");
  if (tb==NULL){
     printf("#ERROR:Can't find table '%s'.\n","atom_site");
     exit(1);
  }

   setup_column_numbers_for_ATOM_HETATM_line(
    tb, &group_PDB,&id, &auth_atom_id,&label_alt_id, &auth_comp_id,&auth_asym_id,&auth_seq_id,&pdbx_PDB_ins_code,
    &Cartn_x, &Cartn_y, &Cartn_z,&occupancy, &B_iso_or_equiv,&type_symbol,&label_asym_id,&pdbx_PDB_model_num);


  head_atom->next = NULL;
  atom = head_atom;
  um = &(P->head_unitmol);
  Natom = 0;
  pre_asym_id[0] = '\0';
  model_num = 0;
  
  Npdbx_PDB_model_num = 0;
  curr_pdbx_PDB_model_num = 0;
  prev_pdbx_PDB_model_num = -1;

  while (um->next != NULL){
    um = um->next;
    rn = um->head_rowline_pointer;

    while ((rn->next != NULL)&&(rn->next != um->tail_rowline_pointer)){
      rn = rn->next;

      if (strcmp(rn->data[label_asym_id],pre_asym_id)!=0){model_num += 1;} /* Different "asym_id" is represented by different "model_num".*/

      accept = 1;
      if (strncmp(rn->data[auth_comp_id],"HOH",3)==0){accept = 0;}
      if (strncmp(rn->data[auth_comp_id],"DDD",3)==0){accept = 0;}

      /* if ((ChainID != '-') && (ChainID != rn->data[auth_asym_id][0])) { accept = 0;} */
      if ( (ChainIDStr[0]!='\0')&&(ChainIDStr[0] != '-') && (strcmp(ChainIDStr,rn->data[auth_asym_id])!=0)){ 
        accept = 0;
      }
      if ( (asym_id_Str[0]!='\0')&&(asym_id_Str[0] != '-') && (strcmp(asym_id_Str,rn->data[label_asym_id])!=0)){ 
        accept = 0;
      }
  
      /** only take the first model  **/ 
      curr_pdbx_PDB_model_num = atoi(rn->data[pdbx_PDB_model_num]);
      if (curr_pdbx_PDB_model_num != prev_pdbx_PDB_model_num){
        Npdbx_PDB_model_num += 1;
      }
      prev_pdbx_PDB_model_num = curr_pdbx_PDB_model_num;
      if (Npdbx_PDB_model_num > 1) { accept = 0;}


      if (AtomSelect == 'R'){
        representative = 0;
        if ((strcmp(rn->data[type_symbol],"C")==0)&& (strcmp(rn->data[auth_atom_id],"CA")==0)){representative = 1;}
        if ((strcmp(rn->data[type_symbol],"P")==0)&& (strcmp(rn->data[auth_atom_id],"P")==0)) {representative = 1;}
        if (representative==0) { accept = 0;}
      }

      


      if (accept==1){
        atom->next = (struct ATOM*)malloc(sizeof(struct ATOM)); 
        Natom += 1;
        atom->next->prev = atom;
        atom = atom->next;
        atom->next = NULL;
        atom->rnext = NULL;
        atom->rprev = NULL;
        atom->num = atom->prev->num + 1;

        insert_ATOM_data(
           atom,
           rn->data[group_PDB], rn->data[id], rn->data[auth_atom_id], rn->data[label_alt_id], rn->data[auth_comp_id], 
           rn->data[auth_asym_id], rn->data[auth_seq_id], rn->data[pdbx_PDB_ins_code],
           atof(rn->data[Cartn_x]), atof(rn->data[Cartn_y]), atof(rn->data[Cartn_z]), 
           rn->data[occupancy], rn->data[B_iso_or_equiv],rn->data[type_symbol]
        );
        atom->model_num = model_num;
 
      /* printf("#%s %s %s\n",atom->Atom,atom->Resi,atom->Rnum); */
      }

      sprintf(pre_asym_id,"%s",rn->data[label_asym_id]);
    }
  }

  printf("#Natom(make_ATOM_list_from_UNITMOLs()) %d\n",Natom);

} /* end of make_ATOM_list_from_UNITMOLs() */



void insert_ATOM_data(
  atom,
  group_PDB, id, auth_atom_id, label_alt_id, auth_comp_id, auth_asym_id, auth_seq_id, pdbx_PDB_ins_code,
  Cartn_x, Cartn_y, Cartn_z, occupancy, B_iso_or_equiv,type_symbol
)
  struct ATOM *atom; 
  char *group_PDB, *id, *auth_atom_id, *label_alt_id, *auth_comp_id, *auth_asym_id,*auth_seq_id, *pdbx_PDB_ins_code;
  float Cartn_x, Cartn_y, Cartn_z;
  char *occupancy, *B_iso_or_equiv, *type_symbol;
{
  int L,Ltype_symbol;

  if (strcmp(group_PDB,"ATOM")==0){
    atom->AHtype = 'A';
  }
  else{
    atom->AHtype = 'H';
  }

  atom->Pos[0] = Cartn_x;
  atom->Pos[1] = Cartn_y;
  atom->Pos[2] = Cartn_z;

/*  [auth_atom_id]  
  
  >> RULE << 
   if (length(auth_atom_id)<4){
     if (length(type_symbol)==2) auth_atom_id starts at 1st position, otherwise starts at 2nd position.  
   }

>> 3icb.pdb <<
          1         2         3         5         5         6         7
01234567890123456789012345678901234567890123456789012345678901234567890123456789
ATOM    586  CA  SER A  74       4.597   6.256   1.954  1.00 37.65           C
ATOM    592  CA  GLN A  75       2.756   4.149   4.526  1.00 45.45           C
HETATM  602 CA    CA A  76      22.527  19.883   5.062  1.00  3.59          CA
HETATM  603 CA    CA A  77      25.690   9.260   6.305  1.00  4.99          CA
  */
  L = strlen(auth_atom_id);
  sprintf(atom->Atom,"    ");
  Ltype_symbol = strlen(type_symbol);
  if (L==4){
    sprintf(atom->Atom,"%s",auth_atom_id);
  }
  else if (L==3){
   if (Ltype_symbol==2){sprintf(atom->Atom,"%s ",auth_atom_id); }
                  else {sprintf(atom->Atom," %s",auth_atom_id); }
  }
  else if (L==2){
   if (Ltype_symbol==2){sprintf(atom->Atom,"%s  ",auth_atom_id); }
                  else {sprintf(atom->Atom," %s ",auth_atom_id); }
  }
  else if (L==1){
   if (Ltype_symbol==2){sprintf(atom->Atom,"%s   ",auth_atom_id); }
                  else {sprintf(atom->Atom," %s  ",auth_atom_id); }
  }

  sprintf(atom->Anum,"%d",atoi(id)%100000);

  L = strlen(auth_comp_id);
       if (L==3){ sprintf(atom->Resi,"%s",auth_comp_id); }
  else if (L==2){ sprintf(atom->Resi," %s",auth_comp_id); }
  else if (L==1){ sprintf(atom->Resi,"  %s",auth_comp_id); }

  sprintf(atom->Rnum,"%4s",auth_seq_id);

/* For example ribosome "4v5b" has auth_asym_id =
   ["A0","A1","A2","A3","A5","AA","AB","AC",..."AZ","BA","BB",...]  

   The last char of auth_asym_id is good for atom->Chain.

 */

  /*
  atom->Chain = auth_asym_id[0];
  */
  L = strlen(auth_asym_id);
  atom->Chain = auth_asym_id[L-1];

  if ((pdbx_PDB_ins_code[0] != '.')&&(pdbx_PDB_ins_code[0] != '?')){
    atom->altLoc = pdbx_PDB_ins_code[0];
  }
  atom->Occup   = atof(occupancy);
  atom->tFactor = atof(B_iso_or_equiv);

} /* end of insert_ATOM_data() */



