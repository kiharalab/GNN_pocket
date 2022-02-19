/*

 <GrPtMergeSort.c>
 
 functions for MergeSort for linked linear list. 
 coded by Takeshi Kawabata.

=========== "ghecom" program ===========

Copyright 2018  Takeshi Kawabata.  All rights reserved.

This software is released under
the GNU Lesser General Public License (LGPL) version 3, see LICENSE.txt.

=========== Installation and Usage ===========

See html file under the "doc/" directory.

=========== Citing "ghecom" program ===========

Please cite:

1) Kawabata T. (2010) Detection of multi-scale pockets on protein surfaces using mathematical morphology. Proteins,78, 1195-1121.




 >> Reference for the sorting algorithm <<:
 Yoshiyuki Kondo. "Algorithms and Data Structures for C programmers".
  Soft Bank Publishing (1998). written in Japanese.  

 <Standard Definition of NODE> 

 struct NODE{
  float value;
  struct NODE *next;
  struct NODE *prev;
 };


 <Standard Usage> 
 
 struct GRIDPNT3D HEAD,*bn;
 
 bn = Merge_Sort_List_GRIDPNT3D(Head.next,'-'); 
 Head.next = bn;
   or 
 Merge_Sort_Double_Linked_List_GRIDPNT3D(head,type)

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globalvar.h" 
#include "pdbstruct.h" 
#include "OpnRotPrb.h" 

/*** FUNCTIONS (GLOBAL) ***/
struct GRIDPNT3D *Merge_Sort_List_GRIDPNT3D();
void Merge_Sort_Double_Linked_List_GRIDPNT3D();

/*** FUNCTIONS (LOCAL) ***/
static void Set_Prev_Pointer();
static struct GRIDPNT3D *Merge_List();


void Merge_Sort_Double_Linked_List_GRIDPNT3D(head,type)
 struct GRIDPNT3D *head;
 char type; /* if 'R', large to small, otherwise small to large */ 
{
 struct GRIDPNT3D *bn;
 bn = Merge_Sort_List_GRIDPNT3D(head->next,type); 
 head->next = bn;
 Set_Prev_Pointer(head);

} /* end of Merge_Sort_Double_Linked_List_GRIDPNT3D() */



struct GRIDPNT3D *Merge_List(a,b,type)
 struct GRIDPNT3D *a,*b;
 char type; /* if 'R', large to small, otherwise small to large */ 
{
 struct GRIDPNT3D head, *p;

 /* Initialize */ 
 p = &head;

 /* Repeat until lists 'a' and 'b' become empty */

 /* The smaller node is attached to the list 'p' */
 if (type != 'R') 
 while ((a!=NULL)&&(b!=NULL)){
  if (a->value <= b->value) { p->next = a; p = a; a = a->next; }
                       else { p->next = b; p = b; b = b->next; }
 }  

 else
 
 /* The larger node is attached to the list 'p' */
 while ((a!=NULL)&&(b!=NULL)){
  if (a->value >= b->value) { p->next = a; p = a; a = a->next; }
                       else { p->next = b; p = b; b = b->next; }
 }  


 /* The rest nodes are attached to the list 'p'. */
       if ((a == NULL) && (b !=NULL)) p->next = b; 
 else if  ((b == NULL) && (a !=NULL)) p->next = a; 
   
 return(head.next);

} /* end of  Merge_List() */






struct GRIDPNT3D *Merge_Sort_List_GRIDPNT3D(x,type)
 struct GRIDPNT3D *x;
 char type;
{
 struct GRIDPNT3D *a, *b, *p;


 /* if x has no nodes, or one nodes */
 if ((x==NULL) || (x->next ==NULL)) return(x); 

 /* a is the 1st node, and 'b' is the 2nd or 3rd node */
 a = x;
 b = x->next;
 if (b->next !=NULL) b = x->next;

 /* 
   Progress a by 1 node, and 'b' by two nodes. 
   if 'b' reaches to the end, 'a' reaches to the midpoint.
 */

 while (b->next !=NULL){
  a = a->next;
  b = b->next;
  if (b->next !=NULL) b = b->next;
  }

 /* Truncate the list on the midpoint */
 p = a->next; a->next = NULL;

 return(Merge_List(Merge_Sort_List_GRIDPNT3D(x,type),Merge_Sort_List_GRIDPNT3D(p,type),type));


} /* end of *Merge_Sort_List_GRIDPNT3D() */




void Set_Prev_Pointer(head)
 struct GRIDPNT3D *head;
{
 struct GRIDPNT3D *bn,*pn;

 bn = head;
 while (bn->next !=NULL){ 
   pn = bn;
   bn = bn->next;
   bn->prev = pn; 
  }

} /* end of Set_Prev_Pointer() */

