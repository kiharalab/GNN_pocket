/*

 <GclsMergeSort.c>
 
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



 Reference for the sorting algorithm:
 Yoshiyuki Kondo. "Algorithms and Data Structures for C programmers".
  Soft Bank Publishing (1998). written in Japanese.  

 <Standard Definition of NODE> 

 struct NODE{
  float value;
  struct NODE *next;
  struct NODE *prev;
 };


 <Standard Usage> 
 
 struct GRID_CLUSTER HEAD,*bn;
 
 bn = Merge_Sort_List_GRID_CLUSTER(Head.next,'-'); 
 Head.next = bn;
   or 
 Merge_Sort_Double_Linked_List_GRID_CLUSTER(head,type)

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globalvar.h" 
#include "pdbstruct.h" 
#include "MscPockClus.h" 

/*** FUNCTIONS (GLOBAL) ***/
struct GRID_CLUSTER *Merge_Sort_List_GRID_CLUSTER();
void Merge_Sort_Double_Linked_List_GRID_CLUSTER();

/*** FUNCTIONS (LOCAL) ***/
static void Set_Prev_Pointer();
static struct GRID_CLUSTER *Merge_List();




void Merge_Sort_Double_Linked_List_GRID_CLUSTER(head,type)
 struct GRID_CLUSTER *head;
 char type; /* if 'R', large to small, otherwise small to large */ 
{
 struct GRID_CLUSTER *bn;
 bn = Merge_Sort_List_GRID_CLUSTER(head->next,type); 
 head->next = bn;
 Set_Prev_Pointer(head);

} /* end of Merge_Sort_Double_Linked_List_GRID_CLUSTER() */



struct GRID_CLUSTER *Merge_List(a,b,type)
 struct GRID_CLUSTER *a,*b;
 char type; /* if 'R', large to small, otherwise small to large */ 
{
 struct GRID_CLUSTER head, *p;

 /* Initialize */ 
 p = &head;

 /* Repeat until lists 'a' and 'b' become empty */

 /* The smaller node is attached to the list 'p' */
 if (type != 'R') 
 while ((a!=NULL)&&(b!=NULL))
 {
  if (a->value <= b->value) { p->next = a; p = a; a = a->next; }
                       else { p->next = b; p = b; b = b->next; }
 }  

 else
 
 /* The larger node is attached to the list 'p' */
 while ((a!=NULL)&&(b!=NULL))
 {
  if (a->value >= b->value) { p->next = a; p = a; a = a->next; }
                       else { p->next = b; p = b; b = b->next; }
 }  


 /* The rest nodes are attached to the list 'p'. */
       if ((a == NULL) && (b !=NULL)) p->next = b; 
 else if  ((b == NULL) && (a !=NULL)) p->next = a; 
   
 return(head.next);

} /* end of  Merge_List() */






struct GRID_CLUSTER *Merge_Sort_List_GRID_CLUSTER(x,type)
 struct GRID_CLUSTER *x;
 char type;
{
 struct GRID_CLUSTER *a, *b, *p;


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

 while (b->next !=NULL)
 {
  a = a->next;
  b = b->next;
  if (b->next !=NULL) b = b->next;
  }

 /* Truncate the list on the midpoint */
 p = a->next; a->next = NULL;

 return(Merge_List(Merge_Sort_List_GRID_CLUSTER(x,type),Merge_Sort_List_GRID_CLUSTER(p,type),type));


} /* end of *Merge_Sort_List_GRID_CLUSTER() */




void Set_Prev_Pointer(head)
 struct GRID_CLUSTER *head;
{
 struct GRID_CLUSTER *bn,*pn;

 bn = head;
 while (bn->next !=NULL) 
 {
   pn = bn;
   bn = bn->next;
   bn->prev = pn; 
  }

} /* end of Set_Prev_Pointer() */

