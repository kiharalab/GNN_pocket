/*

 <C3mapMergeSort.c>

=============================================================
Copyright 2018  Takeshi Kawabata.  All rights reserved.
This software is released under
the GNU Lesser General Public License (LGPL), see LICENSE.txt.
=============================================================

 
 functions for MergeSort for linked linear list. 

 coded by Takeshi Kawabata.

 Reference :
 Yoshiyuki Kondo. "Algorithms and Data Structures for C programmers".
  Soft Bank Publishing (1998). written in Japanese.  

 <Standard Definition of NODE> 

 struct NODE{
  double value;
  int    rank;
  struct NODE *next;
  struct NODE *prev;
 };


 <Standard Usage> 
 
 struct CHAR3DMAP HEAD,*bn;
 
 bn = Merge_Sort_List_CHAR3DMAP(Head.next,'-'); 
 Head.next = bn;
   or 
 Merge_Sort_Double_Linked_List_CHAR3DMAP(head,type)


<option>
 type : if 'D', 'Decreasing':large to small, otherwise small to large  

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globalvar.h" 
#include "Grid3D.h" 

/*** FUNCTIONS (GLOBAL) ***/
struct CHAR3DMAP *Merge_Sort_List_CHAR3DMAP();
void Merge_Sort_Double_Linked_List_CHAR3DMAP();

/*** FUNCTIONS (LOCAL) ***/
static void Set_Prev_Pointer();
static struct CHAR3DMAP *Merge_List();
static void Show_NODE();


void Show_NODE(head)
 struct CHAR3DMAP *head;
{ 
 struct CHAR3DMAP *bn;      

 /* bn = head->next;   */
 bn = head;
 while (bn != NULL){ 
   printf("%f\n",bn->value); 
   bn = bn->next;
 }

} /* end of Show_NODE() */



void Merge_Sort_Double_Linked_List_CHAR3DMAP(head,type)
 struct CHAR3DMAP *head;
 char type; /* if 'D'ecreasing: large to small, otherwise 'I'ncreasing:small to large */ 
{
 struct CHAR3DMAP *bn;
 bn = Merge_Sort_List_CHAR3DMAP(head->next,type); 
 head->next = bn;
 Set_Prev_Pointer(head);
 head->rank = 0;
 bn = head;
 while (bn->next != NULL){
   bn = bn->next;
   bn->rank = bn->prev->rank + 1;
 }

} /* end of Merge_Sort_Double_Linked_List_CHAR3DMAP() */



struct CHAR3DMAP *Merge_List(a,b,type)
 struct CHAR3DMAP *a,*b;
 char type; /* if 'D'ecreasing: large to small, otherwise 'I'ncreasing:small to large */ 
{
 struct CHAR3DMAP head, *p;

 /* Initialize */ 
 p = &head;

 /* Repeat until lists 'a' and 'b' become empty */

 /* The smaller node is attached to the list 'p' */
 if (type != 'D') 
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






struct CHAR3DMAP *Merge_Sort_List_CHAR3DMAP(x,type)
 struct CHAR3DMAP *x;
 char type;
{
 struct CHAR3DMAP *a, *b, *p;


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

 return(Merge_List(Merge_Sort_List_CHAR3DMAP(x,type),Merge_Sort_List_CHAR3DMAP(p,type),type));


} /* end of *Merge_Sort_List_CHAR3DMAP() */




void Set_Prev_Pointer(head)
 struct CHAR3DMAP *head;
{
 struct CHAR3DMAP *bn,*pn;

 bn = head;
 while (bn->next !=NULL){
   pn = bn;
   bn = bn->next;
   bn->prev = pn; 
  }

} /* end of Set_Prev_Pointer() */

