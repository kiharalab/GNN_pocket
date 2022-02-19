/*

 <queue_test.c>

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

static char LastModDate[] = "2018/10/17";


static int *QUEUE_X, *QUEUE_Y, *QUEUE_Z;     /* [0..Q_MAX] (malloc later) */
static int Q_MAX, Q_HEAD, Q_TAIL;

static int PUT_QUEUE_XYZ(x,y,z)
  int x,y,z;
{
  if (Q_TAIL > Q_MAX){
    printf("#ERROR:PUT_QUEUE_XYZ:Q_TAIL %d is > Q_MAX %d .\n",Q_TAIL,Q_MAX);
  }
  QUEUE_X[Q_TAIL] = x;
  QUEUE_Y[Q_TAIL] = y;
  QUEUE_Z[Q_TAIL] = z;
  Q_TAIL += 1;
  if (Q_TAIL > Q_MAX){
     Q_TAIL = 0;
  }
  printf("#PUT(%d %d %d)\n",x,y,z);
  if (Q_TAIL == Q_HEAD){
    printf("#WARNING:QUEUE is FULL !!\n");
    return(0);
  }
  return(1);

} /* end of PUT_QUEUE_XYZ() */


static int GET_QUEUE_XYZ(x, y, z)
  int *x, *y, *z;
{
  if (Q_HEAD != Q_TAIL){
    *x = QUEUE_X[Q_HEAD];
    *y = QUEUE_Y[Q_HEAD];
    *z = QUEUE_Z[Q_HEAD];
    Q_HEAD = Q_HEAD + 1;
    printf("#GET(%d %d %d)\n",*x,*y,*z);
    if (Q_HEAD>Q_MAX){
      Q_HEAD = 0;
    }
    return(1);
  }
  else{
    *x =  *y = *z = 0;
    printf("#GET --failed--\n");
    return(0);
  }
} /* end of PUT_QUEUE_XYZ() */



static void INIT_QUEUE_XYZ()
{
  Q_HEAD = 0;
  Q_TAIL = 0;
} /* end of INIT_QUEUE_XYZ() */

static void QUEUE_XYZ_EMPTY()
{
  Q_HEAD = Q_TAIL = 0;
}  /* end of QUEUE_XYZ_EMPTY() */

static int NUM_QUEUE_XYZ()
{
  int N;

  N = Q_TAIL - Q_HEAD;
  if (N<0){
    N = -N;
  }
  return(N);
} /* end of NUM_QUEUE_XYZ() */


static void SHOW_QUEUE_XYZ()
{
  int i,end,n;

  printf("#SHOW_QUEUE_XYZ() Q_MAX %d Q_HEAD %d Q_TAIL %d  NUM %d \n",Q_MAX,Q_HEAD,Q_TAIL,NUM_QUEUE_XYZ());
  i = Q_HEAD;
  end = 0; 
  n = 0;
  while ((n<=Q_MAX) && (end==0)){
    if (i==Q_TAIL){
      end = 1;
    }
    else{ 
      printf("[%d]%d %d %d\n",i,QUEUE_X[i],QUEUE_Y[i],QUEUE_Z[i]);
      i = (i+1)%(Q_MAX+1); 
      n += 1;
    }
  }
  printf("\n");
} /* end of SHOW_QUEUE_XYZ() */


int main(argc, argv)
  int argc;
  char **argv;
{
  int i,x,y,z;

  Q_MAX = 5;
  printf("#Q_MAX %d\n",Q_MAX);
  QUEUE_X = (int *)malloc(sizeof(int)*(Q_MAX+1));
  QUEUE_Y = (int *)malloc(sizeof(int)*(Q_MAX+1));
  QUEUE_Z = (int *)malloc(sizeof(int)*(Q_MAX+1));

  for (i=0;i<Q_MAX;++i){
    QUEUE_X[i] =  -1;
    QUEUE_Y[i] =  -1;
    QUEUE_Z[i] =  -1;
  } 
  INIT_QUEUE_XYZ();
  GET_QUEUE_XYZ(&x,&y,&z); printf("#get %d %d %d\n",x,y,z);
  SHOW_QUEUE_XYZ();

  for (i=0;i<2*Q_MAX;++i){
    printf("[%d] Q_HEAD %d Q_TAIL %d\n",i,Q_HEAD,Q_TAIL);
    PUT_QUEUE_XYZ(i,i,i);
    SHOW_QUEUE_XYZ();
  }

  for (i=0;i<2*Q_MAX;++i){
    GET_QUEUE_XYZ(&x,&y,&z); printf("#get %d %d %d\n",x,y,z);
    SHOW_QUEUE_XYZ();
  }

  for (i=2*Q_MAX;i<5*Q_MAX;++i){
    printf("[%d] Q_HEAD %d Q_TAIL %d\n",i,Q_HEAD,Q_TAIL);
    PUT_QUEUE_XYZ(i,i,i);
    SHOW_QUEUE_XYZ();
  }

  for (i=0;i<2*Q_MAX;++i){
    GET_QUEUE_XYZ(&x,&y,&z); printf("#get %d %d %d\n",x,y,z);
    SHOW_QUEUE_XYZ();
  }

  free(QUEUE_Z);
  free(QUEUE_Y);
  free(QUEUE_X);
  return(1);
}
