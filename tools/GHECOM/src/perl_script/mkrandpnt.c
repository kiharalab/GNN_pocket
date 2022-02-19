#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_PNT 10000

main(argc,argv)
 int argc;
 char **argv;
{
 int N;
 float R;
 char Type;
 int i,j,ok;
 float x,y,z,len,D,DD,minDD;
 float X[MAX_PNT],Y[MAX_PNT],Z[MAX_PNT],Dnei[MAX_PNT];
 char init;

 if (argc<2){
  printf("mkrandpnt.pl [N] [R]\n");
  exit(1);
 }


N    = atoi(argv[1]);
R    = atof(argv[2]);
Type = argv[3][0];

printf("HEADER   Type %c RANDOM SPHERICAL POINTS\n",Type);
printf("REMARK N %d R %f\n",N,R);
for (i=0;i<N;++i){
 if (Type == 'U'){
  ok = 0;
  while (ok==0){ 
   x = 2.0*((float)rand()/RAND_MAX - 0.5);
   y = 2.0*((float)rand()/RAND_MAX - 0.5);
   z = 2.0*((float)rand()/RAND_MAX - 0.5);
   len = x*x + y*y + z*z;
   if (len>0){len = sqrt(len);}
   if (len<=1.0){ok=1;} 
  }
 }
 else{
  x = 2.0*((float)rand()/RAND_MAX - 0.5);
  y = 2.0*((float)rand()/RAND_MAX - 0.5);
  z = 2.0*((float)rand()/RAND_MAX - 0.5);
  len = x*x + y*y + z*z;
  if (len>0){len = sqrt(len);}
 }
  /* printf("xyz %f %f %f\n",x,y,z); */
  x *= R/len; y *= R/len; z *= R/len;
  X[i] = x;
  Y[i] = y;
  Z[i] = z;
}

for (i=0;i<N;++i){

 init = 1;
 for (j=0;j<N;++j){
  if (i!=j){ 
  DD = 
  (X[i]-X[j])*(X[i]-X[j])
 +(Y[i]-Y[j])*(Y[i]-Y[j])
 +(Z[i]-Z[j])*(Z[i]-Z[j]);
 if ((init==1)||(DD<minDD)) {minDD = DD; init = 0;}
 } 
 }
 D = sqrt(minDD);
 Dnei[i] = D; 
}

for (i=0;i<N;++i){
 printf("HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
       i," C  ","PNT",1, X[i],Y[i],Z[i],Dnei[i],Dnei[i]);

}

printf("HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
       i," O  ","PNT",1,0,0,0,0.0,0.0); ++i;

printf("HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
       i," O  ","PNT",1, R,0,0,0.0,0.0); ++i;
printf("HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
       i," O  ","PNT",1, -R,0,0,0.0,0.0); ++i;

printf("HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
       i," O  ","PNT",1, 0,R,0,0.0,0.0); ++i;
printf("HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
       i," O  ","PNT",1, 0,-R,0,0.0,0.0); ++i;

printf("HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
       i," O  ","PNT",1, 0,0,R,0.0,0.0); ++i;
printf("HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
       i," O  ","PNT",1, 0,0,-R,0.0,0.0); ++i;
}
