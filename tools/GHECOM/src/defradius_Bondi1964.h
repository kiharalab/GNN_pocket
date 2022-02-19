/*
<defradius.h>

## These radius are take from
>> Van der Waals radii taken from Bondi's compilation (1964).[1] <<
Hydrogen(H) 1.20 
Carbon(C) 1.70 
Nitrogen(N) 1.55 
Oxygen(O) 1.52 
Fluorine(F) 1.47 
Phosphorus(P) 1.80 
Sulfur(S) 1.80 
Chlorine(Cl) 1.75 
Copper(Cu) 1.4 

Bondi, A. (1964). "Van der Waals Volumes and Radii". J. Phys. Chem. 68 (3): 441?51. doi:10.1021/j100785a001


*/

#define NRADLINE 16 

static char RADLINE[16][14]=
{
"xHxx xxx 1.20",
"xCxx xxx 1.70",
"xNxx xxx 1.55",
"xOxx xxx 1.52",
"xFxx xxx 1.47",
"xSxx xxx 1.80",
"xPxx xxx 1.80",
"xIxx xxx 1.98",
"xKxx xxx 2.75",
"xCLx xxx 1.75",
"xCUx xxx 1.40",
"xBRx xxx 1.85",
"xZNx xxx 1.39",
"xMGx xxx 1.73",
"xNAx xxx 2.27",
"xxxx xxx 1.70"
};
