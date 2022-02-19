#!/usr/bin/perl

if (scalar(@ARGV)<1)
{
 printf("mkrandpnt.pl [N] [R]\n");
 exit(1);
}

$N    = $ARGV[0];
$R    = $ARGV[1];
$Type = $ARGV[2];

print "HEADER   Type $Type RANDOM SPHERICAL POINTS\n";
print "REMARK N $N R $R\n";
for ($i=0;$i<$N;++$i){
 if ($Type eq 'U'){
  $ok = 0;
  while ($ok==0){ 
   $x = rand(2)-1;
   $y = rand(2)-1;
   $z = rand(2)-1;
   $len = $x*$x + $y*$y + $z*$z;
   if ($len>0){$len = sqrt($len);}
   if ($len<=1.0){$ok=1;} 
   else {print "WOOPS!!\n";}
  }
 }
 else{
  $x = rand(2)-1;
  $y = rand(2)-1;
  $z = rand(2)-1;
  $len = $x*$x + $y*$y + $z*$z;
  if ($len>0){$len = sqrt($len);}
 }
  $x *= $R/$len; $y *= $R/$len; $z *= $R/$len;
  $X[$i] = $x;
  $Y[$i] = $y;
  $Z[$i] = $z;
# printf("$x $y $z\n");
# printf("HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
#       $i," C  ","PNT",1, $x,$y,$z,0.0,0.0);
}

for ($i=0;$i<$N;++$i){

 $init = 1;
 for ($j=0;$j<$N;++$j){
  if ($i!=$j){ 
  $DD = 
  ($X[$i]-$X[$j])*($X[$i]-$X[$j])
 +($Y[$i]-$Y[$j])*($Y[$i]-$Y[$j])
 +($Z[$i]-$Z[$j])*($Z[$i]-$Z[$j]);
 if (($init==1)||($DD<$minDD)) {$minDD = $DD; $init = 0;}
 # print "DD $DD minDD $minDD\n"; 
 } 
 }
 $D = sqrt($minDD);
 $Dnei[$i] = $D; 
 if (($i==0)||($D>$maxD)) {$maxD = $D;}
 #printf("i %d Dnei %f\n",$i,$Dnei[$i]);
}

for ($i=0;$i<$N;++$i){
 printf("HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
       $i," C  ","PNT",1, $X[$i],$Y[$i],$Z[$i],$Dnei[$i],$Dnei[$i]);

}

printf("HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
       $i," O  ","PNT",1,0,0,0,0.0,0.0); ++$i;

printf("HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
       $i," O  ","PNT",1, $R,0,0,0.0,0.0); ++$i;
printf("HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
       $i," O  ","PNT",1, -$R,0,0,0.0,0.0); ++$i;

printf("HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
       $i," O  ","PNT",1, 0,$R,0,0.0,0.0); ++$i;
printf("HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
       $i," O  ","PNT",1, 0,-$R,0,0.0,0.0); ++$i;

printf("HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
       $i," O  ","PNT",1, 0,0,$R,0.0,0.0); ++$i;
printf("HETATM%5d %4s %3s  %5d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
       $i," O  ","PNT",1, 0,0,-$R,0.0,0.0); ++$i;

