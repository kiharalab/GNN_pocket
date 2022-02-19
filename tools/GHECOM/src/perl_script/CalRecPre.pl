#!/usr/bin/perl
#
# <CalRecPre.pl>
#

$LastModDate = "Oct 1, 2009";

if (scalar(@ARGV)<1){
  print "CalRecPre.pl [grid_lig_comp_file]\n";
  print " coded by T.Kawabata. LastModDate:$LastModDate\n";
  print ">>Sample of grid_lig_com_file<<\n";
  print "GRID_LIG_COMP RAYGRID_B7/1jnrA tmp.pdb threA 1 255 NP 61655 NL 1159 NPL 1137 Rec 0.981018 Pre 0.018441 Fme 0.036202\n";
  print "GRID_LIG_COMP RAYGRID_B7/1dctA tmp.pdb threA 1 255 NP 32826 NL  909 NPL  607 Rec 0.667767 Pre 0.018491 Fme 0.035986\n";
  exit(1);
}


$NP = $NL = $NPL = $Nchain = 0;
$ifname = $ARGV[0];
open(F,$ifname)||die "#ERROR:Can't open '$ifname'";

while (<F>){
#GRID_LIG_COMP OPO/1rhsA OLG/1rhsA threA 1 NP 114 NL 0 NPL 0 Rec 0.000000 Pre 0.000000 Fme 0.000000
#GRID_LIG_COMP RAYGRID_B7/1jnrA tmp.pdb threA 1 255 NP 61655 NL 1159 NPL 1137 Rec 0.981018 Pre 0.018441 Fme 0.036202
#GRID_LIG_COMP RAYGRID_B7/1dctA tmp.pdb threA 1 255 NP 32826 NL  909 NPL  607 Rec 0.667767 Pre 0.018491 Fme 0.035986
 if ($_!~/^#/){
   chomp;
   $Nchain += 1;
   @field = split(/\s+/,$_); 
   $thre = $field[4];
   if ($field[5] eq 'NP'){ 
    $np  = $field[6];
    $nl  = $field[8];
    $npl = $field[10];
   }
   if ($field[6] eq 'NP'){ 
    $np  = $field[7];
    $nl  = $field[9];
    $npl = $field[11];
   }
   #print "$np $nl $npl\n";
   $NP  += $np;
   $NL  += $nl;
   $NPL += $npl;
 }
}

$rec = $NPL/$NL;
$pre = $NPL/$NP;

$fme = 2.0/(1.0/$rec + 1.0/$pre);

printf("rec $rec pre $pre fme $fme thre $thre NP $NP NL $NL NPL $NPL Nchain $Nchain NPperch %f %s\n",$NP/$Nchain,$ifname);

