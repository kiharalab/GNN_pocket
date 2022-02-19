#!/usr/bin/perl

$LastModDate = "Dec 12, 2007";

$OPT{'M'} = 'P';
$OPT{'min'} = 0;
$OPT{'oh'}  = 'hist.pdb';

if (scalar(@ARGV)<1)
{
 print "comp_grid.pl [grid_pdbfileA] [grid_pdbfileB] <options>\n";
 printf(" coded by T.Kawabata. LastModDate:%s\n",$LastModDate);
 printf("<options>\n");
 printf("-M   : mode 'P'air, 'L'ist [%s]\n",$OPT{'M'});
 printf("-il  : Input list file [%s]\n",$OPT{'il'});
 printf("-gd  : Grid PDB file directory [%s]\n",$OPT{'gd'});
 printf("-oh  : Output histgram PDB file [%s]\n",$OPT{'oh'});
 printf("-min : Minimum Observation Number for output [%s]\n",$OPT{'min'});
 exit(1);
}

&Read_Options(\@ARGV,\%OPT);

printf("#MODE %s\n",$OPT{'M'});
if ($OPT{'M'} eq 'P')
{
 $gpdbfileA = $ARGV[0]; 
 $gpdbfileB = $ARGV[1]; 
 &Read_Grid_PDB_File($gpdbfileA,\%GridA);
 &Read_Grid_PDB_File($gpdbfileB,\%GridB);
 &Count_Grid_Point(\%Hist,\%tFact,\%GridA);
 &Count_Grid_Point(\%Hist,\%tFact,\%GridB);
 &Write_Grid_Hist_PDB_File($OPT{'oh'},\%Hist,\%tFact);
}

if ($OPT{'M'} eq 'L')
{
 &Read_List_File($OPT{'il'},\@LIST,\%PROP);
 foreach $f (@LIST)
 {
  print "#$f\n"; 
  $fname = sprintf("%s/%s",$OPT{'gd'},$f);
  &Read_Grid_PDB_File($fname,\%Grid);
  &Count_Grid_Point(\%Hist,\%tFact,\%Grid);
 }
 
 &Write_Grid_Hist_PDB_File($OPT{'oh'},\@LIST,\%Hist,\%tFact,$OPT{'min'});
}

###################
#### FUNCTIONS ####
###################
sub Count_Grid_Point{
 my($hist,$tfact,$dat) = @_;
 my($xyz); 

 foreach $xyz (keys(%{$dat}))
 { 
  $hist->{$xyz}  += 1;
  $tfact->{$xyz} += $dat->{$xyz};
 }

} ## end of Count_Grid_Point() ##


sub Read_Grid_PDB_File{
  my($ifname,$dat) = @_;
#REMARK  DATE    Jul 20,2006 21:32:50
#REMARK  grid_size  339  315  308
#REMARK  grid_width    1.000
#REMARK  OrigPos  -34.822  -9.644  -13.773
#REMARK                                              [Radius][MapVal] [i] [j] [k]
#          10        20        30        40        50        60        70 
#01234567890123456789012345678901234567890123456789012345678901234567890123456789
#HETATM    1   C1 GRD     1      41.178 187.356 120.227  0.50  1.00 76 197 134
#HETATM    2   C1 GRD     1      41.178 188.356 119.227  0.50  1.00 76 198 133
#HETATM    3   C1 GRD     1      41.178 188.356 120.227  0.50  1.00 76 198 134
#HETATM    4   C1 GRD     1      41.178 188.356 121.227  0.50  1.00 76 198 135
#HETATM    5   C1 GRD     1      41.178 189.356 118.227  0.50  1.00 76 199 132
#HETATM    6   C1 GRD     1      41.178 189.356 119.227  0.50  1.00 76 199 133
#:
#HETATM36880  C11 GRD   116     223.178 138.356 106.227  0.50116.00 258 148 120
#HETATM36881  C11 GRD   116     223.178 138.356 107.227  0.50116.00 258 148 121
#HETATM36882  C11 GRD   116     223.178 138.356 108.227  0.50116.00 258 148 122
#HETATM36883  C11 GRD   116     223.178 139.356 107.227  0.50116.00 258 149 121
#TER
#END


 my($F); my($head); my($buff1); my($buff2); 
 my($value); 
 my($X); my($Y); my($Z); my($Occup); my($tFactor);
 my($XYZ);

 %{$dat} = ();
 open(F,$ifname)|| die "#ERROR:Can't open \"$ifname\"";;
 while (<F>)
 {
   chomp;
   if (/^HETATM/)
   {
    $val = substr($_,60,6);
    $val =~s/\s+//g;
    $X       = substr($_,30,8);
    $Y       = substr($_,38,8);
    $Z       = substr($_,46,8);
    $XYZ = $X.$Y.$Z;
    $Occup   = substr($_,54,6);
    $tFactor   = substr($_,60,6);
 
    #print "$_\n"; 
    #printf("(i,j,k)=($i,$j,$k) val=($val)\n",$i,$j,$k,$val);  
    $dat->{$XYZ} = $tFactor;
  #  printf("$XYZ tFactor $tFactor\n"); 
 } 
 } 
 close(F);

} ## end of Read_Grid_PDB_File() ##




sub Write_Grid_Hist_PDB_File{
  my($ofname,$List,$Hist,$tFact,$Nobs_min) = @_;
  my($OF); my($xyz);
  my($Nlist) = scalar(@{$List}); 
  my($Natom) = 0; my($p); my($i); my($avetfac);
  printf("#Write_Grid_Hist_PDB_File() --> '%s'\n",$ofname);
  open(OF,">$ofname") || die "#ERROR:Can't write to '$ofname'";
  printf(OF "HEADER  GRID HISTOGRAM FILE\n");
  printf(OF "REMARK  COMMAND %s\n",$OPT{'COMMAND'});
  printf(OF "REMARK  Nobs_min %s\n",$Nobs_min);
  printf(OF "REMARK  Nlist    %d\n",$Nlist);
  $i = 0;
  foreach $p (@{$List})
  {
   printf(OF "REMARK  LIST %3d %s\n",$i,$p);
   ++$i;  
  }
  printf(OF "REMARK                                                [Nobs][Nobs][ave_tFactor]\n");
  foreach $xyz (keys(%{$Hist}))
  {
     if ($Hist->{$xyz}>=$Nobs_min){
     ++$Natom;
     if ($Natom>99999) {$Natom = 1;}
     $avetfac = $tFact->{$xyz}/$Hist->{$xyz}; 
     printf(OF "HETATM%5d %4s %3s %s%5d   %24s%6.2f%6.2f%6.3f\n",
       $Natom," CA ","GRD"," ",$Natom/10, $xyz,$Hist->{$xyz},$Hist->{$xyz},$avetfac);
    }
 }
  close(OF); 
} ## end of Write_Hist_Grid_PDB_File() ##


sub Read_List_File{
 my($fname,$list,$prop) = @_;
 my($F); my($idstr); my($propstr); my($chain);
 my(@T);
 printf("#Read_List_File(\"$fname\")\n");
 @{$list} = ();
 open(F,$fname) || die "#ERROR:Can't open listfile \"$fname\"";
 while (<F>)
 {
  chomp;
  if ($_!~/^#/)
  {
  ($idstr,$propstr) = split(/\s+/,$_);
  # $idstr = substr($idstr,0,4);
  # $chain = substr($idstr,4,1);
   push(@{$list},$idstr);
    $prop->{$idstr} = $chain;
  }
 }
 close(F);
} ## end of Read_List_File() ##
                                                                                                                           
                                                                                                                           


sub Read_Options{
 # $_[0] : ref of \@ARGV
 # $_[1] : ref of \%OPT
 # This script is reading following style options :
 #   psiscan.pl org41list -lib 95pdb01Mar4Mx -tail -I -C
 # In principle, the format is the style like  "[-option] [value]"
 my($x); my($x1); my($i);
 $_[1]->{'COMMAND'} = "$0 ";
 foreach $x (@ARGV) { $_[1]->{'COMMAND'} .= "$x ";}
 $i = 0;
 while ($i<scalar(@ARGV))
 {
  $x  = $_[0]->[$i];
  $x1 = $_[0]->[$i+1];
  if ($x =~/^\-/)
   { $x =~s/^\-//;
     if (($x1 !~ /^\-\w+/)&&(length($x1)>0)) { $_[1]->{$x} = $x1; ++$i; }
   }
  ++$i;
 }
} ## end of Read_Options() ##

