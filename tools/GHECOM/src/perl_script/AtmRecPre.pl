#!/usr/bin/perl
#
# <AtmRecPre.pl>
#

$LastModDate = "Sep 18, 2005";

if (scalar(@ARGV)<1)
{
 print "AtmRecPre.pl <options>\n";
 printf(" coded by T.Kawabata. LastModified:$LastModDate\n");
 printf(" for comparing selected atom regions 'A' and 'B'\n");
 printf("<options>\n");
 printf("*for atoms 'A'\n");
 printf(" -iA : interface_PDBfile[%s]\n",$OPT{'iA'});
 printf(" -cA : Chain ID for iA [%s]\n",$OPT{'cA'});
 printf(" -pA : PDBfile[%s]\n",$OPT{'pA'});
 printf("*for atoms 'B'\n");
 printf(" -pB : PDBfile[%s]\n",$OPT{'pB'});
 printf("*for all the atoms X\n");
 printf(" -aX : all atoms in interface_PDBfile[%s]\n",$OPT{'aX'});
 printf(" -sX : all surface atoms in interface_PDBfile[%s]\n",$OPT{'sX'});
 printf(" -cX : Chain ID for X [%s]\n",$OPT{'cX'});
 printf("*other options\n"); 
 printf(" -or : output rasmol Coloring script file[%s]\n",$OPT{'or'});


 exit(1);
}

&Read_Options(\@ARGV,\%OPT);

### Reading atoms 'A' ###
if ($OPT{'iA'} ne '')
{ &Read_Interface_PDB_File($OPT{'iA'},$OPT{'cA'},\@anumlistA,\%atomA,'I'); $fnameA = $OPT{'iA'};}
if ($OPT{'pA'} ne '')
{ &Read_PDB_File_Surrounding_Pocket($OPT{'pA'},\@anumlistA,\%atomA); $fnameA = $OPT{'pA'};}
$NatomA = scalar(@anumlistA);
print "NatomA $NatomA\n";

### Reading atoms 'B' ###
if ($OPT{'iB'} ne '')
{ &Read_Interface_PDB_File($OPT{'iB'},$OPT{'cB'},\@anumlistB,\%atomB); $fnameB = $OPT{'iB'};}
if ($OPT{'pB'} ne '')
{ &Read_PDB_File_Surrounding_Pocket($OPT{'pB'},\@anumlistB,\%atomB); $fnameB = $OPT{'pB'};}
$NatomB = scalar(@anumlistB);
print "NatomB $NatomB\n";

### Reading all atoms 'X' ###
if ($OPT{'aX'} ne '')
{ &Read_Interface_PDB_File($OPT{'aX'},$OPT{'cX'},\@anumlistX,\%atomX,'A'); }
if ($OPT{'sX'} ne '')
{ &Read_Interface_PDB_File($OPT{'sX'},$OPT{'cX'},\@anumlistX,\%atomX,'S'); }
if (($OPT{'aX'} ne '')||($OPT{'sX'} ne ''))
{
 $NatomX = scalar(@anumlistX);
 print "NatomX $NatomX\n";
 foreach $anum (@anumlistA) 
  { $atomX{$anum}->{'A'} = 1; 
    if ($atomX{$anum}->{'anum'} eq '') 
    {printf("#WARNING:anumA \"$anum\" does not exist in X\n");}
  }

 foreach $anum (@anumlistB) 
  { $atomX{$anum}->{'B'} = 1; 
    ++$NNN;
   if ($atomX{$anum}->{'anum'} eq '') 
    {printf("#WARNING:anumB \"$anum\" does not exist in X\n");}
  }
 print "NNN $NNN\n";
}


## COUNTING INTERSECTION AB ##
$NAB = 0;
foreach $anum (@anumlistA)
{ if ($atomB{$anum}->{'anum'} ne '') { ++$NAB; push(@anumlistAB,$anum); } 
  else { push(@anumlistOnlyA,$anum); } }

foreach $anum (@anumlistB)
{ if ($atomA{$anum}->{'anum'} eq '') { push(@anumlistOnlyB,$anum); } }

$NonlyA = $NatomA - $NAB;
$NonlyB = $NatomB - $NAB;
printf("$fnameA $fnameB NatomA $NatomA NatomB $NatomB NAB $NAB NonlyA $NonlyA NonlyB $NonlyB ");
$Precision  = $NAB/$NatomB;
$Recall     = $NAB/$NatomA;
printf("Pre %.3f Rec %.3f\n",$Precision,$Recall);

if (($OPT{'aX'} ne '')||($OPT{'sX'} ne ''))
{ &Cal_CorrCoeff(\@anumlistX,\%atomX); }

if ($OPT{'or'} ne '')
{
 &Output_Rasmol_Coloring_Script(
  $OPT{'or'},\@anumlistAB,\@anumlistOnlyA,\@anumlistOnlyB,$fnameA,$fnameB);
}

#################
### FUNCTIONS ###
#################

sub Cal_CorrCoeff{
 my($listX,$datX) = @_;
 
 my($anum); my($A); my($B); 
 my($Pre); my($Rec); my($CC);
 my($N00); my($N01); my($N10); my($N11); 
 my($bunshi); my($bunbo);
 my($NA); my($NB); my($N);

 $N = $N00 = $N01 = $N10 = $N11 = $NA = $NB = 0;
 foreach $anum (@{$listX})
 {
  $A = $datX->{$anum}->{'A'};
  $B = $datX->{$anum}->{'B'};
#  printf("$anum A %d B %d\n",$A,$B);
  if (($A==0)&&($B==0)) {++$N00;}
  if (($A==0)&&($B==1)) {++$N01;}
  if (($A==1)&&($B==0)) {++$N10;}
  if (($A==1)&&($B==1)) {++$N11;}
  if ($A==1) {++$NA;}
  if ($B==1) {++$NB;}
  ++$N; 
 }

 $Pre = $N11/($N11+$N01);
 $Rec = $N11/($N11+$N10);

 $bunshi = $N00*$N11 - $N01*$N10;
 $bunbo  = ($N00+$N01)*($N00+$N10)*($N11+$N01)*($N11+$N10);
 if ($bunbo>0.0) {$CC = $bunshi/sqrt($bunbo);}
            else {$CC = 0.0;}


 printf("N $N NA $NA NB $NB Pre %.3f Rec %.3f CC %.3f\n",$Pre,$Rec,$CC);
 printf("N00 $N00 N01 $N01\n");
 printf("N10 $N10 N11 $N11\n");

} ## end of Cal_CorrCoeff() ##



sub Output_Rasmol_Coloring_Script{
 my($ofname,$listAB,$listOnlyA,$listOnlyB,$fnameA,$fnameB) = @_;
 my($OF); my($anum);
 my($colorAB)    = "green";
 my($colorOnlyA) = "blue";
 my($colorOnlyB) = "red";
 
 open(OF,">$ofname") || die "#ERROR:Can't write to \"$ofname\"";
 printf(OF "echo \"## atomA:%s\"\n",$fnameA);
 printf(OF "echo \"## atomB:%s\"\n",$fnameB);
 printf(OF "echo \"## region AandB:color %s\"\n",$colorAB);
 printf(OF "echo \"## region onlyA:color %s\"\n",$colorOnlyA);
 printf(OF "echo \"## region onlyB:color %s\"\n",$colorOnlyB);
 print OF "select all\n";
 print OF "color gray\n";
 $str = &Rasmol_Select_Atomno_String(\@{$listAB},$colorAB);
 print OF "$str\n";
 $str = &Rasmol_Select_Atomno_String(\@{$listOnlyA},$colorOnlyA);
 print OF "$str\n";
 $str = &Rasmol_Select_Atomno_String(\@{$listOnlyB},$colorOnlyB);
 print OF "$str\n";
 print OF "select all\n";
 print OF "spacefill\n";

 close(OF);

} ## end of Output_Rasmol_Coloring_Script() ##



sub Rasmol_Select_Atomno_String{
 my($list,$color) = @_;
 my($str) = '';
 my($i); my($j); my($anum);

 $str = sprintf("select "); 
 $j = 0;
 for ($i=0;$i<scalar(@{$list});++$i)
 { if ((($i%10)==0)&&($i>0)) { $str .= sprintf("\ncolor %s\nselect ",$color); $j=0;}
   if ($j>0) { $str .= sprintf(",");}
   $anum = $list->[$i];
   $anum =~s/\s+//g; 
  $str .= sprintf("atomno=%s",$anum);
  ++$j; }

 $str .= sprintf("\ncolor %s\n",$color);

 return($str);

} ## end of Rasmol_Select_Atomno_String() ##




sub Read_Interface_PDB_File{
 my($ifname,$SpecChainID,$anumlist,$dat,$SelectType) = @_;
 # $SelectType : 'I'nterface, 'S'urface, 'A'll 

 print "#Read_Interface_PDB_File(\"$ifname\" SelectType $SelectType)\n";
 @{$anumlist} = (); %{$dat} = ();

 #<FILE EXAMPLE>
 #REMARK  COMMAND "LeeRich /DB/PDB/er/pdb3er5.ent -pA E:xxx:xxx -pB I:xxx:xxx -M I"
 #REMARK  DATE Sep 17,2005 14:37:23
 #REMARK  RadiusFile ""
 #REMARK  Rsolvent[A] 1.400000 deltaZ[A] 0.100000
 #REMARK  PatternA "E:xxx:xxx" Natom 2389 Nres 330
 #REMARK  PatternB "I:xxx:xxx" Natom   80 Nres   9
 #REMARK REIN : Region('A' or 'B') and Interface (' ' or 'I')
 #REMARK                 RNUM      X       Y       Z     dASA  dACC MonASA CmpASA  ReIn Radius ASP
 #ATOM      1  N   SER E  -2      27.709  33.477  -4.060  0.00   0.0 15.37 15.37 A     1.65 -0.006
 #ATOM      2  CA  SER E  -2      27.934  32.889  -2.722  0.00   0.0 13.98 13.98 A     1.87  0.016
 #ATOM      3  C   SER E  -2      26.749  32.046  -2.297  0.00   0.0  0.00  0.00 A     1.76  0.016
 #:
 #ATOM     95  N   ASP E  12       8.550  21.366  18.387  0.00   0.0  0.00  0.00 A I   1.65 -0.006
 #ATOM     96  CA  ASP E  12       7.715  22.565  18.540  0.00   0.0  0.00  0.00 A I   1.87  0.016
 #ATOM     97  C   ASP E  12       8.221  23.697  17.661  0.04   0.0  0.04  0.00 A I   1.76  0.016
 #ATOM     98  O   ASP E  12       7.506  24.541  17.130  4.12   4.2  4.67  0.54 A I   1.40 -0.006
 #ATOM     99  CB  ASP E  12       7.789  23.035  20.024 10.28   7.7 10.43  0.15 A I   1.87  0.016
 #ATOM    100  CG  ASP E  12       6.666  22.303  20.761  7.70   6.2  7.70  0.00 A I   1.76  0.016
 #ATOM    101  OD1 ASP E  12       5.505  22.268  20.324 15.63  15.9 15.63  0.00 A I   1.40 -0.024
 #ATOM    102  OD2 ASP E  12       6.868  21.717  21.833 10.83  11.0 12.46  1.62 A I   1.40 -0.024
 #ATOM    103  N   ALA E  13       9.526  23.794  17.544  0.00   0.0  0.00  0.00 A I   1.65 -0.006
 #ATOM    104  CA  ALA E  13      10.238  24.774  16.726  0.00   0.0  0.91  0.91 A I   1.87  0.016
 my($IF); my($anum); my($Atom); my($Res); my($ChainID); my($Rnum); 
 my($dASA); my($dACC); my($MonASA); my($CmpASA); my($accept);

 open(IF,$ifname) || die "#ERROR:Can't open interpdbfile \"$ifname\"";
 while (<IF>)
 {
  chomp;
  if (/^ATOM/)
  {
 #Get_Part_Of_Line(Aname,line,12,15);
 #Get_Part_Of_Line(Rname,line,17,19);
 #Get_Part_Of_Line(atminfo,line,12,25);
  $anum    = substr($_,6,5);
  $Atom    = substr($_,12,4);
  $Res     = substr($_,17,3);
  $ChainID = substr($_,21,1);
  $Rnum    = substr($_,22,4);
  $dASA    = substr($_,54,6);
  $dACC    = substr($_,60,6);
  $MonASA  = substr($_,66,6);
  $CmpASA  = substr($_,72,6);
  $accept = 0; 
  if (($SelectType eq 'I') && ($dASA>0.0))   {$accept = 1;}
  if (($SelectType eq 'S') && ($MonASA>0.0))  {$accept = 1;}
  if ($SelectType eq 'A')  {$accept = 1;}

  if (($accept==1) && (($SpecChainID eq '-')||($ChainID eq $SpecChainID)) )
  {
  #  print "$_\n"; 
    push(@{$anumlist},$anum);  
    $dat->{$anum}->{'anum'} = $anum;
    $dat->{$anum}->{'Atom'} = $Atom;
    $dat->{$anum}->{'Res'}  = $Res;
    $dat->{$anum}->{'ChainID'}  = $ChainID;
    $dat->{$anum}->{'Rnum'}  = $Rnum;
    $dat->{$anum}->{'dASA'}  = $dASA;
    $dat->{$anum}->{'dACC'}  = $dACC;
    $dat->{$anum}->{'MonASA'}  = $MonASA;
    $dat->{$anum}->{'CmpASA'}  = $CmpASA;
   }

 } # if ATOM # 

 } # while #

 close(IF);

} ## end of Read_Interface_PDB_File() ##



 
sub Read_PDB_File_Surrounding_Pocket{
 my($ifname,$anumlist,$dat) = @_;

 print "#Read_PDB_File_Surronding_Pocket(\"$ifname\")\n";
 @{$anumlist} = (); %{$dat} = ();

 #<FILE EXAMPLE>
 #HEADER    PROBE CLUSTER ATOMS AND THEIR CONTACTING ATOMS
 #REMARK     FILENAME "3er5E.cpdb"
 #REMARK     COMMAND "phecom 3er5E.pdb -osw 3er5E_phecom6.wrl -rl 6 -rgb 1 0.5 0.5 0.5"
 #REMARK     DATE Sep 16,2005 11:56:24
 #REMARK     Ncluster 12
 #ATOM   1615  C   GLY E 217       6.364  29.279  13.308  1.76  0.50
 #ATOM   1616  O   GLY E 217       6.227  29.640  14.473  1.40  0.84
 #ATOM    254  OD1 ASP E  32       7.942  32.651  11.425  1.40  0.96
 #ATOM   1618  CA  THR E 218       4.304  27.919  13.303  1.87  1.00
 #ATOM   1604  OD1 ASP E 215       3.883  31.341  10.710  1.40  1.51


 my($IF); my($anum); my($Atom); my($Res); my($ChainID); my($Rnum); my($dASA);

 open(IF,$ifname) || die "#ERROR:Can't open interpdbfile \"$ifname\"";
 while (<IF>)
 {
  chomp;
  if (/^ATOM/)
  {
 #Get_Part_Of_Line(Aname,line,12,15);
 #Get_Part_Of_Line(Rname,line,17,19);
 #Get_Part_Of_Line(atminfo,line,12,25);
  $anum    = substr($_,6,5);
  $Atom    = substr($_,12,4);
  $Res     = substr($_,17,3);
  $ChainID = substr($_,21,1);
  $Rnum    = substr($_,22,4);
  if ($dat->{$anum}->{'anum'} eq '')
  { 
  push(@{$anumlist},$anum);  
  $dat->{$anum}->{'anum'} = $anum;
  $dat->{$anum}->{'Atom'} = $Atom;
  $dat->{$anum}->{'Res'}  = $Res;
  $dat->{$anum}->{'ChainID'}  = $ChainID;
  $dat->{$anum}->{'Rnum'}  = $Rnum;
  $dat->{$anum}->{'dASA'}  = $dASA;
  }
 } # if ATOM # 

 } # while #

 close(IF);

} ## end of Read_PDB_File_Surrounding_Pocket() ##


sub Read_Options{
 # $_[0] : ref of \@ARGV
 # $_[1] : ref of \%OPT
 # This script is reading following style options :
 #   psiscan.pl org41list -lib 95pdb01Mar4Mx -tail -I -C
 # In principle, the format is the style like  "[-option] [value]"
 # If [value] is omitted, [option] is set to '1'.
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
     else { $_[1]->{$x} = 1;}
   }
  ++$i;
 }
                                                                                                      
} ## end of Read_Options() ##


