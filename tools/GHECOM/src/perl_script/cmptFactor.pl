#!/usr/bin/perl

$LastModDate = "Nov 4, 2006";

$OPT{'or'} = '';
$OPT{'op'} = '';
$OPT{'od'} = '';
$OPT{'tA'} = '1';
$OPT{'tB'} = '1';
if (scalar(@ARGV)<1)
{
 printf("cmptFactor.pl [pdbfileA] [pdbfileB] <options>\n");
 printf(" for comparison of tFactor values\n");
 printf(" coded by T.Kawabata. LastModified:%s\n",$LastModDate);
 printf("<options>\n");
 printf(" -or : outfile for raw pairwise values[%s]\n",$OPT{'or'});
 printf(" -od : outfile for difference staistics[%s]\n",$OPT{'od'});
 printf(" -os : outfile for pairwise staistics[%s]\n",$OPT{'os'});
 printf(" -oA : outfile for histgram A [%s]\n",$OPT{'oA'});
 printf(" -oB : outfile for histgram B [%s]\n",$OPT{'oB'});
 printf(" -tA : tFactor Column num for pdbfileA (1,2,3) [%s]\n",$OPT{'tA'});
 printf(" -tB : tFactor Column num for pdbfileB (1,2,3) [%s]\n",$OPT{'tB'});
 exit(1);
}

$pdbfileA = $ARGV[0];
$pdbfileB = $ARGV[1];
&Read_Options(\@ARGV,\%OPT);

&Read_PDB_File_Array($pdbfileA,\@atmA,$OPT{'tA'});
printf("#NatomA %d\n",scalar(@atmA));
&Read_PDB_File_Array($pdbfileB,\@atmB,$OPT{'tB'});
printf("#NatomB %d\n",scalar(@atmB));

for ($a=0;$a<scalar(@atmA);++$a)
{
 $ExistA{$atmA[$a]->{'tFactor'}} += 1;
 $ExistB{$atmB[$a]->{'tFactor'}} += 1;
 $diff = $atmA[$a]->{'tFactor'} - $atmB[$a]->{'tFactor'};
 $DiffCount{$diff} += 1;
 $Count{$atmA[$a]->{'tFactor'}}{$atmB[$a]->{'tFactor'}} += 1;
}

### [1] OUTPUT PAIRWISE COMPARISON ###
if ($OPT{'or'} ne '')
{
$ofile = $OPT{'or'};
printf("#Output_Raw_Pairwise_Data() --> '%s'\n",$ofile);
open(OF,">$ofile")||die "#ERROR:Can't write to $ofile";
for ($a=0;$a<scalar(@atmA);++$a)
{
 printf(OF "%f %f %s %s %s %s %s %s %s %s %s %s\n",
$atmA[$a]->{'tFactor'},
$atmB[$a]->{'tFactor'}, 
    $atmA[$a]->{'atmNum'},
    $atmA[$a]->{'resName'},
    $atmA[$a]->{'atmName'},
    $atmA[$a]->{'chainID'},
    $atmA[$a]->{'resNum'},
    $atmB[$a]->{'atmNum'},
    $atmB[$a]->{'resName'},
    $atmB[$a]->{'atmName'},
    $atmB[$a]->{'chainID'},
    $atmB[$a]->{'resNum'});

 }
 close(OF);
}

### [2] OUTPUT DIFFERENCE STATISTICS ###
if ($OPT{'od'} ne '')
{
 $ofile = $OPT{'od'};
 printf("#Output_Difference_Statistics() --> '%s'\n",$ofile);
 open(OF,">$ofile")||die "#ERROR:Can't write to '$ofile'";
 @dlist = keys(%DiffCount);
 @sdlist = sort {$a <=> $b} @dlist;
 foreach $x (@sdlist)
 {
  printf(OF "%s %f %d\n",$x,100.0*$DiffCount{$x}/scalar(@atmA),$DiffCount{$x});
 }
 close(OF);
}

### [3] OUTPUT PAIRWISE STATISTICS ###
if ($OPT{'op'} ne '')
{
 $ofile = $OPT{'op'};
 printf("#Output_Paireise_Statistics() -> '%s'\n",$ofile);
 open(OF,">$ofile")||die "#ERROR:Can't write to '$ofile'";
 @Alist = keys(%ExistA);
 @sAlist = sort {$a <=> $b} @Alist;
 @Blist = keys(%ExistB);
 @sBlist = sort {$a <=> $b} @Blist;
 foreach $a (@sAlist){
  printf(OF "#tFacA %s Count %d\n",$a,$ExistA{$a}); 
  printf(OF "#[tFacA] [tFacB] [Ratio] [Nobs]\n"); 
  foreach $b (@sBlist){
  printf(OF "%s %s %7.3f %d\n",$a,$b,100.0*$Count{$a}{$b}/$ExistA{$a},$Count{$a}{$b});
  }
  printf(OF "\n");
 }
 close(OF);
}

### [4] OUTPUT HISTGRAM FOR EACH tFACTOR ###
if ($OPT{'oA'} ne ''){ &Output_Histgram($OPT{'oA'},\%ExistA,$pdbfileA);}
if ($OPT{'oB'} ne ''){ &Output_Histgram($OPT{'oB'},\%ExistB,$pdbfileB);}

##################
### FUNCTIONS ####
##################

sub Output_Histgram{
 my($ofname,$Count,$comment) = @_;
 my($OF); my($Nall)=0;
 my(@indlist) = keys(%{$Count});
 my(@sindlist) = sort {$a <=> $b} @indlist;
 my($x);
 foreach $x (@indlist) {$Nall += $Count->{$x};}
 
 printf("#Output_Histgram() -->'%s'\n",$ofname);
 open(OF,">$ofname") || die "#ERROR:Can't write to '$ofname'";
 printf(OF "#COMMAND %s\n",$OPT{'COMMAND'});
 printf(OF "#Nall    %d\n",$Nall);
 printf(OF "#COMMENT %s\n",$comment);
 foreach $x (@sindlist)
 { 
  printf(OF "%s %f %d\n",$x,100.0*$Count->{$x}/$Nall,$Count->{$x}); 
 }
 close(OF);

} ## end of Output_Histgram() ##


sub Read_PDB_File_Array{
 my($fname,$dat,$tFactorType) = @_;
 my($IF);
 @{$dat} = ();
 open(IF,$fname) || die "#ERROR:Can't open pdbfile '$fname'\n";
 my($atmName); my($atmNum); my($resName); my($resNum);
 my($chainID); my($X); my($Y); my($Z);  my($Occup); my($tFactor);
 my($tFactor2); my($tFactor3); 
 my($Natom) = 0; 
 while (<IF>)
 {
  chomp;
  if (($_=~/^ATOM/)||($_=~/^HETATM/)){
    $atmNum  = substr($_,6,5);
    $atmName = substr($_,12,4);
    $resName = substr($_,17,3);
    $chainID = substr($_,21,1);
    $resNum  = substr($_,22,4);
    $X       = substr($_,30,8);
    $Y       = substr($_,38,8);
    $Z       = substr($_,46,8);
    $Occup   = substr($_,54,6);
    $tFactor   = substr($_,60,6);
    $tFactor2  = substr($_,66,6);
    $tFactor3  = substr($_,72,6);
    $dat->[$Natom]->{'atmNum'} = $atmNum;
    $dat->[$Natom]->{'atmName'} = $atmName;
    $dat->[$Natom]->{'resName'} = $resName;
    $dat->[$Natom]->{'atmName'} = $atmName;
    $dat->[$Natom]->{'chainID'} = $chainID;
    $dat->[$Natom]->{'resNum'} = $resNum;
    $dat->[$Natom]->{'X'} = $X;
    $dat->[$Natom]->{'Y'} = $Y;
    $dat->[$Natom]->{'Z'} = $Z;
    $dat->[$Natom]->{'Occup'} = $Occup;
    if ($tFactorType eq '1') {$dat->[$Natom]->{'tFactor'}  = $tFactor;}
    if ($tFactorType eq '2') {$dat->[$Natom]->{'tFactor'}  = $tFactor2;}
    if ($tFactorType eq '3') {$dat->[$Natom]->{'tFactor'}  = $tFactor3;}
    ++$Natom;
   }
 }
 close(IF);
} ## end of Read_PDB_File_Array() ##


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

