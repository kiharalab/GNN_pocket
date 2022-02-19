#!/usr/bin/perl
#
# <sum_olg.pl>
#

$LastModDate = "Sep 2, 2008";

$OPT{'ipdb'} = 'ChConPDB';
$OPT{'opp'}  = 'GHECOM';
$OPT{'ilg'}  = 'RotPDB';
$OPT{'olg'}  = 'OLG';
$OPT{'op'}   = '-';
$OPT{'ol'}   = '-';
$OPT{'A'} = 'F';

if (scalar(@ARGV)<1){
 print "sum_olg.pl [listfile] <options>\n";
 print " for making summary of ghecom detected pockets\n";
 printf(" coded by T.Kawabata. LastModified:%s\n",$LastModDate);
 printf(" -olg  : dir for olg  [%s]\n",$OPT{'olg'}); 
 printf(" -lpat : three-letter ligand pattern  [%s]\n",$OPT{'lpat'}); 
 printf(" -op   : output pocket property filename  [%s]\n",$OPT{'op'}); 
 printf(" -ol   : output ligand property filename  [%s]\n",$OPT{'ol'}); 
 printf(" -A    : Action ('T' or 'F') [%s]\n",$OPT{'A'}); 
 exit(1);
}

$ilistfile = $ARGV[0];
&Read_Options(\@ARGV,\%OPT);
&Read_List_File($ilistfile,\@LIST);

$ofname = $OPT{'op'};
open(OP,">$ofname") || die "#ERROR:Can't write to '$ofname'";
printf(OP "#COMMAND %s\n",$OPT{'COMMAND'});
printf(OP "#DATE    %s\n",&Get_Date_String());
printf(OP "#Nlist   %d\n",scalar(@LIST));
close(OP);

$ofname = $OPT{'ol'};
open(OL,">$ofname") || die "#ERROR:Can't write to '$ofname'";
printf(OL "#COMMAND %s\n",$OPT{'COMMAND'});
printf(OL "#DATE    %s\n",&Get_Date_String());
printf(OL "#Nlist   %d\n",scalar(@LIST));

foreach $p (@LIST)
{
 $olgfile = $OPT{'olg'}.'/'.$p;
 &Read_Pocket_Annotated_Ligand_File($olgfile,\@liglist,\%ligdat);

 $Ncluster = scalar(@cluster);
 $Nligand  = scalar(@liglist);
 print OL "#$oppfile Ncluster $Ncluster Nligand $Nligand\n";
}

 &Write_Append_Pocket_Cluster_Information($OPT{'op'},\@cluster);

} ## $p ##

close(OL);

#################
### FUNCTIONS ###
#################



sub Read_Pocket_Annotated_Ligand_File{
 my($fname,$index_list,$lig) = @_;

#>> FILE FORMAT EXAMPLE <<
#REMARK  MULSC_PROBE: NRADIUS 27
#REMARK  MULSC_PROBE: RADIUS   1 th  2.000 Ang  0.500 1/Ang
#REMARK  MULSC_PROBE: RADIUS   2 th  2.500 Ang  0.400 1/Ang
#REMARK  MULSC_PROBE: RADIUS   3 th  3.000 Ang  0.333 1/Ang
#REMARK  MULSC_PROBE: RADIUS   4 th  3.500 Ang  0.286 1/Ang
#:
#REMARK  MULSC_PROBE: RADIUS  25 th 14.000 Ang  0.071 1/Ang
#REMARK  MULSC_PROBE: RADIUS  26 th 14.500 Ang  0.069 1/Ang
#REMARK  MULSC_PROBE: RADIUS  27 th 15.000 Ang  0.067 1/Ang
#REMARK  OUTSIDE_OF_MAX_RADIUS      15.500 Ang  0.065 1/Ang
#REMARK   For ligand atoms,
#REMARK    tFactor  [61-66]:average inverse of Rinaccess in the atom
#REMARK             [67-72]:harmonic means of Rinaccess. (-1 means infinity) in the atom
#REMARK             [73-78]:Shell Accessibility.(Number of noVdW grids within the shell)
#REMARK             [79-  ]:Cluster Num of Pocket included by the ligand atoms
#REMARK                                             [Radius][iRinacc][Rinacc][ShellAcc][ClusterNumbers]
#ATOM    944  OG1 THR A 147      23.071  37.918  63.916  1.40 0.000  0.00 1.000
#ATOM    945  CG2 THR A 147      25.110  37.769  65.177  1.87 0.000  0.00 1.000
#ATOM    946  N   ILE A 148      24.917  37.069  60.768  1.65 0.104  9.60 0.983
#ATOM    947  CA  ILE A 148      24.359  36.656  59.495  1.87 0.132  7.59 0.906
#ATOM    948  C   ILE A 148      24.136  35.153  59.585  1.76 0.158  6.31 0.929
#ATOM    949  O   ILE A 148      24.561  34.520  60.552  1.40 0.154  6.50 0.990
#ATOM   2252  CB  MET A 326       9.240  21.556  86.735  1.87 0.065 15.50 0.921
#ATOM   2253  CG  MET A 326      10.081  22.797  86.902  1.87 0.065 15.37 0.915
#ATOM   2254  SD  MET A 326      11.104  23.113  85.469  1.85 0.065 15.50 0.928
#ATOM   2255  CE  MET A 326      12.676  23.419  86.266  1.87 0.117  8.56 0.767
#HETATM 4281  P   AMP   600      17.232  41.145  52.290  1.90 0.397  2.52 0.563
#HETATM 4282  O1P AMP   600      16.879  42.563  52.262  1.40 0.423  2.36 0.572
#HETATM 4283  O2P AMP   600      18.602  40.844  51.782  1.40 0.400  2.50 0.649
#HETATM 4284  O3P AMP   600      16.956  40.517  53.584  1.40 0.450  2.22 0.520
#HETATM 4285  O5* AMP   600      16.218  40.461  51.324  1.40 0.354  2.83 0.762
#01234567890123456789012345678901234567890123456789012345678901234567890123456789
#          1         2         3         4         5         6         7     
 my($Resi); my($Atom); my($Rnum);  my($str); my(@cnum); my($c);
 my($F); my($index); my(%index_hash) = ();
 my($iRinacc); my($Rinacc); my($shellAcc);
 
 open(F,$fname) || die "#ERROR:Can't open pocket_annotated_ligand_pdbfile '$fname'";
 
 @{$index_list} = (); %{$lig} = ();
 
 while (<F>)
 {
  if ((/^HETATM/)||(/^ATOM/))
  {
   chomp;
   $Atom  = substr($_,13,3);
   $Resi  = substr($_,17,3);
   $Rnum  = substr($_,20,6);
   $index = $Resi.$Rnum;
   $iRinacc   = substr($_,60,6);
   $Rinacc    = substr($_,66,6); 
   $shellACC  = substr($_,72,6); 
   $str = substr($_,79);
   @cnum = split(/\s+/,$str);
   #print "#$_\n";
   foreach $c (@cnum)
   { 
   # print "str \"$str\" c $c\n";
    $lig->{$index}->{'clus_num'}->{$c} += 1; 
  }  

   $lig->{$index}->{'Natom'} += 1;
   $lig->{$index}->{'Rnum'} = $Rnum;
   $lig->{$index}->{'Resi'} = $Resi;
   if ($str ne '')
   {$lig->{$index}->{'Natom_pock'} += 1;}
 
   if ( $index_hash{$index} eq '' )
   {
    push(@{$index_list},$index);
    $index_hash{$index} = 1; 
   }

   ++$n; 
  } 
 }
 close(F);

} ## end of Read_Pocket_Annotated_Ligand_File() ##




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
#  if ($x =~/^\-/)
#   { $x =~s/^\-//;
#     if (($x1 !~ /^\-\w+/)&&(length($x1)>0)) { $_[1]->{$x} = $x1; ++$i; }
#     else { $_[1]->{$x} = 1;}
                                                                                                      
  if ($x =~/^\-/)
   { $x =~s/^\-//;
     if (length($x1)>0) { $_[1]->{$x} = $x1; ++$i; }
     else { $_[1]->{$x} = 1;}
                                                                                                      
   }
  ++$i;
 }
                                                                                                      
} ## end of Read_Options() ##


sub Read_List_File{
 my($fname,$list) = @_;
 my(@D);
                                                                                                      
 ## GENERAL FUNCTION FOR READING LIST FILE
 ## Read only 1st field.
 ## [FILE_EXAMPLE]
 ## 1dlwA  a.1.1.1
 ## 1a6m-  a.1.1.2
 ## 1eca-  a.1.1.2
 ## 1lh1-  a.1.1.2
 ## 1ash-  a.1.1.2
                                                                                                      
open(F,$fname) || die "#ERROR:Can't find listfile \"$fname\"";
 @{$list} = ();
                                                                                                      
 while (<F>)
 {
  if ($_!~/^#/)
  {
   chop;
   @D = split(/\s+/,$_);
   push(@{$list},$D[0]);
  }
 }
 close(F);
} ## end of Read_List_File()
                                                                                                      
                                                                                                      
sub Get_Date_String{
 my(@Month) = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
 my(@date) = localtime(time());
 my($year)  = $date[5]+1900;
 my($month) = $Month[$date[4]];
 my($day)   = $date[3];
 return("$month $day $year, $date[2]:$date[1]:$date[0]");

} ## end of Get_Date_String() ##

