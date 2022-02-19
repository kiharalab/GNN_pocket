#!/usr/bin/perl
#
# <sumolg_DABRL.pl>
#

$LastModDate = "Dec 10, 2008";

$OPT{'olg'}  = 'OLG';
$OPT{'op'}   = '-';
$OPT{'lpat'}  = '.';
$OPT{'olpdb'}  = 'olig.pdb';
$OPT{'idab'} = '';

if (scalar(@ARGV)<1){
 print "sumolg_DABRL.pl [listfile] <options>\n";
 print " for making averaged Rinacc value for each ligand atom type\n";
 printf(" coded by T.Kawabata. LastModified:%s\n",$LastModDate);
 printf(" -olg  : dir for olg  [%s]\n",$OPT{'olg'}); 
 printf(" -lpat : three-letter ligand pattern  [%s]\n",$OPT{'lpat'}); 
 printf(" -olpdb: output ligand property filename in PDB [%s]\n",$OPT{'olpdb'}); 
 printf(" -idab : input dir of PDB with DABRL annotation [%s]\n",$OPT{'idab'}); 
 printf(" -ohead : outputfile header [%s]\n",$OPT{'ohead'}); 
 exit(1);
}

$ilistfile = $ARGV[0];
&Read_Options(\@ARGV,\%OPT);
&Read_List_File($ilistfile,\@LIST);
$ligpat = $OPT{'lpat'};
$ofname = $OPT{'op'};

$ofname = $OPT{'ol'};
$Nligand = 0;

foreach $p (@LIST){
 $olgfile = $OPT{'olg'}.'/'.$p;
 print ">$p $olgfile\n";
 
 ## (1) Read Pocket File ##
 &Read_Pocket_Annotated_Ligand_File($olgfile,\@liglist,\%ligdat);
 %Nlig = ();
 foreach $lig (@liglist){
  if ($lig->{'Resi'} =~/$ligpat/){ 
   printf("#%s\n",$lig->{'Resi'});
   $Nlig{$lig->{'Resi'}} += 1;  
  }
 }

 ## (2) Read DABRL File ##
 @trilist = keys(%Nlig);
 foreach $tri (@trilist){
   $tri =~s/\s+//g;
   printf("'%s'\n",$tri);
   if ($tri!~/[a-z]/){
     &Read_DABRL_Annotated_Ligand_File($OPT{'idab'}.'/'.$tri,\%dabdat);
   }
  }
 
 ## (3) Counting  ##
 foreach $lig (@liglist){
  printf(">%s\n",$lig->{'Resi'});
  if ($lig->{'Resi'} !~/[a-z]/){
  &Count_DABRL_based_Ligand_Rinaccess(\%Count,$lig,\%dabdat);
  ++$Nligand;
  } 
 }

} ## $p ##

$NALL = 0;
@ind = sort {$a <=> $b} keys(%{$Count{'ALL'}});
foreach $i (@ind){ $NALL += $Count{'ALL'}->{$i}; }
foreach $i (@ind) {$Freq_all{$i} = $Count{'ALL'}->{$i}/$NALL;}

print "#Nligand $Nligand\n";
&Output_Count($OPT{'ohead'}.'all.out',\%{$Count{'ALL'}},\%Freq_all);
&Output_Count($OPT{'ohead'}.'D.out',\%{$Count{'D'}},\%Freq_all);
&Output_Count($OPT{'ohead'}.'A.out',\%{$Count{'A'}},\%Freq_all);
&Output_Count($OPT{'ohead'}.'B.out',\%{$Count{'B'}},\%Freq_all);
&Output_Count($OPT{'ohead'}.'R.out',\%{$Count{'R'}},\%Freq_all);
&Output_Count($OPT{'ohead'}.'L.out',\%{$Count{'L'}},\%Freq_all);
exit(1);



#################
### FUNCTIONS ###
#################

sub Output_Count{
 my($ofname,$count,$freq_all) = @_;
 my($OF,$Ntotal);
 printf("#Output_Count()-->'%s'\n",$ofname);
 open(OF,">$ofname");
 @ind = sort {$a <=> $b} keys(%{$count});
 $Ntotal = 0;
 foreach $i (@ind){
   $Ntotal += $count->{$i};
 } 
 printf(OF "#ofname '%s'\n",$ofname);
 printf(OF "#COMMAND  %s\n",$OPT{'COMMAND'});
 printf(OF "#date     %s\n",&Get_Date_String());
 printf(OF "#Ntotal %d\n",$Ntotal);
 printf(OF "#[value] [freq] [freq/freq_all] [count] [freq_all]\n");
 foreach $i (@ind){
   printf(OF "%s %f %f %d %f\n",$i,$count->{$i}/$Ntotal,$count->{$i}/$Ntotal/$freq_all->{$i},$count->{$i},$freq_all->{$i});
 } 
 close(OF);

} ## end of Output_Count() ##



sub Count_DABRL_based_Ligand_Rinaccess {
 my($count,$lig,$dabdat) = @_;

 
 printf("#Count_DABRL_based_Ligand_Rinaccess(Natom %d)\n",scalar(@{$lig->{'ATOM'}}));
 my($a); my($iRinacc,$Rinacc, $bin_width);
 $bin_width = 1.0;

 for ($i=0;$i<scalar(@{$lig->{'ATOM'}});++$i){
  $Resi = $lig->{'ATOM'}->[$i]->{'Resi'};
  $Atom = $lig->{'ATOM'}->[$i]->{'Atom'};
  $d = $dabdat->{$Resi}->{$Atom};
  if ($d ne ''){
  printf("##%s %s DAB %s\n",$lig->{'ATOM'}->[$i]->{'Resi'},$lig->{'ATOM'}->[$i]->{'Atom'},$d);
  $Rinacc = $lig->{'ATOM'}->[$i]->{'Rinacc'};
  $indRinacc = $bin_width * int($Rinacc/$bin_width);
  print "dabrl $d Rinacc $Rinacc $indRinacc\n"; 
  $count->{$d}->{$indRinacc} += 1;
  $count->{'ALL'}->{$indRinacc} += 1;
  }
 }
} ## end of Count_DABRL_based_Ligand_Rinaccess() ##





sub Read_Pocket_Annotated_Ligand_File{
 my($fname,$lig) = @_;

#>> FILE FORMAT EXAMPLE <<
##REMARK  [Ratom](10:55-60):Radius of atom[A]
##REMARK  [tFact](11:61-66):[Riacc] is assigned as tFactor
##REMARK  [shAcc](12:67-72):Shell accessibility. (Ratio of noVdW grids in the shell)[%]
##REMARK  [Riacc](13:73-78):averaged Rinaccess in the 'shell'[A]
##REMARK  [rwAcc](14:79-84):Rinacc-weighted shell accessibility [%]
##REMARK  [pocke](15:85-90):pocketness.ratio of inv Rpocket sum [%]
##REMARK  [ClusterNum](91- ):Cluster Num of Pocket around the atom
##REMARK                                                 Ratom|Riacc|shAcc|Riacc|rwAcc|pocke|ClusterNumbers
##0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
##           1         2         3         4         5         6         7        8         9 
##ATOM   2203  N   ARG B  21      17.709 -28.355  47.219 1.65  10.39  97.8 10.39  98.9  2.88
##ATOM   2204  CA  ARG B  21      16.314 -28.873  47.317 1.87  10.11  93.6 10.11  96.3  3.80
##ATOM   2205  C   ARG B  21      15.348 -27.841  46.739 1.76  10.09  91.5 10.09  96.1  4.90
##ATOM   2206  O   ARG B  21      14.478 -28.171  45.931 1.40   8.20  85.7  8.20  78.1 16.35
##ATOM   2207  CB  ARG B  21      15.961 -29.145  48.782 1.87  10.50  93.5 10.50 100.0  0.00
##ATOM   2208  CG  ARG B  21      17.078 -29.806  49.584 1.87  10.50  94.7 10.50 100.0  0.00
##ATOM   2209  CD  ARG B  21      17.214 -31.294  49.302 1.87  10.30  90.4 10.30  98.1  0.00
##TER
##HETATM 4465 FE   HEM A 401      14.990 -25.959  31.002 1.47   3.43  76.4  3.43  32.7 57.14
##HETATM 4466  CHA HEM A 401      15.090 -22.466  31.321 2.00   3.06  69.4  3.06  29.1 62.18
##HETATM 4467  CHB HEM A 401      11.542 -25.858  31.557 2.00   3.45  74.2  3.45  32.8 57.14
##HETATM 4468  CHC HEM A 401      14.792 -29.516  31.321 2.00   3.09  69.5  3.09  29.5 57.14
##HETATM 4469  CHD HEM A 401      18.321 -26.118  31.119 2.00   2.93  71.1  2.93  27.9 66.52
##HETATM 4470  NA  HEM A 401      13.548 -24.409  31.424 1.65   3.47  76.7  3.47  33.1 57.14
##HETATM 4471  C1A HEM A 401      13.835 -23.082  31.411 1.78   3.33  72.7  3.33  31.7 58.23
##HETATM 5488 HM71 FMN A1401      28.425  11.032  20.682 1.50   2.48  55.9  2.48  23.6 76.54
##HETATM 5489 HM72 FMN A1401      27.577  10.582  19.442 1.50   2.39  56.7  2.39  22.8 70.06
##HETATM 5490 HM73 FMN A1401      28.404   9.524  20.251 1.50   2.31  46.5  2.31  22.0 80.00
##HETATM 5491 HM81 FMN A1401      28.179  12.008  22.200 1.50   2.50  59.8  2.50  23.8 66.94
##HETATM 5492 HM82 FMN A1401      27.895  11.654  23.700 1.50   2.46  58.0  2.46  23.4 50.55
##0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
##           1         2         3         4         5         6         7        8         9 
#
#


my($str); my(@cnum); my($c);
 my($F); my($index); my(%index_hash) = ();
 my($Anum,$Atom, $Resi,$Rnum,$Radius,$iRinac,$Rinacc,$shellACC,$X,$Y,$Z);
 my($N) = 0; my($M) = 0;
 open(F,$fname) || die "#ERROR:Can't open pocket_annotated_ligand_pdbfile '$fname'";

 @{$lig} = ();

 while (<F>){
  if ((/^HETATM/)||(/^ATOM/)){
   chomp;
   $Anum     = substr($_,6,5);
   $Atom     = substr($_,12,4);
   $Resi     = substr($_,17,3);
   $Rnum     = substr($_,20,6);
   $X        = substr($_,30,8);
   $Y        = substr($_,38,8);
   $Z        = substr($_,46,8);
   $Radius   = substr($_,54,6);
   $tFactor  = substr($_,60,6);
   $shellACC = substr($_,66,6);
   $Rinacc   = substr($_,72,6);
   $rwACC    = substr($_,78,6);
   $pocket   = substr($_,84,6);
   $str = substr($_,79);
   if ($Atom !~/^H/){
    $lig->[$M]->{'Resi'} = $Resi;
    $lig->[$M]->{'Type'} = $Resi;
    $lig->[$M]->{'ATOM'}->[$N]->{'Rnum'} = $Rnum;
    $lig->[$M]->{'ATOM'}->[$N]->{'Resi'} = $Resi;
    $lig->[$M]->{'ATOM'}->[$N]->{'Atom'} = $Atom;
    $lig->[$M]->{'ATOM'}->[$N]->{'Anum'} = $Anum;
    $lig->[$M]->{'ATOM'}->[$N]->{'Rnum'} = $Rnum;
    $lig->[$M]->{'ATOM'}->[$N]->{'Resi'} = $Resi;
    $lig->[$M]->{'ATOM'}->[$N]->{'X'} = $X;
    $lig->[$M]->{'ATOM'}->[$N]->{'Y'} = $Y;
    $lig->[$M]->{'ATOM'}->[$N]->{'Z'} = $Z;
    $lig->[$M]->{'ATOM'}->[$N]->{'Rinacc'}   = $Rinacc;
    $lig->[$M]->{'ATOM'}->[$N]->{'shellACC'} = $shellACC;
    $lig->[$M]->{'ATOM'}->[$N]->{'rwACC'}    = $rwACC;
    $lig->[$M]->{'ATOM'}->[$N]->{'pocket'}    = $pocket;
   #print "#$fname Rnum $Rnum Resi $Resi Atom $Atom Anum $Anum shellACC $shellACC Rinacc $Rinacc\n";
   if ($Rinacc==0.0){
    print "#ERROR: $fname Rnum $Rnum Resi $Resi Atom $Atom Anum $Anum shellACC $shellACC Rinacc $Rinacc\n";
    }
    ++$N;
   }
  }
  if (/^TER/) {$N = 0; ++$M;}
 } # while #
 close(F);
 my($Nmax);

} ## end of Read_Pocket_Annotated_Ligand_File() ##




sub Read_DABRL_Annotated_Ligand_File{
 my($fname,$dabdat) = @_;

#HEADER    PLP
#REMARK    Output by write_in_pdb() in mmCIF.py
#REMARK    NAME    '"PYRIDOXAL-5'-PHOSPHATE"'
#REMARK    TYPE    'NON-POLYMER'
#REMARK    FORMULA 'C8 H10 N O6 P'
#REMARK    OCCUP[55-60]: Nneighbor
#REMARK    tFact[61-66]: DABRLnumber(0:H,1:D,2:A,3:B,4:R,5:L,6:X)
#REMARK            [68]: 'D'onor,'A'ccptor,'B'oth, a'R'omatic, a'L'iphatic
#          1         2         3         4         5         6         7        8         9 
#0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
#ATOM      1  N1  PLP     1     -14.333  14.052  -8.386  2.00  2.00 A         N
#ATOM      2  C2  PLP     1     -13.067  14.164  -7.932  3.00  4.00 R         C
#ATOM      3  C2A PLP     1     -12.523  15.554  -7.674  4.00  5.00 L         C
#ATOM      4  C3  PLP     1     -12.313  12.990  -7.718  3.00  5.00 L         C
#ATOM      5  O3  PLP     1     -11.053  13.159  -7.240  2.00  3.00 B         O
#ATOM      6  C4  PLP     1     -12.869  11.745  -7.941  3.00  4.00 R         C
#ATOM      7  C4A PLP     1     -12.000  10.575  -8.113  3.00  5.00 L         C
#ATOM      8  O4A PLP     1     -10.848  10.748  -8.789  1.00  2.00 A         O
#ATOM      9  C5  PLP     1     -14.205  11.723  -8.428  3.00  4.00 R         C
#ATOM     10  C6  PLP     1     -14.909  12.852  -8.636  3.00  4.00 R         C
#ATOM     11  C5A PLP     1     -14.854  10.340  -8.684  4.00  5.00 L         C
#ATOM     12  O4P PLP     1     -14.955   9.680  -7.404  2.00  2.00 A         O
#ATOM     13  P   PLP     1     -15.917   8.413  -7.177  4.00  6.00 X         P
#ATOM     14  O1P PLP     1     -15.467   7.408  -8.239  1.00  2.00 A         O
#ATOM     15  O2P PLP     1     -15.685   7.981  -5.762  2.00  3.00 B         O
#ATOM     16  O3P PLP     1     -17.279   8.894  -7.418  2.00  3.00 B         O
#ATOM     17  H2A1PLP     1     -11.477  15.646  -7.298  1.00  0.00 H         H
#ATOM     18  H2A2PLP     1     -12.634  16.167  -8.598  1.00  0.00 H         H
#ATOM     19  H2A3PLP     1     -13.209  16.091  -6.978  1.00  0.00 H         H
#ATOM     20  HO3 PLP     1     -10.545  12.368  -7.096  1.00  0.00 H         H
#ATOM     21  H4A PLP     1     -12.213   9.561  -7.732  1.00  0.00 H         H
#ATOM     22  H6  PLP     1     -15.946  12.794  -9.007  1.00  0.00 H         H
#ATOM     23  H5A1PLP     1     -15.828  10.405  -9.221  1.00  0.00 H         H
#ATOM     24  H5A2PLP     1     -14.311   9.735  -9.447  1.00  0.00 H         H
#ATOM     25  HOP2PLP     1     -16.253   7.232  -5.627  1.00  0.00 H         H
#ATOM     26  HOP3PLP     1     -17.847   8.145  -7.283  1.00  0.00 H         H
#END
##0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
##           1         2         3         4         5         6         7        8         9 
#
#


 my($str); my(@cnum); my($c);
 my($F); my($index); my(%index_hash) = ();
 my($Anum,$Atom, $Resi,$Rnum,$Radius,$X,$Y,$Z);
 my($N) = 0; my($M) = 0;
 if (-e $fname){
 open(F,$fname) || die "#ERROR:Can't open pocket_annotated_ligand_pdbfile '$fname'";

 @{$lig} = ();

 while (<F>){
  if ((/^HETATM/)||(/^ATOM/)){
   chomp;
   $Anum     = substr($_,6,5);
   $Atom     = substr($_,12,4);
   $Resi     = substr($_,17,3);
   $Rnum     = substr($_,20,6);
   $X        = substr($_,30,8);
   $Y        = substr($_,38,8);
   $Z        = substr($_,46,8);
   $Radius   = substr($_,54,6);
   $tFactor  = substr($_,60,6);
   $DABRL    = substr($_,67,1);
   $str = substr($_,79);
   if ($Atom !~/^H/){
    $dabdat->{$Resi}->{$Atom}=$DABRL;
    print "#$fname Rnum $Rnum Resi $Resi Atom $Atom Anum $Anum DABRL $DABRL\n";
   }
  }
 } # while #
 close(F);
 return(1); 
 }
 else {return(0);}

} ## end of Read_DABRL_Annotated_Ligand_File() ##








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
 while ($i<scalar(@ARGV)){
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

