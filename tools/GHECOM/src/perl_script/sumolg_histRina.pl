#!/usr/bin/perl
#
# <sumolg_histRina.pl>
#

$LastModDate = "Dec 2, 2008";

$OPT{'ipdb'} = 'RotPDB';
$OPT{'opp'}  = 'GHECOM';
$OPT{'ilg'}  = 'RotPDB';
$OPT{'olg'}  = 'OLG';
$OPT{'op'}   = '-';
$OPT{'lpat'}  = '.';
$OPT{'olpdb'}  = 'olig.pdb';
$OPT{'olh'}    = 'hist.dat';
$OPT{'S'} = 'A';

if (scalar(@ARGV)<1){
 print "sumolg_histRina.pl [listfile] <options>\n";
 print " for making histogram Rinacc value for each ligand atom type\n";
 printf(" coded by T.Kawabata. LastModified:%s\n",$LastModDate);
 printf(" -olg  : dir for olg  [%s]\n",$OPT{'olg'}); 
 printf(" -lpat : three-letter ligand pattern  [%s]\n",$OPT{'lpat'}); 
 printf(" -olpdb: output ligand property filename in PDB [%s]\n",$OPT{'olpdb'}); 
 printf(" -olh  : output Rinacc cumulative histogram file[%s]\n",$OPT{'olh'}); 
 printf(" -S    : 'A'tom-based, 'G'roup-based statistics[%s]\n",$OPT{'S'}); 
 exit(1);
}

$ilistfile = $ARGV[0];
&Read_Options(\@ARGV,\%OPT);
&Read_List_File($ilistfile,\@LIST);
$ligpat = $OPT{'lpat'};
$ofname = $OPT{'op'};

$ofname = $OPT{'ol'};

foreach $p (@LIST){
 $olgfile = $OPT{'olg'}.'/'.$p;
 &Read_Pocket_Annotated_Ligand_File($olgfile,\@ligdat);
 #&Set_Type_Of_Ligand_Molecule_Type(\@ligdat);
 printf("$p $olgfile Nlig %d\n",scalar(@ligdat));
 $Nlig = scalar(@ligdat);
 for ($m=0;$m<$Nlig;++$m){
  printf("%s->%s\n",$ligdat[$m]->{'Resi'},$ligdat[$m]->{'Type'});
  if ($ligdat[$m]->{'Type'} =~/$ligpat/){ 
   if ($OPT{'S'} eq 'A'){
     &Count_Atom_Based_Ligand_Rinacc_Histogram(\%CountRinacc,\@ligdat,$m);
   }
   if ($OPT{'S'} eq 'G'){
     $Rinacc = &Count_Group_Based_Ligand_Rinacc_Histogram(\%CountRinacc,\@ligdat,$m);
     printf("#%s %s %s\n",$p,$ligdat[$m]->{'Type'},$Rinacc);
   }
  }
 } ## $m ##


} ## $p ##

&Write_Atom_Based_Ligand_Rinacc_Histogram($OPT{'olh'},\%CountRinacc,$OPT{'lpat'});

#################
### FUNCTIONS ###
#################

sub Count_Atom_Based_Ligand_Rinacc_Histogram{
 my($count,$lig,$Nmol) = @_;

 my($a); my($iRinacc); my($Rinacc); my($ind_Rinacc);
 my($Natom) = scalar(@{$lig->[$Nmol]->{'ATOM'}});

 for ($a=0;$a<$Natom;++$a){
  $Rinacc =  $lig->[$Nmol]->{'ATOM'}->[$a]->{'Rinacc'};
  $shellACC =  $lig->[$Nmol]->{'ATOM'}->[$a]->{'shellACC'};
  if (($shellACC>0.0)&&($shellACC<100.0)){
   $ind_Rinacc = int($Rinacc);  
   #printf("Rinacc $Rinacc ind $ind_Rinacc\n"); 
   $count->{$ind_Rinacc} += 1;
   if ($Rinacc==0){
    print "#ERROR:a $a shellACC $shellACC Rinacc $Rinacc\n";
    exit(1); 
   }
  } 
 }
} ## end of Count_Atom_Based_Ligand_Rinacc_Histogram() ##



sub Count_Group_Based_Ligand_Rinacc_Histogram{
 my($count,$lig,$Nmol) = @_;

 my($a,$iRinacc,$Rinacc,$ind_Rinacc,$sum,$n);
 my($Natom) = scalar(@{$lig->[$Nmol]->{'ATOM'}});

 $sum = 0.0; $n = 0;
 for ($a=0;$a<$Natom;++$a){
  $Rinacc =  $lig->[$Nmol]->{'ATOM'}->[$a]->{'Rinacc'};
  $shellACC =  $lig->[$Nmol]->{'ATOM'}->[$a]->{'shellACC'};
  if (($shellACC>0.0)&&($shellACC<100.0)){
   $ind_Rinacc = int($Rinacc);  
   $sum += $ind_Rinacc;  
   $n += 1;
   #printf("Rinacc $Rinacc ind $ind_Rinacc\n"); 
   if ($Rinacc==0){
    print "#ERROR:a $a shellACC $shellACC Rinacc $Rinacc\n";
    exit(1); 
   }
  } 
 }

 $Rinacc = $sum/$n;
 $ind_Rinacc = int($Rinacc);
 $count->{$ind_Rinacc} += 1;
 return($Rinacc);
 
} ## end of Count_Group_Based_Ligand_Rinacc_Histogram() ##


sub Write_Atom_Based_Ligand_Rinacc_Histogram{
 my($ofname,$count,$lpatstr) = @_;
 my($OF,$r); my(@rlist) = sort {$a <=> $b} keys(%{$count});
 my(@rindlist) = ('-2','3','4','5','6','7-9','10-');
 my(%Cnt); my($Nall) = 0;
 $Cnt{'-2'} = $count->{'2'} + $count{'1'} + $count{'0'}; 
 $Cnt{'3'} = $count->{'3'}; 
 $Cnt{'4'} = $count->{'4'}; 
 $Cnt{'5'} = $count->{'5'}; 
 $Cnt{'6'} = $count->{'6'}; 
 $Cnt{'7-9'} = $count->{'7'}+$count->{'8'}+$count->{'9'}; 
 $Cnt{'10-'} = $count->{'10'}+$count->{'11'}+$count->{'12'} + $count->{'13'}+$count->{'14'}+$count->{'15'}; 
 foreach $r (@rindlist) {$Nall += $Cnt{$r};}
 
 printf("#Write_Atom_Based_Ligand_Rinacc_Histogram() -- >> '%s'\n",$ofname);
 open(OF,">>$ofname") || die "#ERROR:Can't open '$ofname'";
 printf(OF "#COMMAND%s\n",$OPT{'COMMAND'});
 printf(OF "#%s\t",$lpatstr);  
 foreach $r (@rindlist){ printf(OF "%s\t",$r); } printf(OF "\n"); 
 printf(OF "#%s\t",$lpatstr);  
 foreach $r (@rindlist){ printf(OF "%d\t",$Cnt{$r}); } printf(OF "\n"); 
 printf(OF "%s\t",$lpatstr);  
 if ($Nall>0){
 foreach $r (@rindlist){ printf(OF "%f\t",$Cnt{$r}/$Nall); } printf(OF "\n"); 
 }
 else 
 {foreach $r (@rindlist){ printf(OF "%d\t",$Cnt{$r}); } printf(OF "\n"); }
 close(OF);

} ## end of Write_Atom_Based_Ligand_Rinacc_Histogram() ##



sub Read_Pocket_Annotated_Ligand_File{
 my($fname,$lig) = @_;

#>> FILE FORMAT EXAMPLE <<
#REMARK  [Ratom](10:55-60):Radius of atom[A]
#REMARK  [tFact](11:61-66):[Riacc] is assigned as tFactor
#REMARK  [shAcc](12:67-72):Shell accessibility. (Ratio of noVdW grids in the shell)[%]
#REMARK  [Riacc](13:73-78):averaged Rinaccess in the 'shell'[A]
#REMARK  [rwAcc](14:79-84):Rinacc-weighted shell accessibility [%]
#REMARK  [pocke](15:85-90):pocketness.ratio of inv Rpocket sum [%]
#REMARK  [ClusterNum](91- ):Cluster Num of Pocket around the atom
#REMARK                                                 Ratom|Riacc|shAcc|Riacc|rwAcc|pocke|ClusterNumbers
#0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
#           1         2         3         4         5         6         7        8         9 
#ATOM   2203  N   ARG B  21      17.709 -28.355  47.219 1.65  10.39  97.8 10.39  98.9  2.88
#ATOM   2204  CA  ARG B  21      16.314 -28.873  47.317 1.87  10.11  93.6 10.11  96.3  3.80
#ATOM   2205  C   ARG B  21      15.348 -27.841  46.739 1.76  10.09  91.5 10.09  96.1  4.90
#ATOM   2206  O   ARG B  21      14.478 -28.171  45.931 1.40   8.20  85.7  8.20  78.1 16.35
#ATOM   2207  CB  ARG B  21      15.961 -29.145  48.782 1.87  10.50  93.5 10.50 100.0  0.00
#ATOM   2208  CG  ARG B  21      17.078 -29.806  49.584 1.87  10.50  94.7 10.50 100.0  0.00
#ATOM   2209  CD  ARG B  21      17.214 -31.294  49.302 1.87  10.30  90.4 10.30  98.1  0.00
#TER
#HETATM 4465 FE   HEM A 401      14.990 -25.959  31.002 1.47   3.43  76.4  3.43  32.7 57.14
#HETATM 4466  CHA HEM A 401      15.090 -22.466  31.321 2.00   3.06  69.4  3.06  29.1 62.18
#HETATM 4467  CHB HEM A 401      11.542 -25.858  31.557 2.00   3.45  74.2  3.45  32.8 57.14
#HETATM 4468  CHC HEM A 401      14.792 -29.516  31.321 2.00   3.09  69.5  3.09  29.5 57.14
#HETATM 4469  CHD HEM A 401      18.321 -26.118  31.119 2.00   2.93  71.1  2.93  27.9 66.52
#HETATM 4470  NA  HEM A 401      13.548 -24.409  31.424 1.65   3.47  76.7  3.47  33.1 57.14
#HETATM 4471  C1A HEM A 401      13.835 -23.082  31.411 1.78   3.33  72.7  3.33  31.7 58.23


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
   $lig->[$M]->{'ATOM'}->[$N]->{'Rinacc'}   = $Rinacc; 
 #  print "$M $N Rinacc $Rinacc\n";
   $lig->[$M]->{'ATOM'}->[$N]->{'shellACC'} = $shellACC; 
   $lig->[$M]->{'ATOM'}->[$N]->{'rwACC'}    = $rwACC; 
   $lig->[$M]->{'ATOM'}->[$N]->{'pocket'}    = $pocket; 
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


sub Set_Type_Of_Ligand_Molecule_Type{
 my($LIG) = @_;
 my($Nlig) = scalar(@{$LIG}); my($Natom,$m,$a,$atom,%Cnt,$resi);
 for ($m=0;$m<$Nlig;++$m){
  $Natom = scalar(@{$LIG->[$m]->{'ATOM'}}); 
  $resi = $LIG->[$m]->{'Resi'};
  %Cnt = ();
  for ($a=0;$a<$Natom;++$a){  
   $atom = $LIG->[$m]->{'ATOM'}->[$a]->{'Atom'};
  # print "atom '$atom'\n";
   if ($atom eq ' CA ') {$Cnt{' CA '} += 1;}
   if ($atom eq ' N  ') {$Cnt{' N  '} += 1;}
   if ($atom eq ' C  ') {$Cnt{' C  '} += 1;}
   if ($atom eq ' O  ') {$Cnt{' O  '} += 1;}
   if ($atom eq ' P  ') {$Cnt{' P  '} += 1;}
   if ($atom eq ' C1*') {$Cnt{' C1*'} += 1;}
   if ($atom eq ' N1 ') {$Cnt{' N1 '} += 1;}
   if ($atom eq ' O2*') {$Cnt{' O2*'} += 1;}
  } ## $a ##

#   printf("m %d Natom $Natom CA %d P %d\n",$m,$Cnt{' CA '},$Cnt{' P  '});

   if (($Cnt{' CA '}>0)&&($Cnt{' N  '}>0)&&($Cnt{' C  '}>0)&&($Cnt{' O  '}>0)){
     if ($Cnt{' CA '}>10) {$LIG->[$m]->{'Type'} = 'pro';}
     elsif (($resi ne 'SAM')&&($resi ne 'SAH')) {$LIG->[$m]->{'Type'} = 'pep';}
   }

   if (($Cnt{' P  '}>0)&&($Cnt{' C1*'}>0)&&($Cnt{' N1 '}>0)){
        if ($Cnt{' P  '}>=3) 
    {
       if ($Cnt{' O2*'}>0) {$LIG->[$m]->{'Type'} = 'rna';}
     else  {$LIG->[$m]->{'Type'} = 'dna';}
    }   
 }



 } ## $m ##

} ## end of Set_Type_Of_Ligand_Molecule_Type() ##








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
                                                                                                      
 while (<F>){
  if ($_!~/^#/){
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

