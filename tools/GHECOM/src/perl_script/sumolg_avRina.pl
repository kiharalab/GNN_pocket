#!/usr/bin/perl
#
# <sumolg_avRina.pl>
#

$LastModDate = "Feb 12, 2009";

$OPT{'olg'}  = 'OLG';
$OPT{'op'}   = '-';
$OPT{'lpat'}  = '.';
$OPT{'olpdb'}  = 'olig.pdb';

if (scalar(@ARGV)<1){
 print "sumolg_avRina.pl [listfile] <options>\n";
 print " for making averaged Rinacc value for each ligand atom type\n";
 printf(" coded by T.Kawabata. LastModified:%s\n",$LastModDate);
 printf(" -olg  : dir for olg  [%s]\n",$OPT{'olg'}); 
 printf(" -lpat : three-letter ligand pattern  [%s]\n",$OPT{'lpat'}); 
 printf(" -olpdb: output ligand property filename in PDB [%s]\n",$OPT{'olpdb'}); 
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
 print ">$p $olgfile\n";
 &Read_Pocket_Annotated_Ligand_File($olgfile,\@liglist,\%ligdat);

 foreach $lig (@liglist){
  if ($lig->{'Resi'} =~/$ligpat/){ 
   printf("$p Resi '%s' ligpat '%s'\n",$lig->{'Resi'},$ligpat); 
   &Count_Atom_Based_Ligand_Pocketness(\%CountLigAtom,$lig,$p);
   $CountLigAtom{'last_protein'} = $p; 
   $CountLigAtom{'last_index'} = $lig->{'Resi'}; 
   push(@LIGlist,"$p $ind"); 
  } 
 }
} ## $p ##
&Write_Atom_based_Average_invRinacc_Data($OPT{'olpdb'},\%CountLigAtom,\@LIGlist);


#################
### FUNCTIONS ###
#################

sub Write_Atom_based_Average_invRinacc_Data{
 my($ofname,$count,$liglist) = @_; 
 my($OF); my(@atmlist); my($a,$lig,$Nlig,$maxT,$minT,$init);

 printf("#Write_Atom_based_Average_invRinacc_Data()-->'%s'\n",$ofname);
 open(OF,">$ofname")||die "#ERROR:Can't write to '$ofname'";
 printf(OF "HEADER Averaged iRinacc,Rinacc,shellACC for the focused ligand\n");
 printf(OF "REMARK COMMAND %s\n",$OPT{'COMMAND'}); 
 printf(OF "REMARK DATE    %s\n",&Get_Date_String()); 
 $Nlig = scalar(@{$liglist});
 printf(OF "REMARK NLIGAND %d\n",$Nlig);
 foreach $lig (@{$liglist}){
  printf(OF "REMARK LIGAND %s\n",$lig);
 }
 @atmlist = keys(%{$count}); 
 printf(OF "REMARK  [1   ]:Ndata\n");
 printf(OF "REMARK  [2   ]:Average of Rinacc\n");
 printf(OF "REMARK  [3   ]:Average of Rinacc\n");
 printf(OF "REMARK  [4   ]:Average of shellACC\n");
 printf(OF "REMARK  [5   ]:Min of Rinacc\n");
 printf(OF "REMARK  [6   ]:Max of Rinacc\n");
 printf(OF "REMARK  [7   ]:SD of Rinacc\n");
 printf(OF "REMARK  [8   ]:Last protein, Resi and Rnum\n");
 printf(OF "REMARK                                                [1   ][2   ][3   ][4   ][5   ][6   ][7   ]\n");

 $init = 1;
 
 foreach $a (@atmlist){
if ($count->{$a}->{'N'}>0){ 
 
# if (($count->{'last_protein'} eq $count->{$a}->{'last_protein'})&&
#     ($count->{'last_index'}   eq $count->{$a}->{'last_index'}) )  {printf(OF "HETATM");}
#    else  {printf(OF "#ETATM");}
    printf(OF "HETATM");


 $aveRinacc = $count->{$a}->{'Rinacc'}/$count->{$a}->{'N'};
 $SDRinacc = sqrt($count->{$a}->{'RRinacc'}/$count->{$a}->{'N'} - $aveRinacc * $aveRinacc);
 printf(OF "%5s %4s %3s %6s   %8.3f%8.3f%8.3f%6d%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f %s %s\n",
 $count->{$a}->{'Anum'}, $count->{$a}->{'Atom'}, $count->{$a}->{'Resi'}, $count->{$a}->{'Rnum'},
 $count->{$a}->{'X'}, $count->{$a}->{'Y'}, $count->{$a}->{'Z'},
 $count->{$a}->{'N'},
 # ($tfact-$minT)/($maxT-$minT), 
 $count->{$a}->{'Rinacc'}/$count->{$a}->{'N'},
 $count->{$a}->{'Rinacc'}/$count->{$a}->{'N'},
 $count->{$a}->{'shellACC'}/$count->{$a}->{'N'},
 $count->{$a}->{'Rinacc_min'},
 $count->{$a}->{'Rinacc_max'},$SDRinacc,
 $count->{$a}->{'last_protein'},
 $count->{$a}->{'last_index'}
  );
 }

} ## $a ##
 close(OF);
} ## end of Write_Atom_based_Average_invRinacc_Data() ##



sub Count_Atom_Based_Ligand_Pocketness {
 my($count,$lig,$protein) = @_;

 printf("#Count_Atom_Based_Ligand_Pocketness(Natom %d)\n",scalar(@{$lig->{'ATOM'}}));
 my($a); my($iRinacc);

 for ($i=0;$i<scalar(@{$lig->{'ATOM'}});++$i){
  $a = $lig->{'ATOM'}->[$i]->{'Atom'};
  print "atom '$a'\n"; 
  $count->{$a}->{'N'} += 1;
  $count->{$a}->{'shellACC'} += $lig->{'ATOM'}->[$i]->{'shellACC'};
  $iRinacc =  $lig->{'ATOM'}->[$i]->{'iRinacc'};
  $Rinacc  =  $lig->{'ATOM'}->[$i]->{'Rinacc'};
  $count->{'ATOM'}->{$a}->{'iRinacc'}  += $iRinacc;

  if (($count->{$a}->{'iRinacc_max'} eq '')||($iRinacc>$count->{$a}->{'iRinacc_max'}))
   { $count->{$a}->{'iRinacc_max'}  = $iRinacc;}
  if (($count->{$a}->{'iRinacc_min'} eq '')||($iRinacc<$count->{$a}->{'iRinacc_min'}))
   { $count->{$a}->{'iRinacc_min'}  = $iRinacc;}

  if (($count->{$a}->{'Rinacc_max'} eq '')||($Rinacc>$count->{$a}->{'Rinacc_max'}))
   { $count->{$a}->{'Rinacc_max'}  = $Rinacc;}
  if (($count->{$a}->{'Rinacc_min'} eq '')||($Rinacc<$count->{$a}->{'Rinacc_min'}))
   { $count->{$a}->{'Rinacc_min'}  = $Rinacc;}
  
  $count->{$a}->{'Rinacc'} += $Rinacc;
  $count->{$a}->{'RRinacc'} += $Rinacc * $Rinacc;
  $count->{$a}->{'Atom'}  = $lig->{'ATOM'}->[$i]->{'Atom'};
  $count->{$a}->{'Anum'}  = $lig->{'ATOM'}->[$i]->{'Anum'};
  $count->{$a}->{'Rnum'}  = $lig->{'ATOM'}->[$i]->{'Rnum'};
  $count->{$a}->{'Resi'}  = $lig->{'ATOM'}->[$i]->{'Resi'};
  $count->{$a}->{'X'}    = $lig->{'ATOM'}->[$i]->{'X'};
  $count->{$a}->{'Y'}    = $lig->{'ATOM'}->[$i]->{'Y'};
  $count->{$a}->{'Z'}    = $lig->{'ATOM'}->[$i]->{'Z'};
  $count->{$a}->{'last_index'}    = $index;
  $count->{$a}->{'last_protein'}  = $protein;
 }
} ## end of Count_Atom_Based_Ligand_Pocketness() ##





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

