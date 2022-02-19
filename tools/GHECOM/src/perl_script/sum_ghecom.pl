#!/usr/bin/perl
#
# <sum_ghecom.pl>
#

$LastModDate = "Dec 21,2006";

$OPT{'ipdb'} = 'RotPDB';
$OPT{'opp'}  = 'GHECOM';
$OPT{'ilg'}  = 'RotPDB';
$OPT{'olg'}  = 'OLG';
$OPT{'op'}   = '-';
$OPT{'ol'}   = '-';
$OPT{'A'} = 'F';

if (scalar(@ARGV)<1)
{
 print "sum_ghecom.pl [listfile] <options>\n";
 print " for making summary of ghecom detected pockets\n";
 printf(" coded by T.Kawabata. LastModified:%s\n",$LastModDate);
 printf(" -opp  : dir for oppp [%s]\n",$OPT{'opp'}); 
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
 $oppfile = $OPT{'opp'}.'/'.$p;
 $olgfile = $OPT{'olg'}.'/'.$p;
 &Read_Pocket_Grid_PDB_File($oppfile,\@cluster);
 &Read_Pocket_Annotated_Ligand_File($olgfile,\@liglist,\%ligdat);
 $Ncluster = scalar(@cluster);
 $Nligand  = scalar(@liglist);
 print OL "#$oppfile Ncluster $Ncluster Nligand $Nligand\n";
 foreach $ind (@liglist)
 {
  $str = $ind;
  $str =~s/\s/\_/g;
  printf(OL "%s %s Natm %2d Natm_po %2d Nclus %2d",
   $p,$str,$ligdat{$ind}->{'Natom'},$ligdat{$ind}->{'Natom_pock'},
   $ligdat{$ind}->{'Ncluster'});

 if ($ligdat{$ind}->{'Ncluster'}>0)
 {
 $max_c = $ligdat{$ind}->{'clus_num_max'};
  printf(OL " c_max %2d %2d %5.1f",
 $max_c,$ligdat{$ind}->{'clus_num'}->{$max_c},
 100.0*$ligdat{$ind}->{'clus_num'}->{$max_c}/$ligdat{$ind}->{'Natom_pock'});

  $cluster[$max_c-1]->{'ligand_str'} .= ' ';
  $cluster[$max_c-1]->{'ligand_str'} .= $ligdat{$ind}->{'Resi'};

#  printf("%s\n",$ligdat{$ind}->{'clus_num'}->{'1'});
  foreach $c (keys(%{$ligdat{$ind}->{'clus_num'}}))
  { printf(OL " c %2d %2d %5.1f",
    $c,$ligdat{$ind}->{'clus_num'}->{$c},100.0*$ligdat{$ind}->{'clus_num'}->{$c}/$ligdat{$ind}->{'Natom_pock'}); } 
 } 
  printf(OL "\n"); 
}

 &Write_Append_Pocket_Cluster_Information($OPT{'op'},\@cluster);

} ## $p ##

 close(OL);
#################
### FUNCTIONS ###
#################

sub Write_Append_Pocket_Cluster_Information{
 my($ofname,$cluster) = @_;

 my($OP); my($clus);
 open(OP,">>$ofname") || die "#ERROR:Can't write to '$ofname'";
 printf(OP "#[protein] [num] [Ngrid] [Vol] [Rinac_av] [Rinac_mi] [VolRinacW] [ligand_str]\n"); 
 foreach $clus (@{$cluster})
 {
  printf(OP "%s %d %s %s %s %s %s %s\n",
 $p,$clus->{'num'},$clus->{'Ngrid'},$clus->{'Vol'},
 $clus->{'Rinac_av'},$clus->{'Rinac_mi'},$clus->{'VolRinacW'},$clus->{'ligand_str'});
 }
close(OP);

} ## end of Write_Append_Pocket_Cluster_Information() ##


sub Read_Pocket_Grid_PDB_File
{
 my($fname,$pclus) = @_;

#>> FORMAT EXAMPLE <<
#REMARK  CLUSTER_PROPERTY   1 Ngrid   2091 Vol(AAA)  1070.6 Rinac(A) av  4.08 mi  3.00 VolRinacW(AA)   262.6
#REMARK  CLUSTER_PROPERTY   2 Ngrid    979 Vol(AAA)   501.2 Rinac(A) av  4.46 mi  4.00 VolRinacW(AA)   112.4
#REMARK  CLUSTER_PROPERTY   3 Ngrid    352 Vol(AAA)   180.2 Rinac(A) av  4.26 mi  4.00 VolRinacW(AA)    42.3
#REMARK  CLUSTER_PROPERTY   4 Ngrid    339 Vol(AAA)   173.6 Rinac(A) av  2.92 mi  2.50 VolRinacW(AA)    59.5
#REMARK  CLUSTER_PROPERTY   5 Ngrid    308 Vol(AAA)   157.7 Rinac(A) av  3.03 mi  2.50 VolRinacW(AA)    52.0
#REMARK  CLUSTER_PROPERTY   6 Ngrid    286 Vol(AAA)   146.4 Rinac(A) av  2.91 mi  2.50 VolRinacW(AA)    50.3
#REMARK  CLUSTER_PROPERTY   7 Ngrid    176 Vol(AAA)    90.1 Rinac(A) av  2.50 mi  2.50 VolRinacW(AA)    36.0
#REMARK  CLUSTER_PROPERTY   8 Ngrid    142 Vol(AAA)    72.7 Rinac(A) av  4.13 mi  4.00 VolRinacW(AA)    17.6
#REMARK  CLUSTER_PROPERTY   9 Ngrid    138 Vol(AAA)    70.7 Rinac(A) av  4.28 mi  4.00 VolRinacW(AA)    16.5
 my($F); my(@A); my($cnum); my($Ngrid); my($Vol); my($n) = 0; my($end) = 0;
 
 open(F,$fname) || die "#ERROR:Can't open pocket_grid_pdbfile '$fname'";
 @{$pclus} = ();
 
 while (<F>)
 {
  chomp;
  if (/^REMARK  CLUSTER\_PROPERTY /)
  {
   @A = split(/\s+/,$_);
   $cnum  = $A[2];  
   $Ngrid = $A[4];  
   $Vol   = $A[6];  
   $pclus->[$n]->{'num'}   = $cnum; 
   $pclus->[$n]->{'Ngrid'} = $Ngrid; 
   $pclus->[$n]->{'Vol'}   = $Vol; 
   $pclus->[$n]->{'Rinac_av'}   = $A[9]; 
   $pclus->[$n]->{'Rinac_mi'}   = $A[11]; 
   $pclus->[$n]->{'VolRinacW'}   = $A[13]; 
   ++$n;  
  }
  if (/^TER/) {$end = 1;}
  if ($end==1) {close(F);}
 }
 close(F);
} ## end of Read_Pocket_Grid_PDB_File() ##


sub Read_Pocket_Annotated_Ligand_File{
 my($fname,$index_list,$lig) = @_;

#>> FILE FORMAT EXAMPLE <<
# HETATM44451  P   PO4  1980      -1.417 -29.826  30.423  1.90 0.082 12.23 1.000 2
# HETATM44452  O1  PO4  1980      -0.855 -28.482  29.801  1.40 0.048 21.00 1.000 2
# HETATM44453  O2  PO4  1980      -1.406 -29.615  31.956  1.40 0.097 10.35 1.000 2
# HETATM44454  O3  PO4  1980      -0.583 -31.121  30.054  1.40 0.009110.00 1.000 2
# HETATM44455  O4  PO4  1980      -2.948 -29.940  30.001  1.40 0.209  4.78 1.000 2
# HETATM44456  P   PO4  1981     -43.841 -15.659 -26.904  1.90 0.000 -1.00 1.000
# HETATM44457  O1  PO4  1981     -44.547 -14.261 -27.212  1.40 0.000 -1.00 1.000
# HETATM44458  O2  PO4  1981     -42.586 -15.290 -25.975  1.40 0.000 -1.00 1.000
# HETATM44459  O3  PO4  1981     -43.290 -16.373 -28.175  1.40 0.000 -1.00 1.000
# HETATM44460  O4  PO4  1981     -44.901 -16.457 -26.044  1.40 0.000 -1.00 1.000
# HETATM44723 MN    MN  1901      -9.728   3.530  -0.549  1.65 0.097 10.33 0.816 8 32
# HETATM44724 MN    MN  1902      -7.558   0.898   1.139  1.65 0.311  3.21 0.805 8 32
# HETATM44725  K     K  1903      -2.182  -5.905   4.601  1.50 0.000 -1.00 0.926
# HETATM44726  K     K  1904      -9.443   7.271  -0.709  1.50 0.000 -1.00 1.000
# HETATM44727  P   PO4  1906     -10.738   0.500   0.466  1.90 0.371  2.69 1.000 32
# HETATM44728  O1  PO4  1906     -10.804   1.814  -0.230  1.40 0.327  3.06 1.000 8 32
# HETATM44729  O2  PO4  1906      -9.628   0.497   1.477  1.40 0.400  2.50 1.000 32
# HETATM44730  O3  PO4  1906     -10.538  -0.530  -0.555  1.40 0.364  2.75 1.000 32
# HETATM44731  O4  PO4  1906     -11.993   0.361   1.231  1.40 0.220  4.55 1.000 32
# 01234567890123456789012345678901234567890123456789012345678901234567890123456789
 
 my($Resi); my($Rnum);  my($str); my(@cnum); my($c);
 my($F); my($index); my(%index_hash) = ();
 open(F,$fname) || die "#ERROR:Can't open pocket_annotated_ligand_pdbfile '$fname'";
 
 @{$index_list} = (); %{$lig} = ();
 
 while (<F>)
 {
  if (/^HETATM/)
  {
   chomp;
   $Resi  = substr($_,17,3);
   $Rnum = substr($_,20,6);
#   $Rnum =~s/\s+//g;
   $index = $Resi.$Rnum;
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
 my($Nmax);

 ### Set up Maximum Cluster Information ###
 foreach $index (@{$index_list}) 
 {
  $Nmax = 0; $max_c = ''; 
  foreach $c (keys(%{$lig->{$index}->{'clus_num'}})) 
  {
   if ($lig->{$index}->{'clus_num'}->{$c}>$Nmax) 
   {
    $Nmax = $lig->{$index}->{'clus_num'}->{$c}; 
    $max_c = $c;
   }  
  
   $lig->{$index}->{'clus_num_max'} = $max_c;
   $lig->{$index}->{'Ncluster'} = scalar(@{$index_list});
  } 
 }
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

