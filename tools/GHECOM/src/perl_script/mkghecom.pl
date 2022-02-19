#!/usr/bin/perl
#
# <mkghecom.pl>
#

$LastModDate = "June 6, 2010";

$OPT{'ipdb'} = 'ChConPDB';
$OPT{'opp'}  = '';
$OPT{'ilg'}  = 'ChConPDB';
$OPT{'olg'}  = 'OLG';
$OPT{'opdb'} = 'OPDB';
$OPT{'ores'} = 'ORES';
$OPT{'gw'}    = '0.8';
$OPT{'rli'}   = '2';
$OPT{'rlx'}   = '10';
$OPT{'br'}    = '0.5';
$OPT{'A'} = 'F';
$OPT{'M'} = 'M';
$OPT{'div'} = '0/1';
$OPT{'oprb'} = '';
$OPT{'ah'} = 'A';
$OPT{'dwhc'} = '';

if (scalar(@ARGV)<1){
 print "mkghecom.pl [listfile] <options>\n";
 printf(" coded by T.Kawabata. LastModified:%s\n",$LastModDate);
 printf("<options>\n");
 printf(" -M   : MODE 'D'ilation 'E'rosion molceular 'C'losing 'O'pening\n");
 printf("             'P'ocket(masuya_doi) 'p'ocket(kawabata_go) ca'V'ity\n");
 printf("             'M'ultiscale_closing/pocket\n");
 printf("             'I'nterface_bwn_two_chains\n");
 printf("             'G'rid_comparison_binary 'g'rid_comparison_mutiscale [%s]\n",$OPT{'M'});
 printf(" -ipdb : dir for input PDB file [%s]\n",$OPT{'ipdb'}); 
 printf(" -opo  : dir for opo [%s]\n",$OPT{'opo'}); 
 printf(" -opdb : dir for opdb [%s]\n",$OPT{'opdb'}); 
 printf(" -olg  : dir for olg  [%s]\n",$OPT{'olg'}); 
 printf(" -ilg  : dir for ilg  [%s]\n",$OPT{'ilg'}); 
 printf(" -gw   : grid_width  [%s]\n",$OPT{'gw'}); 
 printf(" -mp   : Min small probe number of pocket cluster [%s]\n",$OPT{'mp'});
 printf(" -mg   : Min grid size          of pocket cluster [%s]\n",$OPT{'mg'});
 printf(" -rli  : Radius for min_large probe spheres [%s]\n",$OPT{'rli'});
 printf(" -rlx  : Radius for max_large probe spheres [%s]\n",$OPT{'rlx'});
 printf(" -br   : bin of large probe radius for MODE 'M' [%s]\n",$OPT{'br'});
 printf(" -ah   : Read only ATOM 'A', only HETATM 'H', both 'B' [%s]\n",$OPT{'ah'});
 printf(" -div  : Job division [%s]\n",$OPT{'div'}); 
 printf(" -A    : Action ('T' or 'F') [%s]\n",$OPT{'A'}); 
 printf("<options for spherical probes>\n");
 printf(" -oprb : dir for spherical probes  [%s]\n",$OPT{'oprb'}); 
 printf(" -dwhc : dis thre for weighted hierarchical clustering of spherical probes[%s]\n",$OPT{'dwhc'});
 printf(" -opoc : Output Clusterd Pocket 3D grid PDB files.[%s]\n",$OPT{'opoc'});
 printf(" -rkc  : rank keep for large cluster for MODE 'M' and 'R' [%d]\n",$OPT{'rkc'});
 exit(1);
}

$ilistfile = $ARGV[0];
&Read_Options(\@ARGV,\%OPT);



&Read_List_File($ilistfile,\@LIST);


($bunshi,$bunbo) = split(/\//,$OPT{'div'});
$Nlist = scalar(@LIST);
$Nstart = $bunshi*int($Nlist/$bunbo);
$Nend   = ($bunshi+1)*int($Nlist/$bunbo);
if ($bunshi>=($bunbo-1)) { $Nend = $Nlist;}
print "#Nlist $Nlist bunshi/bunbo $bunshi / $bunbo start $Nstart end $Nend\n";

for ($i=$Nstart;$i<$Nend;++$i){
 $p = $LIST[$i];
 print "#$p\n";
 $chain = substr($p,4);
 $str = sprintf("ghecom %s/%s -M %s -ah %s",$OPT{'ipdb'},$p,$OPT{'M'},$OPT{'ah'});

 if ($chain ne '') { $str .= sprintf(" -ch %s",$chain);}
 if ($OPT{'opo'} ne '')   {$str .= sprintf(" -opo %s/%s",$OPT{'opo'},$p);}
 if ($OPT{'ilg'} ne '')   {$str .= sprintf(" -ilg %s/%s",$OPT{'ilg'},$p);}
 if ($OPT{'opdb'} ne '')  {$str .= sprintf(" -opdb %s/%s",$OPT{'opdb'},$p);}
 if ($OPT{'olg'} ne '')   {$str .= sprintf(" -olg %s/%s",$OPT{'olg'},$p);}
 if ($OPT{'ores'} ne '')  {$str .= sprintf(" -ores %s/%s",$OPT{'ores'},$p);}
 if ($OPT{'gw'}  ne '')   {$str .= sprintf(" -gw %s",$OPT{'gw'});}
 if ($OPT{'rli'}  ne '')  {$str .= sprintf(" -rli %s",$OPT{'rli'});}
 if ($OPT{'rlx'}  ne '')  {$str .= sprintf(" -rlx %s",$OPT{'rlx'});}
 if ($OPT{'br'}  ne '')   {$str .= sprintf(" -br %s",$OPT{'br'});}
 if ($OPT{'oprb'}  ne '')   {$str .= sprintf(" -sprb T -oprb %s/%s",$OPT{'oprb'},$p);}
 if ($OPT{'dwhc'}  ne '')   {$str .= sprintf(" -dwhc %s",$OPT{'dwhc'});}
 if ($OPT{'opoc'} ne '')   {$str .= sprintf(" -opoc %s/%s",$OPT{'opoc'},$p);}
 if ($OPT{'rkc'} ne '')   {$str .= sprintf(" -rkc %s",$OPT{'rkc'});}

 printf("#$str\n");
 if ($OPT{'A'} eq 'T') {system($str);}
}

#################
### FUNCTIONS ###
#################

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
  if ($x =~/^\-/)
   { $x =~s/^\-//;
     if ((length($x1)>0) && ($x1 !~/\-/)) { $_[1]->{$x} = $x1; ++$i; }
     else { $_[1]->{$x} = '';}
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

