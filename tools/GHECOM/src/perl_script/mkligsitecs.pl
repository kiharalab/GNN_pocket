#!/usr/bin/perl
#
# <mkligsitecs.pl>
#

$LastModDate = "Feb 22, 2009";

$OPT{'ipdb'} = 'ChConPDB';

$OPT{'ores'} = 'ORES';
$OPT{'A'} = 'F';
$OPT{'M'} = 'M';
$OPT{'div'} = '0/1';
$OPT{'og'} = 'LIGSITE_GRID';
$OPT{'C'} = 'T';
$OPT{'rk'} = -1;
$OPT{'n'} = 3;
$OPT{'t'} = 6;
$OPT{'tdir'} = '/home/takawaba/tool/pass_2.0.36/pass_tool_kawabata';
$OPT{'sdir'} = '/home/takawaba/work/ghecom/src/perl_script';

$OPT{'tmpf'} = sprintf("tmp%s.pdb",&Time_String());
$OPT{'tmpp'} = sprintf("po%s.pdb",&Time_String());

if (scalar(@ARGV)<1){
 print "mkligsitecs.pl [listfile] <options>\n";
 printf(" coded by T.Kawabata. LastModified:%s\n",$LastModDate);
 printf(" -ipdb : dir for input PDB file [%s]\n",$OPT{'ipdb'}); 
 printf(" -tdir : pass_tool directory [%s]\n",$OPT{'tdir'}); 
 printf(" -sdir : ghecom perl script directory [%s]\n",$OPT{'tdir'}); 
 printf(" -og   : Output dir for grid PDB file [%s]\n",$OPT{'og'}); 
 printf(" -div  : Job division [%s]\n",$OPT{'div'}); 
 printf(" -n    : Number_of_pockets[%s]\n",$OPT{'n'});
 printf(" -t    : SSS_threshold [%s]\n",$OPT{'t'});
 printf(" -tmpf : tempolary input PDB file [%s]\n",$OPT{'tmpf'});
 printf(" -tmpp : tempolary output pocket PDB file [%s]\n",$OPT{'tmpp'});

 printf(" -A    : Action ('T' or 'F') [%s]\n",$OPT{'A'}); 
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
 $chain = substr($p,4,1);

 # Prepare 'correct' PDB file ##
  $ipdbfile = sprintf("%s/%s",$OPT{'ipdb'},$p);
  $str = sprintf("%s/reg_pdb.pl  %s/%s -ch %s -rn T -rr T -lr T -nc T -rm T > %s",
          $OPT{'tdir'},$OPT{'ipdb'},$p,$chain,$OPT{'tmpf'});
  printf("#$str\n");
 if ($OPT{'A'} eq 'T') {system($str);}
 

  $str = sprintf("lcs -i %s -s 0.8 -n %d -t %d -o %s", $OPT{'tmpf'},$OPT{'n'},$OPT{'t'},$OPT{'tmpp'});

  printf("#$str\n");
  if ($OPT{'A'} eq 'T') {system($str);}

  #$str = sprintf("mv pocket_r.pdb  %s/%s", $OPT{'og'},$p);
  $str = sprintf("%s/mod_lcsgrid.pl %s  -of %s/%s", $OPT{'sdir'},$OPT{'tmpp'},$OPT{'og'},$p);
  printf("#$str\n");
  if ($OPT{'A'} eq 'T') {system($str);}

}

#################
### FUNCTIONS ###
#################


sub Read_PDB_and_Cal_Gcenter{
 my($ipdbfile,$Gpos) = @_;
 my($IF);
 my(@Pos);  my($Natom) = 0;

 @{$Gpos} = ();

 printf("#Read_PDB_and_Cal_Gcenter(\"%s\")\n",$ipdbfile);
 open(IF,$ipdbfile) || die "#ERROR:Can't open pdbfile \"$ipdbfile\"";
 while (<IF>){
  if ((/^ATOM/)||(/^HETATM/)){
   chomp;
   $Pos[0] = substr($_,30,8);
   $Pos[1] = substr($_,38,8);
   $Pos[2] = substr($_,46,8);
   $Gpos->[0] += $Pos[0];
   $Gpos->[1] += $Pos[1];
   $Gpos->[2] += $Pos[2];
   ++$Natom;
  }
 } # while <IF> #

 if ($Natom > 0){ 
   $Gpos->[0] /= $Natom;
   $Gpos->[1] /= $Natom;
   $Gpos->[2] /= $Natom; 
 }

 close(IF);

} ## end of Read_PDB_and_Cal_Gcenter() ##



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


sub Time_String{
 my(@Month) = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
 my(@date) = localtime(time());
 my($year)  = $date[5]+1900;
 my($month) = $Month[$date[4]];
 my($day)   = $date[3];
 return("$date[2]_$date[1]_$date[0]");
} ## end of Time_String() ##


