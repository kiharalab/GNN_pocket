#!/usr/bin/perl
#
# <mksurfnet.pl>
#

$LastModDate = "Feb 16, 2009";

$OPT{'ipdb'} = 'ChConPDB';
$OPT{'A'} = 'F';
$OPT{'M'} = 'M';
$OPT{'div'} = '0/1';
$OPT{'sdir'} = '/home/takawaba/tool/surfnet/source';
$OPT{'tdir'} = '/home/takawaba/tool/surfnet/surfnet_tool_kawabata';
$OPT{'ow'} = 'SURFNET_VRML';
$OPT{'og'} = 'SURFNET_GRID';
$OPT{'S'} = 'T';
$OPT{'W'} = 'T';
$OPT{'G'} = 'T';
$OPT{'rk'} = -1;

if (scalar(@ARGV)<1){
 print "mksurfnet.pl [listfile] <options>\n";
 printf(" coded by T.Kawabata. LastModified:%s\n",$LastModDate);
 printf(" -ipdb : dir for input PDB file [%s]\n",$OPT{'ipdb'}); 
 printf(" -sdir : Surfnet source directory [%s]\n",$OPT{'sdir'}); 
 printf(" -tdir : Surfnet tool directory [%s]\n",$OPT{'tdir'}); 
 printf(" -ow   : Output dir for VRML file [%s]\n",$OPT{'ow'}); 
 printf(" -og   : Output dir for grid PDB file [%s]\n",$OPT{'og'}); 
 printf(" -S    : Do surfnet  'T' or 'F' [%s]\n",$OPT{'S'});
 printf(" -W    : Do grf2wrl  'T' or 'F' [%s]\n",$OPT{'W'});
 printf(" -rk   : Max Rank of cluster size for output PDB file[%d]\n",$OPT{'rk'}); 
 printf(" -G    : Do gridcomp 'T' or 'F' [%s]\n",$OPT{'W'});
 printf(" -div  : Job division [%s]\n",$OPT{'div'}); 
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

 ## Do 'S'urfnet ##
 if ($OPT{'S'} eq 'T'){ 
   $str = sprintf("%s/SurfnetCleft.pl %s/%s -ch %s -surfdir %s -Asn T -Asf T",
    $OPT{'tdir'},$OPT{'ipdb'},$p, $chain, $OPT{'sdir'});

   printf("#$str\n");
   if ($OPT{'A'} eq 'T') {system($str);}
 }

 ## Change grf to 'W'rl ##
 if ($OPT{'W'} eq 'T'){ 
   $str = sprintf("%s/grf2wrl.pl gaps.grf -ow %s/%s",$OPT{'tdir'},$OPT{'ow'},$p);
   printf("#$str\n");
   if ($OPT{'A'} eq 'T') {system($str);}
 }

 ## Change wrl to 'G'rid-pdb file##
 if ($OPT{'G'} eq 'T'){ 
   $str = sprintf("gridcomp -M G -iwA %s/%s -og %s/%s -gw 0.8",$OPT{'ow'},$p,$OPT{'og'},$p,$OPT{'rk'});
   printf("#$str\n");
   if ($OPT{'A'} eq 'T') {system($str);}
 }





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
 while ($i<scalar(@ARGV))
 {
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

