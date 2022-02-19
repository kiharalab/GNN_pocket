#!/usr/bin/perl
#
# <mkLigand.pl>
#

$LastModDate = "Sep 9, 2008";
$OPT{'id'}     = 'ChConPDB';
$OPT{'od'} = 'LIG';
$OPT{'O'}  = 'F';
$OPT{'o'}  = 'F';
$OPT{'A'}  = 'F';
$OPT{'T'} = '-';
$OPT{'div'} = '0/1';

if (scalar(@ARGV)<1){
 print "mkLigand.pl [pdblist] <options>\n";
 print " coded by T.Kawabata. LastModified:$LastModDate\n";
 print " <options>\n";
 printf(" -id     : Input Directory of PDB file [%s]\n",$OPT{'id'});
 printf(" -od     : Directory of ligand file [%s]\n",$OPT{'od'});
 printf(" -div    : Division (bunshi/bunbo) [%s]\n",$OPT{'div'});
 printf(" -O      : One line output (T or F) [%s]\n",$OPT{'O'});
 printf(" -o      : oneline_with_molformula (T or F) [%s]\n",$OPT{'o'});
 printf(" -A      : Action (T or F) [%s]\n",$OPT{'A'});
 exit(1);
}

$listfile = $ARGV[0];

&Read_Options(\@ARGV,\%OPT);

&Read_List_File($listfile,\@LIST);
$Nlist = scalar(@LIST);


($bunshi,$bunbo) = split(/\//,$OPT{'div'});
$Nlist = scalar(@LIST);
$Nstart = $bunshi*int($Nlist/$bunbo);
$Nend   = ($bunshi+1)*int($Nlist/$bunbo);
if ($bunshi>=($bunbo-1)) { $Nend = $Nlist;}
                                                                                             
print "#Nlist $Nlist bunshi/bunbo $bunshi / $bunbo start $Nstart end $Nend\n";

for ($i=$Nstart;$i<$Nend;++$i)
{
 $p = $LIST[$i]; 
 $dir   = substr($p,1,2);
 $pdb   = substr($p,0,4);
 $chain = substr($p,4,1);
 $pdbfile = $OPT{'id'}.'/'.$p;
 $outfile = $OPT{'od'}.'/'.$p; 
 $str = sprintf("Ligand $pdbfile -M P -ch $chain -n2h F -tz T -mos T -cos T -cbs T -cbs T -cdrpp T -n2h T > $outfile\n");
 if ($OPT{'O'} eq 'T'){
  $str = sprintf("Ligand $pdbfile -M P -ch $chain -n2h F -tz T -mos T -cos T -cbs T -cbs T -cdrpp T -n2h T -O O\n");
 } 
 if ($OPT{'o'} eq 'T'){
  $str = sprintf("Ligand $pdbfile -M P -ch $chain -n2h F -tz T -mos T -cos T -cbs T -cbs T -cdrpp T -n2h T -O o\n");
 } 
 
 if ($OPT{'A'} eq 'T') { system($str);} else {print "#$str\n";}

} # $i loop #

close(OL);


###################
#### FUNCTIONS ####
###################


sub Get_Date_String{
 my(@Month) = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
 my(@date) = localtime(time());
 my($year)  = $date[5]+1900;
 my($month) = $Month[$date[4]];
 my($day)   = $date[3];
 return("$month $day $year, $date[2]:$date[1]:$date[0]");
} ## end of Get_Date() ##

                                                                                
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


