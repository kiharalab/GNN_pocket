#!/usr/bin/perl
#
# <mkExSqPDB.pl>
#

$LastModDate = "Aug 10, 2007";

$OPT{'ip'} = 'ChConPDB';
$OPT{'op'} = '';
$OPT{'or'} = '';
$OPT{'div'} = '0/1';

if (scalar(@ARGV)<1)
{
 print "mkExSqPDB.pl [listfile] <options>\n";
 printf(" coded by T.Kawabata. LastModifile:%s\n",$LastModDate);
 printf(" -ip : Input Dir. for Input PDB File [%s]\n",$OPT{'ip'});
 printf(" -os : Out Dir. for FASTA sequence files [%s]\n",$OPT{'os'});
 printf(" -ov : Out Dir. for Vertical sequence files [%s]\n",$OPT{'ov'});
 printf(" -A   : Action ('T' or 'F') [%s]\n",$OPT{'A'}); 
 exit(1);
}

$ilistfile = $ARGV[0];
&Read_Options(\@ARGV,\%OPT);

&Read_List_File($ilistfile,\@LIST);

$Nlist = scalar(@LIST);

for ($i=0;$i<$Nlist;++$i)
{
 $p = $LIST[$i];
 print "#$p\n";
 $chain = substr($p,4,1);
 $str = sprintf("ExSqPDB %s/%s -ch %s -os %s/%s -ov %s/%s",
 $OPT{'ip'},$p,$chain,$OPT{'os'},$p,$OPT{'ov'},$p);
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

