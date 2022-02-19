#!/usr/bin/perl
#
# <mkLigRot.pl>
#

$LastModDate = "Dec 16, 2006";
$OPT{'ipdb'} = 'ChConPDB';
$OPT{'opdb'} = 'RotPDB';
$OPT{'ld'} = 'LIG';
$OPT{'A'}  = 'F';
$OPT{'T'} = '-';

if (scalar(@ARGV)<1)
{
 print "mkLigRot.pl [pdblist] <options>\n";
 print " coded by T.Kawabata. LastModified:$LastModDate\n";
 print " <options>\n";
 printf(" -ipdb   : Input Dir for PDB files [%s]\n",$OPT{'ipdb'});
 printf(" -opdb   : Input Dir for Rotated PDB files [%s]\n",$OPT{'opdb'});
 printf("  Pattern: [Atom]|[Res]|[Rnum]  ");
 printf("   ex) ' CA |GLY|543' ' N5 |FAD|*' \n");
 printf(" for PLP : -po ' C2 |PLP|*' -px ' C3 |PLP|*' -py ' N1 |PLP|*'\n");
 printf(" for FAD : -po ' N5 |FAD|*' -px ' C4A|FAD|*' -py ' C5A|FAD|*'\n");  
 printf(" for HEM : -po 'FE  |HEM|*' -px ' N A|HEM|*' -py ' N B|HEM|*'\n");
 printf(" for ADP : -po ' N9 |ADP|*' -px ' C4 |ADP|*' -py ' C8 |ADP|*'\n");
 printf(" -po :Pattern for Origin Atom [%s]\n",$OPT{'po'});
 printf(" -px :Pattern for X-axis Atom [%s]\n",$OPT{'px'});
 printf(" -py :Pattern for Y-axis Atom [%s]\n",$OPT{'py'});
 
 printf(" -A      : Action (T or F) [%s]\n",$OPT{'A'});
 exit(1);
}

$listfile = $ARGV[0];

&Read_Options(\@ARGV,\%OPT);

&Read_List_File($listfile,\@LIST);
$Nlist = scalar(@LIST);

for ($i=0;$i<$Nlist;++$i)
{
 $p = $LIST[$i]; 
 $dir   = substr($p,1,2);
 $pdb   = substr($p,0,4);
 $chain = substr($p,4,1);
 $str = sprintf("LigRot.pl %s/%s -M S -of %s/%s -po '%s' -px '%s' -py '%s'",
  $OPT{'ipdb'},$p,$OPT{'opdb'},$p,$OPT{'po'},$OPT{'px'},$OPT{'py'});
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


