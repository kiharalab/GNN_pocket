#!/usr/bin/perl
#
# <ext_fhijk.pl>
#

$LastModDate = "Sep 15, 2008";

$OPT{'of'} = '-';
$OPT{'ox'} = 'exclude.out';

if (scalar(@ARGV)<1){
 print "ext_fhijk.pl [pdblist] <options>\n";
 print " for removing proteins for class f,h,i,j,k from the SCOP reresentative list\n";
 print " coded by T.Kawabata. LastModified:$LastModDate\n";
 print " <options>\n";
 printf(" -of : output file name [%s]\n",$OPT{'of'});
 printf(" -ox : output file name for removed proteins [%s]\n",$OPT{'ox'});

 exit(1);
}

$listfile = $ARGV[0];

&Read_Options(\@ARGV,\%OPT);

&Read_List_File($listfile,\@LIST,\%CLASS);

$ofname = $OPT{'of'};
$oxname = $OPT{'ox'};
open(OF,">$ofname")||die "#ERROR:can't write to '$ofname'";
open(OX,">$oxname")||die "#ERROR:can't write to '$oxname'";
foreach $pro (@LIST){
 if (($CLASS{$pro} =~/f\./)|| ($CLASS{$pro} =~/h\./)||
     ($CLASS{$pro} =~/i\./)|| ($CLASS{$pro} =~/j\./)||
     ($CLASS{$pro} =~/k\./) ) {
   printf(OX "%s %s\n",$pro,$CLASS{$pro});
 }
 else{
   printf(OF "%s %s\n",$pro,$CLASS{$pro});
 }
}
close(OF);
close(OX);

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
 my($fname,$list,$prop) = @_;
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
 @{$list} = (); %{$prop} = ();

 while (<F>){
  if ($_!~/^#/){
   chomp;
   @D = split(/\s+/,$_);
   $id = $D[0];
   push(@{$list},$id);
   $_ =~s /^$id\s+//;
   $prop->{$id} = $_;
  }
 }
 close(F);
} ## end of Read_List_File()


