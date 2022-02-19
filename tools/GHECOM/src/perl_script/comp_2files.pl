#!/usr/bin/perl

$LastModDate = "July 15, 2007";

$OPT{'cA'} = '';
$OPT{'cB'} = '';

if (scalar(@ARGV)<1){
 print "comp_2files.pl [fileA] [fileB] <options>\n";
 printf(" coded by T.Kawabata. LastModDate:%s\n",$LastModDate);
 printf("<options>\n");
 printf(" -cA : column number for fileA [%s]\n",$OPT{'cA'}); 
 printf(" -cB : column number for fileA [%s]\n",$OPT{'cB'}); 
 exit(1);
}

$datfileA = $ARGV[0];
$datfileB = $ARGV[1];

&Read_Options(\@ARGV,\%OPT);

&Read_File_Lines($datfileA,\@LineA);
&Read_File_Lines($datfileB,\@LineB);

for ($i=0;$i<scalar(@LineA);++$i)
{
 $lineA = $LineA[$i];
 $lineB = $LineB[$i];
 $lineA =~s/^\s+//; 
 $lineB =~s/^\s+//; 
 if (($OPT{'cA'} >0)&&($OPT{'cB'} >0))
 {
  @A = split(/\s+/,$lineA); 
  @B = split(/\s+/,$lineB);
  printf("%s %s ",$A[$OPT{'cA'}-1],$B[$OPT{'cB'}-1]); 
 }
 printf("%s %s\n",$lineA,$lineB);
}


###################
#### FUNCTIONS ####
###################

sub Read_File_Lines{
 my($fname,$line) = @_;
 my($IF);

 @{$line} = ();
 open(IF,$fname)||die "#ERROR:Can't open '$fname'";
 while (<IF>){ 
  chomp;
  if ($_!~/^#/){
   push(@{$line},$_); 
  }
 } 
 close(IF);

} ## end of Read_File_Lines() ##


sub Read_Options{
 # $_[0] : ref of \@ARGV
 # $_[1] : ref of \%OPT
 # This script is reading following style options :
 #   psiscan.pl org41list -lib 95pdb01Mar4Mx -tail -I -C
 # In principle, the format is the style like  "[-option] [value]"
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
     if (($x1 !~ /^\-\w+/)&&(length($x1)>0)) { $_[1]->{$x} = $x1; ++$i; }
   }
  ++$i;
 }
} ## end of Read_Options() ##

