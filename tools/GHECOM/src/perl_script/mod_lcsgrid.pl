#!/usr/bin/perl
#
# <modlcs_grid.pl>
#

$LastModDate = "Feb 22, 2009";

$OPT{'gw'} = 0.8;
$OPT{'of'} = '-';
if (scalar(@ARGV)<1){
 print "modlcs_grid.pl [grid_pdb_file] <options>\n";
 printf(" coded by T.Kawabata. LastModified:%s\n",$LastModDate);
 printf(" -gw   : grid_width [%s]\n",$OPT{'gw'}); 
 printf(" -of   : output filename [%s]\n",$OPT{'of'}); 
 printf(" -A    : Action ('T' or 'F') [%s]\n",$OPT{'A'}); 
 exit(1);
}

$igridfile = $ARGV[0];
&Read_Options(\@ARGV,\%OPT);
&Read_PDB_and_get_Min_Max($igridfile,\@Min,\@Max);
printf("#min %f %f %f max %f %f %f\n",$Min[0],$Min[1],$Min[2],$Max[0],$Max[1],$Max[2]);

for ($i=0;$i<3;++$i){
 $OrigPos[$i] = $Min[$i] - $OPT{'gw'};
 $Ngrid[$i] = ($Max[$i] - $OrigPos[$i])/$OPT{'gw'} + 1;
}
printf("#OrigPos %f %f %f Ngrid %d %d %d\n",$OrigPos[0],$OrigPos[1],$OrigPos[2],$Ngrid[0],$Ngrid[1],$Ngrid[2]);


&Read_PDB_and_Write_PDB_with_Grid_Info($igridfile,$OPT{'of'},\@OrigPos,\@Ngrid,$OPT{'gw'});


#################
### FUNCTIONS ###
#################

sub Read_PDB_and_Write_PDB_with_Grid_Info{
 my($ipdbfile,$opdbfile,$OrigPos,$Ngrid,$grid_width) = @_;
 my($IF,$OF,$i,$init);
 my(@Pos);  my($Natom) = 0;
 my(@min,@max);
 
 @{$Gpos} = ();

 printf("#Read_PDB_and_Write_PDB_with_Grid_Info(\"%s\") --> '%s'\n",$ipdbfile,$opdbfile);
 open(IF,$ipdbfile) || die "#ERROR:Can't open pdbfile \"$ipdbfile\"";
 open(OF,">$opdbfile") || die "#ERROR:Can't write to pdbfile \"$opdbfile\"";
 $init = 1;
 printf(OF "REMARK  COMMAND %s\n",$OPT{'COMMAND'});
 printf(OF "REMARK  DATE %s\n",&Get_Date_String());
 printf(OF "REMARK  grid_size  %4d %4d %4d\n",$Ngrid[0],$Ngrid[1],$Ngrid[2]);
 printf(OF "REMARK  grid_width   %8.3f\n",$OPT{'gw'}); 
 printf(OF "REMARK  OrigPos  %8.3f %.3f %8.3f\n",$OrigPos[0],$OrigPos[1],$OrigPos[2]); 
 
 while (<IF>){
  if ((/^ATOM/)||(/^HETATM/)){
   chomp;
   $head = substr($_,0,56);
   substr($head,0,6)  = "HETATM"; 
   substr($head,17,3) = "GRD"; 
   substr($head,13,2) = "C "; 
   $Pos[0] = substr($_,30,8);
   $Pos[1] = substr($_,38,8);
   $Pos[2] = substr($_,46,8);
   $rnum   = substr($_,23,4);
   $x = ($Pos[0]-$OrigPos[0])/$grid_width; 
   $y = ($Pos[1]-$OrigPos[1])/$grid_width; 
   $z = ($Pos[2]-$OrigPos[2])/$grid_width; 
   printf(OF "%s%6.2f%6.2f %d %d %d\n",$head,$grid_width,$rnum,$x,$y,$z);
   ++$Natom;
  }
  else {printf(OF "%s",$_);}
 } # while <IF> #
 close(IF);
 close(OF);
}



sub Read_PDB_and_get_Min_Max{
 my($ipdbfile,$Min,$Max) = @_;
 my($IF,$i,$init);
 my(@Pos);  my($Natom) = 0;
 my(@min,@max);
 
 @{$Gpos} = ();

 printf("#Read_PDB_and_get_Min_Max(\"%s\")\n",$ipdbfile);
 open(IF,$ipdbfile) || die "#ERROR:Can't open pdbfile \"$ipdbfile\"";
 $init = 1;
 while (<IF>){
  if ((/^ATOM/)||(/^HETATM/)){
   chomp;
   $Pos[0] = substr($_,30,8);
   $Pos[1] = substr($_,38,8);
   $Pos[2] = substr($_,46,8);
   if ($init==1){
    for ($i=0;$i<3;++$i){ $min[$i] = $max[$i] = $Pos[$i];}
    $init = 0;
    }
   
  for ($i=0;$i<3;++$i){ 
   if ($Pos[$i]<$min[$i]) {$min[$i] = $Pos[$i];}
   if ($Pos[$i]>$max[$i]) {$max[$i] = $Pos[$i];}
  }
   ++$Natom;
  }
 } # while <IF> #

 close(IF);
 for ($i=0;$i<3;++$i){
  $Min->[$i] = $min[$i];
  $Max->[$i] = $max[$i];
 }

} ## end of Read_PDB_and_get_Min_Max() ##



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


sub Get_Date_String{
 my(@Month) = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
 my(@date) = localtime(time());
 my($year)  = $date[5]+1900;
 my($month) = $Month[$date[4]];
 my($day)   = $date[3];
 return("$month $day $year, $date[2]:$date[1]:$date[0]");

} ## end of Get_Date_String() ##

