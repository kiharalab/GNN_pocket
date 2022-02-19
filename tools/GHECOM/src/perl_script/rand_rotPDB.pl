#!/usr/bin/perl
#
# <rand_rotPDB.pl>
# for rotating and translating PDB 
#

$LastModDate = "Nov 2, 2006";
$OPT{'ang'} = '0:0:0';
$OPT{'tvc'} = '0:0:0';
$OPT{'ip'} = '';
$OPT{'op'} = '-';
$OPT{'ig'} = '';
$OPT{'og'} = '-';
$OPT{'G'}  = 'F';
$OPT{'R'}  = 'T';
$OPT{'sr'} = 0;
if (scalar(@ARGV)<1)
{
 print "rand_rotPDB.pl [ipdbfile] <options>\n";
 print " coded by T.Kawabata. LastModified:$LastModDate\n"; 
 print " transform : X' = Rmat * X + tvec\n";
 print "<options>\n";
 printf(" -ang: phi:theta:psi (degree) [%s]\n",$OPT{'ang'}); 
 printf(" -tvc: Tx:Ty:Tz    (angstrom) [%s]\n",$OPT{'tvc'}); 
 printf(" -ip : Input PDB file [%s]\n",$OPT{'ip'}); 
 printf(" -op : Output transformed PDB file [%s]\n",$OPT{'op'}); 
 printf(" -G  : set Gcenter to origin (T or F) [%s]\n",$OPT{'G'});
 printf(" -R  : Random Rotation (T or F) [%s]\n",$OPT{'R'});
 printf(" -sr : seed for random number [%s]\n",$OPT{'sr'});
 exit(1);
}

$ipdbfile = $ARGV[0];

&Read_Options(\@ARGV,\%OPT);
if ($OPT{'R'}  eq 'T') {srand($OPT{'sr'});}

($phi,$theta,$psi) = split(/:/,$OPT{'ang'});
@Tvec = split(/:/,$OPT{'tvc'});
$deg_to_rad = 3.14159/180;
$phi   *= $deg_to_rad;
$theta *= $deg_to_rad;
$psi   *= $deg_to_rad;

### (1) Make Rotation Matrix ###
if ($OPT{'R'} eq 'T')
{
 $phi   = rand(360.0);
 $theta = rand(360.0);
 $psi   = rand(360.0);
}
&Make_Rotation_Matrix_By_Euler_Angle(\@Rmat,$phi,$theta,$psi);

### (2) Set Gcenter to the Origin ###

if (($OPT{'G'} eq 'T')||($OPT{'R'} eq 'T'))
{
  &Read_PDB_and_Cal_Gcenter($ipdbfile,\@Gcen); 
  printf("#Gcenter %f %f %f\n",$Gcen[0],$Gcen[1],$Gcen[2]);
  $Tvec[0] = -$Gcen[0]; $Tvec[1] = -$Gcen[1]; $Tvec[2] = -$Gcen[2]; 
}


printf("#Tvec   %f %f %f\n",$Tvec[0],$Tvec[1],$Tvec[2]);
printf("#Rmat0  %f %f %f\n",$Rmat[0]->[0],$Rmat[0]->[1],$Rmat[0]->[2]);
printf("#Rmat1  %f %f %f\n",$Rmat[1]->[0],$Rmat[1]->[1],$Rmat[1]->[2]);
printf("#Rmat2  %f %f %f\n",$Rmat[2]->[0],$Rmat[2]->[1],$Rmat[2]->[2]);



### (3) Read PDB files and output transformed one ###
&Read_PDB_and_Output_Transformed_PDB($ipdbfile,\@Rmat,\@Tvec,$OPT{'op'}); 


#################
### FUNCTIONS ###
#################


sub Make_Rotation_Matrix_By_Euler_Angle{
 my($R,$phi,$the,$psi) = @_;

 @{$R} = ();

 print "#phi $phi theta $the psi $psi\n";
 my($Cphi) = cos($phi);
 my($Sphi) = sin($phi);
 my($Cthe) = cos($the);
 my($Sthe) = sin($the);
 my($Cpsi) = cos($psi);
 my($Spsi) = sin($psi);

 print "#Cphi $Cphi Sphi $Sphi\n";
 print "#Cthe $Cthe Sthe $Sthe\n";
 print "#Cpsi $Cpsi Sphi $Spsi\n";

 $R->[0]->[0] =  $Cpsi*$Cthe*$Cphi - $Spsi*$Sphi;
 $R->[1]->[0] = -$Spsi*$Cthe*$Cphi - $Cpsi*$Sphi;
 $R->[2]->[0] =  $Sthe*$Cphi;
 
 $R->[0]->[1] =  $Cpsi*$Cthe*$Sphi + $Spsi*$Cphi;
 $R->[1]->[1] = -$Spsi*$Cthe*$Sphi + $Cpsi*$Cphi;
 $R->[2]->[1] =  $Sthe*$Sphi;

 $R->[0]->[2] = -$Cpsi*$Sthe;
 $R->[1]->[2] =  $Spsi*$Sthe;
 $R->[2]->[2] =  $Cthe;


} ## end of Make_Rotation_Matrix_By_Euler_Angle() ##







sub Read_PDB_and_Output_Transformed_PDB{
 my($ipdbfile,$Rmat,$Tvec,$opdbfile) = @_;
 my($IF); my($OF);
 my(@Pos); my(@TPos); my($head); my($tail);
 my($i); my($j);

 if ($opdbfile ne '-')
 { printf("#Read_PDB_and_Output_Transformed_PDB()-->\"%s\"\n",$opdbfile); } 
 open(IF,$ipdbfile) || die "#ERROR:Can't open pdbfile \"$ipdbfile\"";
 open(OF,">$opdbfile") || die "#ERROR:Can't write to pdbfile \"$opdbfile\"";
 printf(OF "REMARK  COMMAND %s\n",$OPT{'COMMAND'});
 while (<IF>)
 { 
  if ((/^ATOM/)||(/^HETATM/))
  {
   chomp;
   $Pos[0] = substr($_,30,8); 
   $Pos[1] = substr($_,38,8); 
   $Pos[2] = substr($_,46,8);
   $head = substr($_,0,30);
   $tail = substr($_,54);
   $ch = substr($head,21,1);
   &Transform_Rx_plus_tvec(\@TPos,\@Pos,$Rmat,$Tvec);
   printf(OF "%30s%8.3f%8.3f%8.3f%s\n",$head,$TPos[0],$TPos[1],$TPos[2],$tail);
  } 
  else {print OF "$_";}
 
 } # while <IF> #
 close(OF);
 close(IF);

} ## end of Read_PDB_and_Output_Transformed_PDB() ##




sub Read_PDB_and_Cal_Gcenter{
 my($ipdbfile,$Gpos) = @_;
 my($IF); 
 my(@Pos);  my($Natom) = 0;

 @{$Gpos} = ();

 printf("#Read_PDB_and_Cal_Gcenter(\"%s\")\n",$ipdbfile);
 open(IF,$ipdbfile) || die "#ERROR:Can't open pdbfile \"$ipdbfile\"";
 while (<IF>)
 { 
  if ((/^ATOM/)||(/^HETATM/))
  {
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

 if ($Natom > 0)
 { $Gpos->[0] /= $Natom;
   $Gpos->[1] /= $Natom;
   $Gpos->[2] /= $Natom; }
 close(IF);

} ## end of Read_PDB_and_Cal_Gcenter() ##










sub Transform_Rx_plus_tvec{
 my($tpos,$pos,$Rmat,$Tvec) = @_;

 my($i); my($j);

 for ($i=0;$i<3;++$i)
 {
  $tpos->[$i] = $Tvec->[$i];
  for ($j=0;$j<3;++$j) { $tpos->[$i] += $Rmat->[$i]->[$j] * $pos->[$j];}
 }

} ## end of Transform_Rx_plus_tvec() ## 



sub Transform_RAtR_Matrix3D{
 my($C,$A,$R) = @_; ## C = R*A*tR ## 
 my($i); my($j); my($k);  my(@buff);

 ## (1) buff = R * A ## 
 for ($i=0;$i<3;++$i)
 {
  for ($j=0;$j<3;++$j)
  {
    $buff[$i]->[$j] = 0.0;
    for ($k=0;$k<3;++$k) { $buff[$i]->[$j] += $R->[$i]->[$k]* $A->[$k]->[$j]; }
  }
 }                                                                                                       
 
 ## (2) C = buff * tR ##
 for ($i=0;$i<3;++$i)
 { 
  for ($j=0;$j<3;++$j)
  {
    $C->[$i]->[$j] = 0.0;
    for ($k=0;$k<3;++$k) { $C->[$i]->[$j] += $buff[$i]->[$k] * $R->[$j]->[$k];}
  }
 }                                                                                                        

} ## end of Transform_RAtR_Matrix3D() ## 






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
     if (($x1 !~ /^\-\w+/)&&(length($x1)>0)) { $_[1]->{$x} = $x1; ++$i; }
     else { $_[1]->{$x} = 1;}
   }
  ++$i;
 }
                                                                                               
} ## end of Read_Options() ##

