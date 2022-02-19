#!/usr/bin/perl

$LastModDate = "Sep 18, 2005";

$OPT{'ow'}   = '';
$OPT{'op'}   = '';
$OPT{'rgbt'} = '0:1:0:0';

if (scalar(@ARGV)<1)
{
 print "grf2wrl.pl [input_grffile] <options>[output_wrlfile]\n";
 printf(" for changing SURFNET *.grf file to VRML file format\n");
 printf(" coded by T.Kawabata. LastModified:$LastModDate\n");
 printf("<options>\n");
 printf(" -ow   : output VRML(*.wrl) file [%s]\n",$OPT{'ow'});
 printf(" -rgbt : RGBcolor (R:G:B:Trans)[%s]\n",$OPT{'rgbt'});
 printf(" -op   : output PDB file [%s]\n",$OPT{'op'});
 exit(1); 
}

$igrffile = $ARGV[0];
&Read_Options(\@ARGV,\%OPT);

&Read_Grf_File($igrffile,\@Polygon,\@Vertex);

if ($OPT{'ow'} ne '') {&Write_VRML_File($OPT{'ow'},\@Polygon,\@Vertex,$OPT{'rgbt'});}
if ($OPT{'op'} ne '') {&Write_PDB_Vertex_File($OPT{'op'},\@Polygon,\@Vertex);}


exit(1);

$Npolygon = scalar(@Polygon)-1;
$Nvertex  = scalar(@Vertex)-1;
print "Npolygon $Npolygon Nvertex $Nvertex\n";
for ($i=1;$i<=$Npolygon;++$i)
{
 printf("poly %d %s %s %s\n",
$i,
$Polygon[$i]->[0],
$Polygon[$i]->[1],
$Polygon[$i]->[2]);
}



for ($i=1;$i<=$Nvertex;++$i)
{
 printf("%d %s %s %s\n",
$i,
$Vertex[$i]->[0],
$Vertex[$i]->[1],
$Vertex[$i]->[2]);
}


#################
### FUNCTIONS ###
#################


sub Read_Grf_File{
 my($fname,$Polygon,$Vertex) = @_;
 # 
 # $Polygon[$N]->[0] : Number of 1st Vertex
 # $Polygon[$N]->[1] : Number of 2nd Vertex
 # $Polygon[$N]->[2] : Number of 3rd Vertex
 #   ($N : 1 .. $Npolygon)
 # 
 # $Vertex[$M]->[0] : X-cooridinate of Vertex
 # $Vertex[$M]->[1] : Y-cooridinate of Vertex
 # $Vertex[$M]->[2] : Z-cooridinate of Vertex
 #   ($M : 1 .. $Nvertex)


 my($F); my($state) = ''; my($Nvertex) = 1; my($Npolygon) = 1;
 my($i); my($val);  my($xyz) = 0; my($vernum) = 0;
 @{$Polygon} = ();
 @{$Vertex} = ();
 open(F,$fname) || die "#ERROR:Can't open grffile \"$fname\"\n";
 while (<F>)
 {
  chomp;
  #printf("%s:%s\n",$state,$_);

  if ($state eq 'Line')
  {
  #Polygon starts:
  #       1       4       7      10      13      16      19      22      25      28
  #      31      34      37      40      43      46      49      52      55      58

  for ($i=0;$i<10;++$i) 
  { 
    $val = substr($_,8*$i,8);
    $val =~s/\s+//g;
    if ($val =~/[0-9]/)
    {
     $Polygon->[$Npolygon]->[$vernum] = $val;
     if ($vernum>=2) {$vernum = 0; ++$Npolygon;} else {++$vernum;}
    }
  }

  }

 
  if ($state eq 'Vertex')
  {
  #Vertex coordinates:
  #    -6.400    19.200   -11.224    -5.600    19.200   -11.374    -6.400    20.000
  for ($i=0;$i<8;++$i)   
   { 
     $val = substr($_,10*$i,10);
     $val =~s/\s+//g;
     if ($val ne '')
     {
      $Vertex->[$Nvertex]->[$xyz] = $val;
      if ($xyz>=2) {$xyz = 0; ++$Nvertex;} else {++$xyz;}
     } 
    }  
  } 


 
  if (/^Polygon starts/)     {$state = "Polygon";}
  if (/^Line pointers/)      {$state = "Line";}
  if (/^Paired polygon/)     {$state = "Paired";}
  if (/^Vertex coordinates/) {$state = "Vertex";}

 } # while <F> #
 close(F);

} ## end of Read_Grf_File() 



sub Write_VRML_File{
 my($ofname,$Polygon,$Vertex,$RGBTstring) = @_;
 my($Npolygon) = scalar(@{$Polygon}) - 1;
 my($Nvertex)  = scalar(@{$Vertex})  - 1;
 my($OF); my($i);
 my($R,$G,$B,$T) = split(/:/,$RGBTstring);

 printf("#Write_VRML_File(Npolygon %d Nvertex %d RGBT %s) -> \"%s\"\n",
  $Npolygon,$Nvertex,$RGBTstring,$ofname);

 open(OF,">$ofname") || die "#ERROR:Can't write to \"$ofname\""; 


 printf(OF "#VRML V2.0 utf8\n");
 printf(OF "# OFILE   %s\n",$ofname);
 printf(OF "# COMMAND %s\n",$OPT{'COMMAND'});
 printf(OF "Shape{\n");
 printf(OF " appearance Appearance{\n");
 printf(OF " material Material{\n");
 printf(OF "  diffuseColor %.2f %.2f %.2f\n",$R,$G,$B);
 printf(OF "  transparency %.2f\n",$T);
 printf(OF " }\n");
 printf(OF "}\n");
 printf(OF "geometry IndexedFaceSet{\n");
 printf(OF "coord Coordinate{\n");
 printf(OF "point[\n");
 for ($i=1;$i<=$Nvertex;++$i)
 {
  printf(OF "%s %s %s",
  $Vertex->[$i]->[0], $Vertex->[$i]->[1], $Vertex->[$i]->[2]);
  if ($i<$Nvertex) {printf(OF ",");}
  printf(OF "\n");  
 }

 printf(OF "]}\n");
 printf(OF " coordIndex[\n");
 for ($i=1;$i<=$Npolygon;++$i)
 {
  printf(OF "%d,%d,%d,-1",
  $Polygon->[$i]->[0]-1, $Polygon->[$i]->[1]-1, $Polygon->[$i]->[2]-1);
  if ($i<$Npolygon) {printf(OF ",");}
  printf(OF "\n");  
 } 
 
 printf(OF "]\n");
 printf(OF "}}\n");

 close(OF);
 
} ## end of Write_VRML_File() ##



sub Write_PDB_Vertex_File{
 my($ofname,$Polygon,$Vertex) = @_;
 my($Npolygon) = scalar(@{$Polygon}) - 1;
 my($Nvertex)  = scalar(@{$Vertex})  - 1;
 my($OF); my($i);

 printf("#Write_PDB_Vertex_File(Npolygon %d Nvertex %d RGBT %s) -> \"%s\"\n",
  $Npolygon,$Nvertex,$RGBTstring,$ofname);

 open(OF,">$ofname") || die "#ERROR:Can't write to \"$ofname\""; 

 for ($i=1;$i<=$Nvertex;++$i)
 {
  printf(OF "HETATM%5d %4s %3s %s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
        $i,"MAP","MAP",' ','1', 
        $Vertex->[$i]->[0], $Vertex->[$i]->[1], $Vertex->[$i]->[2]);
 }

 close(OF);
 
} ## end of Write_PDB_Vertex_File() ##





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

