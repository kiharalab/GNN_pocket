#!/usr/bin/perl
#
# <cmpPoGrdLig.pl>
#

$LastModDate = "June 13, 2010";

$OPT{'ipdb'} = 'ChConPDB';
$OPT{'igA'}  = 'OPO';
$OPT{'ilg'}  = 'OLG';
$OPT{'tgAi'} = 1;
$OPT{'tgAx'} = 254;
$OPT{'A'} = 'F';
$OPT{'M'} = 'M';
$OPT{'div'} = '0/1';
$OPT{'nh'} = 7;
$OPT{'nc'} = 3;
$OPT{'sele'} = 'T';
$OPT{'thc'} = 0;
$OPT{'thci'} = 0;
$OPT{'gw'} = 0.8;
$OPT{'nn'} = 6;
$OPT{'cgA'} = '';

if (scalar(@ARGV)<1){
 print "cmpPoGrdLig.pl [listfile] <options>\n";
 print "  for comparison between pocket grids and ligand atoms\n";
 printf(" coded by T.Kawabata. LastModified:%s\n",$LastModDate);
 printf(" -igA  : dir for opo [%s]\n",$OPT{'igA'}); 
 printf(" -ilg  : dir for ilg  [%s]\n",$OPT{'ilg'}); 
 printf(" -tgAi : min threshold value for grid A [%s]\n",$OPT{'tgAi'}); 
 printf(" -tgAx : max threshold value for grid A [%s]\n",$OPT{'tgAx'}); 
 printf(" -gw   : grid_width  [%s]\n",$OPT{'gw'}); 
 printf(" -sele : Do selection ligands 'T' or 'F' [%s]\n",$OPT{'sele'});
 printf(" -nh   : minimum number of heavy atoms for ligands [%d]\n",$OPT{'nh'});
 printf(" -nc   : minimum number of carbons for ligands [%d]\n",$OPT{'nc'});
 printf(" -thc  : max thre. char value for clustering ('255' for all)(normaly same as -tgAx) [0]\n",$OPT{'thc'});
 printf(" -thci : min thre. char value for clustering [0]\n",$OPT{'thci'});
 printf(" -rkc  : rank keep for large cluster. '-1' or '255' for everything  [-1]\n",$OPT{'rkc'});
 printf(" -nn   : Neighbor Number for clustering (6,18,or 26) [%d]\n",$OPT{'nn'});
 printf(" -cgA : Cluster_num(chainID) select for gridA. '1':1,'2':1,2,'3':1,2,3 '-':all. [-]\n",$OPT{'cgA'});
 printf(" -A    : Action ('T' or 'F') [%s]\n",$OPT{'A'}); 
 exit(1);
}

$ilistfile = $ARGV[0];
&Read_Options(\@ARGV,\%OPT);

printf("#%s\n",$OPT{'COMMAND'});
&Read_List_File($ilistfile,\@LIST);
$Nlist = scalar(@LIST);

$PID = $$;
$tmpligfile = "tmp$PID.pdb";

for ($i=0;$i<$Nlist;++$i){
 $p = $LIST[$i];
 $chain = substr($p,4,1);
 if ($OPT{'sele'} eq 'T'){
   $ipdbfile = sprintf("%s/%s",$OPT{'ilg'},$p);
   if ($OPT{'A'} eq 'T'){ 
     &Read_PDB_Molecule_file($ipdbfile,\@mol);
     #printf("%s Nmol %d\n",$ipdbfile,scalar(@mol));
     &Output_Proper_Ligand_PDB_File($tmpligfile,\@mol,$OPT{'nh'},$OPT{'nc'});
   }
 }
 else{
   $str = sprintf("cp %s/%s $tmpligfile",$OPT{'ilg'},$p);
   if ($OPT{'A'} eq 'T') {system($str);}
   else {printf("#$str\n");}
 }


 $str = sprintf("ghecom -M L -igA %s/%s -ilg $tmpligfile -gw %s -tgAi %s -tgAx %s -nn %s",
   $OPT{'igA'},$p,$OPT{'gw'},$OPT{'tgAi'},$OPT{'tgAx'},$OPT{'nn'});
  
 if ($OPT{'thc'}>0) { $str .= sprintf(" -thc %d",$OPT{'thc'}); }
 if ($OPT{'thci'}>0) { $str .= sprintf(" -thci %d",$OPT{'thci'}); }
 if ($OPT{'rkc'}>0) { $str .= sprintf(" -rkc %d",$OPT{'rkc'}); }
 if ($OPT{'cgA'}>0) { $str .= sprintf(" -cgA %s",$OPT{'cgA'}); }

 $str .= '|grep GRID_LIG_COMP'; 
 if ($OPT{'A'} eq 'T') {$ret = system($str);}
 else {printf("#$str\n");}
}

#################
### FUNCTIONS ###
#################

sub Read_PDB_Molecule_file()
{
 my($ifname,$mol) = @_;
 open(IF,$ifname)||die "#ERROR:Can't open pdbfile '$ifname'";
 @{$mol} = (); 
 $Nmol = 0;
 while (<IF>){
   chomp;
  if ((/^ATOM/)||(/^HETATM/)){
    push(@{$mol->[$Nmol]},$_); 
  }
  elsif ((/^TER/)||(/^END/)){
    $Nmol += 1; 
  }
 }
 close(IF);

} # end of Read_PDB_Molecule_File() #


sub Output_Proper_Ligand_PDB_File(){
 my($ofname,$mol,$NheavyMin,$NcarbonMin) = @_;
 my($m,$a,$out,$Nheavy,$Ncarbon);

 @notallow = ("pro","dna","rna","GOL","MPD","TRS","MES","BOG","EPE","DTT","MRD","PG4");
 foreach $a (@notallow) {$NotAllow{$a} = 1;}
 
 #printf("#Output_Proper_Ligand_PDB_File() --> '%s'\n",$ofname);
 open(OF,">$ofname")||die "#ERROR:Can't write to '$ofname'";
 for ($m=0;$m<scalar(@{$mol});++$m){
   ### (1) Check Molecule ##
   $accept = 1;
   $Nheavy = $Ncarbon = 0;
   for ($a=0;$a<scalar(@{$mol->[$m]});++$a){
     $res  = substr($mol->[$m]->[$a],17,3);
     if ($NotAllow{$res}==1){$accept = 0;}
     $atom = substr($mol->[$m]->[$a],12,4);
     #print "res '$res' atom '$atom'\n";
     if ($atom!~/^H/){
       $ele = substr($atom,1,1);
       if ($ele eq 'C') {$Ncarbon += 1;}
       $Nheavy += 1;
      } 
   } 
   if (($Nheavy < $NheavyMin)||($Ncarbon < $NcarbonMin)) {$accept = 0;}
   ### (2) Output Molecule ##
   
   #printf("#Nheavy %d Ncarbon %d accept %d\n",$Nheavy,$Ncarbon,$accept);
  if ($accept == 1){ 
      for ($a=0;$a<scalar(@{$mol->[$m]});++$a){ printf(OF "%s\n",$mol->[$m]->[$a]); } 
    } 
    else {
     for ($a=0;$a<scalar(@{$mol->[$m]});++$a){ printf(OF "#%s\n",$mol->[$m]->[$a]); } 
    }
    printf(OF "TER\n");
  } # $m # 

 close(OF);
} ## end of Output_Proper_Ligand_PDB_File() ##


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

