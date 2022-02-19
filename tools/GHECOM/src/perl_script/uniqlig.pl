#!/usr/bin/perl
#
#

if (scalar(@ARGV)<1){
 print "uniqlig.pl  [one-line ligand file]\n";
 exit(1);
}

$ifname = $ARGV[0];

open(F,$ifname)||die "#ERROR:Can't open '$ifname'";
#>> FILE FORMAT EXAMPLE <<
#1kvdB   77 pro A SO4 - SO4 - SO4 -
#1upsA  402 pro B  CA A
#2a5sA  278 GLU -
#1j26A  112
#1qp6A   35 pro B
#2czcA  334 pro D PO4 - PO4 - NAD -
#2iuwA  205  FE A  FE A  FE A AKG A BME A BME A BME A BME A
#1pyuA   28 pro B
#1q08A   94 pro B PO4 -  ZN A  ZN A

while (<F>){
 chomp;
 if ($_!~/^#/){
  @col = split(/\s+/,$_);
  $i = 2;
  while ($i<scalar(@col)){
    $Nlig{$col[$i]} += 1;
    $i += 2;
  }
 } 
}
close(F);


@list = sort {$Nlig{$b} <=> $Nlig{$a}} keys(%Nlig);

foreach $lig (@list){
 if ($lig=~/[a-z]/){ printf("#");}
 printf("%s %d\n",$lig,$Nlig{$lig});
}

