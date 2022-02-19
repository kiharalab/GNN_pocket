##
## <res_base_pock.py>
##  dealing with residue-based pocket file 
##  
##>> FILE FORMAT EXAMPLE <<
#>> FILE EXAMPLE <<
##COLUMN  1|RNUM              |Residue Number
##COLUMN  2|CHAIN             |Chain Identifier
##COLUMN  3|RES               |Three-letter residue name
##COLUMN  4|shellAcc          |shell accessibility (%)
##COLUMN  5|Rinacc            |averaged Rinaccess (A)
##COLUMN  6|Natom             |Number of atoms
##COLUMN  7|Natom_contact     |Number of contacting atoms with ligand
##COLUMN  8|pocketness        |sum of 1/[Rpocket] /(1/[Rmin]*[vol of shell]) (%)
##COLUMN  9|pocketness_clus 1 |sum of 1/[Rpocket_for_pocketcluster 1] /(1/[Rmin]*[vol of shell]) (%)
##COLUMN 10|pocketness_clus 2 |sum of 1/[Rpocket_for_pocketcluster 2] /(1/[Rmin]*[vol of shell]) (%)
##COLUMN 11|pocketness_clus 3 |sum of 1/[Rpocket_for_pocketcluster 3] /(1/[Rmin]*[vol of shell]) (%)
##COLUMN 12|pocketness_clus 4 |sum of 1/[Rpocket_for_pocketcluster 4] /(1/[Rmin]*[vol of shell]) (%)
##COLUMN 13|pocketness_clus 5 |sum of 1/[Rpocket_for_pocketcluster 5] /(1/[Rmin]*[vol of shell]) (%)
#   2  A ILE  77.59  4.818  8  0   8.85   5.08   0.00   0.00   0.00   0.00
#   3  A SER  55.38  2.584  6  0   4.12   2.17   0.00   0.00   0.00   0.00
#   4  A PRO  73.31  5.473  7  0   1.43   0.92   0.00   0.00   0.00   0.00
#   5  A ILE  57.08  2.662  8  0   4.62   4.00   0.00   0.00   0.00   0.00
#:
# 554  A ALA  74.79  6.153  5  0   0.40   0.00   0.00   0.00   0.00   0.00
# 555  A GLY  70.89  5.673  4  0   0.23   0.00   0.00   0.00   0.00   0.00
# 556  A ILE  73.38  5.889  9  0   0.95   0.00   0.00   0.00   0.00   0.00
#   2  B ILE  91.79  8.544  8  0   3.19   3.19   0.00   0.00   0.00   0.00
#   3  B SER  73.38  4.938  6  0   6.91   6.91   0.00   0.00   0.00   0.00
#   4  B PRO  76.58  6.020  7  0   3.23   3.23   0.00   0.00   0.00   0.00
#   5  B ILE  80.33  6.758  8  0   1.84   1.84   0.00   0.00   0.00   0.00
#   6  B GLU  73.02  5.093  9  0   5.12   5.12   0.00   0.00   0.00   0.00
#:
# 426  B TRP  72.99  2.386 14  0  26.65  26.65   0.00   0.00   0.00   0.00
# 427  B TYR  77.60  4.309 12  0  19.58  19.58   0.00   0.00   0.00   0.00
# 428  B GLN  80.82  5.999 10  0  10.67  10.67   0.00   0.00   0.00   0.00



import sys
import os 

LastModDate = 'Sep 9, 2008'
   
class Pocket:
  
  def __init__(self):
      self.filename = ''   
      self.rnumchlist    = []
      self.rnum        = {} 
      self.chain       = {} 
      self.res         = {} 
      self.shellAcc    = {} 
      self.Rinacc      = {} 
      self.Natom       = {} 
      self.Natom_contact  = {} 
      self.pocketness  = {} 
      self.pocketness_clus1  = {} 
      self.pocketness_clus2  = {} 
      self.pocketness_clus3  = {} 
      self.pocketness_clus4  = {} 
      self.pocketness_clus5  = {} 
    
  def read(self,fname):
      print "#respock.Pocket.read(\"%s\") "%(fname)
      if not os.access(fname,os.R_OK):
        print "#ERROR:Can't open res-based pocketness file '%s'" % fname
        sys.exit(1)

      self.filename = fname
      f = open(fname)
      start_align = 0 
      for line in f: 
        if line:
          line = line.rstrip('\n')
        if (line.startswith("#")==0):
#265  B ASN  56.50  1.697  8  0  14.59   0.92   0.00   0.00  13.66   0.00
##COLUMN  1|RNUM              |Residue Number
##COLUMN  2|CHAIN             |Chain Identifier
##COLUMN  3|RES               |Three-letter residue name
##COLUMN  4|shellAcc          |shell accessibility (%)
##COLUMN  5|Rinacc            |averaged Rinaccess (A)
##COLUMN  6|Natom             |Number of atoms
##COLUMN  7|Natom_contact     |Number of contacting atoms with ligand
##COLUMN  8|pocketness        |sum of 1/[Rpocket] /(1/[Rmin]*[vol of shell]) (%)
##COLUMN  9|pocketness_clus 1 |sum of 1/[Rpocket_for_pocketcluster 1] /(1/[Rmin]*[vol of shell]) (%)
##COLUMN 10|pocketness_clus 2 |sum of 1/[Rpocket_for_pocketcluster 2] /(1/[Rmin]*[vol of shell]) (%)
##COLUMN 11|pocketness_clus 3 |sum of 1/[Rpocket_for_pocketcluster 3] /(1/[Rmin]*[vol of shell]) (%)
##COLUMN 12|pocketness_clus 4 |sum of 1/[Rpocket_for_pocketcluster 4] /(1/[Rmin]*[vol of shell]) (%)
##COLUMN 13|pocketness_clus 5 |sum of 1/[Rpocket_for_pocketcluster 5] /(1/[Rmin]*[vol of shell]) (%)
          field = line.split()
          rnum = field[0]
          chain = field[1]
          rnumch = rnum + ' ' + chain
          self.rnumchlist.append(rnumch)
          self.rnum[rnumch]          = field[0]
          self.chain[rnumch]         = field[1]
          self.res[rnumch]           = field[2]
          self.shellAcc[rnumch]      = float(field[3])
          self.Rinacc[rnumch]        = float(field[4])
          self.Natom[rnumch]         = int(field[5])
          self.Natom_contact[rnumch] = int(field[6])
          self.pocketness[rnumch]    = float(field[7])
          self.pocketness_clus1[rnumch]    = float(field[8])
          self.pocketness_clus2[rnumch]    = float(field[9])
          self.pocketness_clus3[rnumch]    = float(field[10])
          self.pocketness_clus4[rnumch]    = float(field[11])
          self.pocketness_clus5[rnumch]    = float(field[12])
          #print rnum,line
      f.close

  def __str__(self):
      s = "#res_base_pock.Pocket:'%s' "%(self.filename)
      return s  
  
 
def _main():
    if (len(sys.argv)<2):
      print "#ERROR:Insufficient arguments"
      sys.exit()
    P = Pocket()   
    P.read(sys.argv[1]) 
    print P 
    for rnum in (P.rnumchlist):
      print "'%s' shellAcc %s pocket %s"%(rnum, P.shellAcc[rnum],P.pocketness[rnum])
if __name__ == '__main__':_main()
