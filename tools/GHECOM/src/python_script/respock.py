##
## <respock.py>
##  dealing with residue-based pocket file 
##  
##>> FILE FORMAT EXAMPLE <<
#>> FILE EXAMPLE <<
# #COLUMN 1|RNUM      |Residue Number
# #COLUMN 2|CHAIN     |Chain Identifier
# #COLUMN 3|RES       |Three-letter residue name
# #COLUMN 4|shellAcc  |shell accessibility (%)
# #COLUMN 5|invRinacc |Average of inverse of Rinacc (1/A)
# #COLUMN 6|pocketness|[invRinacc] */(1/[RprobeS]) (%)
# #COLUMN 7|Rinvacc   |inverse of [invRinacc] (A)
#    1  A ALA  39.22 0.011712   2.19 85.380
#    2  A VAL  52.52 0.051019   9.54 19.601
#    3  A ASN  77.56 0.041348   7.73 24.185
#    4  A GLY  55.22 0.055561  10.39 17.998
#    5  A LYS  65.11 0.041982   7.85 23.820
#    6  A GLY  54.51 0.042088   7.87 23.760
#    7  A MET  44.75 0.013610   2.55 73.473
#    8  A ASN  52.22 0.030910   5.78 32.351
# :
# 1206  A ARG  85.19 0.040583   7.59 24.641
# 1207  A ASN  87.74 0.033620   6.29 29.744
# 1208  A ASP  90.74 0.021993   4.11 45.470
##COLUMN 1|RNUM        |Residue Number
##COLUMN 2|CHAIN       |Chain Identifier
##COLUMN 3|RES         |Three-letter residue name
##COLUMN 4|shellAcc    |shell accessibility (%)
##COLUMN 5|Rinacc      |averaged Rinaccess (A)
##COLUMN 6|rw_shellAcc |Rinacc-weighted shell accessibility (%)
##COLUMN 7|pocketness  |sum of 1/[Rpocket] /(1/[Rmin]*[vol of shell]) (%)
#   2  A THR  87.29  9.145  87.10   0.38
#   3  A PHE  88.18  9.246  88.06   0.20
#   4  A ASN  81.13  8.293  78.98   2.68
#   :
#1068  A MET  39.04  4.010  38.19   0.89
#1069  A HIS  51.79  5.056  48.15   3.54
#1070  A ALA  65.93  6.922  65.93   0.00
#1071  A GLN  63.68  6.500  61.90   1.89
#1072  A ILE  53.75  5.516  52.53   1.20
#1073  A LYS  88.37  9.112  86.78   1.60

import sys
import os 

LastModDate = 'Sep 9, 2008'
   
class Pocket:
  
  def __init__(self):
      self.filename = ''   
      self.rnumlist    = []
      self.chain       = {} 
      self.res         = {} 
      self.shellAcc    = {} 
      self.Rinacc      = {} 
      self.rw_shellAcc = {} 
      self.pocketness  = {} 
    
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
          field = line.split()
          rnum = field[0] + ' ' + field[1]
          self.rnumlist.append(rnum)
          self.chain[rnum]        = field[1]
          self.res[rnum]          = field[2]
          self.shellAcc[rnum]     = field[3]
          self.Rinacc[rnum]       = field[4]
          self.rw_shellAcc[rnum]  = field[5]
          self.pocketness[rnum]   = field[6]
          #print rnum,line
      f.close

  def __str__(self):
      s = "#Pocket:'%s' "%(self.filename)
      return s  
  
 
def _main():
    if (len(sys.argv)<2):
      print "#ERROR:Insufficient arguments"
      sys.exit()
    P = Pocket()   
    P.read(sys.argv[1]) 
    print P 
    for rnum in (P.rnumlist):
      print "'%s' shellAcc %s pocket %s"%(rnum, P.shellAcc[rnum],P.pocketness[rnum])
if __name__ == '__main__':_main()
