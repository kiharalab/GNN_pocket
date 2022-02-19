##
## <pairali.py>
##  dealing with pairwise alignment
##  
##>> FILE FORMAT EXAMPLE <<
#NPRO 2
#PRO1 1mbdA.bssp
#PRO2 4hhbA.bssp
#COMMENT Naa1  153 Naa2  141
#COMMENT Ncomp 141 SqID 27.0 RMS 1.560 DRMS 1.398
#COMMENT ScDis 149312.6 Rdis 70.7
#PARAM_FOR_SUPERIMPOSING
##Afit=R*(A-Ga)+Gb
#Ga 15.52624 20.25461  5.17376
#Gb 13.49291  8.08652 -8.30567
#R0  0.82302 -0.33401  0.45944
#R1 -0.30589  0.42092  0.85397
#R2 -0.47862 -0.84336  0.24426
#ALIGNMENT
#V 1 V 1
#L 2 L 2
#S 3 S 3
#E 4 P 4
#G 5 A 5
#E 6 D 6
#:
#END

import sys
import os 
import math 

LastModDate = 'Sep 9, 2008'
   
class Alignment:
  
  def __init__(self):
      self.filename = ''   
      self.proname1 = ''   
      self.proname2 = '' 
      self.ali_resnum1 = []
      self.ali_resnum2 = []
      self.seq1 = {} 
      self.seq2 = {} 
      self.Nalign  = 0
      self.header   = ''   
    
  def read_vertical(self,fname):
      print "#pairali.Alignment.read_vertical_style(\"%s\") "%(fname)
      f = open(fname)
      start_align = 0 
      for line in f: 
        if line:
          line = line.rstrip('\n')
        field = line.split()
        if (line.startswith("PRO1 ")):
          self.proname1 = field[1]
        if (line.startswith("PRO2 ")):
          self.proname2 = field[1]
        if (line.startswith("END")):
          start_align = 0
        if (start_align==1):
          self.ali_resnum1.append(field[1])
          self.ali_resnum2.append(field[3])
          if (field[1] != "-1"):
            self.seq1[field[1]] = field[0] 
          if (field[3] != "-1"):
            self.seq2[field[3]] = field[2] 

          #print line 
          pass
        if (line.startswith("ALIGNMENT")):
          start_align = 1 
      f.close
      self.Nalign = len(self.ali_resnum1)

  def __str__(self):
      s = "#Alignment: pro1 '%s' pro2 '%s' Nalign %d"%(self.proname1,self.proname2,self.Nalign)
      return s  
  
 
def _main():
    if (len(sys.argv)<2):
      print "#ERROR:Insufficient arguments"
      sys.exit()
    A = Alignment()   
    A.read_vertical(sys.argv[1]) 
    print A
    for a in range(A.Nalign):
      rnum1 = A.ali_resnum1[a]
      rnum2 = A.ali_resnum2[a]
      aa1 = '-'
      aa2 = '-'
      if (rnum1 != "-1"):
        aa1 = A.seq1[rnum1]
      if (rnum2 != "-1"):
        aa2 = A.seq1[rnum1]
      print "%d %s %s %s %s"%(a, aa1,rnum1,aa2,rnum2)
if __name__ == '__main__':_main()
