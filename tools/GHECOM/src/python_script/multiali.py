##
## <multiali.py>
##  dealing with pairwise alignment
##  
##>> FILE FORMAT EXAMPLE <<
#NPRO 4
#PRO1 4hhbA
#PRO2 4hhbB
#PRO3 1mbdA
#PRO4 1ecdA
#COMMENT DATE Sep 19, 2008 TIME 14:46:48
#ALIGNMENT
#-   -1 V    1 -   -1 -   -1
#V    1 H    2 V    1 -   -1
#L    2 L    3 L    2 L    1
#S    3 T    4 S    3 S    2
#P    4 P    5 E    4 A    3
#A    5 E    6 G    5 D    4
#D    6 E    7 E    6 Q    5
#K    7 K    8 W    7 I    6
#T    8 S    9 Q    8 S    7
#N    9 A   10 L    9 T    8
#V   10 V   11 V   10 V    9
#K   11 T   12 L   11 Q   10
#A   12 A   13 H   12 A   11
#A   13 L   14 V   13 S   12
#W   14 W   15 W   14 F   13
#G   15 G   16 A   15 D   14
#K   16 K   17 K   16 K   15
#V   17 V   18 V   17 V   16
#-   -1 -   -1 E   18 K   17
#-   -1 -   -1 A   19 G   18
#G   18 -   -1 -   -1 -   -1
#A   19 -   -1 -   -1 -   -1
#H   20 N   19 D   20 -   -1
#A   21 V   20 V   21 -   -1
#:
#:
#-   -1 -   -1 Q  152 -   -1
#-   -1 -   -1 G  153 -   -1
#END


import sys
import os 
import math 

LastModDate = 'Sep 19, 2008'
   
class MultiAlign:
  
  def __init__(self):
      self.filename = ''   
      self.Nprotein = 0
      self.proname = []  # list of protein name   
      self.Nalign  = 0
      self.ali_resnum = [] # residue_number[nprotein][nalign] 
      self.seq = []
      self.comment   = ''   
    
  def read_vertical(self,fname):
      print "#multiali.MulgiAlign.read_vertical_style(\"%s\") "%(fname)
      f = open(fname)
      start_align = 0 
      for line in f: 
        if line:
          line = line.rstrip('\n')
        #print line
        field = line.split()
        if (line.startswith("NPRO")):
          print field 
          self.Nprotein = int(field[1])
          print "#Nprotein %d"%(self.Nprotein)
          self.proname    = [' ' for i in range(self.Nprotein)]
          self.ali_resnum = [[]  for i in range(self.Nprotein)]
          self.seq        = [{}  for i in range(self.Nprotein)]
          pass 
        if (line.startswith("PRO")):
           npro = int(field[0][3:])
           print "'%s'-->'%d'"%(line,npro)
           self.proname[npro-1] = field[1]
           #print "npro %d"%(npro)
        if (line.startswith("END")):
          start_align = 0

        if (start_align==1):
          for i in range(self.Nprotein):
            rnum = field[2*i+1]
            self.ali_resnum[i].append(rnum)
            if (rnum != "-1"):
              self.seq[i][rnum] = field[2*i]

        if (line.startswith("ALIGNMENT")):
           start_align = 1 
      f.close
      self.Nalign = len(self.ali_resnum[0])

  def __str__(self):
      s = '' 
      s = "#MultiAlign  Nprotein %d Nalign %d"%(self.Nprotein,self.Nalign)
      return s  
  
 
def _main():
    if (len(sys.argv)<2):
      print "#ERROR:Insufficient arguments"
      sys.exit()
    M = MultiAlign()   
    M.read_vertical(sys.argv[1]) 
    print M
    for a in range(M.Nalign):
      sys.stdout.write("%3d"%(a))
      for p in range(M.Nprotein):
        rnum = M.ali_resnum[p][a]
        if (rnum=="-1"):
          aa = '-'
        else:
          aa = M.seq[p][rnum] 
        sys.stdout.write(" %s %3s"%(aa, M.ali_resnum[p][a]))
      sys.stdout.write("\n")

if __name__ == '__main__':_main()
