##
## <pdb.py>
##
## 

import sys
import os 
import math 

LastModDate = 'Oct 1, 2009'
   
class Molecule:
  
  def __init__(self):
      self.filename = ''   
      self.header = ''
      self.title = ''
      self.Natom  = 0
      self.Nheavy = 0
      self.Nres  = 0
      self.line   = [] # [0..Natom-1] 'entire line for pdb'
      self.anum  = []  # [0..Natom-1]
      self.atom  = []  # [0..Natom-1]
      self.resi  = []  # [0..Natom-1]
      self.rnum  = []  # [0..Natom-1]
      self.posX  = []  # [0..Natom-1]
      self.posY  = []  # [0..Natom-1]
      self.posZ  = []  # [0..Natom-1]
      self.chain = []  # [0..Natom-1]
      self.AHtype = [] # [0..Natom-1]
      self.res_num = [] # [0..Natom-1] residue number 0... Nres
      self.anum2int = {} 
      self.cAt1 = [] 
      self.cAt2 = [] 
      self.cAt3 = [] 
      self.conservation = []
 
  def read(self,filename,AHtype="B",Chain="-"):
  
#REMARK  [Ratom](10: 55 - 60):Radius of atom[A]
#REMARK  [tFact](11: 61 - 66):[Riacc] is assigned as tFactor
#REMARK  [shAcc](12: 67 - 72):Shell accessibility. (Ratio of noVdW grids in the shell)[%]
#REMARK  [Riacc](13: 73 - 78):averaged Rinaccess in the 'shell'[A]
#REMARK  [co]   (14: 79 - 81):Contact with other molecules (0 or 1)
#REMARK  [pocke](15: 82 - 87):pocketness. sum of 1/[Rpocket] /(1/[Rmin]*[vol of shell]) (%)
#REMARK  [pock1](16: 88 - 93):pocketness_clus[1] sum of 1/[Rpocket_for_pocke_tcluster1] /(1/[Rmin]*[vol of shell]) (%)
#REMARK  [pock2](17: 94 - 99):pocketness_clus[2] sum of 1/[Rpocket_for_pocke_tcluster2] /(1/[Rmin]*[vol of shell]) (%)
#REMARK  [pock3](18:100 -105):pocketness_clus[3] sum of 1/[Rpocket_for_pocke_tcluster3] /(1/[Rmin]*[vol of shell]) (%)
#REMARK  [pock4](19:106 -111):pocketness_clus[4] sum of 1/[Rpocket_for_pocke_tcluster4] /(1/[Rmin]*[vol of shell]) (%)
#REMARK  [pock5](20:112 -117):pocketness_clus[5] sum of 1/[Rpocket_for_pocke_tcluster5] /(1/[Rmin]*[vol of shell]) (%)
#REMARK  [cAt1] (21:118 -122):Atom number of contact protein atom 1
#REMARK  [cAt2] (22:123 -127):Atom number of contact protein atom 2
#REMARK  [cAt3] (23:128 -133):Atom number of contact protein atom 3
#REMARK  [cons] (24:133 -139):Conservation score
#REMARK                                                 Ratom|cons |shAcc|Riacc|co|pocke|pock1|pock2|pock3|pock4|pock5|cAt1|cAt2|cAt3|cons|
#HETATM    1  CA  PRB 1    1     83.668  42.417  49.795 1.87  100.0  52.6  2.74  0 66.67 66.67  0.00  0.00  0.00  0.00   68   81   85 100.0
#HETATM    2  CA  PRB 1    1     83.446  42.166  50.873 1.87  100.0  51.5  2.63  0 64.10 64.10  0.00  0.00  0.00  0.00   68  478  481 100.0
#HETATM    3  CA  PRB 1    1     83.728  42.891  51.646 1.87  100.0  49.0  2.59  0 65.38 65.38  0.00  0.00  0.00  0.00   68  709  710 100.0
#HETATM    4  CA  PRB 1    1     83.902  43.410  50.549 1.87  100.0  50.9  2.75  0 66.67 66.67  0.00  0.00  0.00  0.00   68  710  728 100.0
#HETATM    5  CA  PRB 1    1     81.754  42.489  47.401 1.87  100.0  68.1  3.15  0 57.28 57.28  0.00  0.00  0.00  0.00   81   84   89 100.0
#HETATM    6  CA  PRB 1    1     82.697  42.807  48.663 1.87  100.0  61.4  2.95  0 63.17 63.17  0.00  0.00  0.00  0.00   81   84  728 100.0
#HETATM    7  CA  PRB 1    1     81.229  41.916  47.862 1.87  100.0  67.8  3.05  0 57.49 57.49  0.00  0.00  0.00  0.00   81   89  481 100.0


 
  
      if not os.access(filename,os.R_OK):
        print "#ERROR:Can't open '%s'" % filename
        sys.exit()  
     
      f = open(filename)
      self.filename = filename
  
      #print f 
      resi0  = ''
      rnum0  = ''
      chain0 = ''
      self.Natom = 0 
      self.Nres  = 0 
      for line in f: 
        line = line.rstrip('\n')
    
    #          1         2         3         4         5         6         7
    #01234567890123456789012345678901234567890123456789012345678901234567890123456789
    #ATOM    676  CB  GLU A  85      10.440  29.552  12.788  6.00 16.96           C
    #ATOM    680  OE2 GLU A  85      10.230  30.451  16.374  8.00 41.03           O
    #ATOM    682  CA  LEU A  86       7.618  29.487   9.238  6.00 12.23           C
    #HETATM 1236  O4  SO4   154      33.810  28.815  -4.624  8.00 14.90           O
    #HETATM 1237 FE   HEM   155      15.271  27.962   0.622 24.00  7.86          FE
     
        if line:
          if (line[0:6]=='HEADER'):
            self.header = line
          if (line[0:6]=='TITLE '):
            self.title = line
          read_it = 0
          if (line[0:6]=='ATOM  ') and (AHtype=='A'):
            read_it = 1
          if (line[0:6]=='HETATM') and (AHtype=='H'):
            read_it = 1
          #if (((line[0:6]=='ATOM  ') or (line[0:6]=='HETATM')) and (AHtype=='B')):
          #  read_it = 1
          if (line[0:6]=='ATOM  ') and (AHtype=='B'):
            read_it = 1
          if (line[0:6]=='HETATM') and (AHtype=='B'):
            read_it = 1
          if (read_it==1):
            chain = line[21:22] 
            if (Chain=="-") or (chain == Chain):
              anum  = line[6:11] 
              atom  = line[12:16] 
              resi  = line[17:20] 
              chain = line[21:22] 
              if (chain == ' '):
                chain = '-'
              rnum  = line[22:27] 
              x     = line[31:38] 
              y     = line[38:46] 
              z     = line[46:54] 
              if (line[0:6]=='ATOM  '):
                self.AHtype.append('A')
              if (line[0:6]=='HETATM'):
                self.AHtype.append('H')
            #print "'%s' '%s' '%s' '%s' '%s' %s %s %s\n" %(anum,atom,resi,chain,rnum,x,y,z)
              if (resi != resi0) or (rnum != rnum0) or (chain != chain0):
                self.Nres += 1 
              self.line.append(line)
              self.anum.append(anum)
              self.anum2int[int(anum)] = self.Natom
              self.atom.append(atom)
              self.resi.append(resi)
              self.rnum.append(rnum)
              self.chain.append(chain)
              self.posX.append(float(x))
              self.posY.append(float(y))
              self.posZ.append(float(z))

              if (len(line)>=132):
                cAt1 = int(line[117:122])
                cAt2 = int(line[122:127])
                cAt3 = int(line[127:133])
              else:
                cAt1 = cAt2 = cAt3 = 0
              self.cAt1.append(cAt1)
              self.cAt2.append(cAt2)
              self.cAt3.append(cAt3)

              if (len(line)>=138):
                cons = float(line[132:139]) 
              else:
                cons = 0.0
              self.conservation.append(cons)

              self.res_num.append(self.Nres)
              n = self.Natom
           

              #sys.stdout.write("[%d] %s %s %s %d %d %d\n"%(n,self.atom[n],self.resi[n],self.rnum[n],self.cAt1[n],self.cAt2[n],self.cAt3[n]))
              resi0 = resi
              rnum0 = rnum 
              chain0 = chain 
              if (atom[1] != 'H'):
                self.Nheavy += 1
              self.Natom += 1
      f.close
      self.Natom = len(self.anum)
      print "#pdbmol.read(%s AHtype %s) --> Natom %d"%(filename,AHtype,self.Natom)

  def __str__(self):
      s = '' 
      s = "#pdb.Molecule() Natom %d"%(self.Natom)
      return s  
  
 
def _main():
    print "woops!!"
    if (len(sys.argv)<2):
      print "#ERROR:Insufficient arguments"
      sys.exit()
    p = Molecule()
    print p
    p.read(sys.argv[1])
    print p
if __name__ == '__main__':_main()
