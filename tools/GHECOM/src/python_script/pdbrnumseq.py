##
## <pdbrnumseq.py>
##  
##>> FILE FORMAT EXAMPLE <<
#
#>1YDRI Exp X Res 2.2 Naa   20 Mch 100 Unk   0 Mut - Tax - Nreg 0
##[num] [AA] [RnumPDB] [ChainID] [TriAA]
#1    T    5  I THR
#2    T    6  I THR
#3    Y    7  I TYR
#4    A    8  I ALA
#5    D    9  I ASP
#6    F   10  I PHE
#7    I   11  I ILE
#8    A   12  I ALA
#9    S   13  I SER
#10   G   14  I GLY
#11   R   15  I ARG
#12   T   16  I THR
#13   G   17  I GLY
#14   R   18  I ARG
#15   R   19  I ARG
#16   N   20  I ASN
#17   A   21  I ALA
#18   I   22  I ILE
#19   H   23  I HIS
#20   D   24  I ASP
#//



import sys
import os 

LastModDate = 'Sep 16, 2009'
   
class PdbRnumSeq:
  
  def __init__(self):
      self.filename = ''
      self.length   = 0
      self.rnum     = [0]
      self.seq      = ['-']   ## [1..length]
      self.rnumpdb  = ['-']   ## [1..length]
      self.chain    = ['-']   ## [1..length]
      self.triseq   = ['---'] ## [1..length]
      self.rnumP2S  = {}
 
  def read(self,fname):
      print "#pdbrnumseq.PdbRnumSeq.read(\"%s\") "%(fname)
      if not os.access(fname,os.R_OK):
        print "#ERROR:Can't open PdbRnumSeq '%s'" % fname
        sys.exit(1)
      f = open(fname)
      nsite = 0
      for line in f:
        if ((len(line)>5) and (line.startswith('>')==0) and (line.startswith('#')==0) and (line.startswith('/')==0)):
          nsite += 1
          field = line.split()
          self.rnum.append(int(field[0]))
          self.seq.append(field[1])
          self.rnumpdb.append(field[2])
          self.chain.append(field[3])
          self.triseq.append(field[4])
          self.rnumP2S[field[2]] = int(field[0])
      self.length = nsite
      f.close()

  def __str__(self):
      s = "#PdbRnumSeq:'%s' length %d"%(self.filename,self.length)
      return s


def _main():
    if (len(sys.argv)<2):
      print "#ERROR:Insufficient arguments"
      sys.exit()
    P = PdbRnumSeq()
    P.read(sys.argv[1])
    print P
if __name__ == '__main__':_main()


