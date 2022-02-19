#!/usr/bin/python

import sys

f = open(sys.argv[1])
N = {}

for line in f:
  line = line.rstrip('\n')
#  print line 
  if ((len(line)>0) and (line[0]!='#')):
    field = line.split() 
    a = int(field[0]) 
    b = int(field[1]) 
    c = int(field[2]) 
    if (N.has_key(a)==0):
      N[a] = {}
    N[a][b] = c
f.close()

for k in range (1,8):
  Nall = 0  
  for i in range(0,18):
    Nall += N[k][i]
  sys.stdout.write("%d\t"%(k))
  sys.stdout.write("%f\t"%((N[k][1]+N[k][2]+N[k][3])/float(Nall)))
  sys.stdout.write("%f\t"%((N[k][4]+N[k][5])/float(Nall)))
  sys.stdout.write("%f\t"%((N[k][6]+N[k][7])/float(Nall)))
  sys.stdout.write("%f\t"%((N[k][8]+N[k][9])/float(Nall)))
  sys.stdout.write("%f\t"%((N[k][10]+N[k][11])/float(Nall)))
  sys.stdout.write("%f\t"%((N[k][12]+N[k][13])/float(Nall)))
  sys.stdout.write("%f\t"%((N[k][14]+N[k][15]+N[k][16]+N[k][17]+N[k][0])/float(Nall)))
  sys.stdout.write("\n")


