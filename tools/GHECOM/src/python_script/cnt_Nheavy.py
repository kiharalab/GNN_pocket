#!/usr/bin/python


import os
import sys
import pdb

LastModDate = "Oct 1, 2009"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]


def read_list_file(ifname,list,prop):
 ## [FILE_EXAMPLE]
 ## 1dlwA  a.1.1.1
 ## 1a6m-  a.1.1.2
 ## 1eca-  a.1.1.2
 ## 1lh1-  a.1.1.2
 ## 1ash-  a.1.1.2
 if not os.access(ifname,os.R_OK):
   print "#ERROR:Can't open listfile '%s'" % ifname
   return(1)

# list = []
# prop = {}
 f = open(ifname)
 for line in f:
   line = line.rstrip('\n')
   if (line.startswith('#')==0) and (len(line)>2):
     field = line.split()
     list.append(field[0])
     prop[field[0]] = line
 print "#Nlist %d"%(len(list)) 
 f.close()

##############
#### MAIN ####
##############

OPT = {}
OPT['ipdb'] = 'ChConPDB'
OPT['A']  = 'F'
OPT['of'] = '-'
OPT['ol'] = '-'
OPT['nh'] = 10000

if (len(sys.argv)<2):
  print "cnt_Nheavy.py [listfile] <options>"
  print " counting Number of heavy atoms for the protein chain list"
  print "<options>"
  print " -ipdb : Input  dir. for selected PDB chains [%s]"%(OPT['ipdb'])
  print " -of   : output list of Nheavy atoms  [%s]"%(OPT['of'])
  print " -ol   : output list with less Nheavy atoms  [%s]"%(OPT['ol'])
  print " -nh   : maximum Nheavy atoms  [%d]"%(OPT['nh'])
  sys.exit(1)

read_option(sys.argv,OPT)
ilistfile = sys.argv[1]
LIST = []
PROP = {}
read_list_file(ilistfile,LIST,PROP)

if (OPT['ol'] == '-'):
  fol = sys.stdout
else:
  fol = open(OPT['ol'],'w')

if (OPT['of'] == '-'):
  fof = sys.stdout
else:
  fof = open(OPT['of'],'w')

for i in range(len(LIST)):
  pdbch = LIST[i] 
  chain = pdbch[4:]
  P = pdb.Molecule()
  P.read(OPT['ipdb']+'/'+pdbch,AHtype='A',Chain=chain)
  fof.write("%s Natom %d Nheavy %d\n"%(pdbch,P.Natom,P.Nheavy))
  if (P.Nheavy < int(OPT['nh'])):
    fol.write("%s\n"%(PROP[pdbch])) 
  else:
    fof.write("#%s is not chosen(Nheavy %d)\n"%(pdbch,P.Nheavy)) 


if (OPT['ol'] != '-'):
  fol.close()

if (OPT['of'] != '-'):
  fof.close()
