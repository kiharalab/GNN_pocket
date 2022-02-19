#!/usr/bin/python

import os
import sys


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
     if (len(field)>2):
       prop[field[0]] = field[1] 
 print "#Nlist %d"%(len(list)) 
 f.close()

##############
#### MAIN ####
##############

OPT = {}
OPT['of'] = ''
OPT['od'] = 'ChConPDB'
OPT['A']  = 'F'

if (len(sys.argv)<2):
  print "mkExChPDB.py [listfile] <options>"
  print " to extracting specfied chains and contacting ligands"
  print "<options>"
  print " -od     : output directory for selected chains [%s]"%(OPT['od'])
  print " -of     : output file [%s]"%(OPT['of'])
  print " -A      : Action ('T' or 'F') [%s]"%(OPT['A'])
  sys.exit(1)

read_option(sys.argv,OPT)
ilistfile = sys.argv[1]
LIST = []
PROP = {}
read_list_file(ilistfile,LIST,PROP)
for X in (LIST):
  pdbid  = X[0:4]
  chains = X[4:]
  print "%s %s"%(pdbid,chains)
  command = "ExChPDB.py %s -ch %s -of %s/%s"%(pdbid,chains,OPT['od'],X)
  if (OPT['A'] == 'T'):
    os.system(command)
  else:
    print "#%s"%(command)

