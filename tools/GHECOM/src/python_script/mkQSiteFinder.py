#!/usr/bin/python

import os
import sys

LastModDate = "Aug 5, 2009"

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
OPT['ipdb'] = 'ChConPDB'
OPT['iqsf'] = 'QSiteOUTPUT'
OPT['og']   = 'QSiteGRID'
OPT['G']  = 'F'
OPT['A']  = 'F'
OPT['gw'] = 0.8
OPT['rk'] = 'A';

if (len(sys.argv)<2):
  print "mkQSiteFinder.py [listfile] <options>"
  print " for changing QSiteFinder output PDB file to grids. "
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -ipdb  : Input  dir. for selected PDB chains [%s]"%(OPT['ipdb'])
  print " -iqsf  : Input  dir. for output PDB file for QSiteFinder [%s]"%(OPT['iqsf'])
  print " -og    : Output dir. for grid PDB file [%s]"%(OPT['og'])
  print " -gw    : grid width [%s]"%(OPT['gw'])
  print " -rk    : Cluster rank 'A'll,only '1'st, '2'nd, '3'rd [%s]"%(OPT['rk'])
  print " -G     : Do gridcomp  ('T' or 'F') [%s]"%(OPT['G'])
  print " -A     : Action ('T' or 'F') [%s]"%(OPT['A'])
  sys.exit(1)

read_option(sys.argv,OPT)
ilistfile = sys.argv[1]
LIST = []
PROP = {}
read_list_file(ilistfile,LIST,PROP)
for pdbch in (LIST):
  command = "gridcomp -M G -ipA %s/%s -gw %s -rA 1.7 -aA A -og %s/%s"%(OPT['iqsf'],pdbch,OPT['gw'],OPT['og'],pdbch)
  if (OPT['rk']=='A'):
    command = command + " -resA Yxx "
  elif (OPT['rk']=='1'):
    command = command + " -resA 'YAA' "
  elif (OPT['rk']=='2'):
    command = command + " -resA 'YAA|YAB' "
  elif (OPT['rk']=='3'):
    command = command + " -resA 'YAA|YAB|YAC' "
 
  if (OPT['G'] == 'T'):
    os.system(command)
  else:
    print "#%s"%(command)

