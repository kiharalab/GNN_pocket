#!/usr/bin/python

import sys
import os
import math

LastModDate = "Sep 25, 2008"

def read_option(argv,opt_dic):
  for i in range(1,len(argv)):
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]


def read_protein_list(fname,prolist):
  f = open(fname)
  for line in f:
    if line:
      line = line.rstrip('\n')
    if (line.startswith("#")==0):
      print line
      prolist.append(line)
      print prolist
  f.close()

def action(command,action):
  print "#%s : %s"%(command,action)
  if (action == 'T'):
    os.system(command)

OPT={}
OPT['A']   = 'F'
OPT['T']   = 'P'
OPT['odc'] = 'ChConPDB'
OPT['odp'] = 'POCKET'
OPT['oda'] = 'ACC'
OPT['odl'] = 'LIGAND'

if (len(sys.argv)<2):
  print "doghecom.py  [prolistfile] <options>"
  print " do ghecom/LeeRich/Ligand for protein lists"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate) 
  print "<options>"
  print " -T   : type 'C'hain-extracted pdb, 'P'ocket 'A'cc 'L'igand [%s]"%(OPT['T']) 
  print " -odc : output dir for chain splited PDB file [%s]"%(OPT['odc']) 
  print " -odp : output dir for pocket [%s]"%(OPT['odp']) 
  print " -oda : output dir for acc    [%s]"%(OPT['oda']) 
  print " -odl : output dir for ligand [%s]"%(OPT['odl']) 
  print " -A : Action (T or F)[%s]"%(OPT['A']) 
  sys.exit(1)

read_option(sys.argv,OPT)
prolistfile = sys.argv[1]
prolist = []
read_protein_list(prolistfile,prolist)
print prolist
for pro in (prolist):
  chain  = pro[-1] 
  pdb = pro[0:4] 
  print ">%s %s %s"%(pro,pdb,chain)
  if (OPT['T'] == 'C'):
     str = "Ligand /DB/PDBv3/%s/pdb%s.ent -M - -ch %s -tz T -mos T -cos T -cbs T -cdrpp T -n2h T -omp %s/%s"%(pdb[1:3],pdb,chain,OPT['odc'],pro);
     action(str,OPT['A'])

  if (OPT['T'] == 'P'):
    str = "ghecom %s/%s -ch %s -M M -ores %s/%s"%(OPT['odc'],pro,chain,OPT['odp'],pro)
    action(str,OPT['A'])

  if (OPT['T'] == 'A'):
    str = "LeeRich %s/%s -ch %s -or %s/%s"%(OPT['odc'],pro,chain,OPT['oda'],pro)
    action(str,OPT['A'])

  if (OPT['T'] == 'L'):
    str = "Ligand %s/%s -ch %s -M P > %s/%s"%(OPT['odc'],pro,chain,OPT['odl'],pro)
    action(str,OPT['A'])
