#!/usr/bin/python

import sys
import os
import math
import glob
from datetime import datetime

import res_base_pock
import psiprof
import pdbrnumseq
import pdb

LastModDate = "Sep 30, 2009"


def read_option(argv,opt_dic):
  now = datetime.now()
  opt_dic['DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
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
      #print line
      prolist.append(line)
      #print prolist
  f.close()

def action(command,action):
  print "#%s : %s"%(command,action)
  if (action == 'T'):
    os.system(command)


def cal_lig_prb_overlap(lig,prb,Dthre):
  lig.mark = [0 for i in range(lig.Natom)]
  prb.mark = [0 for i in range(prb.Natom)]
  DDthre = Dthre * Dthre
  Nprb_select = 0
  
  for i in range(prb.Natom):
    if (prb.select[i]==1):
      for j in range(lig.Natom):
        dx = prb.posX[i] - lig.posX[j]
        dy = prb.posY[i] - lig.posY[j]
        dz = prb.posZ[i] - lig.posZ[j]
        DD = dx * dx + dy * dy + dz * dz
        if (DD<DDthre):
          prb.mark[i] = 1 
          lig.mark[j] = 1 
      Nprb_select += 1

  Nlig_over = Nprb_over = 0 

  for i in range(prb.Natom):
    Nprb_over += prb.mark[i]
  for j in range(lig.Natom):
    Nlig_over += lig.mark[j]

  return([Nlig_over,Nprb_over,lig.Natom,Nprb_select])


def trim_top_cluster(prb,toptype):
  for i in range(prb.Natom):
    if (toptype=='1') and (prb.chain[i]!='1'):
      prb.select[i] = 0
    if (toptype=='2') and (prb.chain[i]!='1') and (prb.chain[i]!='2'):
      prb.select[i] = 0
    if (toptype=='3') and (prb.chain[i]!='1') and (prb.chain[i]!='2') and (prb.chain[i]!='3'):
      prb.select[i] = 0

def trim_bad_conservation_probes(prb,cons_thre):
  for i in range(prb.Natom):
    if (prb.conservation[i] < cons_thre):
      prb.select[i] = 0




##########################
#########  MAIN ##########
##########################

OPT = {}
OPT['iprb'] = 'OPRB'
OPT['ilig'] = 'LigPDB'
OPT['dth']  = 3.74
OPT['top']  = 'x' 
OPT['con']  = 0.0
if (len(sys.argv)<2):
  print "prbRecPre.py  [prolistfile] <options>"
  print " for calcualting recall-precision of probes for binding liggnds"  
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate) 
  print "<options>"
  print " -iprb : input dir for probe PDB file [%s]"%(OPT['iprb']) 
  print " -ilig : input dir for ligand PDB file [%s]"%(OPT['ilig']) 
  print " -dth  : threshold of matching distance [%f]"%(OPT['dth']) 
  print " -top  : prb_select only top '1','2','3',...'x' [%s]"%(OPT['top']) 
  print " -con  : conservation thre. if (<con) remove that probe.  [%f]"%(OPT['con']) 
  sys.exit(1)

read_option(sys.argv,OPT)
prolistfile = sys.argv[1]
prolist = []
read_protein_list(prolistfile,prolist)
print prolist

Nlig_over = 0
Nprb_over = 0
Nlig = 0
Nprb = 0
for id in (prolist):

  lig = pdb.Molecule()
  lig.read(OPT['ilig']+'/'+id)

  prb = pdb.Molecule()
  prb.read(OPT['iprb']+'/'+id)
  prb.select = [1 for i in range(prb.Natom)]

  trim_top_cluster(prb,OPT['top'])
  trim_bad_conservation_probes(prb,float(OPT['con']))

  [nlig_over, nprb_over,nlig,nprb] = cal_lig_prb_overlap(lig,prb,float(OPT['dth']))

  recall = precision = 0.0
  if (nlig_over>0) and (nlig>0):
    recall   = float(nlig_over)/nlig
  if (nprb_over>0) and (nprb>0):
    precision = float(nprb_over)/nprb

  print "%s nlig_over %d nprb_over %d Rec %f Pre %f"%(id,nlig_over,nprb_over,recall,precision)
  Nlig_over += nlig_over
  Nprb_over += nprb_over
  Nlig += nlig 
  Nprb += nprb 

Recall    = float(Nlig_over)/float(Nlig)
Precision = float(Nprb_over)/float(Nprb)
Fmeasure  = 2.0*Recall*Precision/(Recall+Precision)
print "top %s con %f Recall %f Precision %f Fmeasure %f"%(OPT['top'],float(OPT['con']),Recall,Precision,Fmeasure)

