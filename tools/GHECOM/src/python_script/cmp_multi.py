#!/usr/bin/python

import sys
import os
import math
import pairali 
import multiali 
import respock 
import resacc 
import ligand 

LastModDate = "Sep 22, 2008"

def read_option(argv,opt_dic):
  for i in range(1,len(argv)):
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]


def write_multi_aligned_prorperties(ofname,mali,pock,acc,type):
  print "#write_multi_aligned_prorperties(type %s) -->'%s'"%(type,ofname)
  of = open(ofname,"w")
  for p in range(mali.Nprotein):
    of.write("#protein %3d %s\n"%(p, mali.proname[p]))
  if (type=='P'):
    of.write("#VALUE pocketness\n")
  if (type=='A'):
    of.write("#VALUE relative accessibility\n")

  of.write("#COLUMN   1|nsite |site number\n"); 
  of.write("#COLUMN   2|Npro  |number of aligned protein for the site\n"); 
  of.write("#COLUMN   3|ave   |average value\n"); 
  of.write("#COLUMN   4|SD    |standard deviation\n"); 
  of.write("#COLUMN   5|min   |minimum value\n"); 
  of.write("#COLUMN   6|max   |maximum value\n"); 
  for p in range(mali.Nprotein):
    of.write("#COLUMN %3d|val%-3d|value for protein %d '%s'\n"%(p+7,p+1,p+1,mali.proname[p])); 

  of.write("#[nsite][Npro][ave][SD][min][max]"); 
  for p in range(mali.Nprotein):
    of.write("  [val%d]"%(p+1)); 
  of.write("\n"); 

  for a in range(mali.Nalign):
    nalipro = 0
    S  = 0.0
    SS = 0.0
    min = 0.0
    max = 0.0
    mean = 0.0
    SD = 0.0
    # calculate mean, SD, min max #
    for p in range(mali.Nprotein):
      rnum = mali.ali_resnum[p][a]
      if (rnum != "-1"):
        nalipro += 1
 
        value = 0.0
        if (type=='P'):
          value = float(pock[p].pocketness[rnum])
        if (type=='A'):
          value = float(acc[p].RACC[rnum])

        S  += value
        SS += value*value
        if (p==0):
          min = value
          max = value
        else:
          if (value<min):
            min = value 
          if (value>max):
            max = value 
      if (nalipro>0):
        mean = S/nalipro
        SD = SS/nalipro - mean*mean
        if (SD>0.0):
          SD = math.sqrt(SD)
    if (nalipro<mali.Nprotein):
      of.write("#")
    of.write("%4d %3d %6.2f %6.2f %6.2f %6.2f"%(a,nalipro,mean,SD,min,max))

    # output properties for proteins #
    for p in range(mali.Nprotein):
      rnum = mali.ali_resnum[p][a]
      if (rnum != "-1"):
        if (type=='P'):
          value = float(pock[p].pocketness[rnum])
        if (type=='A'):
          value = float(acc[p].RACC[rnum])

        of.write(" %6.2f"%(value))
      else:
        of.write(" ------")

    of.write("\n") 
  of.close() 
  pass


def write_multi_aligned_contact_ligands(ofname,mali,lig):
  print "#write_multi_aligned_contact_ligands() -->'%s'"%(ofname)
  of = open(ofname,"w")
  for p in range(mali.Nprotein):
    of.write("#protein %3d %s\n"%(p, mali.proname[p]))
  if (type=='P'):
    of.write("#VALUE pocketness\n")
  if (type=='A'):
    of.write("#VALUE relative accessibility\n")

  of.write("#COLUMN   1|nsite |site number\n"); 
  of.write("#COLUMN   2|Npro  |number of aligned protein for the site\n"); 
  of.write("#COLUMN   3|ave   |average value\n"); 
  for p in range(mali.Nprotein):
    of.write("#COLUMN %3d|val%-3d|value for protein %d '%s'\n"%(p+4,p+1,p+1,mali.proname[p])); 

  of.write("#[nsite][Npro]"); 
  for p in range(mali.Nprotein):
    of.write("  [val%d]"%(p+1)); 
  of.write("\n"); 

  for a in range(mali.Nalign):
    nalipro = 0
    nligpro = 0
    # calculate mean, SD, min max #
    for p in range(mali.Nprotein):
      rnum = mali.ali_resnum[p][a]
      if (rnum != "-1"):
        nalipro += 1
        if (lig[p].seq.has_key(rnum)):
          nligpro += 1

    if (nalipro<mali.Nprotein):
      of.write("#")
    of.write("%4d %3d %3d"%(a,nalipro,nligpro))
    
    # output properties for proteins #
    for p in range(mali.Nprotein):
      rnum = mali.ali_resnum[p][a]
      of.write(" %s"%(rnum))
      if (rnum != "-1"):
        if (lig[p].seq.has_key(rnum)):
          of.write(" %s"%(lig[p].conligseq[rnum][0]))
        else:
          of.write(" ===")
      else:
        of.write(" ---")
    of.write("\n")
  of.close() 











#################
##### MAIN ######
#################

OPT={}
OPT['ima']    = ''
OPT['dpo']   = ''
OPT['epo']   = ''
OPT['dac']   = ''
OPT['eac']   = ''
OPT['dli']   = ''
OPT['eli']   = ''

OPT['ofp']    = 'pmul'
OPT['ofa']    = 'amul'
OPT['ofl']    = 'lmul'

if (len(sys.argv)<2):
  print "cmp_multi.py  <options>"
  print " compared pairwise-aligned two protein structure properties"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate) 
  print "<options>"
  print " -ima  : input (vertical) multiple alignment file [%s]"%(OPT['ima']) 
  print " -dpo  : directory for pocket file [%s]"%(OPT['dpo']) 
  print " -epo  : file extention for pocket file [%s]"%(OPT['epo']) 
  print " -dac  : directory for accessibility file [%s]"%(OPT['dac']) 
  print " -eac  : file extention for accessibility file [%s]"%(OPT['eac']) 
  print " -dli  : directory for ligand file [%s]"%(OPT['dac']) 
  print " -eli  : file extention for ligand file [%s]"%(OPT['eac']) 
  print " -ofp  : output file for pocketness[%s]"%(OPT['ofp']) 
  print " -ofa  : output file for accessiblity[%s]"%(OPT['ofa']) 
  print " -ofl  : output file for ligand[%s]"%(OPT['ofl']) 
  sys.exit(1)

read_option(sys.argv,OPT)

if (OPT['ima'] != ''):
  mali = multiali.MultiAlign()
  mali.read_vertical(OPT['ima'])
  print mali
  for p in range(mali.Nprotein):
   print "%d %s"%(p,mali.proname[p])
else:
  print "#ERROR:option -ima is obligation."
  sys.exit()

pock  = [respock.Pocket() for i in range(mali.Nprotein)] 
acc   = [resacc.SolventAcc() for i in range(mali.Nprotein)] 
lig   = [ligand.ResLigand() for i in range(mali.Nprotein)] 

if (OPT['dpo'] != ''):
  for p in range(mali.Nprotein):
    fname = OPT['dpo'] + '/' + mali.proname[p]
    if (OPT['epo'] != ''):
      fname = fname + '.' + OPT['epo'] 
    print "#fname_pock:%s"%(fname)
    pock[p].read(fname)
  write_multi_aligned_prorperties(OPT['ofp'],mali,pock,acc,'P')


if (OPT['dac'] != ''):
  for p in range(mali.Nprotein):
    fname = OPT['dac'] + '/' + mali.proname[p]
    if (OPT['eac'] != ''):
      fname = fname + '.' + OPT['eac'] 
    print "#fname_acc:%s"%(fname)
    acc[p].read(fname)
  write_multi_aligned_prorperties(OPT['ofa'],mali,pock,acc,'A')


if (OPT['dli'] != ''):
  for p in range(mali.Nprotein):
    fname = OPT['dli'] + '/' + mali.proname[p]
    if (OPT['eli'] != ''):
      fname = fname + '.' + OPT['eli'] 
    print "#fname_acc:%s"%(fname)
    lig[p].read(fname)
  write_multi_aligned_contact_ligands(OPT['ofl'], mali,lig)

sys.exit(1)

