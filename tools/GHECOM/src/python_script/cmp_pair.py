#!/usr/bin/python

import sys
import os
import math
import pairali 
import respock 
import resacc 
import ligand 

LastModDate = "Nov 12, 2008"

def read_option(argv,opt_dic):
  for i in range(1,len(argv)):
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]

OPT={}
OPT['iali']  = ''
OPT['ipoA']  = ''
OPT['ipoB']  = ''
OPT['iacA']  = ''
OPT['iacB']  = ''
OPT['ilgA']  = ''
OPT['ilgB']  = ''


if (len(sys.argv)<2):
  print "cmppair.py  <options>"
  print " compared pairwise-aligned two protein structure properties"
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate) 
  print "<options>"
  print " -iali : alignment file [%s]"%(OPT['iali']) 
  print " -ipoA : pocket file for protein A [%s]"%(OPT['ipoA']) 
  print " -ipoB : pocket file for protein B [%s]"%(OPT['ipoB']) 
  print " -iacA : Solvent Acc file for protein A [%s]"%(OPT['iacA']) 
  print " -iacB : Solvent Acc file for protein B [%s]"%(OPT['iacB']) 
  print " -ilgA : Ligand file for protein A [%s]"%(OPT['ilgA']) 
  print " -ilgB : Ligand file for protein B [%s]"%(OPT['ilgB']) 
  sys.exit(1)


read_option(sys.argv,OPT)

if (OPT['iali'] != ''):
  ali = pairali.Alignment()
  ali.read_vertical(OPT['iali'])
  print ali
else:
  print "#ERROR:option -iali is obligation."
  sys.exit()


if (OPT['ipoA'] != ''):
  poA = respock.Pocket()  
  poA.read(OPT['ipoA'])

if (OPT['ipoB'] != ''):
  poB = respock.Pocket()  
  poB.read(OPT['ipoB'])


if (OPT['iacA'] != ''):
  acA = resacc.SolventAcc()  
  acA.read(OPT['iacA'])

if (OPT['iacB'] != ''):
  acB = resacc.SolventAcc()  
  acB.read(OPT['iacB'])

if (OPT['ilgA'] != ''):
  lgA = ligand.ResLigand()  
  lgA.read(OPT['ilgA'])

if (OPT['ilgB'] != ''):
  lgB = ligand.ResLigand()  
  lgB.read(OPT['ilgB'])

print "#[alinum(1)] [rnumA(2)] [aaA(3)] [rnumB(4)] [aaB(5)]"
print "#[AccA(6)] [AccB(7)] [RaccA(8)] [RaccB(9)]"
print "#[shellAccA(10)] [shellAccB(11)] [pocketA(12)] [pocketB(13)]"
print "#[NligA(14)] [NligB(15)]"
for a in range(ali.Nalign):
  rnumA = ali.ali_resnum1[a]
  rnumB = ali.ali_resnum2[a]
  shellAccA  = shellAccB = 0.0
  pocketnessA = pocketnessB = 0.0
  AccA = AccB = 0.0
  RaccA = RaccB = 0.0
  NligA = NligB = 0.0
  if (rnumA != "-1") and (rnumB != "-1"): 
    aaA = ali.seq1[rnumA]
    aaB = ali.seq2[rnumB]
    if (OPT['ipoA'] != '') and (OPT['ipoB'] != ''):
      shellAccA    = float(poA.shellAcc[rnumA])
      shellAccB    = float(poB.shellAcc[rnumB])
      pocketnessA = float(poA.pocketness[rnumA])
      pocketnessB = float(poB.pocketness[rnumB])
    if (OPT['iacA'] != '') and (OPT['iacB'] != ''):
      AccA    = float(acA.ACC[rnumA])
      AccB    = float(acB.ACC[rnumB])
      RaccA   = float(acA.RACC[rnumA])
      RaccB   = float(acB.RACC[rnumB])
    if (OPT['ilgA'] != '') and (OPT['ilgB'] != ''):
      if (lgA.conligseq.has_key(rnumA)):
        NligA    = len(lgA.conligseq[rnumA])
      if (lgB.conligseq.has_key(rnumB)):
        NligB    = len(lgB.conligseq[rnumB])

    sys.stdout.write("%-4d %4s %3s %4s %3s"%(a,rnumA,aaA,rnumB,aaB))
    sys.stdout.write(" %6.2f %5.1f %6.2f %5.1f"%(AccA,AccB,RaccA,RaccB))
    sys.stdout.write(" %6.2f %6.2f %6.2f %6.2f"%(shellAccA,shellAccB,pocketnessA,pocketnessB))
    sys.stdout.write(" %2d %2d"%(NligA,NligB))
    sys.stdout.write("\n")
