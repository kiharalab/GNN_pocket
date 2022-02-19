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


def merge_res_pocket_and_profile(rpock,prof,vaa,chain,contype):
  print "def merge_res_pocket_and_profile(rpock,prof,vaa,chain,contype):"
  for i in range(1,vaa.length):
    if (contype == 'F'):
      cons = prof.freq[i][prof.qseq[i]]
    elif (contype == 'I'):
      cons  = prof.info[i]
    rnumpdb = vaa.rnumpdb[i]
    index = rnumpdb + ' ' + chain
    #print "i %d qseq %s index '%s' cons %f"%(i,prof.qseq[i],index,cons)
    rpock.conservation[index] = cons
  pass





def write_res_base_pocket_with_profile(ofname,rpock,proid):
  print "#write_res_base_pocket_with_profile() --> '%s'"%(ofname)
  of = open(ofname,"w")
  of.write("#>>RESIDUE-BASED POCKET FILE WITH PROFILE SCORE<<\n")
  of.write("#COLUMN  1|RNUM              |Residue Number\n")
  of.write("#COLUMN  2|CHAIN             |Chain Identifier\n")
  of.write("#COLUMN  3|RES               |Three-letter residue name\n")
  of.write("#COLUMN  4|shellAcc          |shell accessibility (%)\n")
  of.write("#COLUMN  5|Rinacc            |averaged Rinaccess (A)\n")
  of.write("#COLUMN  6|Natom             |Number of atoms\n")
  of.write("#COLUMN  7|Natom_contact     |Number of contacting atoms with ligand\n")
  of.write("#COLUMN  8|pocketness        |sum of 1/[Rpocket] /(1/[Rmin]*[vol of shell]) (%)\n")
  of.write("#COLUMN  9|pocketness_clus 1 |sum of 1/[Rpocket_for_pocketcluster 1] /(1/[Rmin]*[vol of shell]) (%)\n")
  of.write("#COLUMN 10|pocketness_clus 2 |sum of 1/[Rpocket_for_pocketcluster 2] /(1/[Rmin]*[vol of shell]) (%)\n")
  of.write("#COLUMN 11|pocketness_clus 3 |sum of 1/[Rpocket_for_pocketcluster 3] /(1/[Rmin]*[vol of shell]) (%)\n")
  of.write("#COLUMN 12|pocketness_clus 4 |sum of 1/[Rpocket_for_pocketcluster 4] /(1/[Rmin]*[vol of shell]) (%)\n")
  of.write("#COLUMN 13|pocketness_clus 5 |sum of 1/[Rpocket_for_pocketcluster 5] /(1/[Rmin]*[vol of shell]) (%)\n")
  of.write("#COLUMN 14|conservation      |freq_aa of query amino acids calcualted by PSI-BLAST (%)\n")
  for r in (rpock.rnumchlist):
    of.write("%4s  %s %s %6.2f %6.3f %2d %2d %6.2f"%(rpock.rnum[r],rpock.chain[r],rpock.res[r],rpock.shellAcc[r],rpock.Rinacc[r],rpock.Natom[r],rpock.Natom_contact[r],rpock.pocketness[r]))
    of.write(" %6.2f %6.2f %6.2f %6.2f %6.2f"%(rpock.pocketness_clus1[r], rpock.pocketness_clus2[r], rpock.pocketness_clus3[r], rpock.pocketness_clus4[r], rpock.pocketness_clus5[r]))
    if (rpock.conservation.has_key(r)):
      of.write(" %6.1f"%(rpock.conservation[r]))
    else: 
      of.write(" %6.1f"%(-1.0))
    of.write(" %s\n"%(proid))
  of.close() 


def read_pdb_file_and_write_it_with_profile(ifname,ofname,rpock):
  print "#read_pdb_file_and_write_it_with_profile('%s')-->'%s'"%(ifname,ofname)
  fi = open(ifname,"r")
  fo = open(ofname,"w")
  for line in fi:
    if (line.startswith('ATOM')) or (line.startswith('HETATM')): 
      line = line.rstrip('\n')
      chain   = line[21:22]
      if (chain == ' '):
        chain = '-'
      rnumstr = line[22:27].replace(' ','')
      index = rnumstr + ' ' + chain
      tFactor = 0.0
      if (rpock.conservation.has_key(index)):
        tFactor = float(rpock.conservation[index])
      fo.write("%s%6.1f\n"%(line[0:60],tFactor))
    else:
      fo.write("%s"%(line))
  fi.close()
  fo.close() 


def set_conservation_to_probes(prb,pdb,rpock):
  for i in range(prb.Natom):
    a1 = a2 = a3 = -1 
    ind1 = ind2 = ind3 = ''
    cons1 = cons2 = cons3 = -1.0 
    N = 0
    sum = 0.0
    if (pdb.anum2int.has_key(prb.cAt1[i])):
     a1 = pdb.anum2int[prb.cAt1[i]]
     ind1 = pdb.rnum[a1].replace(' ','') + ' ' + pdb.chain[a1]
     if (rpock.conservation.has_key(ind1)):
       cons1 = rpock.conservation[ind1]
       sum += cons1
       N   += 1 
    if (pdb.anum2int.has_key(prb.cAt2[i])):
     a2 = pdb.anum2int[prb.cAt2[i]]
     ind2 = pdb.rnum[a1].replace(' ','') + ' ' + pdb.chain[a1]
     if (rpock.conservation.has_key(ind2)):
       cons2 = rpock.conservation[ind2]
       sum += cons2
       N   += 1 
    if (pdb.anum2int.has_key(prb.cAt3[i])):
     a3 = pdb.anum2int[prb.cAt3[i]]
     ind3 = pdb.rnum[a1].replace(' ','') + ' ' + pdb.chain[a1]
     if (rpock.conservation.has_key(ind3)):
       cons3 = rpock.conservation[ind3]
       sum += cons3
       N   += 1 

    #print "%s %d %d %d -> %d %d %d -> '%s' '%s' '%s' -> %f %f %f"%(pdb.filename,prb.cAt1[i],prb.cAt2[i],prb.cAt3[i],a1,a2,a3,ind1,ind2,ind3,cons1,cons2,cons3)
    if (N>0): 
      prb.conservation[i] = float(sum)/N


def write_probe_PDB_file_with_conservation(ofname,prb):
  print "#write_probe_PDB_file_with_conservation()-->'%s'"%(ofname)
  of = open(ofname,"w")
  if (prb.header != ''):
    of.write("%s\n"%(prb.header))
  if (prb.title != ''):
    of.write("%s\n"%(prb.title))
  of.write("REMARK  COMMAND %s\n"%(OPT['COMMAND']))
  of.write("REMARK  DATE %s\n"%(OPT['DATE']))
  of.write("REMARK  [Ratom](10: 55 - 60):Radius of atom[A]\n")
  of.write("REMARK  [tFact](11: 61 - 66):[Riacc] is assigned as tFactor\n")
  of.write("REMARK  [shAcc](12: 67 - 72):Shell accessibility. (Ratio of noVdW grids in the shell)[%]\n")
  of.write("REMARK  [Riacc](13: 73 - 78):averaged Rinaccess in the 'shell'[A]\n")
  of.write("REMARK  [co]   (14: 79 - 81):Contact with other molecules (0 or 1)\n")
  of.write("REMARK  [pocke](15: 82 - 87):pocketness. sum of 1/[Rpocket] /(1/[Rmin]*[vol of shell]) (%)\n")
  of.write("REMARK  [pock1](16: 88 - 93):pocketness_clus[1] sum of 1/[Rpocket_for_pocke_tcluster1] /(1/[Rmin]*[vol of shell]) (%)\n")
  of.write("REMARK  [pock2](17: 94 - 99):pocketness_clus[2] sum of 1/[Rpocket_for_pocke_tcluster2] /(1/[Rmin]*[vol of shell]) (%)\n")
  of.write("REMARK  [pock3](18:100 -105):pocketness_clus[3] sum of 1/[Rpocket_for_pocke_tcluster3] /(1/[Rmin]*[vol of shell]) (%)\n")
  of.write("REMARK  [pock4](19:106 -111):pocketness_clus[4] sum of 1/[Rpocket_for_pocke_tcluster4] /(1/[Rmin]*[vol of shell]) (%)\n")
  of.write("REMARK  [pock5](20:112 -117):pocketness_clus[5] sum of 1/[Rpocket_for_pocke_tcluster5] /(1/[Rmin]*[vol of shell]) (%)\n")
  of.write("REMARK  [cAt1] (21:118 -122):Atom number of contact protein atom 1\n")
  of.write("REMARK  [cAt2] (22:123 -127):Atom number of contact protein atom 2\n")
  of.write("REMARK  [cAt3] (23:128 -133):Atom number of contact protein atom 3\n")
  of.write("REMARK  [cons] (24:133 -139):Conservation score\n")

  of.write("REMARK                                                 Ratom|cons |shAcc|Riacc|co|pocke|pock1|pock2|pock3|pock4|pock5|cAt1|cAt2|cAt3|cons|\n");
  for i in range(prb.Natom):
    line = prb.line[i]
    head = line[0:60]
    tail = line[66:]
    of.write("%s%6.1f%s%6.1f\n"%(head,prb.conservation[i],tail,prb.conservation[i]))
  pass

##########################
#########  MAIN ##########
##########################

OPT = {}
OPT['ires'] = 'ORES'
OPT['iprf'] = 'PROF'
OPT['ivaa'] = 'VAAseq'
OPT['oprf'] = 'ORES_PROF'
OPT['ipdb'] = 'RecPDB'
OPT['opdb'] = ''
OPT['con']  = 'F'
OPT['iprb'] = ''
OPT['oprb'] = ''

if (len(sys.argv)<2):
  print "mrgProf2Res.py  [prolistfile] <options>"
  print " merge PSI-BLAST profile files into ghecom residue-based pocketness files"  
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate) 
  print "<options>"
  print "  -con : type of conservation score 'F'req_aa of query 'I'nfomation [%s]"%(OPT['con'])
  print " -iprf : input dir for PSI-BLAST ascii profile [%s]"%(OPT['iprf']) 
  print " -ivaa : input dir for vertical sequences [%s]"%(OPT['ivaa']) 
  print "<options for res-based pocketness files>" 
  print " -ires : input  dir for res-based pocketness files [%s]"%(OPT['ires']) 
  print " -oprf : output dir for res-based pocketness with profile[%s]"%(OPT['oprf']) 
  print "<options for receptor PDB files>" 
  print " -ipdb : input  dir for receptor PDB file [%s]"%(OPT['ipdb']) 
  print " -opdb : output dir for receptor PDB file with conservation score[%s]"%(OPT['opdb']) 
  print "<options for probe PDB files>" 
  print " -iprb : input  dir for probe PDB file [%s]"%(OPT['iprb']) 
  print " -oprb : output dir for probe PDB file with conservation score[%s]"%(OPT['oprb']) 
  print "  (The '-ipdb' option is also necessary.)"
  sys.exit(1)

read_option(sys.argv,OPT)
prolistfile = sys.argv[1]
prolist = []
read_protein_list(prolistfile,prolist)
print prolist

for id in (prolist):
  #pdb = id[0:4] 

  if (OPT['ires'] != ''):
    Pck = res_base_pock.Pocket()
    Pck.read(OPT['ires']+'/'+id)
    Pck.conservation = {}
  
  seqlist = []
  seqlist = glob.glob(OPT['iprf']+'/'+id+'*')
  for seqfile in (seqlist):
    print seqfile
    Prf = psiprof.Profile()
    Prf.read(seqfile) 
    [dir,seq] = seqfile.split('/')
    chain  = seq[-1] 
    Vaa = pdbrnumseq.PdbRnumSeq()
    Vaa.read(OPT['ivaa']+'/'+seq) 
    if (OPT['ires'] != ''):
      merge_res_pocket_and_profile(Pck,Prf,Vaa,chain,OPT['con'])

  write_res_base_pocket_with_profile(OPT['oprf']+'/'+id,Pck,id)

  if (OPT['ipdb'] != '') and (OPT['opdb'] != ''):
    read_pdb_file_and_write_it_with_profile(OPT['ipdb']+'/'+id,OPT['opdb']+'/'+id,Pck)
    pass    

  if (OPT['iprb'] != '') and (OPT['oprb']!='') and (OPT['ipdb'] != ''):
    Pdb = pdb.Molecule()
    Pdb.read(OPT['ipdb']+'/'+id)
    Prb = pdb.Molecule()
    Prb.read(OPT['iprb']+'/'+id)
    set_conservation_to_probes(Prb,Pdb,Pck)
    write_probe_PDB_file_with_conservation(OPT['oprb']+'/'+id,Prb)
    pass
