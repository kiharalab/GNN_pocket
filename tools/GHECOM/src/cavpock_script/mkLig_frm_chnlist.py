#!/usr/bin/env python

##
## <mkLig_frm_chnlist.py>
##
## 

import sys
import os 
import math 
from datetime import datetime
import glob
import mmCIF

LastModDate = '2019/07/10'

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  now = datetime.now()
  opt_dic['START_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        opt_dic[argv[i][1:]] = argv[i+1]


def read_list_file(ifname):
  if (os.access(ifname,os.R_OK)==0):
    print "#ERROR:Can't open list file '%s'"%(ifname)
    sys.exit(1)
  f = open(ifname)
  items = []
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
      columns = line.split()
      items.append(columns[0])
  return(items)


def count_Natom_in_asym_id_molecule(M,asym_id):
  A = M['atom_site']
  Natom = 0
  Nheavyatom = 0
  Nelement = {}
  for i in range(len(A['label_asym_id'])):
    if (A['label_asym_id'][i] == asym_id):
      element = A['type_symbol'][i]
      Nelement[element] = Nelement.get(element,0) + 1
      Natom += 1
      if (element != 'H'):
        Nheavyatom += 1

  return(Natom,Nheavyatom,Nelement)


def cal_contact_bwn_two_asym_ids(M,asym_idA, asym_idB):
  A = M['atom_site']
  posAlist = []
  posBlist = []
  Dthre = 4.0
  DDthre = Dthre*Dthre
  for i in range(len(A['label_asym_id'])):
    if (A['label_asym_id'][i] == asym_idA):
      pos = [0.0,0.0,0.0]
      pos[0] = float(A['Cartn_x'][i])
      pos[1] = float(A['Cartn_y'][i])
      pos[2] = float(A['Cartn_z'][i])
      posAlist.append(pos)
    if (A['label_asym_id'][i] == asym_idB):
      pos = [0.0,0.0,0.0]
      pos[0] = float(A['Cartn_x'][i])
      pos[1] = float(A['Cartn_y'][i])
      pos[2] = float(A['Cartn_z'][i])
      posBlist.append(pos)

  for posA in (posAlist):
    for posB in (posBlist):
      DD =  (posA[0]-posB[0])*(posA[0]-posB[0]) + (posA[1]-posB[1])*(posA[1]-posB[1]) + (posA[2]-posB[2])*(posA[2]-posB[2])
      if (DD <= DDthre):
        return(1)
  return(0)



################
##### MAIN #####
################

OPT = {}
OPT['ghecom'] = 'ghecom'
OPT['odir'] = 'tmpout/'
OPT['A'] = 'F'
OPT['olog'] = 'mkLig_frm_chnlist.log'
OPT['chaingrid'] = '-'
OPT['icifdir'] = '/db1/mmCIF'
OPT['minNheavyatom'] = '7'
OPT['minNcarbon']    = '3'
OPT['minNnuc']    = '3'
OPT['minNaa']     = '10'
OPT['avoid_comp_ids'] = 'BOG:DTT:EPE:GOL:MES:MPD:MRD:PG4:TRS'
OPT['olist'] = 'outlig.list'
OPT['oligdir'] = ''

if (len(sys.argv)<2):
  print "mkLig_frm_chnlist.py <options>"
  print "  for comparing grid pdb files for list of mmCIFs."
  print "  LastModDate:%s"%(LastModDate)
  print "<options>"
  print " -ilist  : input list file with [pdb_id][chain]  [%s]"%(OPT.get('ilist',''))
  print " -icifdir : Input mmCIF directory [%s]"%(OPT.get('icifdir',''))
  print " -minNheavyatom : min Nheavyatom  [%s]"%(OPT['minNheavyatom'])
  print " -minNcarbonn   : min Ncarbon     [%s]"%(OPT['minNcarbon'])
  print " -minNnuc       : min Nnucleotide [%s]"%(OPT['minNnuc'])
  print " -minNaa        : min Naa         [%s]"%(OPT['minNaa'])
  print " -avoid_comp_ids  : avoiding comp_ids  [%s]"%(OPT['avoid_comp_ids'])
  print " -oligdir      : output ligand dir [%s]"%(OPT['oligdir'])
  print " -olog   : output log file [%s]"%(OPT['olog'])
  print " -olist  : output list file [%s]"%(OPT['olist'])
  print " -A : Action ('T' or 'F')[%s]"%(OPT['A'])
  sys.exit()

read_option(sys.argv, OPT)

minNheavyatom = int(OPT['minNheavyatom'])
minNcarbon    = int(OPT['minNcarbon'])
minNnuc       = int(OPT['minNnuc'])
minNaa        = int(OPT['minNaa'])
avoid_comp_ids = OPT['avoid_comp_ids'].split(':')

olist = open(OPT['olist'],'w')
print "#write_list() --> '%s'"%(OPT['olist'])

PDBCHlist = read_list_file(OPT.get('ilist',''))

Npdbch = len(PDBCHlist)
print "#Npdbch %s"%(Npdbch)
olist.write("#COMMAND %s\n"%(OPT['COMMAND']))
olist.write("#DATE    %s\n"%(OPT['START_DATE']))
olist.write("# -ilist  : input list file with [pdb_id][chain]  [%s]\n"%(OPT.get('ilist','')))
olist.write("# -icifdir : Input mmCIF directory [%s]\n"%(OPT.get('icifdir','')))
olist.write("# -minNheavyatom : min Nheavyatom  [%s]\n"%(OPT['minNheavyatom']))
olist.write("# -minNcarbonn   : min Ncarbon     [%s]\n"%(OPT['minNcarbon']))
olist.write("# -minNnuc       : min Nnucleotide [%s]\n"%(OPT['minNnuc']))
olist.write("# -minNaa        : min Naa         [%s]\n"%(OPT['minNaa']))
olist.write("# -avoid_comp_ids  : avoiding comp_ids  [%s]\n"%(OPT['avoid_comp_ids']))
olist.write("# -olog   : output log file [%s]\n"%(OPT['olog']))
olist.write("#Npdbch %s\n"%(Npdbch))
olist.write("#[pdb_id:1] [asym_id:2] [chain_id(auth_asym_id):3] [Nligand:4]\n")
olist.write("#[ligand asym_id:5+2*i]:[ligand comp_id/poly_type:6+2*i]\n")

if (OPT.get('oligdir','') != ''):
  ologfile = OPT['oligdir'] + '/' + OPT['olog']
  olog = open(ologfile,'w')
  olog.write("#COMMAND %s\n"%(OPT['COMMAND']))
  olog.write("#DATE    %s\n"%(OPT['START_DATE']))
  olog.write("# -ilist  : input list file with [pdb_id][chain]  [%s]\n"%(OPT.get('ilist','')))
  olog.write("# -icifdir : Input mmCIF directory [%s]\n"%(OPT.get('icifdir','')))
  olog.write("# -minNheavyatom : min Nheavyatom  [%s]\n"%(OPT['minNheavyatom']))
  olog.write("# -minNcarbonn   : min Ncarbon     [%s]\n"%(OPT['minNcarbon']))
  olog.write("# -minNnuc       : min Nnucleotide [%s]\n"%(OPT['minNnuc']))
  olog.write("# -minNaa        : min Naa         [%s]\n"%(OPT['minNaa']))
  olog.write("# -avoid_comp_ids  : avoiding comp_ids  [%s]\n"%(OPT['avoid_comp_ids']))
  olog.write("# -olog   : output log file [%s]\n"%(OPT['olog']))
  olog.write("#Npdbch %s\n"%(Npdbch))



for pdbch in (PDBCHlist):
  pdb_id   = pdbch[0:4]
  chain_id = pdbch[4:]
  iciffile = OPT['icifdir'] + '/' +  pdb_id[1:3] + '/' + pdb_id + ".cif.gz"

  print "%s %s %s"%(pdb_id, chain_id,iciffile)
  M = {}
  mmCIF.read_mmCIF_file(iciffile,M)
  ENTITY = []
  UNITMOL = []
  mmCIF.make_entity_from_mmCIFdic(M,ENTITY,'T')
  mmCIF.make_unitmol_from_mmCIFdic(M,UNITMOL,ENTITY)

  for u in (UNITMOL):
    print "#>%s"%(u.get('asym_id',''))
    for key in (u.keys()):
      print "'%s': '%s'"%(key,u[key])


  receptor = ''
  ligands = []
  for u in (UNITMOL):
    print "#>asym_id %s auth_asym_id %s nresidue %d chain_id %s"%(u.get('asym_id',''),u.get('auth_asym_id',''),int(u.get('nresidue',0)),chain_id)
    if (u['auth_asym_id'] == chain_id) and (int(u.get('nresidue',0))>=40):
      receptor = u

  if (receptor == ''):
    print "#ERROR: no receptor for %s %s.\n"%(pdb_id,chain_id)
    olist.write("#ERROR: no receptor for %s %s.\n"%(pdb_id,chain_id))
    receptor = {}
  else:
    for u in (UNITMOL):
      if (u['asym_id'] != receptor['asym_id']):
        (natom,nheavyatom,nelement) = count_Natom_in_asym_id_molecule(M,u['asym_id'])
        ncarbon = nelement.get('C',0)
        nresidue = int(u.get('nresidue','0'))
        accept = 1
        if (nheavyatom < minNheavyatom):
          accept = 0 
        if (ncarbon < minNcarbon):
          accept = 0 
        if (u['comp_id'] in avoid_comp_ids):
          accept = 0
        if (u['poly_type'] == 'polypeptide(L)') and (nresidue >= minNaa):
          accept = 0 
        if (u['poly_type'] == 'polydeoxyribonucleotide') and (int(u['nresidue']) >= minNnuc):
          accept = 0 
        if (u['poly_type'] == 'polyribonucleotide') and (int(u['nresidue']) >= minNnuc):
          accept = 0 
        print "#>%s comp_id '%s' nheavyatom %d ncarbon %d nresidue %d poly_type '%s' accept %d"%(u['asym_id'],u.get('comp_id',''),nheavyatom,ncarbon,nresidue,u.get('poly_type',''),accept) 
        if (accept == 1):
          contact = cal_contact_bwn_two_asym_ids(M,receptor['asym_id'], u['asym_id'])
          if (contact == 1):
            ligands.append(u) 
  
   #print "receptor",receptor
  olist.write("%s %s %s %d"%(pdb_id,receptor.get('asym_id','NONE'),receptor.get('auth_asym_id','NONE'),len(ligands)))
  for u in (ligands):
    olist.write(" %s:%s%s"%(u['asym_id'],u['comp_id'],u['poly_type']))
  olist.write("\n")

  if (OPT.get('oligdir','') != ''):
    olog.write(">%s %s %s\n"%(pdb_id,u['asym_id'],u['auth_asym_id'])) 
    for u in (ligands):
      opdbfile = OPT['oligdir'] + '/' + pdb_id + '_' + u.get('asym_id','') + '.pdb'
      mmCIF.write_asym_id_molecule_in_PDB(opdbfile,OPT,M,u.get('asym_id',''))
      olog.write("%s\n"%(opdbfile)) 
      pass


print "#write_list() --> '%s'"%(OPT['olist'])
olist.close()

if (OPT.get('oligdir','') != ''):
  olog.close()
