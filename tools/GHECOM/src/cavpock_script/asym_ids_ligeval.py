#!/usr/bin/env python

##
## <asym_ids_ligeval.py>
##
## 

import sys
import os 
import math 
from datetime import datetime

LastModDate = '2019/07/11'


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




def read_pdb_asym_list_file(ifname):
  if (os.access(ifname,os.R_OK)==0):
    print "#ERROR:Can't open table file '%s'"%(ifname)
    sys.exit(1)
  f = open(ifname)
  PDB_ASYMS = []
  LIG_ASYMS = []
  LINES     = []

# #[pdb_id:1] [asym_id:2] [chain_id(auth_asym_id):3] [Nligand:4]
# #[ligand asym_id:5+2*i]:[ligand comp_id/poly_type:6+2*i]
# 12as A A 2 C:ASN D:AMP
# 1a0i A A 1 B:ATP
# 1a34 A A 1 E:U
# 1a8p A A 1 B:FAD
# 1a8r A A 2 P:GTP T:GTP
# 1a9x A A 4 EA:ORN FA:NET CA:ADP DA:ADP
# 1ab8 A A 2 C:FOK D:FOK
# 1ac0 A A 12 C:BGC B:GLC D:GLC G:GLO F:BGC I:GLC H:GLC J:GLC M:GLC L:GLC O:GLC N:GLC
# 1ddg A A 1 C:FAD
# 1ddl C A 1 B:polyribonucleotide
# 1de4 C C 2 J:NAG L:NAG
# 1xg0 A A 6 G:DBV F:DBV I:PEB K:PEB J:PEB M:PEB
# 1xg0 B B 5 G:DBV F:DBV M:PEB O:PEB N:PEB
# 1xg0 C C 6 G:DBV F:DBV I:PEB K:PEB J:PEB M:PEB
# 1xg5 A A 1 E:NAP

  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
      columns = line.split()
      pdb_id  = columns[0]
      asym_id = columns[1]
      pdb_asym = pdb_id + '_' + asym_id
      lig_asyms = []
      for i in range(4,len(columns)):
        asym_comp_id  = columns[i]
        (asym,comp_id) = asym_comp_id.split(':')
        lig_asyms.append(asym)

      PDB_ASYMS.append(pdb_asym)
      LIG_ASYMS.append(lig_asyms)
      LINES.append(line)

  print "#read_pdb_asym_list_file('%s') Nline %d"%(ifname,len(LINES))
  return(PDB_ASYMS,LIG_ASYMS,LINES)



def perform_ghecom_and_get_key_values(command):
#>>Example of -KV output of ghecom<<
# #COMMAND ghecom -M GcmpS -igridpdb tmpout/P_g2.0_l20.0_s10.0.pdb -isphepdb consphe.pdb
# #grid_width 2.000000
# #OrigPos 66.260002 68.241997 -36.641998
# #Ngrid_read 9612
# #Output_CHAR3DMAP_PDB_format(mode w min_val 1) -->"pock_grid.pdb" #M->N 111 105 213
# #Natom 7893
# #Compare_Two_CHAR3DMAP(A,minP 1 maxP 255 B minL 1 maxL 255)
# NP 9612
# NL 7893
# NPL 6404
# RECALL 0.811352
# PRECISION 0.666251
# FMESURE 0.731677
# TANIMOTO 0.576885
  valkey = {}
  f = os.popen(command)
  for line in f:
    line = line.rstrip('\n')
    field = line.split()
    if (line.startswith('#')==0) and (len(field)>=2):
      key   = field[0]
      value = field[1]
      valkey[key] = value
    pass
  f.close()
  return(valkey)


def merge_muptiple_files(ifilenames,omergefile,opt):
  print "#merge_muptiple_files(%s) --> '%s'"%(ifilenames, omergefile)
  of = open(omergefile,'w') 
  of.write("REMARK  COMMAND %s\n"%(opt['COMMAND']))
  of.write("REMARK  DATE    %s\n"%(opt['START_DATE']))
  for ifile in (ifilenames):
    if (os.access(ifile,os.R_OK)==0):
      print "#ERROR:Can't open '%s'."%(ifile)
      sys.exit(1)
    f = open(ifile)
    for line in f:
      of.write("%s"%(line))
  of.close()
  pass

def merge_muptiple_asym_ids_PDB_files(pdbdir,pdb_id,asym_ids,omergefile,opt):
  print "#merge_muptiple_files(%s %s %s) --> '%s'"%(pdbdir,pdb_id,asym_ids, omergefile)
  of = open(omergefile,'w') 
  of.write("REMARK  COMMAND %s\n"%(opt['COMMAND']))
  of.write("REMARK  DATE    %s\n"%(opt['START_DATE']))
  for asym_id in (asym_ids):
    ifile = "%s/%s_%s.pdb"%(pdbdir,pdb_id,asym_id)
    if (os.access(ifile,os.R_OK)==0):
      print "#ERROR:Can't open '%s'."%(ifile)
      sys.exit(1)
    f = open(ifile)
    for line in f:
      of.write("%s"%(line))
  of.close()
  pass


################
##### MAIN #####
################

OPT = {}
OPT['ghecom'] = 'ghecom'
OPT['odir'] = 'tmpout/'
OPT['A'] = 'F'
OPT['chaingrid'] = '-'
OPT['osum'] = 'summary.out'
OPT['ilist'] = ''
OPT['odir'] = 'tmpout/'
OPT['opdb'] = 'PDB/'
OPT['ligdir'] = 'LIGAND_PDB'
OPT['model_num'] = '-1'
OPT['tmpmerge'] = 'tmpmerge_for_asym_ids_ligeval.pdb'

if (len(sys.argv)<2):
  print "asym_ids_ligeval.py <options>"
  print "  for making pockets for list of mmCIFs."
  print "  LastModDate:%s"%(LastModDate)
  print "<options>"
  print " -ilist   : input list file with (pdb_id) (asym_id) [%s]"%(OPT['ilist'])
  print " -odir    : Output pocket directory [%s]"%(OPT['odir'])
  print " -ligdir  : ligand PDB file directory [%s]"%(OPT['ligdir'])
  print " -ghecom  : ghecom program [%s]"%(OPT['ghecom'])
  print " -osum    : output summary file [%s]"%(OPT['osum'])
  print "-model_num: select model number  '-1':no-select '1':model=1 'le2' model<=2 [%s]"%(OPT['model_num'])
  print "-tmpmerge : temporary merged ligand PDB file [%s]"%(OPT['tmpmerge'])
  print " -A       : Action ('T' or 'F')[%s]"%(OPT['A'])
  sys.exit()

read_option(sys.argv, OPT)


PID =os.getpid()

tmp_lig_pdbfile = "%s/%s.%d"%(OPT['odir'],OPT['tmpmerge'],PID)

if (OPT['A']=='T'):
  osum = open(OPT['osum'],'w')
else:
  osum = sys.stdout

osum.write("#COMMAND    %s\n"%(OPT['COMMAND']))
osum.write("#START_DATE %s\n"%(OPT['START_DATE']))
osum.write("# -ilist  : input list file with (pdb_id) (asym_id) [%s]\n"%(OPT['ilist']))
osum.write("# -odir    : Output pocket directory [%s]\n"%(OPT['odir']))
osum.write("# -ligdir  : ligand PDB file directory [%s]\n"%(OPT['ligdir']))
osum.write("# -osum    : output summary file [%s]\n"%(OPT['osum']))

(PDB_ASYMS, LIG_ASYMS, LINES) = read_pdb_asym_list_file(OPT['ilist'])

NPtotal  = 0
NLtotal  = 0
NPLtotal = 0
osum.write("#[pdb_asym:1] [NATOM_LIG:2] [NP:3] [NL:4] [NPL:5]\n")
osum.write("#[RECALL:6] [PRECISION:7] [FMEASURE:8] [TANIMOTO:9]\n")

for i in range(len(LINES)):
  pdb_asym  = PDB_ASYMS[i] 
  lig_asyms = LIG_ASYMS[i] 
  (pdb_id, asym_id) = pdb_asym.split('_')
  #print pdb_id, asym_id
  opocpdbfile = "%s/%s_%s_poc.pdb"%(OPT['odir'], pdb_id,asym_id)
  #ligpdbfiles = glob.glob(OPT['ligdir'] + '/' + pdb_id + '_*.pdb')
  #merge_muptiple_files(ligpdbfiles,tmp_lig_pdbfile, OPT)
  merge_muptiple_asym_ids_PDB_files(OPT['ligdir'],pdb_id,lig_asyms,tmp_lig_pdbfile,OPT)

  opdbfile    = "%s/%s_%s.pdb"%(OPT['opdb'], pdb_id,asym_id)
  command = "%s GcmpL -igridpdb %s -iligpdb %s"%(OPT['ghecom'],opocpdbfile,tmp_lig_pdbfile)
  if (OPT['model_num'] != '-1'):
    command += " -model_num %s"%(OPT['model_num'])
  print "#COMMAND %s"%(command)
  #command = "%s -M %s -ipdb %s -model M -gw %s -clus T -rl %f -rs %f -opocpdb %s"%(OPT['ghecom'],OPT['M'],OPT.get('ipdb',''),OPT['gw'],rl,rs,opocpdb)
  if (OPT['A'] == 'T'):
    D = perform_ghecom_and_get_key_values(command)
    #osum.write("#COMMAND %s\n"%(command))
    NATOM_LIG  = int(D.get('NATOM_LIG','0'))   
    NP  = int(D.get('NP','0'))   
    NL  = int(D.get('NL','0'))   
    NPL = int(D.get('NPL','0'))   
    NPtotal  += NP
    NLtotal  += NL
    NPLtotal += NPL
    RECALL    = float(D.get('RECALL','0.0'))   
    PRECISION = float(D.get('PRECISION','0.0'))   
    FMEASURE  = float(D.get('FMEASURE','0.0'))   
    TANIMOTO  = float(D.get('TANIMOTO','0.0'))   
    osum.write("%s %4d %5d %5d %5d %f %f %f %f %s\n"%(pdb_asym,NATOM_LIG,NP,NL,NPL,RECALL,PRECISION,FMEASURE,TANIMOTO,LINES[i]))


osum.write("#NPtotal  %d\n"%(NPtotal))
osum.write("#NLtotal  %d\n"%(NLtotal))
osum.write("#NPLtotal %d\n"%(NPLtotal))
RECALLtotal    = 0.0
PRECISIONtotal = 0.0
FMEASUREtotal  = 0.0
TANIMOTOtotal  = 0.0

if (NLtotal>0):
  RECALLtotal = float(NPLtotal)/float(NLtotal)

if (NPtotal>0):
  PRECISIONtotal = float(NPLtotal)/float(NPtotal)

if ((NPtotal+NLtotal-NPtotal)>0):
  TANIMOTOtotal = float(NPLtotal)/float(NPtotal + NLtotal - NPLtotal)
if ((RECALLtotal+PRECISIONtotal)>0.0):
  FMEASUREtotal = (RECALLtotal*PRECISIONtotal)/(RECALLtotal+PRECISIONtotal)

osum.write("#RECALLtotal    %f\n"%(RECALLtotal))
osum.write("#PRECISIONtotal %f\n"%(PRECISIONtotal))
osum.write("#FMEASUREtotal  %f\n"%(FMEASUREtotal))
osum.write("#TANIMOTOtotal  %f\n"%(TANIMOTOtotal))
osum.write("#%s %s REC_PREC_FMEAS_TANI %f %f %f %f\n"%(OPT['odir'],OPT['model_num'],RECALLtotal,PRECISIONtotal,FMEASUREtotal,TANIMOTOtotal))

if (OPT['A'] == 'T'):
  now = datetime.now()
  osum.write("#END_DATE %s\n"%(now.strftime("%Y/%m/%d %H:%M:%S")))
  osum.close()
  print "#output_summary() --> '%s'"%(OPT['osum'])
sys.exit(1)

