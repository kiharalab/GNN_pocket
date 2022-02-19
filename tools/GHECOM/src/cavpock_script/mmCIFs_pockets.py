#!/usr/bin/env python

##
## <mmCIFs_pockets.py>
##
## 

import sys
import os 
import math 
from datetime import datetime
import glob

LastModDate = '2019/06/29'


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


def read_tab_splited_table_file(ifname):
  if (os.access(ifname,os.R_OK)==0):
    print "#ERROR:Can't open table file '%s'"%(ifname)
    sys.exit(1)
  f = open(ifname)
  lines = []
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
      columns = line.split('\t')
      lines.append(columns)
  print "#read_tab_splited_table_file('%s') Nline %d"%(ifname,len(lines))
  return(lines)




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


################
##### MAIN #####
################

OPT = {}
OPT['M'] = 'P'
OPT['ghecom'] = 'ghecom'
OPT['gw'] = '2.0'
OPT['rl'] = '10'
OPT['rs'] = '5'
OPT['ghecom'] = 'ghecom'
OPT['odir'] = 'tmpout/'
OPT['A'] = 'F'
OPT['olog'] = 'mmCIFs_pockets.log'
OPT['chaingrid'] = '-'
OPT['osum'] = 'summary.out'
OPT['itab'] = ''
OPT['icifdir'] = '/db1/mmCIF'
OPT['odir'] = 'tmpout/'
OPT['opdb'] = 'PDB/'
OPT['div'] = '0/1'

if (len(sys.argv)<2):
  print "mmCIFs_pockets.py <options>"
  print "  for making pockets for list of mmCIFs."
  print "  LastModDate:%s"%(LastModDate)
  print "<options>"
  print " -itab   : input table file for mmCIF assembly [%s]"%(OPT['itab'])
  print " -icifdir: Input mmCIF directory [%s]"%(OPT.get('icifdir',''))
  print " -M      : 'P'ocket, 'CP':Cave-pocket 'V':ca'V'ity, [%s]"%(OPT.get('M',''))
  print " -odir   : Output directory [%s]"%(OPT['odir'])
  print " -opdb   : Output PDB directory [%s]"%(OPT['opdb'])
  print " -gw     : grid_width [%s]"%(OPT.get('gw',''))
  print " -rl     : Rlarge [%s]"%(OPT.get('rl','10:20:1'))
  print " -rs     : Rlarge [%s]"%(OPT.get('rl','5:10:1'))
  print " -ghecom : ghecom program [%s]"%(OPT['ghecom'])
  print " -olog   : output log file [%s]"%(OPT['olog'])
  print " -div    : Job division [%s]"%(OPT['div'])
  print " -A      : Action ('T' or 'F')[%s]"%(OPT['A'])
  sys.exit()

read_option(sys.argv, OPT)

[bunshi,bunbo] = OPT['div'].split('/')
bunshi = int(bunshi)
bunbo  = int(bunbo)


if (OPT['A']=='T'):
  ologfile = "%s/%s.%d"%(OPT['odir'],OPT['olog'],bunshi)
  olog = open(ologfile,'w')
else:
  olog = sys.stdout

olog.write("#COMMAND    %s\n"%(OPT['COMMAND']))
olog.write("#START_DATE %s\n"%(OPT['START_DATE']))
olog.write("# -itab  : input table file for mmCIF assembly [%s]"%(OPT['itab']))
olog.write("# -icifdir : Input mmCIF directory [%s]\n"%(OPT.get('icifdir','')))
olog.write("# -M: 'P'ocket, 'CP':Cave-pocket 'V':ca'V'ity, [%s]\n"%(OPT['M']))
olog.write("# -odir    : Output directory [%s]\n"%(OPT['odir']))
olog.write("# -gw  : grid_width [%s]\n"%(OPT['gw']))
olog.write("# -rl  : Rlarge [%s]\n"%(OPT['rl']))
olog.write("# -rs  : Rsmall [%s]\n"%(OPT['rs']))

lines = read_tab_splited_table_file(OPT['itab'])
# #[pdb_id:1] [assembly_id:2] [oligometic_count:3] [Nasym_id:4]
# #[asym_id:5] [auth_asym_id:6] [nresidue:7] [resolution_high:8] [expltl_method:9] [uniprod_id:10]
# 12as    1       2       6       B       B       327     2.2     X-RAY DIFFRACTION       ASNA_ECOLI
# 148l    1       2       5       B       S       4       1.9     X-RAY DIFFRACTION       
# 16vp    1       1       3       A       A       311     2.1     X-RAY DIFFRACTION       ATIN_HHV11
# 173d    1       4       4       D       D       6       3.0     X-RAY DIFFRACTION       

rl = float(OPT['rl'])
rs = float(OPT['rs'])

Nline = len(lines);
Nstart   = bunshi*int(Nline/bunbo);
Nend     = (bunshi+1)*int(Nline/bunbo);
if (bunshi>=(bunbo-1)):
  Nend = Nline
print "#Nline %d bunshi/bunbo %d/%d start %d end %d"%(Nline,bunshi,bunbo,Nstart,Nend)


#for line in (lines):
for i in range(Nstart,Nend):
  line = lines[i]
  pdb_id           = line[0]
  assembly_id      = line[1]
  oligomeric_count = line[2]
  print pdb_id, assembly_id, oligomeric_count
  iciffile = OPT['icifdir'] + '/' + pdb_id[1:3] + '/' + pdb_id + '.cif.gz'
  #opocpdbfile = OPT['odir'] + '/' +  pdb_id + '_' + assembly_id + '.pdb'
  #opocpdbfile = "%s/%s_%s_L%.0f_S%.0f.pdb"%(OPT['odir'], pdb_id,assembly_id,rl,rs)
  opocpdbfile = "%s/%s_%s_poc.pdb"%(OPT['odir'], pdb_id,assembly_id)
  opdbfile    = "%s/%s_%s.pdb"%(OPT['opdb'], pdb_id,assembly_id)
  command = "%s -M %s -icif %s -assembly %s -gw %s -clus T -rl %f -rs %f -opocpdb %s"%(OPT['ghecom'],OPT['M'],iciffile,assembly_id,OPT['gw'],rl,rs,opocpdbfile)
  command += " -opdb %s"%(opdbfile)
  print "#COMMAND %s"%(command)
  #command = "%s -M %s -ipdb %s -model M -gw %s -clus T -rl %f -rs %f -opocpdb %s"%(OPT['ghecom'],OPT['M'],OPT.get('ipdb',''),OPT['gw'],rl,rs,opocpdb)
  if (OPT['A'] == 'T'):
    olog.write("#COMMAND %s\n"%(command))
    os.system(command)

if (OPT['A'] == 'T'):
  now = datetime.now()
  olog.write("#END_DATE %s\n"%(now.strftime("%Y/%m/%d %H:%M:%S")))
  olog.close()

sys.exit(1)

