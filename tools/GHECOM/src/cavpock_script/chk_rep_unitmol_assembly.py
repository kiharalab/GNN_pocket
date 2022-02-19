#!/usr/bin/env python

import sys
import os
import math
import re
import pgdb
from datetime import datetime
import re

LastModDate = "2019/06/28"


def read_option(argv,opt_dic):
  now = datetime.now()
  opt_dic['START_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
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
  return(lines)
# #COLUMN_1 unitmol.pdb_id
# #COLUMN_2 unitmol.asym_id
# #COLUMN_3 unitmol.auth_asym_id
# #COLUMN_4 unitmol.nresidue
# #COLUMN_5 unitmol.nheavyatom
# #COLUMN_6 unitmol.entity_id
# #COLUMN_7 unitmol.cluster95
# #COLUMN_8 unitmol.cluster40
# #COLUMN_9 unitmol.clusterE4
# #COLUMN_10 unitmol.uniprot_id
# #COLUMN_11 pdb.resolution_high
# #COLUMN_12 pdb.exptl_method
# #COLUMN_13 Nmember
# 6cnw    A       A       116     883     1       68434   26      26              0.92    X-RAY DIFFRACTION               15491
# 2vb1    A       A       129     1001    1       None    None    None    LYSC_CHICK      0.65    X-RAY DIFFRACTION               11140
# 5le5    Z       Z       213     1641    12      8650    6271    1104    PSB1_HUMAN      1.8     X-RAY DIFFRACTION               10676
# 6hmq    A       A       327     2759    1       3267    1512    156     CSK22_HUMAN     0.97    X-RAY DIFFRACTION               7232
# 3f1l    B       B       244     1897    1       9764    6979    171     YCIK_ECOLI      0.95    X-RAY DIFFRACTION               3399


#############
### MAIN ####
#############


OPT = {}
OPT['dbname']    = 'pdb_mmcif_sheep'
OPT['user']      = 'takawaba'
#OPT['dbhost'] = 'ipproo'
OPT['dbhost'] = 'localhost'
OPT['A'] = 'F'
OPT['irep'] = ''
OPT['oassembly'] = 'out_assembly.list'
OPT['only1'] = 'T'
OPT['exnmr'] = 'T'
OPT['maxreso'] = '4.0'

if (len(sys.argv)<2):
  print "chk_rep_unitmol_assembly.py <options>"
  print " to cheak 'assembly' for representative unitmol list."
  print "  coded by T.Kawabata. LastModDate:%s"%(LastModDate)

  print "<option>"
  print " -irep      : input representative unitmol file [%s]"%(OPT['irep']) 
  print " -dbname    : dbname of SQL db [%s]"%(OPT['dbname']) 
  print " -dbhost    : hostname of SQL db [%s]"%(OPT['dbhost'])
  print " -user      : username of SQL db [%s]"%(OPT['user']) 
  print " -oassembly : output assembly list [%s]"%(OPT['oassembly'])
  print " -only1     : restrict assembly=1 ('T' or 'F')[%s]"%(OPT['only1'])
  print " -exnmr     : exclude NMR ('T' or 'F')[%s]"%(OPT['exnmr'])
  print " -maxreso   : maximum allowed resolution [%s]"%(OPT['maxreso'])
  sys.exit(1)

read_option(sys.argv,OPT)

LINES = read_tab_splited_table_file(OPT['irep'])
RepUnitmolList = []
for cols in (LINES):
  dat = {}
  dat['pdb_id']       = cols[0]
  dat['asym_id']      = cols[1]
  dat['auth_asym_id']  = cols[2]
  dat['nresidue']     = cols[3]
  dat['uniprot_id']     = cols[9]
  dat['resolution_high']     = cols[10]
  dat['exptl_method']     = cols[11]
# #COLUMN_10 unitmol.uniprot_id
# #COLUMN_11 pdb.resolution_high
# #COLUMN_12 pdb.exptl_method
  RepUnitmolList.append(dat)


db = pgdb.connect(database=OPT['dbname'])
cur = db.cursor()

AssemblyDic = {} 
pdb_id_ass_dic = {}
for d in (RepUnitmolList):
  print "#>%s %s"%(d['pdb_id'],d['asym_id'])
  #command = "SELECT pdb_id,assembly_id,oligomeric_count,asym_ids FROM assembly WHERE pdb_id='%s';"%(d['pdb_id'])
  command = "SELECT pdb_id,assembly_id,oligomeric_count,array_length(asym_ids,1),asym_ids FROM assembly WHERE pdb_id='%s' AND '%s' = ANY(asym_ids);"%(d['pdb_id'],d['asym_id'])
  print "#COMMAND '%s'"%(command)
  cur.execute(command)
  db.commit() 
  HitLines = cur.fetchall()
  print HitLines 
  for line in (HitLines):
    a = {}
    a['pdb_id']           = line[0] 
    a['asym_id']          =  d['asym_id'] 
    a['auth_asym_id']     =  d['auth_asym_id'] 
    a['nresidue']         =  d['nresidue'] 
    a['uniprot_id'] = d['uniprot_id']
    a['resolution_high'] = d['resolution_high']
    a['exptl_method']  = d['exptl_method']
    a['assembly_id']      = line[1] 
    a['oligomeric_count'] = line[2] 
    a['Nasym_id']         = line[3] 
    index = a['pdb_id'] + '_' + a['assembly_id']
    accept = 1
    if (OPT['only1'] == 'T') and (a['assembly_id'] != '1'):
      accept = 0

    if (OPT['exnmr'] == 'T') and (a['exptl_method'] == 'SOLUTION NMR'):
      accept = 0
    if (OPT['exnmr'] == 'T') and (a['exptl_method'] == 'SOLID-STATE NMR'):
      accept = 0

    if (a['resolution_high'] == 'None'):
      accept = 0
    elif (float(a['resolution_high'])>float(OPT['maxreso'])):
      accept = 0

    if (index in pdb_id_ass_dic):
      accept = 0    

    pdb_id_ass_dic[index] = pdb_id_ass_dic.get(index,0) + 1   
 
    if (accept == 1): 
      AssemblyDic[index] = a

db.close()

pdb_ass_list = sorted(AssemblyDic.keys(),lambda x,y:cmp(x,y))

of = open(OPT['oassembly'],'w')
print "#write_assembly_list() --> '%s'"%(OPT['oassembly'])
of.write("#COMMAND %s\n"%(OPT['COMMAND']))
of.write("#DATE    %s\n"%(OPT['START_DATE']))
of.write("# -irep      : input representative unitmol file [%s]\n"%(OPT['irep'])) 
of.write("# -dbname    : dbname of SQL db [%s]\n"%(OPT['dbname'])) 
of.write("# -dbhost    : hostname of SQL db [%s]\n"%(OPT['dbhost']))
of.write("# -user      : username of SQL db [%s]\n"%(OPT['user'])) 
of.write("# -oassembly : output assembly list [%s]\n"%(OPT['oassembly']))
of.write("# -only1     : restrict assembly=1 ('T' or 'F')[%s]\n"%(OPT['only1']))
of.write("# -exnmr     : exclude NMR ('T' or 'F')[%s]\n"%(OPT['exnmr']))
of.write("# -maxreso   : maximum allowed resolution [%s]\n"%(OPT['maxreso']))
of.write("#Nrep_unitmol %d\n"%(len(RepUnitmolList)))
of.write("#Nassembly    %d\n"%(len(pdb_ass_list)))
of.write("#[pdb_id:1] [assembly_id:2] [oligometic_count:3] [Nasym_id:4]\n")
of.write("#[asym_id:5] [auth_asym_id:6] [nresidue:7] [resolution_high:8] [expltl_method:9] [uniprod_id:10]\n")


for index in (pdb_ass_list):
  a = AssemblyDic[index]
  of.write("%s\t%s\t%s\t%s\t%s\t%s\t%s"%(a['pdb_id'],a['assembly_id'],a['oligomeric_count'],a['Nasym_id'],a['asym_id'],a['auth_asym_id'],a['nresidue']))
  of.write("\t%s\t%s\t%s\n"%(a['resolution_high'],a['exptl_method'],a['uniprot_id']))

print "#write_assembly_list() --> '%s'"%(OPT['oassembly'])
of.close()


