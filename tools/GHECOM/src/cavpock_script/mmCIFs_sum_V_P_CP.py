#!/usr/bin/env python

##
## <mmCIFs_sum_V_P_CP.py>
##
## 

import sys
import os 
import math 
from datetime import datetime
import glob

LastModDate = '2019/07/04'

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
OPT['ghecom'] = 'ghecom'
OPT['gw'] = '2.0'
OPT['rl'] = '10'
OPT['rs'] = '5'
OPT['ghecom'] = 'ghecom'
OPT['odir'] = 'tmpout/'
OPT['A'] = 'F'
OPT['olog'] = 'mmCIFs_pockets.log'
OPT['chaingrid'] = '-'
OPT['itab'] = ''
OPT['icifdir'] = '/db1/mmCIF'
OPT['osum'] = 'summary.out'
OPT['gw'] = '2.0'
OPT['rl'] = '10.0'

if (len(sys.argv)<2):
  print "mmCIFs_sum_V_P_CP.py <options>"
  print "  for comparing grid pdb files for list of mmCIFs."
  print "  LastModDate:%s"%(LastModDate)
  print "<options>"
  print " -itab  : input table file for mmCIF assembly [%s]"%(OPT['itab'])
  print " -icifdir : Input mmCIF directory [%s]"%(OPT.get('icifdir',''))
  print " -odirV   : Output dirV  [%s]"%(OPT.get('odirV',''))
  print " -odirP   : Output dirP  [%s]"%(OPT.get('odirP',''))
  print " -odirCP  : Output dirCP [%s]"%(OPT.get('odirCP',''))
  print " -gw      : grid width [%s]"%(OPT.get('gw',''))
  print " -rl      : Rlarge [%s]"%(OPT.get('rl',''))
  print " -ghecom : ghecom program [%s]"%(OPT['ghecom'])
  print " -osum  : output summary file [%s]"%(OPT['osum'])
  print " -olog : output log file [%s]"%(OPT['olog'])
  print " -A : Action ('T' or 'F')[%s]"%(OPT['A'])
  sys.exit()

read_option(sys.argv, OPT)

if (OPT['A']=='T'):
  #olog = open(OPT['odir']+'/'+ OPT['olog'],'w')
  osum = open(OPT['osum'],'w')
  print "#write_summary() --> '%s'"%(OPT['osum'])
else:
  osum = sys.stdout

osum.write("#COMMAND    %s\n"%(OPT['COMMAND']))
osum.write("#START_DATE %s\n"%(OPT['START_DATE']))
osum.write("# -itab  : input table file for mmCIF assembly [%s]\n"%(OPT['itab']))
osum.write("# -odirV   : Output dirV  [%s]\n"%(OPT.get('odirV','')))
osum.write("# -odirP   : Output dirP  [%s]\n"%(OPT.get('odirP','')))
osum.write("# -odirCP  : Output dirCP [%s]\n"%(OPT.get('odirCP','')))
osum.write("# -odirA   : Output directoryA [%s]\n"%(OPT.get('odirA','')))
osum.write("# -odirB   : Output directoryB [%s]\n"%(OPT.get('odirB','')))
osum.write("# -gw      : grid width [%s]\n"%(OPT.get('gw','')))
osum.write("# -rl      : Rlarge [%s]\n"%(OPT.get('rl','')))
osum.write("# -icifdir : Input mmCIF directory [%s]\n"%(OPT.get('icifdir','')))

lines = read_tab_splited_table_file(OPT['itab'])
# #[pdb_id:1] [assembly_id:2] [oligometic_count:3] [Nasym_id:4]
# #[asym_id:5] [auth_asym_id:6] [nresidue:7] [resolution_high:8] [expltl_method:9] [uniprod_id:10]
# 12as    1       2       6       B       B       327     2.2     X-RAY DIFFRACTION       ASNA_ECOLI
# 148l    1       2       5       B       S       4       1.9     X-RAY DIFFRACTION       
# 16vp    1       1       3       A       A       311     2.1     X-RAY DIFFRACTION       ATIN_HHV11
# 173d    1       4       4       D       D       6       3.0     X-RAY DIFFRACTION       

osum.write("#[pdb_id:1] [assembly_id:2] [oligomeric_count:3]\n")
osum.write("#[N_CHAIN:4] [N_RESIDUE:5] [VOL_VDW:6]\n")
osum.write("#[volV:7] [volP:8] [volCP:9]\n")
osum.write("#[volV_P:10] [volV_CP:11] [volP_CP:12]\n")


for line in (lines):
  pdb_id           = line[0]
  assembly_id      = line[1]
  oligomeric_count = line[2]
  print pdb_id, assembly_id, oligomeric_count
  iciffile = OPT['icifdir'] + '/' + pdb_id[1:3] + '/' + pdb_id + '.cif.gz'
  commandX  = "%s X -icif %s -assembly %s -gw %s -rl %s"%(OPT['ghecom'],iciffile,assembly_id,OPT.get('gw','2.0'),OPT.get('rl','10.0'))
  #opocpdbfile = OPT['odir'] + '/' +  pdb_id + '_' + assembly_id + '.pdb'
  #opocpdbfile = "%s/%s_%s_L%.0f_S%.0f.pdb"%(OPT['odir'], pdb_id,assembly_id,rl,rs)
  opocpdbfileV  = "%s/%s_%s_poc.pdb"%(OPT.get('odirV',''),  pdb_id,assembly_id)
  opocpdbfileP  = "%s/%s_%s_poc.pdb"%(OPT.get('odirP',''),  pdb_id,assembly_id)
  opocpdbfileCP = "%s/%s_%s_poc.pdb"%(OPT.get('odirCP',''), pdb_id,assembly_id)
  commandV_P  = "%s GcmpG -igridpdbA %s -igridpdbB %s"%(OPT['ghecom'],opocpdbfileV,opocpdbfileP)
  commandV_CP = "%s GcmpG -igridpdbA %s -igridpdbB %s"%(OPT['ghecom'],opocpdbfileV,opocpdbfileCP)
  commandP_CP = "%s GcmpG -igridpdbA %s -igridpdbB %s"%(OPT['ghecom'],opocpdbfileP,opocpdbfileCP)
  print "#commandX %s"%(commandX)
  print "#commandV_P %s"%(commandV_P)
  print "#commandV_CP %s"%(commandV_CP)
  print "#commandP_CP %s"%(commandP_CP)
  #command = "%s -M %s -ipdb %s -model M -gw %s -clus T -rl %f -rs %f -opocpdb %s"%(OPT['ghecom'],OPT['M'],OPT.get('ipdb',''),OPT['gw'],rl,rs,opocpdb)
  if (OPT['A'] == 'T'):
    print "#commandX %s"%(commandX)
    X = perform_ghecom_and_get_key_values(commandX)

    #osum.write("#COMMAND %s\n"%(command))
    print "#commandV_P %s"%(commandV_P)
    V_P  = perform_ghecom_and_get_key_values(commandV_P)
    print "#commandV_CP %s"%(commandV_CP)
    V_CP = perform_ghecom_and_get_key_values(commandV_CP)
    print "#commandP_CP %s"%(commandP_CP)
    P_CP = perform_ghecom_and_get_key_values(commandP_CP)
    osum.write("%s %s %3d"%(pdb_id,assembly_id,int(oligomeric_count)))
    osum.write(" %3d %7d %10.1f"%(int(X.get('N_CHAIN','0')),int(X.get('N_RESIDUE','0')), float(X.get('VOL_VDW','0.0')) ) ) 
    #osum.write(" %10.1f %10.1f %10.1f"%(float(D.get('VP','0')),float(D.get('VL','0')),float(D.get('VPL','0')) ))
    #osum.write("  %.3f  %.3f  %.3f"%(float(D.get('RECALL','0')),float(D.get('PRECISION','0')),float(D.get('TANIMOTO','0')) ))
    osum.write(" %10.1f %10.1f %10.1f"%(float(V_P.get('VP','0.0')),float(V_P.get('VL','0')), float(V_CP.get('VL','0.0')) ) ) 
    osum.write(" %10.1f %10.1f %10.1f"%(float(V_P.get('VPL','0.0')),float(P_CP.get('VPL','0')), float(P_CP.get('VPL','0.0')) ) ) 
    osum.write("\n")

if (OPT['A'] == 'T'):
  now = datetime.now()
  osum.write("#END_DATE %s\n"%(now.strftime("%Y/%m/%d %H:%M:%S")))
  print "#write_summary() --> '%s'"%(OPT['osum'])
  osum.close()

sys.exit(1)

