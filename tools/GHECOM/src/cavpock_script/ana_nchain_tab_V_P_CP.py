#!/usr/bin/env python

##
## <ana_nchain_tab_V_P_CP.py>
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


def read_space_splited_table_file(ifname):
  if (os.access(ifname,os.R_OK)==0):
    print "#ERROR:Can't open table file '%s'"%(ifname)
    sys.exit(1)
  f = open(ifname)
  lines = []
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
      #columns = line.split('\t')
      columns = line.split()
      lines.append(columns)
  print "#read_tab_splited_table_file('%s') Nline %d"%(ifname,len(lines))
  return(lines)

def nchain_class(nchain,boundary):
  Nboundary = len(boundary)
  for i in range(len(boundary)):
    if (nchain>=boundary[i]):
      return(i)


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
  print "ana_nchain_tab_V_P_CP.py <options>"
  print "  for analyzing tables V_P_CP."
  print "  LastModDate:%s"%(LastModDate)
  print "<options>"
  print " -itab  : input table file for mmCIF assembly [%s]"%(OPT['itab'])
  print " -A : Action ('T' or 'F')[%s]"%(OPT['A'])
  sys.exit()

read_option(sys.argv, OPT)

if (OPT['A']=='T'):
  #olog = open(OPT['odir']+'/'+ OPT['olog'],'w')
  osum = open(OPT['osum'],'w')
  print "#write_summary() --> '%s'"%(OPT['osum'])
else:
  osum = sys.stdout

lines = read_space_splited_table_file(OPT['itab'])
# #[pdb_id:1] [assembly_id:2] [oligomeric_count:3]
# #[N_CHAIN:4] [N_RESIDUE:5] [VOL_VDW:6]
# #[volV:7] [volP:8] [volCP:9]
# #[volV_P:10] [volV_CP:11] [volP_CP:12]
# 1a8l 1   2   2     452    42960.0        0.0        0.0        0.0        0.0        0.0        0.0
# 1ae9 1   6   6    1011    93016.0     4120.0     2408.0     2408.0        0.0     2408.0     2408.0
# 1aht 1   3   3     252    27496.0        0.0        0.0        0.0        0.0        0.0        0.0
# 1bpo 1   3   3    1467   135440.0        0.0    13136.0    13136.0        0.0    13136.0    13136.0
# 1c4v 1   3   3     254    28224.0        0.0        0.0        0.0        0.0        0.0        0.0
# 1c8n 1 180 180   35820  3204360.0  6310288.0        0.0  6460208.0        0.0        0.0        0.0

#print lines
BOUNDARY_NCHAIN = [60,30,12,7,6,5,4,3,2,1]


cnt_NchCls   = {}
cnt_NchCls_V = {}
cnt_NchCls_P = {}
cnt_NchCls_CP = {}

for line in (lines):
  N_CHAIN = int(line[3])
  cls = nchain_class(N_CHAIN,BOUNDARY_NCHAIN)
  #print "N_CHAIN %d class_N_CHAIN %d"%(N_CHAIN,cls)
  volV  = float(line[6])
  volP  = float(line[7])
  volCP = float(line[8])
  V  = 0
  P  = 0
  CP = 0
  if (volV>0.0):
    V = 1
  if (volP>0.0):
    P = 1
  if (volCP>0.0):
    CP = 1

  cnt_NchCls[cls] = cnt_NchCls.get(cls,0) + 1

  if (cnt_NchCls_V.get(cls,'')==''):
    cnt_NchCls_V[cls] = {}
  cnt_NchCls_V[cls][V] = cnt_NchCls_V[cls].get(V,0) + 1

  if (cnt_NchCls_P.get(cls,'')==''):
    cnt_NchCls_P[cls] = {}
  cnt_NchCls_P[cls][P] = cnt_NchCls_P[cls].get(P,0) + 1

  if (cnt_NchCls_CP.get(cls,'')==''):
    cnt_NchCls_CP[cls] = {}
  cnt_NchCls_CP[cls][CP] = cnt_NchCls_CP[cls].get(CP,0) + 1

CHAINlist = sorted(cnt_NchCls.keys(),lambda x,y:cmp(int(x),int(y)))

print "#COMMAND %s"%(OPT['COMMAND'])
print "#DATE    %s"%(OPT['START_DATE'])
for i in range(len(BOUNDARY_NCHAIN)):
  cls = len(BOUNDARY_NCHAIN) - i -1 
  #print "cls %s cnt %d V %d"%(cls,cnt_NchCls[cls], cnt_NchCls_V[cls].get(1,0))
  sys.stdout.write("%2d -- %2d "%(BOUNDARY_NCHAIN[cls],BOUNDARY_NCHAIN[cls-1]-1))
  sys.stdout.write("cnt %3d"%(cnt_NchCls[cls]))
  sys.stdout.write(" V %3d"%(cnt_NchCls_V[cls].get(1,0)))
  sys.stdout.write(" P %3d"%(cnt_NchCls_P[cls].get(1,0)))
  sys.stdout.write(" CP %3d"%(cnt_NchCls_CP[cls].get(1,0)))

  sys.stdout.write(" V %5.1f"%( 100.0*cnt_NchCls_V[cls].get(1,0)/float(cnt_NchCls[cls])))
  sys.stdout.write(" P %5.1f"%( 100.0*cnt_NchCls_P[cls].get(1,0)/float(cnt_NchCls[cls])))
  sys.stdout.write(" CP %5.1f"%(100.0*cnt_NchCls_CP[cls].get(1,0)/float(cnt_NchCls[cls])))
  sys.stdout.write("\n")
