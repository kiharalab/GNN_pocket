#!/usr/bin/env python

##
## <mk_pocket_rl_rs.py>
##
## 

import sys
import os 
import math 
from datetime import datetime
import glob

LastModDate = '2019/07/03'


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
OPT['rlmin'] = '10'
OPT['rlmax'] = '20'
OPT['rlbin'] = '1'

OPT['rsmin'] = '5'
OPT['rsmax'] = '10'
OPT['rsbin'] = '1'
OPT['ghecom'] = 'ghecom'
OPT['odir'] = 'tmpout/'
OPT['gw'] = '2.0'
OPT['A'] = 'F'
OPT['olog'] = 'mk_pockets_rl_rs.log'
OPT['mkeval'] = 'M'
OPT['chaingrid'] = '-'
OPT['osum'] = 'summary.out'
OPT['div'] = '0/1'


if (len(sys.argv)<2):
  print "mk_pockets_rl_rs.py <options>"
  print "  to make pockets with different rl and rs."
  print "  LastModDate:%s"%(LastModDate)
  print "<options>"
  print " -mkeval: 'M'ake pocket 'E'val pocket [%s]"%(OPT['mkeval'])
  print " -M: 'P'ocket(masuya_doi), 'CP':Cave-pocket 'V':ca'V'ity, [%s]"%(OPT.get('M',''))
  print " -ipdb : Input pdbfile [%s]"%(OPT.get('ipdb',''))
  print " -gw  : grid_width [%s]"%(OPT.get('gw',''))
  print " -rl  : Rlarge [min]:[max]:[bin] [%s]"%(OPT.get('rl','10:20:1'))
  print " -rs  : Rlarge [min]:[max]:[bin] [%s]"%(OPT.get('rl','5:10:1'))
#  print " -rlmin: min of rl [%s]"%(OPT.get('rlmin',''))
#  print " -rlmin: min of rl [%s]"%(OPT.get('rlmin',''))
#  print " -rlmax: max of rl [%s]"%(OPT.get('rlmax',''))
#  print " -rlbin: div of rl [%s]"%(OPT.get('rlbin',''))
#  print " -rsmin: min of rl [%s]"%(OPT.get('rsmin',''))
#  print " -rsmax: max of rl [%s]"%(OPT.get('rsmax',''))
#  print " -rsbin: div of rl [%s]"%(OPT.get('rsbin',''))
  print " -odir : output directory [%s]"%(OPT.get('odir',''))
  print " -ghecom : ghecom program [%s]"%(OPT['ghecom'])
  print " -olog : output log file [%s]"%(OPT['olog'])
  print " -A : Action ('T' or 'F')[%s]"%(OPT['A'])
  print " -div    : Job division (only for -mkeval M) [%s]"%(OPT['div'])
  print " <options only for '-mkeval E'>"
  print "  -isphepdb  :Input reference sphere file in PDB [%s]"%(OPT.get('isphepdb',''))
  print "  -iligpdb   :Input reference ligand file in PDB [%s]"%(OPT.get('iligpdb',''))
  print "  -chaingrid : Cluster_num(chainID) select for grid. '1':1,'2':1,2,'3':1,2,3 '-':all. [%s]"%(OPT.get('chaingrid','-'))
  print "  -osum      : output summary file [%s]"%(OPT['osum'])
  sys.exit()

read_option(sys.argv, OPT)

[bunshi,bunbo] = OPT['div'].split('/')
bunshi = int(bunshi)
bunbo  = int(bunbo)

if (os.path.isdir(OPT['odir'])==0):
  command = "mkdir %s"%(OPT['odir'])
  print "#%s\n"%(command)
  if (OPT['A']=='T'):
    os.system(command) 

gw    = float(OPT['gw'])
(rlmin,rlmax,rlbin) = OPT['rl'].split(':')
rlmin = float(rlmin)
rlmax = float(rlmax)
rlbin = float(rlbin)
rlNdiv = int((rlmax-rlmin)/rlbin)

(rsmin,rsmax,rsbin) = OPT['rs'].split(':')
rsmin = float(rsmin)
rsmax = float(rsmax)
rsbin = float(rsbin)
rsNdiv = int((rsmax-rsmin)/rsbin)

rsNdiv_start = 0
rsNdiv_end   = rsNdiv + 1 
if (OPT['mkeval'] == 'M'):
  rsNdiv_start   = bunshi*int((rsNdiv+1)/bunbo);
  rsNdiv_end     = (bunshi+1)*int((rsNdiv+1)/bunbo);
  if (bunshi>=(bunbo-1)):
    Nend = rsNdiv  + 1
  print "#rsNdiv %d bunshi/bunbo %d/%d start %d end %d"%(rsNdiv,bunshi,bunbo,rsNdiv_start,rsNdiv_end)


Nall = (rlNdiv+1)*(rsNdiv+1)


if (OPT.get('mkeval','') == 'E') and (OPT['A']=='T'):
  OPT['olog'] = 'eval_' + OPT['olog'] 

if (OPT['A']=='T'):
  ologfile = "%s/%s.%d"%(OPT['odir'],OPT['olog'],bunshi)
  olog = open(ologfile,'w')
else:
  olog = sys.stdout


olog.write("#COMMAND    %s\n"%(OPT['COMMAND']))
olog.write("#START_DATE %s\n"%(OPT['START_DATE']))
olog.write("#M %s\n"%(OPT['M']))

print "#rlmin %f rlmax %f rldiv %f rlNdiv %d"%(rlmin, rlmax, rlbin, rlNdiv)
print "#rsmin %f rsmax %f rsdiv %f rsNdiv %d"%(rsmin, rsmax, rsbin, rsNdiv)
olog.write("#gw %f\n"%(gw))
olog.write("#rlmin %f rlmax %f rldiv %f rlNdiv %d\n"%(rlmin, rlmax, rlbin, rlNdiv))
olog.write("#rsmin %f rsmax %f rsdiv %f rsNdiv %d\n"%(rsmin, rsmax, rsbin, rsNdiv))
olog.write("#Nall %d\n"%(Nall))
olog.write("#A %s\n"%(OPT['A']))

if (OPT.get('mkeval','') == 'E') and (OPT['A']=='T'):
  osum = open(OPT['odir']+'/'+OPT['osum'],'w')
  print "#output_summary() --> '%s'"%(OPT['osum'])
  osum.write("#COMMAND    %s\n"%(OPT['COMMAND']))
  osum.write("#START_DATE %s\n"%(OPT['START_DATE']))
  osum.write("#[Rlarge:1] [Rsmall:2] [Tanimoto:3] [Recall:4] [Precision:5] [Vpocket:6] [Vligand:7]\n")
for i in range(rlNdiv+1):
  rl = rlmin + i*rlbin
  #for j in range(rsNdiv+1):
  for j in range(rsNdiv_start,rsNdiv_end):
    rs = rsmin + j*rsbin
    opocpdb = OPT['odir'] + '/' + "%s_L%.0f_S%.0f.pdb"%(OPT['M'],rl,rs)
    print "#[%d/%d %d/%d] rl %f rs %f"%(i,rlNdiv,j,rsNdiv,rl,rs)
    if (OPT.get('mkeval','') == 'M'):
      command = "%s %s -ipdb %s -model M -gw %s -clus T -rl %f -rs %f -opocpdb %s"%(OPT['ghecom'],OPT['M'],OPT.get('ipdb',''),OPT['gw'],rl,rs,opocpdb)
      print "#COMMAND %s"%(command)
      olog.write("#COMMAND %s\n"%(command))
      if (OPT['A'] == 'T'):
        os.system(command)

    if (OPT.get('mkeval','') == 'E'):
      if (OPT.get('isphepdb','') != ''):
        command = "%s -M GcmpS -igridpdb %s -isphepdb %s -chaingrid %s"%(OPT['ghecom'],opocpdb,OPT.get('isphepdb',''),OPT.get('chaingrid','-'))
      if (OPT.get('iligpdb','') != ''):
        command = "%s -M GcmpL -igridpdb %s -iligpdb %s -chaingrid %s"%(OPT['ghecom'],opocpdb,OPT.get('iligpdb',''),OPT.get('chaingrid','-'))
      print "#COMMAND %s"%(command)
      olog.write("#COMMAND %s\n"%(command))
      if (OPT['A'] == 'T'):
        #os.system(command)
        dat = perform_ghecom_and_get_key_values(command)
        osum .write("%f %f %s %s %s %s %s\n"%(rl,rs,dat.get('TANIMOTO','0.0'),dat.get('RECALL','0.0'),dat.get('PRECISION','0.0'),dat.get('VP','0.0'),dat.get('VL','0.0')))

  if (OPT.get('mkeval','') == 'E') and (OPT['A']=='T'):
    osum.write("\n")

if (OPT.get('mkeval','') == 'M') and (OPT['A']=='T'):
  now = datetime.now()
  olog.write("#END_DATE %s\n"%(now.strftime("%Y/%m/%d %H:%M:%S")))
  olog.close()
    
if (OPT.get('mkeval','') == 'E') and (OPT['A']=='T'):

  print "#output_summary() --> '%s'"%(OPT['odir']+'/'+OPT['osum'])
  osum.close()



