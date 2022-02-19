#!/usr/bin/python


import os
import sys

LastModDate = "Aug 5, 2009"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]


def read_list_file(ifname,list,prop):
 ## [FILE_EXAMPLE]
 ## 1dlwA  a.1.1.1
 ## 1a6m-  a.1.1.2
 ## 1eca-  a.1.1.2
 ## 1lh1-  a.1.1.2
 ## 1ash-  a.1.1.2
 if not os.access(ifname,os.R_OK):
   print "#ERROR:Can't open listfile '%s'" % ifname
   return(1)

# list = []
# prop = {}
 f = open(ifname)
 for line in f:
   line = line.rstrip('\n')
   if (line.startswith('#')==0) and (len(line)>2):
     field = line.split()
     list.append(field[0])
     if (len(field)>2):
       prop[field[0]] = field[1] 
 print "#Nlist %d"%(len(list)) 
 f.close()

##############
#### MAIN ####
##############

OPT = {}
OPT['of'] = ''
OPT['ipdb'] = 'ChConPDB'
OPT['opo']  = 'PHECOM_PRB'
OPT['oR']   = 'PHECOM_RAC'
OPT['A']  = 'F'
OPT['pl']  = 'M'
OPT['rl']  = float(6.0) 
OPT['nrl'] = 8
OPT['div'] = '0/1'
OPT['G'] = 'F'
OPT['og'] = 'PHECOM_GRID'
OPT['xtA'] = -1.0

if (len(sys.argv)<2):
  print "mkphecom.py [listfile] <options>"
  print " to perform 'phecom' calculation for the protein chain list"
  print "<options>"
  print " -ipdb : Input  dir. for selected PDB chains [%s]"%(OPT['ipdb'])
  print " -pl   : Type of Large Probe Generation"
  print "       : 'I'cosa '3'-con '2'-con 'X':3-con-fast 'P'lane 'S'pheHull"
  print "       : 'M'ultiple-size probes 'm':write multiple size probes [%s]"%(OPT['pl'])
  print " -opo  : Output dir for  pocket probe PDB File (only probes) [%s]"%(OPT['opo'])
  print " -oR  : Output Pocket Residue Raccess File [%s]"%(OPT['oR'])
  print " -nrl  : Number of trying Rlarge for '-pl M' option.[%d]"%(OPT['nrl'])
  print " -og   : Output dir. for grid PDB file [%s]"%(OPT['og'])
  print " -of   : Output file [%s]"%(OPT['of'])
  print " -div  : division (bunshi/bunsbo)[%s]"%(OPT['div'])
  print " -xtA  : selection for probe Rinaccess (3,4,5,...10)[%f]"%(OPT['xtA'])
  print " -A    : Action ('T' or 'F') [%s]"%(OPT['A'])
  print " -G    : probes to Grid Calculation ('T' or 'F') [%s]"%(OPT['G'])
  sys.exit(1)

read_option(sys.argv,OPT)
ilistfile = sys.argv[1]
OPT['xtA'] = float(OPT['xtA'])
LIST = []
PROP = {}
read_list_file(ilistfile,LIST,PROP)


[bunshi,bunbo] = OPT['div'].split('/');
bunshi = int(bunshi)
bunbo = int(bunbo)
Nlist = len(LIST);
Nstart = bunshi*int(Nlist/bunbo);
Nend   = (bunshi+1)*int(Nlist/bunbo);
if (bunshi>=(bunbo-1)):
  Nend = Nlist
print "#Nlist %d bunshi/bunbo %d/%d start %d end %d"%(Nlist,bunshi,bunbo,Nstart,Nend)

for i in range(Nstart,Nend):
  pdbch = LIST[i] 
  chain = pdbch[4:]
  if (OPT['G'] == 'T'):
    command = "gridcomp -M G -ipA %s/%s -aA H -resA PRB -gA T -xtA %f -og %s/%s"%(OPT['opo'],pdbch,OPT['xtA'],OPT['og'],pdbch)
  else:
    command = "phecom %s/%s -ch %s -pl %s -nrl %d -opo %s/%s -oR %s/%s -op tmp.mpdb"%(OPT['ipdb'],pdbch,chain,OPT['pl'],int(OPT['nrl']),OPT['opo'],pdbch,OPT['oR'],pdbch)

  if (OPT['A'] == 'T'):
    os.system(command)
  else:
    print "#%s"%(command)

