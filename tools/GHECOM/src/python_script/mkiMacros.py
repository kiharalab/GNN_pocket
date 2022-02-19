#!/usr/bin/python

import os
import sys

LastModDate = "Oct 1, 2009"

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
OPT['of'] = 'out.iim'
OPT['dp'] = 'ChConPDB'
OPT['dq'] = 'QSiteOUTPUT'
OPT['base'] = '/home/takawaba/etc/ghecom_scop1.73'
OPT['nolig'] = 'F'

if (len(sys.argv)<2):
  print "mkiMacros.py [listfile] <options>"
  print " for making firefox iMacros file for Q-SiteFinder"
  print " coded by T. Kawabata. LastModDate:%s"%(LastModDate)
  print "<options>"
  print " -base   : base direcotory [%s]"%(OPT['base'])
  print " -dp     : inputdirectory for PDBfiles [%s]"%(OPT['dp'])
  print " -dq     : output directory for selected chains [%s]"%(OPT['dq'])
  print " -nolig  : no ligand in PDB file ('T' or 'F') [%s]"%(OPT['nolig'])
  print " -of     : output file [%s]"%(OPT['of'])
  print "         : put *.iim file in '~/iMacros/Macros'."
  sys.exit(1)

read_option(sys.argv,OPT)
ilistfile = sys.argv[1]
LIST = []
PROP = {}
read_list_file(ilistfile,LIST,PROP)

if (OPT['of'] == '-'):
  of  = sys.stdout
else:
  of = open(OPT['of'],'w')
  print "#write to -->'%s'"%(OPT['of'])

for pdbch in (LIST):
  of.write("VERSION BUILD=6240709 RECORDER=FX\n")
  of.write("TAB T=1\n")
  of.write("URL GOTO=http://www.modelling.leeds.ac.uk/qsitefinder/\n")
  of.write("TAG POS=1 TYPE=INPUT:FILE FORM=ACTION:http://www.modelling.leeds.ac.uk/cgi-bin/qsitefinder/qsf_first.cgi ATTR=ID:protein CONTENT=%s/%s/%s\n"%(OPT['base'],OPT['dp'],pdbch))
  of.write("TAG POS=1 TYPE=INPUT:SUBMIT FORM=ACTION:http://www.modelling.leeds.ac.uk/cgi-bin/qsitefinder/qsf_first.cgi ATTR=VALUE:Submit\n")
  if (OPT['nolig'] != 'T'):
    of.write("TAG POS=1 TYPE=INPUT:SUBMIT FORM=ACTION:http://www.modelling.leeds.ac.uk/cgi-bin/qsitefinder/qsf_pleasewait.cgi ATTR=VALUE:Submit\n")
  of.write("TAG POS=1 TYPE=A ATTR=TXT:Download\n")
  of.write("ONDOWNLOAD FOLDER=%s/%s FILE=%s WAIT=YES\n"%(OPT['base'],OPT['dq'],pdbch))

if (OPT['of'] != '-'):
  of.close()

