#!/usr/bin/python

import os
import sys

LastModDate = "Mar 11, 2010"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]


##############
#### MAIN ####
##############

OPT = {}
OPT['of'] = 'out.ras'
OPT['nc'] = 24
OPT['min'] = 0 
OPT['max'] = 1000 
OPT['col'] = 'BGR' 
OPT['tgt'] = 'protein' 

if (len(sys.argv)<2):
  print "mkrasscr.py <options>"
  print " for making rasmol coloring script file"
  print " coded by T. Kawabata. LastModDate:%s"%(LastModDate)
  print "<options>"
  print " -nc     : number of color [%d]"%(OPT['nc'])
  print " -min    : minimum temprature [%d]"%(OPT['min'])
  print " -max    : maximum temprature [%d]"%(OPT['max'])
  print " -of     : output file  [%s]"%(OPT['of'])
  print " -col    : color scheme [%s]"%(OPT['col'])
  print " -tgt    : target molecule 'protein', 'GRD', 'ligand' [%s]"%(OPT['tgt'])
  sys.exit(1)

read_option(sys.argv,OPT)
OPT['min'] = int(OPT['min'])
OPT['max'] = int(OPT['max'])

print "#write_rasmol_script -->'%s'"%(OPT['of'])
of = open(OPT['of'],'w')


of.write("echo \"%s\"\n"%(OPT['COMMAND']))
of.write("background white\n")
of.write("set specular true\n")

for n in range (0,int(OPT['nc'])+1):
  vL  = float(n-1)/(float(OPT['nc']))
  vH  = float(n)/(float(OPT['nc']))
  lowT  = OPT['min'] + vL*(OPT['max'] - OPT['min']) 
  highT = OPT['min'] + vH*(OPT['max'] - OPT['min']) 

  v = vH
  print n,v,lowT,highT

  if (OPT['col']=='WR'):
    r = int(255*1.0) 
    g = int(255*(1.0-v)) 
    b = int(255*(1.0-v))

  if (OPT['col']=='RW'):
    r = int(255*1.0) 
    g = int(255*v) 
    b = int(255*v)

  if (OPT['col']=='BGR'):
    if (0.0<=v) and (v<0.25):
      r = int(255*0.0) 
      g = int(255*4*v) 
      b = int(255*1.0)
    elif (0.25<=v) and (v<0.50):
      r = int(255*0.0) 
      g = int(255*1.0) 
      b = int(255*(1.0-4*(v-0.25)))
    elif (0.50<=v) and (v<0.75):
      r = int(255*(4*(v-0.5)))
      g = int(255*1.0) 
      b = int(255*0.0)
    elif (0.75<=v) and (v<=1.00):
      r = int(255*1.0) 
      g = int(255*(1.0-4*(v-0.75)))
      b = int(255*0.0) 

  if (OPT['col']=='WBGR'):
    if (0.0<=v) and (v<2.0/6.0):
      r = int(255*(1.0-3*v))
      g = int(255*(1.0-3*v))
      b = int(255*1.0) 
    elif (2.0/6.0<=v) and (v<3.0/6.0):
      r = int(255*0)
      g = int(255*(6*(v-2.0/6.0)))
      b = int(255*1.0) 
    elif (3.0/6.0<=v) and (v<4.0/6.0):
      r = int(255*0)
      g = int(255*1)
      b = int(255*(1.0-6*(v-3.0/6.0))) 
    elif (4.0/6.0<=v) and (v<5.0/6.0):
      r = int(255*(6*(v-4.0/6.0)))
      g = int(255*1)
      b = int(255*0) 
    elif (5.0/6.0<=v) and (v<=6.0/6.0):
      r = int(255*1)
      g = int(255*(1.0-6*(v-5.0/6.0)))
      b = int(255*0) 

  if (OPT['col']=='WGR'):
    if (0.0<=v) and (v<0.25):
      r = int(255*(1.0-2*v)) 
      g = int(255*1.0) 
      b = int(255*(1.0-2*v))
    elif (0.25<=v) and (v<0.50):
      r = int(255*(0.5-2*(v-0.25))) 
      g = int(255*1.0) 
      b = int(255*(0.5-2*(v-0.25)))
    elif (0.50<=v) and (v<0.75):
      r = int(255*(4*(v-0.5)))
      g = int(255*1.0) 
      b = int(255*0.0)
    elif (0.75<=v) and (v<=1.00):
      r = int(255*1.0) 
      g = int(255*(1.0-4*(v-0.75)))
      b = int(255*0.0) 

  if (OPT['col']=='WBR'):
    if (0.0<=v) and (v<0.25):
      r = int(255*(1.0-2*v)) 
      g = int(255*(1.0-2*v))
      b = int(255*1.0) 
    elif (0.25<=v) and (v<0.50):
      r = int(255*(0.5-2*(v-0.25))) 
      g = int(255*(0.5-2*(v-0.25)))
      b = int(255*1.0) 
    elif (0.50<=v) and (v<0.75):
      r = int(255*(4*(v-0.5)))
      g = int(255*0.0)
      b = int(255*1.0) 
    elif (0.75<=v) and (v<=1.00):
      r = int(255*1.0) 
      g = int(255*0.0) 
      b = int(255*(1.0-4*(v-0.75)))




  if (n==0):
    of.write("select %s && temperature<=%.0f\n"%(OPT['tgt'],highT)) 
  elif (n==int(OPT['nc'])):
    of.write("select %s && temperature>=%.0f\n"%(OPT['tgt'],lowT)) 
  else: 
    of.write("select %s && temperature>%.0f && temperature<=%.0f\n"%(OPT['tgt'],lowT,highT)) 
  of.write("color [%d,%d,%d]\n"%(r,g,b))
of.write("select %s\n"%(OPT['tgt']))
of.write("spacefill\n")
of.write("select all\n")

of.close()
