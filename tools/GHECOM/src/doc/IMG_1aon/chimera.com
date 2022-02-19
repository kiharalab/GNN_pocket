#2020/07/21
#chimera pdb1aon.ent 1aon_V.pdb 1aon_V.map 1aon_P.pdb 1aon_P.map 1aon_CP.pdb 1aon_CP.map
# #0 pdb1aon.ent 
# #1 1aon_V.pdb 
# #2 1aon_V.map 
# #3 1aon_P.pdb 
# #4 1aon_P.map
# #5 1aon_CP.pdb 
# #6 1aon_CP.map
reset
windowsize 800 1000
turn x 90
set bgcolor white
delete :HOH
~ribbon
~display
ribbon #0
color white #0
display :GRD
repr sphere :GRD
vdwdef 2.0 :GRD
color blue #1
color green #3
color red   #5
volume #2 style surface level 0.5 color blue 
volume #4 style surface level 0.5 color green 
volume #6 style surface level 0.5 color red

