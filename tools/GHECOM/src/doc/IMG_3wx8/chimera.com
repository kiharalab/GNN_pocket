#2020/07/21
#chimera pdb3wz8.ent 3wz8_P.pdb 3wz8_P_1.map 3wz8_P_2.map 3wz8_P_3.map 3wz8_M_rec.pdb 3wz8_M.pdb 3wz8_M.map 
# #0:pdb3wz8.ent 
# #1:3wz8_P.pdb 
# #2:3wz8_P_1.map 
# #3:3wz8_P_2.map 
# #4:3wz8_P_3.map 
# #5:3wz8_M_rec.pdb 
# #6:3wz8_M.pdb 
# #7:3wz8_M.map
reset
windowsize 1000 800
set bgcolor white
delete :HOH
~ribbon
~display
display protein
repr stick #0
repr sphere #5
display :GRD
repr sphere :GRD
vdwdef -1.2 :GRD
repr sphere :GRD & #6
rangecol bfactor,a 0 blue 5 green 10 red :#5 & protein
rangecol bfactor,a 3 red 5 green 10 blue :#6 & :GRD
color gray #0
rainbow model :GRD & #1
display :IXV
color byelement :IXV
repr sphere :IXV
volume #2 style surface level 0.5 color red 
volume #3 style surface level 0.5 color green 
volume #4 style surface level 0.5 color blue 
volume #7 style surface level 0.01 color green 


