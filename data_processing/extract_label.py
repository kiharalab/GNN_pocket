import os
import numpy as np

def extract_label(label_dir):
    dir = "/net/kihara/home/ykagaya/Share/20220119-SHREC2022/training"
    ids = os.listdir(dir)
    radius = []
    for id in ids:
        label = []
        #print(id)
        #exit(0)
        path = os.path.join(dir+'/'+id,"structure.pqr")
        with open(path,'r') as file:
            lines = file.readlines()
        atom_coordinates = []
        for line in lines:
            line = line.strip('\n')
            line = line.split(' ')
            line= [i for i in line if(len(str(i))!=0)]
            #print(line)
            #print(line[5],line[6],line[7])
            label.append(1 if float(line[8])>0 else 0)
        np.save(os.path.join(label_dir,id+'.npy'),np.array(label))