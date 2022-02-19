import os
import numpy as np
from scipy.spatial import distance

def extract_edges(single_dir,edges_dir,water_bond):
    dir = "/net/kihara/home/ykagaya/Share/20220119-SHREC2022/training"
    if water_bond:
        edges_dir +='_water'
        print(edges_dir)
    if not os.path.exists(edges_dir):
        os.mkdir(edges_dir)
    ids = os.listdir(dir)
    import pickle
    # with open("/net/kihara/home/zhang038/Shrec2022/remain.pkl",'rb') as file:
    #     ids = list(pickle.load(file))
    for id in ids:
        #print(id)
        radius = []
        path = os.path.join(dir+'/'+id,"structure.pqr")
        radius = []
        #path = "/Users/yuanyuanzhang/Documents/kihara_lab/SHREC2022/dataset/2/structure.pqr"
        with open(path,'r') as file:
            lines = file.readlines()
        atom_coordinates = []
        for line in lines:
            line = line.strip('\n')
            line = line.split(' ')
            line= [i for i in line if(len(str(i))!=0)]
            #print(line)
            #print(line[5],line[6],line[7])
            x = float(line[5])
            y = float(line[6])
            z = float(line[7])
            radius.append(float(line[-1]))
            atom_coordinates.append(np.array([x,y,z]))
        atom_coordinates = np.array(atom_coordinates)
        atom_adjacent = distance.cdist(atom_coordinates,atom_coordinates,'euclidean')
        for i in range(len(atom_adjacent)):
            for j in range(len(atom_adjacent[0])):
                if water_bond:
                    atom_adjacent[i][j] -= (radius[i]+ radius[j])/2+2.8
                else:
                    atom_adjacent[i][j] -= (radius[i]+ radius[j])/2
                if atom_adjacent[i][j]>0:
                    atom_adjacent[i][j]=0
                else:
                    atom_adjacent[i][j]=1
        np.save(os.path.join(edges_dir,id+".npy"),atom_adjacent)
        
        #np.save(os.path.join(single_dir,"atom_feature/"+id+".npy"),np.sum(atom_adjacent,axis=1))
        # print(np.sum(atom_adjacent,axis=1))
        # print(np.sum(atom_adjacent,axis=0))
#extract_edges("../features/", "../edges", True)
