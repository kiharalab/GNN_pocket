import os
import numpy as np
from scipy.spatial import distance

def extract_atoms(single_dir):
    dir = "/net/kihara/home/ykagaya/Share/20220119-SHREC2022/training"
    ids = os.listdir(dir)
    import pickle
    # with open("/net/kihara/home/zhang038/Shrec2022/remain.pkl",'rb') as file:
    #     ids = list(pickle.load(file))
    with open("/net/kihara/home/zhang038/Shrec2022/atom_types.pkl",'rb') as file:
        all_types = pickle.load(file)
    for id in ids:
        #print(id)
        path = os.path.join(dir+'/'+id,"structure.pqr")
        with open(path,'r') as file:
            lines = file.readlines()
        atom_types_0 = []
        atom_types_1 = []
        for line in lines:
            atom_type_1 = np.zeros(13)
            atom_type_0 = np.zeros(1)
            line = line.strip('\n')
            line = line.split(' ')
            line= [i for i in line if(len(str(i))!=0)]

            atom_type_0[0] = all_types.index(line[-1])
            atom_types_0.append(atom_type_0)

            atom_type_1[all_types.index(line[-1])] = 1
            atom_types_1.append(atom_type_1)
        atom_types_0 = np.array(atom_types_0)
        atom_types_1 = np.array(atom_types_1)

        odir_0 = os.path.join(single_dir,"atom_types_0")
        odir_1 = os.path.join(single_dir,"atom_types_1")
        if not os.path.exists(odir_0):
            os.mkdir(odir_0)
        if not os.path.exists(odir_1):
            os.mkdir(odir_1)
        np.save(os.path.join(odir_0,id+".npy"),atom_types_0)
        np.save(os.path.join(odir_1,id+".npy"),atom_types_1)
        #print(np.sum(atom_adjacent,axis=1))
        print(atom_types_0.shape)
        print(atom_types_1.shape)
        #visualize_atom.visualize_graph("../visulaization",str(id),atom_coordinates,pairs)
