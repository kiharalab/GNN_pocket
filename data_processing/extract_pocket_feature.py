import os
import pickle
import numpy as np
from scipy.spatial import distance

def extract_pocket_feature(single_dir):
    visgrid_path = "/fast-scratch/zhang038/Shrec2022/VisGrid_tool/VisGrid"
    dir = "/fast-scratch/zhang038/Shrec2022/dataset/train_processed/pdb"
    ids = os.listdir(dir)
    for id in ids:
        # if os.path.exists(os.path.join(single_dir,"visGrid_feature_8/"+id[:-4]+".npy")):
        #     continue
        path = os.path.join(dir,id)
        output_lines = os.popen(visgrid_path+ ' '+ '-v'+ ' ' +path).readlines()
        #print(output_lines)
        #exit(0)
        voxels = []
        for i in range(2,len(output_lines)):
            if "Voxels:" in output_lines[i]:
                continue
            line = output_lines[i].replace('\n','').replace('\t','')
            if line=='':
                continue
            line = np.array([float(i) for i in line.split(' ')])
            voxels.append(line)
        atoms = []
        with open(path,'r') as file:
            lines = file.readlines()
        for line in lines:
            line = line.strip('\n')
            line = line.replace('-',' ')
            line = line.split(' ')
            line= [i for i in line if(len(str(i))!=0)]
            x = float(line[5])
            #print(id)
            y = float(line[6])
            
            z = float(line[7])
            atoms.append(np.array([x,y,z]))
        # print(len(voxels[0]),len(atoms[0]))
        # exit(0)
        voxels = np.array(voxels)
        atoms = np.array(atoms)
        distances = distance.cdist(atoms,voxels)
        # print(distances.shape)
        # exit(0)
        distances = np.where(distances>8,0,1)

        # voxel_map = np.zeros([len(distances),len(distances[0])])
        # print(voxel_map.shape)
        # for index,i in enumerate(np.argmin(distances,axis=0)):
        #     voxel_map[i][index]=1

        voxel_feature = np.sum(distances, axis=1)
        max_value = np.sum(voxel_feature)
        voxel_feature = voxel_feature/(max_value+1e-9)
        # for i in voxel_feature:
        #     if i>0:
        #         print(i)
        #print(voxel_feature,len(voxel_feature))
        np.save(os.path.join(single_dir,"nor_visGrid_feature_8/"+id[:-4]+".npy"),voxel_feature)

extract_pocket_feature("/fast-scratch/zhang038/Shrec2022/GCN_pocket/single_feats/")
    
        

