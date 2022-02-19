import os
import numpy as np
from scipy.spatial import distance

def extract_vertex_feature(single_dir):
    dir = "/net/kihara/home/ykagaya/Share/20220119-SHREC2022/training"
    ids = os.listdir(dir)
    radius = []
    for id in ids:
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
            x = float(line[5])
            y = float(line[6])
            z = float(line[7])
            radius.append(float(line[-1]))
            atom_coordinates.append(np.array([x,y,z]))
        atom_coordinates = np.array(atom_coordinates)
        atom_distance = distance.cdist(atom_coordinates,atom_coordinates,'euclidean')
        #print(len(atom_distance))
        path = os.path.join(dir+'/'+id,"triangulatedSurf.off")
        vertex_distance = []
        vertex_coordinates = []
        with open(path,'r') as file:
            lines = file.readlines()
        for i in range(4,len(lines)):
            line = lines[i].strip('\n')
            line = line.split(' ')
            line = [i for i in line if(len(str(i))!=0)]
            if len(line)>3:
                break
            else:
                vertex_coordinates.append(np.array([float(line[0]),float(line[1]),float(line[2])]))
        vertex_coordinates = np.array(vertex_coordinates)
        #print(len(vertex_coordinates))
        vertex_distance = distance.cdist(atom_coordinates,vertex_coordinates,'euclidean')
        vertex_map = np.zeros([len(vertex_distance),len(vertex_distance[0])])
        print(vertex_map.shape)
        for index,i in enumerate(np.argmin(vertex_distance,axis=0)):
            vertex_map[i][index]=1
        vertex_feature = np.sum(vertex_map, axis=1)
        print(vertex_feature,len(vertex_feature))

        #np.save("../atom_distance/"+id+".npy", atom_distance)
        np.save(os.path.join(single_dir,"vertex_feature/"+id+".npy"),vertex_feature)
        #print(atom_distance,len(atom_distance[0]),len(atom_distance[1]))
        #print(len(vertex_distance),len(vertex_distance[0]))
