import os
import numpy as np
from scipy.spatial import distance

def extract_face_feature(single_dir):
    dir = "/net/kihara/home/ykagaya/Share/20220119-SHREC2022/training"
    ids = os.listdir(dir)
    radius = []
    for id in ids:
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
                face_start = i
                break
            else:
                vertex_coordinates.append(np.array([float(line[0]),float(line[1]),float(line[2])]))
        faces = []
        for i in range(face_start,len(lines)):
            line = lines[i].strip('\n')
            line = line.split(' ')
            line = [i for i in line if(len(str(i))!=0)]
            faces.append(np.array([int(line[1]),int(line[2]),int(line[3])]))
        faces = np.array(faces)
        vertex_coordinates = np.array(vertex_coordinates)
        #print(len(vertex_coordinates))
        vertex_distance = distance.cdist(atom_coordinates,vertex_coordinates,'euclidean')
        vertex_map = np.zeros([len(vertex_distance),len(vertex_distance[0])])
        face_map = {}
        for index,i in enumerate(np.argmin(vertex_distance,axis=0)):
            vertex_map[i][index]=1
            # if index in face_map.keys():
            #     value = face_map[index]
            #     value.append(i)
            #     face_map[index] = value
            # else:
            face_map[index] = i
        #print(face_map)
        #exit(0)
        #vertex_feature = np.sum(vertex_map, axis=1)
        #print(vertex_feature,len(vertex_feature))
        #print(face_map)
        face_feature_dic = {}
        for face in faces:
            if face_map[face[0]]==face_map[face[1]] and face_map[face[0]]==face_map[face[2]] and face_map[face[2]]==face_map[face[1]]:
                if face_map[face[0]] in face_feature_dic.keys():
                    face_feature_dic[face_map[face[0]]] += 1
                else:
                    face_feature_dic[face_map[face[0]]] = 1
            # for k,v in face_map.items():
            #     if len(set(face) | set(v)) == len(set(v)):
            #         face_feature_dic[i] += 1s
        #print(face_feature_dic)

        print(len(list(face_feature_dic.keys())))
        face_feature = np.zeros(len(atom_coordinates))
        for key in face_feature_dic.keys():
            face_feature[key] = face_feature_dic[key]
        #print(np.sum(np.array(list(face_feature_dic.values()))))
        # print(face_feature)
        # print(len(face_feature))
        np.save(os.path.join(single_dir,"face_feature/"+id+".npy"), face_feature)
        #np.save("../atom_distance/"+id+".npy", atom_distance)
        #np.save("../vertex_feature/"+id+".npy",vertex_feature)
        #print(atom_distance,len(atom_distance[0]),len(atom_distance[1]))
        #print(len(vertex_distance),len(vertex_distance[0]))
