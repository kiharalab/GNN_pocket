import os
import numpy as np
from scipy.spatial.distance import cdist
from utils.pqr2pdb import pqr2pdb
from scipy.spatial import cKDTree

def extract_dataset(input_dir,output_dir):
    label_dir = os.path.join(output_dir,'label')
    feat_dir = os.path.join(output_dir,'feat')
    nor_vis_8_dir = os.path.join(output_dir,'nor_vis8')
    edge_dir = os.path.join(output_dir,'edges')
    edge_water_dir = os.path.join(output_dir,'edges_water')
    if not os.path.exists(label_dir):
        os.mkdir(label_dir)
    if not os.path.exists(feat_dir):
        os.mkdir(feat_dir)
    if not os.path.exists(edge_dir):
        os.mkdir(edge_dir)
    if not os.path.exists(edge_water_dir):
        os.mkdir(edge_water_dir) 
    if not os.path.exists(nor_vis_8_dir):
        os.mkdir(nor_vis_8_dir)
    ids = os.listdir(input_dir)

    for id in ids:
#         if os.path.exissts(os.path.join(feat_dir,str(id)+'.npy')):
#             continue
        path = os.path.join(input_dir,str(id)+'/structure.pqr')
        with open(path,'r') as file:
            lines = file.readlines()
        coords = np.zeros((len(lines),3))
        label = np.zeros(len(lines))
        radius = np.zeros(len(lines))
        for idx,line in enumerate(lines):
            line = line.strip('\n')
            line = line.split(' ')
            line = [i for i in line if len(i)!=0]
            coords[idx][0] = float(line[5])
            coords[idx][1] = float(line[6])
            coords[idx][2] = float(line[7])
            label[idx] = float(line[8])
            radius[idx] = float(line[9])
        np.save(os.path.join(label_dir,str(id)+'.npy'),label)
        #edges = np.zeros((len(lines),len(lines)))
        radius = radius.reshape(-1,1)
        #extract edges, with different mode
        distance = cdist(coords,coords)
        radius_sum = cdist(radius,radius,lambda u, v: (u+v)/2)
        edges = radius_sum - distance
        edges = np.where(edges<0,0,1)
        np.save(os.path.join(edge_dir,str(id)+'.npy'),edges)

        #if params['water'] == 1:
        radius_sum = cdist(radius,radius,lambda u, v: (u+v)/2+2.8)
        edges = radius_sum - distance
        edges = np.where(edges<0,0,1)
        np.save(os.path.join(edge_water_dir,str(id)+'.npy'),edges)
        
        # convert pqr to pdb
        pdb_write_dir = os.path.join(output_dir,'pdb')
        if not os.path.exists(pdb_write_dir):
            os.mkdir(pdb_write_dir)
        pdb_write_path = os.path.join(pdb_write_dir,str(id)+'.pdb')
        pqr2pdb(path,pdb_write_path)

        # running visgrid and prepare visgrid features
        tools_dir = './tools'
        vis_exe = os.path.join(tools_dir,'VisGrid/visgrid-alltaggedatoms')
        os.system("chmod 777 "+ vis_exe)
        vis_pred = os.popen(vis_exe+ ' '+ pdb_write_path)
        index = []
        for line in vis_pred:
            line = line.split(' ')
            line = list(filter(None, line))
            line = [i for i in line if len(i)!=0]
            index.append(int(line[1])-1)
        vis_atom_feat = np.zeros(len(coords))
        for i in index:
            vis_atom_feat[i] = 1.0
        vis_atom_feat=vis_atom_feat.reshape(-1,1)

        #running ghecom and prepare ghecom features
        ghecom_exe = os.path.join(tools_dir,'GHECOM/ghecom')
        os.system("chmod 777 "+ ghecom_exe)
        ghecom_odir = os.path.join(output_dir,'ghecom_output')
        if not os.path.exists(ghecom_odir):
            os.mkdir(ghecom_odir)
        ghecom_opath = os.path.join(ghecom_odir, str(id)+'.pdb')
        command = ghecom_exe +" -M P -ipdb "+ pdb_write_path +" -opocpdb " + ghecom_opath
        os.system(command)
        ghecom_coords = []
        with open(ghecom_opath,'r') as file:
            lines = file.readlines()
        for line in lines:
            if "HETATM" in line:
                line = line.strip('\n').split(' ')
                line = [i for i in line if len(i)>0]
                ghecom_coords.append(tuple(float(line[i]) for i in [6,7,8]))
        ghecom_feat = np.zeros(len(coords))
        if len(ghecom_coords)>0:
            atoms_arr = tuple(a for a in coords)
            atoms_tree = cKDTree(atoms_arr)
            (dists, indices) = atoms_tree.query(ghecom_coords)
            count = {}
            for i in indices:
                if i in count.keys():
                    count[i]+=1
                else:
                    count[i]=1
            for key in count.keys():
                ghecom_feat[key] = count[key]
            max_value = ghecom_feat.max()
            ghecom_feat = ghecom_feat* 1/max_value
        ghecom_feat = ghecom_feat.reshape(-1,1)

        #prepare visgrid 8A visibility features

        #save vis_8 feature
        tools_dir = './tools'
        vis_pred_exe = os.path.join(tools_dir,'VisGrid/VisGrid')
        os.system("chmod 777 "+ vis_pred_exe)
        output_lines = os.popen(vis_pred_exe+ ' '+ '-v'+ ' ' +path).readlines()
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

        voxels = np.array(voxels)
        atoms = np.array(atoms)
        distances = cdist(atoms,voxels)
 
        distances = np.where(distances>8,0,1)
        voxel_feature = np.sum(distances, axis=1)
        max_value = np.sum(voxel_feature)
        voxel_feature = voxel_feature/(max_value+1e-9)

        visgird_8A_feat = np.sum(distances, axis=1)
        visgird_8A_feat = visgird_8A_feat.reshape(-1,1)

        #combine all features
        feat = np.concatenate((vis_atom_feat,ghecom_feat,visgird_8A_feat),axis=1).reshape(-1,3)
        #print(feat.shape)
        np.save(os.path.join(feat_dir,str(id)+'.npy'),feat)
        np.save(os.path.join(nor_vis_8_dir,str(id)+'.npy'),voxel_feature)