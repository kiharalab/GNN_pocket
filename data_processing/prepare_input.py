import os
import numpy as np

def prepare_input(single_feat_dir, feat_dir,atom_type):
    if atom_type == 0:
        atom_dir  = os.path.join(single_feat_dir,"atom_types_0")
    else:
        atom_dir  = os.path.join(single_feat_dir,"atom_types_1")
    vertex_dir = os.path.join(single_feat_dir,"vertex_feature")
    face_dir  =  os.path.join(single_feat_dir,"face_feature")
    voxel_dir = os.path.join(single_feat_dir,"visGrid_feature")
    visgrid_dir = os.path.join(single_feat_dir,"visgrid")
    ghecom_dir = os.path.join(single_feat_dir,"ghecom_feat")
    ids = os.listdir(atom_dir)
    for id in ids:
        feats = []
        atom_feat_path = os.path.join(atom_dir,id)
        vertex_feat_path = os.path.join(vertex_dir,id)
        face_feat_path = os.path.join(face_dir,id)
        voxel_feat_path = os.path.join(voxel_dir,id)
        vis_feat_path = os.path.join(visgrid_dir,id)
        ghecom_feat_path = os.path.join(ghecom_dir,id)
        with open(atom_feat_path,'rb') as file:
            atom_feat  = np.load(file)
        with open(vertex_feat_path,'rb') as file:
            vertex_feat  = np.load(file)
        with open(face_feat_path,'rb') as file:
            face_feat  = np.load(file)
        with open(voxel_feat_path,'rb') as file:
            voxel_feat  = np.load(file) 
        with open(vis_feat_path,'rb') as file:
            visgrid_feat = np.load(file)
        with open(ghecom_feat_path,'rb') as file:
            ghecom_feat = np.load(file)
        for i in range(len(atom_feat)):
            f = list(atom_feat[i])
            f.append(vertex_feat[i])
            f.append(face_feat[i])
            f.append(voxel_feat[i])
            f.append(visgrid_feat[i])
            f.append(ghecom_feat[i])
            feat = np.array(f)
            #feat.append(atom_feat[i],vertex_feat[i],face_feat[i],voxel_feat[i])
            feats.append(feat)
        output_path = os.path.join(feat_dir,id)
        if not os.path.exists(feat_dir):
            os.mkdir(feat_dir)
        np.save(output_path,feats)


