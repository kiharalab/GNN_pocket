
from matplotlib.pyplot import axis
from torch.utils.data import Dataset, DataLoader
from utils.normalize import normalize
import numpy as np
import os
 
class MyDataset(Dataset):
    def __init__(self, vis,feat_num,label_dir, feature_dir, edge_dir):
        self.label_dir = os.path.join(os.getcwd(),label_dir)
        self.feat_num = feat_num
        self.vis_mode = vis
        self.feature_dir = feature_dir
        self.edge_dir = edge_dir
        self.dirs = self.get_nums()
 
    def __getitem__(self, i):
        index = self.dirs[i]
        #print("i={},index={}".format(i, index))
        feat_path = os.path.join(self.feature_dir,index)
        edge_path = os.path.join(self.edge_dir,index)
        label_path = os.path.join(self.label_dir,index)
        with open(feat_path,'rb') as file:
            feat = np.load(file)
        if self.vis_mode == 8:
            vis_path = os.path.join('./single_feats/visGrid_feature_8',index)
            vis_8_feat = np.load(vis_path).reshape(-1,1)
            feat = np.concatenate((feat[:,-2:].reshape(-1,2),vis_8_feat),axis = 1).reshape(-1,3)
        else:
            feat = feat[:,0-self.feat_num:].reshape(-1,self.feat_num)
        #feat = feat[:,-2].reshape(-1,self.feat_num)
        #feat = np.where(feat>0,1,0)
            # print(feat.shape)
            # exit(0)
        if self.feat_num>1 or self.vis_mode:
            feat = normalize(feat)
        with open(edge_path,'rb') as file:
            adj = np.load(file)
            adj = normalize(adj)
        with open(label_path,'rb') as file:
            label = np.load(file)
        
        return feat,adj,label
 
    def __len__(self):
        return len(self.dirs)
    
    def get_nums(self):
        dirs = os.listdir(self.label_dir)
        dirs.sort()
        return dirs