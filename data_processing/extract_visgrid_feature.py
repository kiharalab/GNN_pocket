import os
import pickle
import numpy as np
from scipy.spatial import distance

def extract_visgrid_feature(single_dir):
    dir = "/net/kihara/home/zhang038/Shrec2022/dataset/pdb"
    bin_path = "/net/kihara/home/zhang038/Shrec2022/visgrid-alltaggedatoms"
    ids = os.listdir(dir)
    for id in ids:
        path = os.path.join(dir,id)
        with open(path,'r') as file:
            rows = len(file.readlines())
        vector = np.zeros(rows)
        output_lines = os.popen(bin_path + " "+path)
        opath = os.path.join(single_dir,'visgrid')
        if not os.path.exists(opath):
            os.mkdir(opath)
        opath = os.path.join(opath,id[:-4]+'.npy')
        for line in output_lines:
            line = line.split(' ')
            line = list(filter(None, line))
            vector[int(line[1])-1] = 1
        np.save(opath,vector)

extract_visgrid_feature("../single_feats")