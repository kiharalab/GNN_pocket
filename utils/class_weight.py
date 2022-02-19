import os
import numpy as np

def class_weight(params):
    weights = []
    label_dir = os.path.abspath(params['label'])
    dirs = os.listdir(label_dir)
    for dir in dirs:
        class_1 = 0
        count = 0
        path = os.path.join(label_dir,dir)
        with open(path,'rb') as file:
            label = np.load(file)
            count = len(label)
            class_1 = np.sum(label)/count
        weights.append(class_1)
    return weights
