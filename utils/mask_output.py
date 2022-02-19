import torch
import numpy as np
import random
def masking(output,label):#size:batch * atom_number
    output = output.detach().cpu().numpy()
    label = label.detach().cpu().numpy()
    count_pos = len(np.nonzero(label))
    a = np.zeros_like(label)
    mask = a+label
    pos_zeros =[list(i) for i in np.argwhere(label == 0)]  #index
    #print(pos_zeros)
    remain_zeros = random.sample(pos_zeros,count_pos) #index
    for index in remain_zeros:
        mask[index[0],index[1]] = 1
    return torch.FloatTensor(mask)
    




