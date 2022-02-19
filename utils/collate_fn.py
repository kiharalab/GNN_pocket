import numpy as np
import torch

def collate_fn(batch):
    max_natoms = max([len(item[0]) for item in batch if item is not None])
    nfeat = len(batch[0][0][0])
    F = np.zeros((len(batch), max_natoms, nfeat))
    A = np.zeros((len(batch), max_natoms, max_natoms))
    # A2 = np.zeros((len(batch), max_natoms, max_natoms))
    V = np.zeros((len(batch), max_natoms))
    #keys = []
    #print(F.shape, A.shape,V.shape)
    Atoms_Number=[]
    for i in range(len(batch)):
        natom = len(batch[i][0])
        F[i, :natom] = batch[i][0]
        A[i, :natom, :natom] = batch[i][1]
        #A2[i, :natom, :natom] = batch[i]['A2']
        V[i, :natom] = batch[i][2]
        #keys.append(batch[i]['key'])
        Atoms_Number.append(natom)
    F = torch.from_numpy(F).float()
    A = torch.from_numpy(A).float()
    #A2 = torch.from_numpy(A2).float()
    V = torch.from_numpy(V).float()
    Atoms_Number = np.array(Atoms_Number)
    Atoms_Number=torch.from_numpy(Atoms_Number).int()

    return F, A, V,Atoms_Number#, keys