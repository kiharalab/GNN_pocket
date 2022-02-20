import torch
import numpy as np
import os
from scipy.sparse import diags
from utils.normalize import normalize
from model.model import GCN

def predict(params):
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    for i in [3,5,7,10]:
        model_path = "./ensemble_models/model_"+str(i)+".ckpt"
        if i==5:
            nfeat = 1
        elif i==7:
            nfeat=2
        else:
            nfeat=3
        model = GCN(nfeat,dropout=0.0)
        model.load_state_dict(torch.load(model_path))
        model = model.to(device)
        model.eval()
        if i == 5:
            edge_path = os.path.join(params['test_odir'],'edges_water')
        else:
            edge_path = os.path.join(params['test_odir'],'edges')
        feat_path = os.path.join(params['test_odir'],'feat')
        vis_feat_path = os.path.join(params['test_odir'],'nor_vis8')
    
        dirs = os.listdir(feat_path)

        for id in dirs:
            adj=np.load(os.path.join(edge_path,id))
            #print(adj)
            adj = normalize(adj)
            feat = np.load(os.path.join(feat_path,id))
            if i==3:
                vis_feat = np.load(os.path.join(vis_feat_path,id)).reshape(-1,1)
                new_feat = np.concatenate((feat[:,:2],vis_feat),axis=1).reshape(-1,3)
                feat = normalize(new_feat)
            if i == 5:
                feat = feat[:,0].reshape(-1,1)
            elif i==7:
                feat = normalize(feat[:,:2])
            else:
                feat = normalize(feat)
        
            feat = torch.FloatTensor(feat).unsqueeze(0).to(device)
            adj = torch.FloatTensor(adj).unsqueeze(0).to(device)

            output = model(feat,adj)
            output = output.detach().cpu().numpy().reshape(-1)

            pdb = './dataset/test_processed/pdb/' + str(id[:-4])+'.pdb'
            pred_lines = '' 
            with open(pdb,'r') as file:
                lines = file.readlines()
            for idx in range(len(output)):
                if output[idx]>=0.5:
                    pred_lines +=lines[idx].strip('\n') + ' ' + str(output[idx]) + '\n'
            odir = "./test_prediction/pdb_model_"+str(i)
            odir = os.path.abspath(odir)
            if not os.path.exists(odir):
                os.makedirs(odir)
            ofile = os.path.join(odir,id[:-4]+'.pdb')
            with open(ofile,'w') as file:
                file.write(pred_lines)
            
            odir = "./test_prediction/model_"+str(i)
            odir = os.path.abspath(odir)
            if not os.path.exists(odir):
                os.makedirs(odir)
            ofile = os.path.join(odir,id)
            np.save(ofile,output)


