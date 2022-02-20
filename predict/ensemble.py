import os
import numpy as np


def ensemble_results(id):

    model_3_dir = './test_prediction/model_3'
    model_5_dir = './test_prediction/model_5'
    model_7_dir = './test_prediction/model_7'
    model_10_dir = './test_prediction/model_10'

    model_3_pred =  np.load(os.path.join(model_3_dir,id)).reshape(-1,1)
    model_5_pred =  np.load(os.path.join(model_5_dir,id)).reshape(-1,1)
    model_7_pred =  np.load(os.path.join(model_7_dir,id)).reshape(-1,1)
    model_10_pred =  np.load(os.path.join(model_10_dir,id)).reshape(-1,1)
    
    ensemble_pred = np.concatenate((model_3_pred,model_5_pred,model_7_pred,model_10_pred),axis=1)

    return np.sum(ensemble_pred,axis=1)/4

def ensemble(params):
    feat_path = os.path.join(params['test_odir'],'feat')
    dirs = os.listdir(feat_path)

    avg_dir = "./test_prediction/avg_pred"
    avg_pdb_dir = "./test_prediction/avg_pdb_pred"
    if not os.path.exists(avg_dir):
        os.mkdir(avg_dir)
    if not os.path.exists(avg_pdb_dir):
        os.mkdir(avg_pdb_dir)
    for id in dirs:
        pdb = './dataset/test_processed/pdb/' + str(id[:-4])+'.pdb'
        pred_lines = '' 
        with open(pdb,'r') as file:
            lines = file.readlines()
        avg_pred = ensemble_results(id)
        avg_path = os.path.join(avg_dir,id)
        np.save(avg_path,avg_pred)
        pred_lines = '' 
        for idx in range(len(avg_pred)):
            if avg_pred[idx]>=0.5:
                pred_lines +=lines[idx].strip('\n') + ' ' + str(avg_pred[idx]) + '\n'
        ofile = os.path.join(avg_pdb_dir,id[:-4]+'.pdb')
        with open(ofile,'w') as file:
            file.write(pred_lines)
           



        

