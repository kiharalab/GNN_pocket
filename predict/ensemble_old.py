from audioop import avg
import os
import pickle
from matplotlib.pyplot import axis
import numpy as np
from sklearn.metrics import accuracy_score,precision_score,recall_score,f1_score

def ensemble_results(id):
    # label_dir  = './label'
    model_3_dir = '../prediction/model_3'
    model_5_dir = '../prediction/model_5'
    model_7_dir = '../prediction/model_7'
    model_10_dir = '../prediction/model_10'
    # Acc = 0
    # Precision = 0
    # Recall = 0
    # F1 = 0
    # for id in ids:
        #label = np.load(os.path.join(label_dir, id+'.npy')).reshape(-1,1)
    model_3_pred =  np.load(os.path.join(model_3_dir,id+'.npy')).reshape(-1,1)
    model_5_pred =  np.load(os.path.join(model_5_dir,id+'.npy')).reshape(-1,1)
    model_7_pred =  np.load(os.path.join(model_7_dir,id+'.npy')).reshape(-1,1)
    model_10_pred =  np.load(os.path.join(model_10_dir,id+'.npy')).reshape(-1,1)
    
    ensemble_pred = np.concatenate((model_3_pred,model_5_pred,model_7_pred,model_10_pred),axis=1)

    #if params['ensemble_mode']==0: #union
    #union_pred = np.sum(ensemble_pred,axis=1)
    #union_pred = np.where(union_pred>=0.5,union_pred,0)
    #elif params['ensemble_mode']==1: #average:
    avg_pred = np.sum(ensemble_pred,axis=1)/4
    #avg_pred = np.where(avg_pred>=0.5,avg_pred,0)
    #elif params['ensemble_mode']==2: #voting
    #vote_pred = np.where(ensemble_pred>=0.5,1,0)
    #vote_pred = np.sum(ensemble_pred,axis=1)
    #vote_pred = np.where(vote_pred>=2.0,vote_pred,0)
    #else:
    #weighted_pred = ensemble_pred.dot([0.6,0.4,0])
    #weighted_pred = np.sum(weighted_pred,axis=1)
    #weighted_pred = np.where(weighted_pred>=0.5,weighted_pred,0)
    return avg_pred
        #print(ensemble_pred.shape,label.shape)
    #     Acc += accuracy_score(label,ensemble_pred)
    #     Precision +=precision_score(label,ensemble_pred)
    #     Recall += recall_score(label,ensemble_pred)
    #     F1 +=f1_score(label,ensemble_pred)
    # print(Acc/len(ids),Precision/len(ids),Recall/len(ids),F1/len(ids))
            
        #else: #reweight
def ensemble_old(params):
    feat_path = os.path.join(params['train_odir'],'feat')
    #dirs = os.listdir(feat_path)
    with open('valid_index.pkl','rb') as file:
        ids = pickle.load(file)

    avg_dir = "../prediction/avg_pred"
    avg_pdb_dir = "../prediction/avg_pdb_pred"
    if not os.path.exists(avg_dir):
        os.mkdir(avg_dir)
    if not os.path.exists(avg_pdb_dir):
        os.mkdir(avg_pdb_dir)
    for id in ids:
        pdb = '../dataset/train_processed/pdb/' + id+'.pdb'
        pred_lines = '' 
        with open(pdb,'r') as file:
            lines = file.readlines()
        avg_pred = ensemble_results(id)
        avg_path = os.path.join(avg_dir,id+'.npy')
        np.save(avg_path,avg_pred)
        pred_lines = '' 
        for idx in range(len(avg_pred)):
            if avg_pred[idx]>=0.3:
                pred_lines +=lines[idx].strip('\n') + ' ' + str(avg_pred[idx]) + '\n'
        ofile = os.path.join(avg_pdb_dir,id+'.pdb')
        with open(ofile,'w') as file:
            file.write(pred_lines)




# for id in index:
#     union_pred,avg_pred,vote_pred,weighted_pred = ensemble(id)
#     union_path = os.path.join(union_dir,id+'.npy')
#     union_pdb_path = os.path.join(union_pdb_dir,id+'.pdb')
#     np.save(union_path,union_pred)
#     pdb = '../dataset/pdb/' + str(id)+'.pdb'
#     pred_lines = '' 
#     with open(pdb,'r') as file:
#         lines = file.readlines()
#     for idx in range(len(union_pred)):
#         if union_pred[idx]>=0.5:
#             pred_lines +=lines[idx].strip('\n') + ' ' + str(union_pred[idx]) + '\n'
#     ofile = os.path.join(union_pdb_dir,id+'.pdb')
#     with open(ofile,'w') as file:
#         file.write(pred_lines)
    
#     avg_path = os.path.join(avg_dir,id+'.npy')
#     np.save(avg_path,avg_pred)
    

#     pred_lines = '' 
#     for idx in range(len(avg_pred)):
#         if avg_pred[idx]>=0.5:
#             pred_lines +=lines[idx].strip('\n') + ' ' + str(avg_pred[idx]) + '\n'
#     ofile = os.path.join(avg_pdb_dir,id+'.pdb')
#     with open(ofile,'w') as file:
#         file.write(pred_lines)
    
#     vote_path = os.path.join(vote_dir,id+'.npy')
#     np.save(vote_path,vote_pred)
#     pred_lines = '' 
#     for idx in range(len(vote_pred)):
#         if vote_pred[idx]>=0.5:
#             pred_lines +=lines[idx].strip('\n') + ' ' + str(vote_pred[idx]) + '\n'
#     ofile = os.path.join(vote_pdb_dir,id+'.pdb')
#     with open(ofile,'w') as file:
#         file.write(pred_lines)
    
#     weight_path = os.path.join(weighted_dir,id+'.npy')
#     np.save(weight_path,weighted_pred)
#     pred_lines = '' 
#     for idx in range(len(weighted_pred)):
#         if weighted_pred[idx]>=0.5:
#             pred_lines +=lines[idx].strip('\n') + ' ' + str(weighted_pred[idx]) + '\n'
#     ofile = os.path.join(weighted_pdb_dir,id+'.pdb')
#     with open(ofile,'w') as file:
#         file.write(pred_lines)
    



        

