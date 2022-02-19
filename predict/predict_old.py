import torch
import numpy as np
import os
import numpy  as np
from scipy.sparse import diags
from utils.normalize import normalize
from model.model import GCN
from torch.utils.data import DataLoader,random_split
from data_processing.dataset import MyDataset
from utils.collate_fn import collate_fn
from model.dice_loss import DiceLoss
import time
import pickle

def predict_old(params):
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    #nfeat = params['feat_num']
    with open("./valid_index.pkl",'rb') as file:
    	ids = pickle.load(file)

    for i in [3,5,7,10]:
        model_path = "./ensemble_models/model_"+str(i)+".ckpt"
        if i==5:
            nfeat = 1
        elif i==7:
            nfeat=2
        else:
            nfeat=3
        
        model = GCN(nfeat,dropout=0.3)
        model.load_state_dict(torch.load(model_path))
        model = model.to(device)
        model.eval()
        if i == 5:
            edge_path = os.path.join(params['train_odir'],'edges_water')
        else:
            edge_path = os.path.join(params['train_odir'],'edges')
        feat_path = os.path.join(params['train_odir'],'feat')
        vis_feat_path = os.path.join('/fast-scratch/zhang038/Shrec2022/GCN_pocket/single_feats/','nor_visGrid_feature_8')
    
        dirs = os.listdir(feat_path)

        for id in ids:
            adj=np.load(os.path.join(edge_path,id+'.npy'))
            #print(adj)
            adj = normalize(adj)
            feat = np.load(os.path.join(feat_path,id+'.npy'))
            if i==3:
                vis_feat = np.load(os.path.join(vis_feat_path,id+'.npy')).reshape(-1,1)
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
        #print(output)
        #pred = np.array([1 if i >=0.5 else 0 for i in output]).astype(int)
            pdb = '../dataset/train_processed/pdb/' + id+'.pdb'
            pred_lines = '' 
            with open(pdb,'r') as file:
                lines = file.readlines()
            for idx in range(len(output)):
                if output[idx]>=0.5:
                    pred_lines +=lines[idx].strip('\n') + ' ' + str(output[idx]) + '\n'
            odir = "../prediction/pdb_model_"+str(i)
            odir = os.path.abspath(odir)
            if not os.path.exists(odir):
                os.makedirs(odir)
            ofile = os.path.join(odir,id+'.pdb')
            with open(ofile,'w') as file:
                file.write(pred_lines)
            
            odir = "../prediction/model_"+str(i)
            odir = os.path.abspath(odir)
            if not os.path.exists(odir):
                os.makedirs(odir)
            ofile = os.path.join(odir,id)
            np.save(ofile,output)
            
def metrics(output, labels):
    # TP = 0
    # acc_correct = 0
    # all = 0
    # all_pred = 0
    # all_label = 0
    
    acc_correct = 0
    pred = np.array([1 if i >=0.5 else 0 for i in output]).astype(int)
    label = np.array(labels).astype(int)
        #acc_correct = np.sum(np.equal(pred,label))
    for i in range(len(label)):
        if label[i]==pred[i]:
            acc_correct +=1
    all = len(label)
    TP = np.sum((pred & label))
    all_pred = np.sum(pred)
    all_label = np.sum(label)
        #print(TP,all_pred,all_label)
        #preds.append(torch.LongTensor([1 if i >=0.5 else 0 for i in batch]))
    acc = acc_correct/(all+1e-9)
    pre = TP/(all_pred +1e-9)
    rec = TP/(all_label+1e-9)
    f1 = 2*pre*rec/(pre+rec)
    return acc,pre,rec,f1

# device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
# #model_path = "./model_test_dict/88_1.1235955056179774e-05.ckpt"
# model_path = 
# with open('valid_index.pkl','rb') as file:
#     dirs = pickle.load(file)
# # dir = "./features"
# # dirs = os.listdir(dir)
# model = GCN(feat_num=1,atom_type=0,dropout=0.3)
# model.load_state_dict(torch.load(model_path))
# model = model.to(device)
# model.eval()
# atoms = 0
# Acc =0
# Pre =0
# Rec =0
# F1 =0
# for i in dirs:
#     label_path = './label/' + i+'.npy'
#     feat_path = './features/'+i+'.npy'
#     edge_path = './edges_water/'+i+'.npy'
#     import numpy as np
#     label = np.load(label_path).astype(int)
#     feat = np.load(feat_path)[:,-2:-1]
#     adj = normalize(np.load(edge_path))

#     pdb = '../dataset/pdb/' + str(i)+'.pdb'
#     odir = "./test_pocket_score"
#     if not os.path.exists(odir):
#         os.mkdir(odir)
#     ofile = './test_pocket_score/pred_' + str(i)+'.pdb'
    

#     feat = torch.FloatTensor(feat).unsqueeze(0)
#     feat = feat.to(device)
#     adj = torch.FloatTensor(adj).unsqueeze(0)
#     adj = adj.to(device)
    

#     #output = torch.Tensor([0.])
#     output = model(feat,adj)
#     output = output.detach().cpu().numpy().reshape(-1)
#     print(output)
#     pred = np.array([1 if i >=0.5 else 0 for i in output]).astype(int)

#     pred_lines = '' 
#     with open(pdb,'r') as file:
#         lines = file.readlines()
#     for idx in range(len(pred)):
#         if pred[idx]>0:
#             print(idx)
#             atoms+=1
#             pred_lines +=lines[idx].strip('\n') + ' ' + str(output[idx]) + '\n'
#     # acc,pre,rec,f1 = metrics(pred,label)
#     # Acc += acc
#     # Pre += pre
#     # Rec += rec
#     # F1 += f1
#     # pred_lines += str(acc)+','+str(pre)+','+str(rec)+','+str(f1)
#     with open(ofile,'w') as file:
#         file.write(pred_lines)
# print(atoms)
# print(Acc / len(dirs),Pre / len(dirs),Rec / len(dirs),F1 / len(dirs))

# def predict(params):
#     torch.manual_seed(34)
#     device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
#     feat_dir = os.path.abspath(params['feat'])
#     label_dir = os.path.abspath(params['label'])
#     edge_dir = os.path.abspath(params['edge'])

#     train_data = MyDataset(params['vis'],params['feat_num'],label_dir,feat_dir,edge_dir)

#     train_size = int(len(train_data) * 0.8)
#     valid_size = int(len(train_data) * 0.2)
#     train_dataset, valid_dataset= random_split(train_data, [train_size, valid_size])
#     #train_loader = DataLoader(train_dataset, batch_size=params['batch_size'], shuffle=False, num_workers=1,collate_fn=collate_fn)
#     validate_loader = DataLoader(valid_dataset, batch_size=params['batch_size'], shuffle=False, num_workers=1,collate_fn=collate_fn)


#     model = GCN(feat_num=params['feat_num'], atom_type=params['atom_type'],dropout=params['dropout'])
#     model.load_state_dict(torch.load('88_1.1235955056179774e-05.ckpt'))
#     model = model.to(device=device)

#     for batch_idx,sample in enumerate(validate_loader):
#         feat,adj,label,natoms= sample

#         dice_loss = DiceLoss(natoms=natoms)
#         t = time.time()
#         adj = torch.FloatTensor(adj).to(device)
#         feat = torch.FloatTensor(feat).to(device)
#         label = torch.FloatTensor(label).to(device)
#         size = label.shape[0]

#         with torch.no_grad():
#             model.eval()
#             output = model(feat, adj)
#             output = output.reshape(size,-1)
#             label = label.reshape(size,-1)
#             #acc = accuracy(output, label)
#             acc,precision,recall,f1 = metrics(output,label)

#             batch_loss.append(loss.item()),batch_acc.append(acc),batch_precision.append(precision),batch_recall.append(recall),batch_f1.append(f1)
#         print('Val Loss: {:.3f}, Val Acc: {:.3f},,Val Precision: {:.3f}, Val Recall: {:.3f}, Val F1: {:.3f}'.format(loss.item(), acc,precision,recall,f1))
#     val_loss.append(np.sum(np.array(batch_loss))/len(validate_loader))
#     val_acc.append(np.sum(np.array(batch_acc))/len(validate_loader))
#     val_precision.append(np.sum(np.array(batch_precision))/len(validate_loader))
#     val_recall.append(np.sum(np.array(batch_recall))/len(validate_loader))
#     val_f1.append(np.sum(np.array(batch_f1))/len(validate_loader))
#     print('Epochs: {}, Train Loss: {:.3f}, Train Acc: {:.3f},Train Precision: {:.3f},Train Recall: {:.3f},Train F1: {:.3f}, Validation Loss: {:.3f}, \
#     Validation Acc: {:.3f},Validation Precision: {:.3f},Validation Recall: {:.3f},Validation F1: {:.3f}'\
#         .format(i, train_loss[-1], train_acc[-1],train_precision[-1],train_recall[-1],train_f1[-1], \
#             val_loss[-1], val_acc[-1],val_precision[-1],val_recall[-1],val_f1[-1]))



