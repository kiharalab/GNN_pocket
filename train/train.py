from cmath import isnan
import torch
import time
import torch.optim as optim
from torch.optim.lr_scheduler import LambdaLR as LambdaLR
from torch.optim.lr_scheduler import MultiStepLR as MultiStepLR
from torch.optim.lr_scheduler import StepLR as StepLR
from model.model import GCN
from data_processing.prepare_input import prepare_input
from data_processing.dataset import MyDataset
from utils.normalize import normalize
import os
from torch.utils.data import DataLoader,random_split,WeightedRandomSampler
import torch.nn.functional as F
import torch.nn as nn
from utils.class_weight import class_weight
from utils.collate_fn import collate_fn
from  utils.mask_output import masking

from model.focal_loss import BCEFocalLoss
from model.dice_loss import DiceLoss, DiceLoss_atom
import numpy as np

from torch.utils.tensorboard import SummaryWriter
summaryWriter = SummaryWriter("../tensorboard_runs/")

# model = GCN(nfeat=n_features,
#             nhid=20, #hidden = 16
#             nclass=n_labels,
#             dropout=0.5) #dropout = 0.5
# optimizer = optim.Adam(model.parameters(),
#                        lr=0.001, weight_decay=5e-4)


def train(params):
    torch.manual_seed(34)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    feat_dir = os.path.join(params['train_odir'],'feat')
    label_dir = os.path.join(params['train_odir'],'label')
    if params["water"] == 1:
        edge_dir = os.path.join(params['train_odir'],'edges_water')
    else:
        edge_dir = os.path.join(params['train_odir'],'edges')

    train_data = MyDataset(params['nfeat'],label_dir,feat_dir,edge_dir)

    train_size = int(len(train_data) * 0.8)
    valid_size = int(len(train_data) * 0.2)
    train_dataset, valid_dataset= random_split(train_data, [train_size, valid_size])
    #indexes,_,_,_ = valid_dataset
    # indexes = []
    # for i in valid_dataset:
    #     index,_,_,_ = i
    #     indexes.append(index[:-4])
    # print(indexes)
    # exit(0)
    # weights = class_weight(params)
    # sampler = WeightedRandomSampler(weights,len(weights))

    train_loader = DataLoader(train_dataset, batch_size=params['batch_size'], shuffle=False, num_workers=1,collate_fn=collate_fn)
    validate_loader = DataLoader(valid_dataset, batch_size=params['batch_size'], shuffle=False, num_workers=1,collate_fn=collate_fn)


    model = GCN(feat_num=params['nfeat'],dropout=params['dropout'])
    model = model.to(device=device)
    optimizer= optim.Adam(model.parameters(),lr=params['lr'], weight_decay=params['weight_decay'])
    scheduler = LambdaLR(optimizer, lr_lambda=lambda epoch:1 / (epoch+1))
    #scheduler = torch.optim.lr_scheduler.MultiStepLR(optimizer, milestones=[10,20,30,40,50,60], gamma=0.5, last_epoch=-1)

    #focal_loss = BCEFocalLoss()
    
    train_loss, train_acc,train_precision,train_recall,train_f1 = [], [],[],[],[]
    val_loss, val_acc,val_precision,val_recall,val_f1 = [], [],[],[],[]
    

    for i in range(1,params['epochs']+1):
        batch_loss = []
        batch_acc = []
        batch_precision = []
        batch_recall = []
        batch_f1 = []
        for batch_index,sample in enumerate(train_loader):
            feat,adj,label,natoms= sample
            # print(natoms[0])
            # print(feat.shape)
            # exit(0)
            if params['loss_type'] ==0:
                dice_loss = DiceLoss(natoms=natoms)
            else:
                dice_loss = DiceLoss_atom(natoms=natoms)
            adj = torch.FloatTensor(adj)
            feat = torch.FloatTensor(feat)
            label = torch.FloatTensor(label)

            size = label.shape[0]
            adj = adj.to(device)
            feat = feat.to(device)
            label = label.to(device)


            model.train()
            optimizer.zero_grad()
            output = model(feat, adj)
            output = output.reshape(size,-1)
            label = label.reshape(size,-1)
            #mask = masking(output,label).to(device)
           # output_masked = output.masked_fill(mask == 0, 0).to(device)

            #loss = focal_loss(output,label)
            loss = dice_loss(output,label)
            loss.backward()
            optimizer.step()

            if params['loss_type']==0:
                acc, precision,recall,f1 = metrics(output,label)
            else:
                acc, precision,recall,f1 = metrics_atom(output,label)
            batch_loss.append(loss.item()),batch_acc.append(acc),batch_precision.append(precision),batch_recall.append(recall),batch_f1.append(f1)
            print('Train Loss: {:.3f}, Train Acc: {:.3f}, Train Precision: {:.3f}, Train Recall: {:.3f}, Train F1: {:.3f}'.format(loss.item(), acc,precision,recall,f1))
        scheduler.step()
        print("lr at epoch "+ str(i) +":",optimizer.state_dict()['param_groups'][0]['lr'])
        train_loss.append(np.sum(np.array(batch_loss))/len(train_loader))
        train_acc.append(np.sum(np.array(batch_acc))/len(train_loader))
        train_precision.append(np.sum(np.array(batch_precision))/len(train_loader))
        train_recall.append(np.sum(np.array(batch_recall))/len(train_loader))
        train_f1.append(np.sum(np.array(batch_f1))/len(train_loader))


        #evaluate part
        batch_loss = []
        batch_acc = []
        batch_precision = []
        batch_recall = []
        batch_f1 = []
        for batch_idx,sample in enumerate(validate_loader):
            feat,adj,label,natoms= sample

            if params['loss_type'] ==0:
                dice_loss = DiceLoss(natoms=natoms)
            else:
                dice_loss = DiceLoss_atom(natoms=natoms)
            t = time.time()
            adj = torch.FloatTensor(adj)
            # feat = torch.DoubleTensor(feat)
            # feat = feat.float()
            feat = torch.FloatTensor(feat)
            # label = torch.LongTensor(label)
            # label = label.float()
            label = torch.FloatTensor(label)

            size = label.shape[0]

            adj = adj.to(device)
            feat = feat.to(device)
            label = label.to(device)
            with torch.no_grad():
                model.eval()
                output = model(feat, adj)
                output = output.reshape(size,-1)
                label = label.reshape(size,-1)
                loss = dice_loss(output,label)
                #acc = accuracy(output, label)
                if params['loss_type']==0:
                    acc, precision,recall,f1 = metrics(output,label)
                else:
                    acc, precision,recall,f1 = metrics_atom(output,label)

                # print("epoch",i, "batch:",batch_index,acc, precision,recall,f1)
                # x=torch.tensor(1.).to(device)
                # y=torch.tensor(0.).to(device)
                # print(torch.sum(torch.where(output>=0.5,x,y))/torch.sum(natoms))

                batch_loss.append(loss.item()),batch_acc.append(acc),batch_precision.append(precision),batch_recall.append(recall),batch_f1.append(f1)
            print('Val Loss: {:.3f}, Val Acc: {:.3f},,Val Precision: {:.3f}, Val Recall: {:.3f}, Val F1: {:.3f}'.format(loss.item(), acc,precision,recall,f1))
        val_loss.append(np.sum(np.array(batch_loss))/len(validate_loader))
        val_acc.append(np.sum(np.array(batch_acc))/len(validate_loader))
        val_precision.append(np.sum(np.array(batch_precision))/len(validate_loader))
        val_recall.append(np.sum(np.array(batch_recall))/len(validate_loader))
        val_f1.append(np.sum(np.array(batch_f1))/len(validate_loader))
        print('Epochs: {}, Train Loss: {:.3f}, Train Acc: {:.3f},Train Precision: {:.3f},Train Recall: {:.3f},Train F1: {:.3f}, Validation Loss: {:.3f}, \
        Validation Acc: {:.3f},Validation Precision: {:.3f},Validation Recall: {:.3f},Validation F1: {:.3f}'\
            .format(i, train_loss[-1], train_acc[-1],train_precision[-1],train_recall[-1],train_f1[-1], \
                val_loss[-1], val_acc[-1],val_precision[-1],val_recall[-1],val_f1[-1]))
        if not os.path.exists(params['save_dir']):
            os.mkdir(params['save_dir'])
        torch.save(model.state_dict(),params['save_dir']+'/'+str(i)+'_'+str(optimizer.state_dict()['param_groups'][0]['lr'])+'.ckpt')


from sklearn.metrics import f1_score
from sklearn.metrics import recall_score
from sklearn.metrics import precision_score
from sklearn.metrics import accuracy_score

def metrics(output, labels):
    acc = 0
    precision = 0
    recall = 0
    F1 = 0
    output = output.detach().cpu().numpy()
    labels = labels.detach().cpu().numpy()
    for batch_idx in range(len(output)):
        acc_correct = 0
        pred = np.array([1 if i >=0.5 else 0 for i in output[batch_idx]]).astype(int)
        label = np.array(labels[batch_idx]).astype(int)
        acc +=accuracy_score(label,pred)
        precision +=precision_score(label,pred,zero_division=1)
        recall +=recall_score(label,pred)
        f1 = f1_score(label,pred)
        if isnan(f1):
            f1 = np.nan_to_num(f1)
        F1 += f1
    acc = acc/len(output)
    precision = precision/len(output)
    recall = recall/len(output)
    F1 = F1/len(output)
    return acc,precision,recall,F1

def metrics_atom(output, labels):
    acc = 0
    precision = 0
    recall = 0
    F1 = 0
    output = output.detach().cpu().numpy()
    labels = labels.detach().cpu().numpy()
    pred = []
    label = []
    for batch_idx in range(len(output)):
        pred.append(np.array([1 if i >=0.5 else 0 for i in output[batch_idx]]).astype(int))
        label.append(np.array(labels[batch_idx]).astype(int))
    pred = np.array(pred).reshape(-1)
    label = np.array(label).reshape(-1)
    acc =accuracy_score(label,pred)
    precision =precision_score(label,pred,zero_division=1)
    recall =recall_score(label,pred)
    f1 = f1_score(label,pred)
    if isnan(f1):
        f1 = np.nan_to_num(f1)
    F1 = f1
    # acc_correct = np.sum(np.equal(pred,label))
    # for i in range(len(label)):
    #     if label[i]==pred[i]:
    #         acc_correct +=1
    # all = len(label)
    # TP = np.sum((pred & label))
    # all_pred = np.sum(pred)
    # all_label = np.sum(label)

    # acc = acc_correct/(all+1e-9)
    # precision = TP/(all_pred +1e-9)
    # recall = TP/(all_label+1e-9)
    # F1 = 2*precision*recall/(precision+recall)
    return acc,precision,recall,F1




