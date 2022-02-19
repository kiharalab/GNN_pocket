import parser
import argparse

from click import argument

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode',type=int,required=True,help='0: training 1: predicting')
    parser.add_argument('--train_dir',type=str,default="../dataset/train")
    parser.add_argument('--train_odir',type=str,default="../dataset/train_processed")
    parser.add_argument('--test_dir',type=str,default="../dataset/test")
    parser.add_argument('--test_odir',type=str,default="../dataset/test_processed")
    parser.add_argument('--water',type=int,default=0,help='0:use edges without water, 1:use edges with water')
    parser.add_argument('--nfeat',type=int,default=3,help='how many features needed')
    parser.add_argument('--gpu',type=str,default='0',help='Choose gpu id, example: \'1,2\'(specify use gpu 1 and 2)')
    parser.add_argument("--batch_size", help="batch_size", type=int, default=16)
    parser.add_argument("--epochs", help="epochs", type=int, default=100)
    parser.add_argument("--lr",help="learning_rate", type=float, default=0.001)
    parser.add_argument("--dropout", help="dropout_rate", type=float, default=0.3)
    parser.add_argument('--loss_type',type=int,default=0,help='0: graph loss 1: atom loss')

    parser.add_argument('--ensemble',type=int,default=1)
    #parser.add_argument("--ensemble_mode",type=int,default=0,help="0 union,1 avg 2 voting")

    # parser.add_argument('--single',type=str, default="./single_feats",help='single feature save dir')#File path for decoy dir
    # parser.add_argument('--feat',type=str, default="./features",help='combined feature save dir')
    # parser.add_argument('--edge',type=str, default="./edges",help='combined feature save dir')
    # parser.add_argument('--label',type=str, default="./label",help='label save dir')
    
    
    # parser.add_argument('--vis',type=int,default=0,help='The mode for visgrid voxel feature,0:closest,8:8A')
    # parser.add_argument('--atom_type',type=int,required=True,help='0: single, 1: one hot')
    
    
    parser.add_argument('--model_type',type=int,default=1,help='which model to predict')
    parser.add_argument('--save_dir',type=str,default='./saved_model',help='Model saved dir')
   
    parser.add_argument("--tensorboard", help="tensorboard_dir", type=str, default="../tensorboard_runs/")
    parser.add_argument("--weight_decay",help="weight_decay",type=float,default=5e-4)
    parser.add_argument('--seed',type=int,default=888,help='random seed for shuffling')
    parser.add_argument('--ensemble_mode',type=int,default=0,help='model ensembling,0 for avg,1 for voting')
    #parser.add_argument('--fold',required=True,help='specify fold model for prediction',type=int,default=-1)
    args = parser.parse_args()
    params = vars(args)
    return params
