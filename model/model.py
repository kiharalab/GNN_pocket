from __future__ import print_function
from __future__ import division

import torch
import torch.nn as nn
import torch.nn.functional as F

from .layers import GraphConvolution



class GCN(nn.Module):
    def __init__(self,nfeat,nhid1=128, nhid2=64,nhid3=32,nhid4=32,nhid5=16,nout=1, dropout=0.3):
        super(GCN, self).__init__()
        self.fc1 = nn.Linear(nfeat,nhid1)

        self.gc1 = GraphConvolution(nhid1, nhid2)
        self.gc2 = GraphConvolution(nhid2, nhid3)
        self.gc3 = GraphConvolution(nhid3, nhid4)
        self.gc4 = GraphConvolution(nhid4, nhid5)

        self.fc2 = nn.Linear(nhid5,nout)
        self.dropout = dropout

    def forward(self, x, adj):

        x = self.fc1(x)
        x = F.dropout(x, self.dropout)

        x = F.relu(self.gc1(x, adj))
        x = F.dropout(x, self.dropout)

        x = F.relu(self.gc2(x, adj))
        x = F.dropout(x, self.dropout)

        x = F.relu(self.gc3(x, adj))
        x = F.dropout(x, self.dropout)

        x = F.relu(self.gc4(x, adj))
        x = F.dropout(x, self.dropout)

        x = self.fc2(x)
        return torch.sigmoid(x)
        #return x