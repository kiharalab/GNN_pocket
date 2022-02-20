from __future__ import print_function
from __future__ import division
import torch
import math
from torch.nn.parameter import Parameter
from torch.nn.modules.module import Module




class GraphConvolution(Module):
    
    def __init__(self, in_features, out_features, bias=True):
        super(GraphConvolution, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.weight = Parameter(torch.FloatTensor(in_features, out_features))
        if bias:
            self.bias = Parameter(torch.FloatTensor(out_features))
        else:
            self.register_parameter('bias', None)
        self.reset_parameters()

    def reset_parameters(self):
        stdv = 1. / math.sqrt(self.weight.size(1))
        self.weight.data.uniform_(-stdv, stdv)
        if self.bias is not None:
            self.bias.data.uniform_(-stdv, stdv)

    def forward(self, input, adj):
        #device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        support = torch.einsum('ijk, kh -> ijh', input, self.weight)
        #support = support.to(device)
        #support = torch.bmm(input, self.weight)
        output = torch.einsum('ijk, ikh -> ijh',adj, support)
        #output = output.to(device)
        #output = torch.bmm(adj, support)
        if self.bias is not None:
            return output + self.bias.T
        else:
            return output

    def __repr__(self):
        return self.__class__.__name__ + ' (' \
               + str(self.in_features) + ' -> ' \
               + str(self.out_features) + ') '