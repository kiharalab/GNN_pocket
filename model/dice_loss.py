import torch.nn as nn
import torch
import numpy as np

torch.autograd.set_detect_anomaly(True)

class DiceLoss(nn.Module):
    def __init__(self, natoms,reduction='elementwise_mean'):
        super().__init__()
        self.natoms = natoms
        self.reduction = reduction
 
    def	forward(self, input, target):
        # print(input.shape, target.shape)
        #print(self.natoms)
        smooth= 1
        max_atoms = input.shape[1]
        pred = []
        y = []
        input_1 =input
        target_1 = target
        final_input = []
        final_target = []
        losses = 0
        for i in range(len(input_1)):
            final_input = []
            final_target = []
            final_input.append(input_1[i][:int(self.natoms[i])].unsqueeze(0))
            final_target.append(target_1[i][:int(self.natoms[i])].unsqueeze(0))
            final_input = torch.cat(final_input,dim=1)
            final_target = torch.cat(final_target,dim=1)
            intersection = final_input * final_target
            loss = 2 * (intersection.sum(1) + smooth) / (final_input.sum(1) + final_target.sum(1) + smooth)
            if self.reduction == 'elementwise_mean':
                loss = 1 - loss.sum()
            losses +=loss
        N = target_1.size(0)

        # final_input = torch.cat(final_input,dim=1)
        # final_target = torch.cat(final_target,dim=1)

        # N = target_1.size(0)
        # input_flat = final_input
        # target_flat = final_target
        
        # #input_flat = input_1.view(N, -1)
        # #target_flat = target_1.view(N, -1)

        # intersection = input_flat * target_flat

        # loss = 2 * (intersection.sum(1) + smooth) / (input_flat.sum(1) + target_flat.sum(1) + smooth)
        # if self.reduction == 'elementwise_mean':
        #     loss = 1 - loss.sum() / N
        losses = (losses/N).to(torch.device("cuda:0" if torch.cuda.is_available() else "cpu"))
        return losses

class DiceLoss_atom(nn.Module):
    def __init__(self, natoms,reduction='elementwise_mean'):
        super().__init__()
        self.natoms = natoms
        self.reduction = reduction
 
    def	forward(self, input, target):
        # print(input.shape, target.shape)
        #print(self.natoms)
        smooth= 1
        input_1 =input
        target_1 = target
        final_input = []
        final_target = []
        for i in range(len(input_1)):
            final_input = []
            final_target = []
            final_input.append(input_1[i][:int(self.natoms[i])].unsqueeze(0))
            final_target.append(target_1[i][:int(self.natoms[i])].unsqueeze(0))
            final_input = torch.cat(final_input,dim=1)
            final_target = torch.cat(final_target,dim=1)
            # intersection = final_input * final_target
            # loss = 2 * (intersection.sum(1) + smooth) / (final_input.sum(1) + final_target.sum(1) + smooth)
            # if self.reduction == 'elementwise_mean':
            #     loss = 1 - loss.sum()
            # losses +=loss
        N = target_1.size(0)

        input_flat = final_input
        target_flat = final_target
        
        input_flat = input_1.view(N, -1)
        target_flat = target_1.view(N, -1)

        intersection = input_flat * target_flat

        loss = 2 * (intersection.sum(1) + smooth) / (input_flat.sum(1) + target_flat.sum(1) + smooth)
        if self.reduction == 'elementwise_mean':
            loss = 1 - loss.sum() / N
        #losses = (losses/N).to(torch.device("cuda:0" if torch.cuda.is_available() else "cpu"))
        return loss