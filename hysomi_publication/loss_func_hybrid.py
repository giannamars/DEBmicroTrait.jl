import torch
import torch.nn as nn
import random
import numpy as np

#%%

def softargmax2d(input, beta=10e6):
    *_, h, w = input.shape

    input = input.reshape(*_, h * w)
    input = nn.functional.softmax(beta * input, dim=-1)

    indices_c, indices_r = np.meshgrid(
        np.linspace(0, 1, w),
        np.linspace(0, 1, h),
        indexing='xy'
    )

    indices_r = torch.tensor(np.reshape(indices_r, (-1, h * w)))
    indices_c = torch.tensor(np.reshape(indices_c, (-1, h * w)))

    result_r = torch.sum((h - 1) * input * indices_r, dim=-1)
    result_c = torch.sum((w - 1) * input * indices_c, dim=-1)

    result = torch.stack([result_r, result_c], dim=-1)

    return result


def softargmax1d(input, beta=1e4):
    *_, n = input.shape
    input = nn.functional.softmax(beta * input, dim=-1)
    indices = torch.linspace(0, 1, n)
    result = torch.sum((n - 1) * input * indices, dim=-1)
    return result

def peak_width_V2(input,th,k=1e3):
    return torch.sum(torch.sigmoid((input - th)*k),dim=1)
        
class AutomaticWeightedLossL2(nn.Module):
    """automatically weighted multi-task loss
    """
    def __init__(self, num=15):
        super(AutomaticWeightedLossL2, self).__init__()
        params = torch.ones(num, requires_grad=True)
        self.params = torch.nn.Parameter(params)

    def forward(self, x):
        loss_sum = 0
        for i, loss in enumerate(x):
            loss_sum += 0.5 / (self.params[i] ** 2) * loss + torch.log(1 + self.params[i] ** 2)
        return loss_sum

class no_constraint_module(nn.Module):
    """use for unconstrained runs
    """
    def __init__(self, num=2):
        super(no_constraint_module, self).__init__()
        params = torch.ones(num, requires_grad=True)
        self.params = torch.nn.Parameter(params)

    def forward(self, x):

        return x[0]

class loss_constraint_autoW3(nn.Module):
    def __init__(self,data_sel,k):
        super(loss_constraint_autoW3,self).__init__()
        self.data_loss = torch.nn.MSELoss()
        self.RELU=torch.nn.ReLU()
        self.L1 = torch.nn.L1Loss()
        self.L1_no_reduction = torch.nn.L1Loss(reduction='none')

        self.data_sel=data_sel #index on the state variables to train on. for CO2 only, use [5]. for all states : [0,1,2,3,4,5]
        self.test_data_slice = [0,1,2,3,4] #every state variables but CO2
        self.loss_labels=['total_loss','MSE', 'L1' ,'par_L1','unmasked_MSE','unmasked_L1','CO2MSE','B_aMSE','B_iMSE','C_LMSE','C_HMSE','C_l_LMSE']
        self.n_loss=12
        self.k=k #sharpness parameter for peak_width_V2() and softargmax1d()


    def constraint(self,a,b):#A<B
        return self.RELU((a-b)/a).float().mean().reshape(1)
    

    def forward(self,pred,data,par,true_par,par_trans_list):

        #PARAMETER CONSTRAINTS RELATED LOSS########################
        C1=self.constraint(par[:,10],par[:,0])
        C2=self.constraint(par[:,0],par[:,5])
        C3=self.constraint(par[:,5],par[:,4])

        par_constraint_tensor=torch.cat((C1,C2,C3))
        ###########################################################

        #PROCESS CONSTRAINTS RELATED LOSS##########################
        Btot=pred[:,:,0]+pred[:,:,1]
        Ctot=pred[:,:,2]+pred[:,:,3]+pred[:,:,4]
        C_H=pred[:,:,3]
        B_i=pred[:,:,1]
        B_a=pred[:,:,0]
        B_a_max=torch.max(pred[:,:,0],dim=1)[0]
        Ba_ratio = pred[:,:,0]/B_a_max[:,None]
        B_a_whidth = peak_width_V2(Ba_ratio,0.5,k=self.k)

        PC1 = self.constraint(Btot,3.180)
        PC2 = self.constraint(0.05,Btot)
        PC3 = self.constraint(Ctot,65)
        PC4 = self.constraint(1.5,Ctot)
        PC5 = self.constraint(0.5,torch.max(Btot,dim=1)[0])
        PC6 = self.constraint(0.5*B_a_max,B_i[:,-1])
        PC7 = self.constraint(softargmax1d(B_a,beta=self.k),48)
        PC8 = self.constraint(15,softargmax1d(B_a,beta=self.k))
        PC9 = self.constraint(B_a_whidth,48)
        PC10 = self.constraint(5,B_a_whidth)
        PC11 = self.constraint(0.9*C_H[:,0],C_H[:,-1])

        process_constraint_tensor = torch.cat((PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11))
        ###########################################################
        constraint_tensor = torch.cat((par_constraint_tensor,process_constraint_tensor))
        #DATA RELATED LOSS#########################################

        loss_pred=self.data_loss(pred[:,:,self.data_sel],data[:,:,self.data_sel]).reshape(1)
        
        #loss for analysis
        loss_pred_L1=self.L1(pred[:,:,self.data_sel],data[:,:,self.data_sel]).reshape(1)
        loss_pred_unmasked = self.data_loss(pred[:,:,self.test_data_slice],data[:,:,self.test_data_slice]).reshape(1)#MSE_Hidden
        loss_pred_L1_unmasked = self.L1(pred[:,:,self.test_data_slice],data[:,:,self.test_data_slice]).reshape(1)

        par_unranged=(par-par_trans_list[1])/par_trans_list[0]
        true_par_unranged=(true_par-par_trans_list[1])/par_trans_list[0]
        hidden_par_loss = self.L1(par_unranged,true_par_unranged).reshape(1)#L1_theta
        
        loss_pred_CO2=self.data_loss(pred[:,:,5],data[:,:,5]).reshape(1)
        loss_pred_B_a=self.data_loss(pred[:,:,0],data[:,:,0]).reshape(1)
        loss_pred_B_i=self.data_loss(pred[:,:,1],data[:,:,1]).reshape(1)
        loss_pred_C_L=self.data_loss(pred[:,:,2],data[:,:,2]).reshape(1)
        loss_pred_C_H=self.data_loss(pred[:,:,3],data[:,:,3]).reshape(1)
        loss_pred_C_l_L=self.data_loss(pred[:,:,4],data[:,:,4]).reshape(1)

        theta_L1_tensor = self.L1_no_reduction(par_unranged,true_par_unranged).mean(dim=0)
        
        ana_loss = torch.cat((loss_pred_L1,hidden_par_loss,loss_pred_unmasked,loss_pred_L1_unmasked,loss_pred_CO2,loss_pred_B_a,loss_pred_B_i,loss_pred_C_L,loss_pred_C_H,loss_pred_C_l_L))
        ###########################################################
        loss_tensor=torch.cat((loss_pred,constraint_tensor))

        
        return loss_tensor,loss_pred,constraint_tensor,ana_loss,theta_L1_tensor
