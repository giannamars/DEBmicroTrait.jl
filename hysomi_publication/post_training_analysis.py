import torch
import numpy as np
import json
import os
from pathlib import Path
from torchdiffeq import odeint

import synthetic_data as sd
import hybrid_training as ht
import hybrid_PBM as PBM
#%%
##### UTILS ###########
def flatten_hyperdict(hyperdict):
    items=[]
    keys=[]
    for key,value in hyperdict.items():
        for subkey,subvalue in value.items():
            items.append(subvalue)
            keys.append(subkey)

    return dict(zip(keys,items))

def classify_run(hyper_dict):
    if 'Auto' in hyper_dict['awl']:
            exp='HySoMi'#HySoMi
    else :
        if hyper_dict['data_sel'][0]==5:
            exp='Unconstrained'
        if hyper_dict['data_sel'][0]==0:
            exp='All states'
    return exp

def load_config_dict(path):
    list_files=os.listdir(path)
    json_path = path+'/'+list(filter(lambda x: 'hyper_dict.json' in x,list_files))[0]
    with open(json_path) as json_file:
        hyper_dict = json.load(json_file)

    hyper_dict = flatten_hyperdict(hyper_dict)
    del hyper_dict['model_repr']
    return hyper_dict

def load_model(config_dict):
    if 'state' in config_dict['model']:
        module = __import__('hybrid_PBM')
        NN = getattr(module, config_dict['model'])
        mod=NN(config_dict['input_dim'],
                          config_dict['output_dim'],
                          config_dict['hidden_dim'],
                          config_dict['depth'],
                          config_dict['dropout'])
    else:  
        module = __import__('hybrid_PBM')
        NN = getattr(module, config_dict['model'])
        mod=NN(torch.tensor(config_dict['par_range']),
                            config_dict['input_dim'],
                            config_dict['output_dim'],
                            config_dict['hidden_dim'],
                            config_dict['depth'],
                            config_dict['dropout'])
    return mod

def add_sum_constraints(train_dicts,val_dicts):
    constraint_sel=('C1','C2','C3','P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11')
    labels=list(train_dicts.keys())

    for k in labels:
        constraints=np.array(list({k: v for k, v in train_dicts[k].items() if k in constraint_sel}.values()))
        constraint_sum=np.sum(constraints,axis=0)
        train_dicts[k]['constraint_sum']=constraint_sum

    for k in labels:
        constraints=np.array(list({k: v for k, v in val_dicts[k].items() if k in constraint_sel}.values()))
        constraint_sum=np.sum(constraints,axis=0)
        val_dicts[k]['constraint_sum']=constraint_sum

    return train_dicts,val_dicts

################################################

def retrieve_runs(output_library_path,hash_dict): 
    #FUNCTION LOADING TRAINING OUTPUT. input : output library path, dictionary of hash codes. output : dictionaries of training, validation loss and training NN models
    #list every directories
    listing=[x[0] for x in os.walk(output_library_path)]
    run_sel=[]
    #filter with hash id
    for id in list(hash_dict.values()):
        run_sel.append(list(filter(lambda x: id in x,listing))[0])
    
    #extract NN hyperparameters
    model_list=[]
    hyper_dict_list= []
    for run in run_sel:
        hyper_dict_flat = load_config_dict(run)
        hyper_dict_list.append(hyper_dict_flat)
        
        model_list.append(load_model(hyper_dict_flat))
        
    #LOAD LOSS DICT AND TRAINED MODEL DICT
    loss_sel=[]
    model_dict_sel=[]
    vloss_sel=[]
    for run in run_sel:
        list_files=os.listdir(run)
        loss_sel.append(torch.load(run+'/'+list(filter(lambda x: x.startswith('loss'),list_files))[0],weights_only=False))
        model_dict_sel.append(torch.load(run+'/'+list(filter(lambda x: x.startswith('trained_NN'),list_files))[0],weights_only=True))
        vloss_sel.append(torch.load(run+'/'+list(filter(lambda x: x.startswith('val'),list_files))[0],weights_only=False))
    #combine models with their trained parameters
    trained_model_list=[]
    for i in range(len(model_list)):
        model = model_list[i]
        model.load_state_dict(model_dict_sel[i])
        trained_model_list.append(model)

    loss_sel_dict=dict(zip(list(hash_dict.keys()),loss_sel))
    vloss_sel_dict=dict(zip(list(hash_dict.keys()),vloss_sel))
    trained_model_dict = dict(zip(list(hash_dict.keys()),trained_model_list))
    config_dict=dict(zip(hash_dict.keys(),hyper_dict_list))
    return loss_sel_dict,vloss_sel_dict,trained_model_dict,config_dict

def load_dataset_dict(dataset_path):
    listing=os.listdir(dataset_path)
    name_data=['trait','X','outputs','ini','soil']
    dataset_dict={}
    for i in range(len(name_data)):
        data=torch.load(dataset_path / list(filter(lambda x: x.startswith(name_data[i]),listing))[0])
        dataset_dict[name_data[i]]=data

    return dataset_dict

######################################################

def inference(trained_model,hyperparameter_dict,test_dataset_dict):

    #DATA PREPARATION
    days=hyperparameter_dict['days']
    time_data = torch.linspace(0.0,days,int(days*24))
    ODE_model=PBM.traitSOM_single_pool_BATCH_DEV_PROT_V3()

    #LOSS FUNCTION INITALISATION
    module = __import__('loss_func_hybrid')
    func = getattr(module, hyperparameter_dict['loss_fn'])
    loss_function= func(hyperparameter_dict['data_sel'],hyperparameter_dict['k'])

    par_trans_list=[trained_model.range_mult,trained_model.range_low]
    cons_lab=['C1','C2','C3','P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11']
    loss_labels=loss_function.loss_labels+cons_lab+['constraint_sum']


    #INFERENCE
    theta = trained_model(test_dataset_dict['trait'])
    pred=odeint(lambda t,y: ODE_model(t,y,theta,test_dataset_dict['soil']),test_dataset_dict['ini'],time_data)

    pred_r=torch.permute(pred,(1,0,2))
    out_data_r=torch.permute(test_dataset_dict['outputs'],(1,0,2))

    loss,loss_MSE,constraints,ana_loss,theta_L1 = loss_function(pred_r,out_data_r,theta,test_dataset_dict['X'],par_trans_list)
    constraint_sum=torch.sum(constraints).reshape(1)
    loss_sum=torch.sum(loss).reshape(1)
    loss_long=torch.cat((loss_sum,loss_MSE,ana_loss,constraints,theta_L1,constraint_sum)).detach().numpy()
    loss_dict=dict(zip(loss_labels,loss_long))

    return loss_dict