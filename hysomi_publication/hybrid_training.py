import sys
import os
from pathlib import Path

current_directory = Path(os.path.dirname(os.path.abspath(__file__)))

import torch
import torch.nn as nn
from torchdiffeq import odeint, odeint_adjoint
from torch.utils.data import Dataset, DataLoader
import hybrid_PBM as PBM
import synthetic_data as sd
import numpy as np
import random
from tqdm import tqdm
import json
import hashlib
import timeit
#%%
torch.set_default_dtype(torch.float64)
torch.set_printoptions(sci_mode=False)
import datetime
#%%

#Training functions

def hyperpar_to_json(NN_hyperpar,loss_hyperpar,training_hyperpar,training_outputs_library_path,length_ID=12):
    #CREATE HYPERPARAMETER DICT
    hyper_dict={'NN_hyperpars':NN_hyperpar,'loss_hyperpar':loss_hyperpar,'training_hyperpar':training_hyperpar}

    #CREATE UNIQUE ID AND ASSOCIATED FOLDER
    raw = hashlib.sha256(json.dumps(hyper_dict, sort_keys=True).encode()).hexdigest()
    run_hash_id=raw[:min(length_ID,64)]
    out_dir = Path(training_outputs_library_path+'/train_'+run_hash_id)
    out_dir.mkdir(parents=True, exist_ok=True)

    #SAVE JSON
    path_to_file=str(out_dir)+'/hyper_dict.json'

    with open(path_to_file, 'w') as file:
        json_string = json.dumps(hyper_dict, default=lambda o: o.__dict__, sort_keys=True, indent=2)
        file.write(json_string)
    print(run_hash_id)
    return hyper_dict,run_hash_id
    
#####################
#AUTOMATIC WEIGHTED LOSS TRAINING##
####################

def train_hybrid_model_autoW(NN_model,ODE_model,data,vdata,time_data,minmax_tuple,EPOCHS,loss_func,batchsize,lr,awl,step_gamma,gamma_lr,save_model_dict_path):
    #MAIN TRAINING FUNCTION FOR HYBRID MODELS
    #input : initialised NN model, initialised PMB, dataset object, validation dataset object, time data, minmax_tuple, epoch number, 
    # initialised loss_function, batchsize, lr, initialised awl module, steps for lr scheduler, multiplier for lr scheduler, output_path
    #TO BE USED BY train_from_json_autoW() FUNCTION
    start = timeit.default_timer()
    trainloader = DataLoader(data, batch_size=batchsize, shuffle=True)#define the data loader
    vtrainloader = DataLoader(vdata, batch_size=batchsize, shuffle=True)#define the validation data loader
    
    optimizer = torch.optim.Adam([{'params': NN_model.parameters()},{'params': awl.parameters(), 'weight_decay': 0}], lr=lr)#define the optimizer : Adam method
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=step_gamma, gamma=gamma_lr) 
    loss_fn = loss_func #define the loss function
    
    N_cons=14 #number of constraints
    N_theta=13 #number of PBM parameters
    loss_part_number=loss_fn.n_loss+N_cons+N_theta

    loss_array=np.zeros((EPOCHS,loss_part_number)) #initializing loss array
    running_loss=np.zeros(loss_part_number)
    
    vloss_array=np.zeros_like(loss_array) #initializing validation loss array
    vMSE_array=np.zeros(EPOCHS)
    lambda_array=[]
    theta_array=[]
    vtheta_array=[]

    best_vloss = 1_000_000.
    
    train_timestamp=str(datetime.datetime.now())
    
    range_mult=NN_model.range_mult
    range_low=NN_model.range_low
    par_trans_list=[range_mult,range_low]

    for epoch in range(EPOCHS):#MAIN LOOP
        NN_model.train(True)
        batch_index=0
        print(f'EPOCH: {epoch}')
        theta_batch=[]
        for batchdata in tqdm(trainloader,ncols=100):#EPOCH LOOP
            traits_batch,process,soil_batch,ini_batch,true_par= batchdata #extract data
            batch_index+=1

            optimizer.zero_grad()

            theta = NN_model(traits_batch) #run NN to output parameter set

            pred=odeint(lambda t,y: ODE_model(t,y,theta,soil_batch),ini_batch,time_data) #run PBM forward

            pred_r=torch.permute(pred,(1,0,2)) #get PBM output in right shape
            pred_r_scaled=(pred_r-minmax_tuple[0])/(minmax_tuple[1]-minmax_tuple[0])#

            loss,loss_MSE,constraints,ana_loss,theta_L1 = loss_fn(pred_r_scaled.double(),process,theta,true_par,par_trans_list) #loss
            wloss=awl(loss) #loss weighting

            print(awl.params.clone().detach().numpy())
            
            wloss.backward() #backward run

            optimizer.step()
            lambda_array.append(awl.params.clone().detach().numpy())
            loss_long=torch.cat((wloss.reshape(1),loss_MSE,ana_loss,constraints,theta_L1)).detach().numpy()
            running_loss += loss_long

            theta_batch.append(theta)

        #LOGGING LOSS
            
        loss_epoch=np.array(running_loss / len(trainloader))
        loss_array[epoch]=loss_epoch

        theta_epoch=torch.vstack(theta_batch)
        theta_array.append(theta_epoch)

        running_loss = running_loss*0.0

        scheduler.step()

        running_loss=np.zeros(loss_part_number)
        
        running_vMSE = 0.0
        running_vloss=np.zeros(loss_part_number)

        with torch.no_grad(): #VALIDATION LOOP
            NN_model.eval()
            print('RUNNING VALIDATION LOOP')
            theta_vepoch=[]
            for vdata in tqdm(vtrainloader,ncols=100):
                vtraits, vprocess,vsoil_batch,vini_batch,vtrue_par = vdata
                vtheta = NN_model(vtraits)
                vpred=odeint(lambda t,y: ODE_model(t,y,vtheta,vsoil_batch),vini_batch,time_data) #run PBM forward

                vpred_r=torch.permute(vpred,(1,0,2))
                vpred_r_scaled=(vpred_r-minmax_tuple[0])/(minmax_tuple[1]-minmax_tuple[0])
                
                vloss,vloss_MSE,vconstraints,vana,vtheta_L1 = loss_fn(vpred_r_scaled, vprocess,vtheta,vtrue_par,par_trans_list)
                vlossw=awl(vloss)#loss weighting
                
                vloss_long=torch.cat((vlossw.reshape(1),vloss_MSE,vana,vconstraints,vtheta_L1)).detach().numpy()
                running_vMSE += vloss_MSE
                running_vloss += vloss_long
                theta_vepoch.append(vtheta)
            
            vtheta_array.append(torch.vstack(theta_vepoch))
            avg_vMSE = running_vMSE / len(vtrainloader)
            avg_vloss = running_vloss / len(vtrainloader)
            vMSE_array[epoch]=avg_vMSE
            vloss_array[epoch]=avg_vloss
            
            if avg_vMSE < best_vloss:
                best_vloss = avg_vMSE
                torch.save(NN_model.state_dict(),save_model_dict_path/'trained_NN_model.pt')

    vtheta_log=torch.stack(vtheta_array)
    theta_log=torch.stack(theta_array)

    #CREATE LOG FILES

    cons_lab=['C1','C2','C3','P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11']
    theta_lab= ['mu_max','K_L','Y_mu','S_tr','r_B','d_B','K_H','vmax','fm','Y_m','m_max','beta','tau']
    loss_labels=['total_loss','MSE', 'L1' ,'par_L1','unmasked_MSE','unmasked_L1','CO2MSE','B_aMSE','B_iMSE','C_LMSE','C_HMSE','C_l_LMSE']+cons_lab+theta_lab
    log_lambda=dict(zip(['MSE']+cons_lab,np.array(lambda_array).transpose()))
    log_dict=dict(zip(loss_labels,list(loss_array.transpose())))
    log_vdict = dict(zip(loss_labels,list(vloss_array.transpose())))

    stop = timeit.default_timer()

    print('Time: ', stop - start)  

    return log_dict,log_vdict,log_lambda,theta_log,vtheta_log

def train_from_json_autoW(training_outputs_library_path,hash_id,save_theta=False):
    #From a training output library path and a hash_id, start a training run. hyperpar_to_json() to create the hash, folder and json configuration file.
    #LOAD JSON HYPERPARAMETER FILE####################################
    json_path=Path(training_outputs_library_path)/str('train_'+hash_id+'/hyper_dict.json')
    out_dir=Path(training_outputs_library_path)/str('train_'+hash_id)
    print(json_path)
    print(out_dir)

    with open(json_path) as json_file:
        hyper_dict = json.load(json_file)

    #DEFINE MODELS####################################################
    NN_hyperdict=hyper_dict['NN_hyperpars']
    seed_run = hyper_dict['training_hyperpar']['seed']
    torch.manual_seed(seed_run)
    np.random.seed(seed_run)
    random.seed(seed_run)

    module = __import__('hybrid_PBM')
    NN = getattr(module, NN_hyperdict['model'])
    NN_model= NN(torch.tensor(NN_hyperdict['par_range']),
                          NN_hyperdict['input_dim'],
                          NN_hyperdict['output_dim'],
                          NN_hyperdict['hidden_dim'],
                          NN_hyperdict['depth'],
                          NN_hyperdict['dropout'])

    ODE_model = PBM.traitSOM_single_pool_BATCH_DEV_PROT_V3()

    print(NN_model)

    #DEFINE LOSS FUNCTION##############################################
    loss_hyperdict = hyper_dict['loss_hyperpar']

    module = __import__('loss_func_hybrid')
    func = getattr(module, loss_hyperdict['loss_fn'])
    awl_mod = getattr(module, loss_hyperdict['awl'])
    loss_function= func(loss_hyperdict['data_sel'],loss_hyperdict['k'])
    print('AUTO WEIGHTING LAMBDA')
    print(f'data_sel : {loss_hyperdict['data_sel']}')
    
    awl=awl_mod(15)
    #LOAD DATASET######################################################
    training_hyperdict = hyper_dict['training_hyperpar']
    dataset_path = training_hyperdict['dataset_path']+'/'

    listing=os.listdir(dataset_path)
    name_data=['trait','X','outputs','ini','soil','min','max']
    dataset_dict={}
    for i in range(len(name_data)):
        dataset_dict[name_data[i]]=torch.tensor(torch.load(
            dataset_path+list(filter(lambda x: x.startswith(name_data[i]),listing))[0]))
    
    min_max_tuple=(dataset_dict['min'],dataset_dict['max'])
    train_validation_ratio=0.8
    test_db_partition=int(len(dataset_dict['trait'])*train_validation_ratio)

    tdata_cut=training_hyperdict['data_cut']

    test_training_set = sd.sim_trait_process(dataset_dict['trait'][:test_db_partition-tdata_cut],
                                                 dataset_dict['outputs'][:,:test_db_partition-tdata_cut],
                                                 dataset_dict['soil'][:test_db_partition-tdata_cut],
                                                 dataset_dict['ini'][:test_db_partition-tdata_cut],
                                                 dataset_dict['X'][:test_db_partition-tdata_cut])
    
    test_validation_set = sd.sim_trait_process(dataset_dict['trait'][test_db_partition:],
                                                 dataset_dict['outputs'][:,test_db_partition:],
                                                 dataset_dict['soil'][test_db_partition:],
                                                 dataset_dict['ini'][test_db_partition:],
                                                 dataset_dict['X'][test_db_partition:])
    
    time_data = torch.linspace(0.0,training_hyperdict['days'],training_hyperdict['days']*24)
    
    log_dict,log_vdict,log_lambda,log_theta,log_vtheta = train_hybrid_model_autoW(NN_model,ODE_model,
                                                            test_training_set,test_validation_set,
                                                            time_data,
                                                            min_max_tuple,
                                                            training_hyperdict['n_epochs'],
                                                            loss_function,
                                                            training_hyperdict['batch_size'],
                                                            training_hyperdict['lr'],
                                                            awl,
                                                            training_hyperdict['scheduler_step'],
                                                            training_hyperdict['scheduler_gamma'],
                                                            out_dir)

    torch.save(log_dict,out_dir / 'loss_dict.pt')
    torch.save(log_vdict,out_dir / 'val_loss_dict.pt')
    torch.save(log_lambda,out_dir /'lambda_dict.pt')
    if save_theta is True:
        torch.save(log_theta,out_dir /'log_theta.pt')
        torch.save(log_vtheta,out_dir /'log_val_theta.pt')

    return log_dict