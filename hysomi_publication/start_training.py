#%%
from pathlib import Path
import os
import sys
current_directory = Path(os.path.dirname(os.path.abspath(__file__)))

import torch
import numpy as np
import random
import hybrid_training as ht
import post_training_analysis as pa


run_seed=100 #run seed

torch.set_default_dtype(torch.float64)
torch.manual_seed(run_seed)
random.seed(run_seed)
np.random.seed(run_seed)
torch.set_printoptions(sci_mode=False)



#TRAINING FOR HYBRID MODEL#
#SETTING UP MODELS
par_range = torch.tensor([[0.01,10],#mu_B_O
                                [1,500],#K_B_O
                                [0.20,0.995],#Y_B_O
                                [1e-7,0.3],#s_tr_O
                                [0.01,100],#r_B_O
                                [0.01,100],#d_B_O
                                [10e-6,10],#K_C_H
                                [1e-5,1],#vmax
                                [0.2,0.9],#f_m
                                [0.25,0.995],#Y_m
                                [1e-4,2],#m_max
                                [1e-4,1],#beta
                                [1e-4,1]])#tau



# %%
data_path = current_directory / 'data/synthetic_data/' #current_directory / 'data/synthetic_training_set/'
print(data_path)

NN_hyper={'model':'MLP_V6',
          'par_range':par_range.tolist(),
          'input_dim':38,
          'output_dim':13,
          'hidden_dim':38,
          'depth':5,
          'dropout':0}
         
          
module = __import__('hybrid_PBM')
NN = getattr(module, NN_hyper['model'])
NN_model= NN(torch.tensor(NN_hyper['par_range']),
                    input_dim=NN_hyper['input_dim'],
                    output_dim=NN_hyper['output_dim'],
                    hidden_dim=NN_hyper['hidden_dim'],
                    depth=NN_hyper['depth'],
                    dropout = NN_hyper['dropout'])
                          
repr=NN_model.__repr__()
NN_hyper['model_repr']=repr

loss_hyper={'data_sel':[5],
            'awl': 'AutomaticWeightedLossL2',# HySoMi : 'AutomaticWeightedLossL2'     #unconstrained or all states : 'no_constraint_module'
            'loss_fn':'loss_constraint_autoW3',#loss_constraint_autoW3
            'k':1e8} #sharpness parameter for differentiable argmax and peak width
            

training_hyper={'dataset_path':str(data_path),
                'days':7, #simulation lenght. default is 7 days, must be same as in data 
                'n_epochs':3, #number of epochs
                'batch_size':512,
                'lr':0.5e-2, #learning rate
                'scheduler_step':10,
                'scheduler_gamma':0.9,
                'seed':run_seed,
                'data_cut':0} #number of sample to be removed from the synthetic dataset. keep it to 0


training_outputs_library_path = current_directory / 'data/outputs/'
print(training_outputs_library_path)
#%%
print('running hybrid model training')
hybrid_dict,hash_id=ht.hyperpar_to_json(NN_hyper,loss_hyper,training_hyper,str(training_outputs_library_path)) #create hyperparameter json file and associated folder and hash ID
print(hash_id)

print(hybrid_dict)

log_dict=ht.train_from_json_autoW(training_outputs_library_path,hash_id,save_theta=False) #training function

print(f'results saved at{str(training_outputs_library_path)}/train_{hash_id}')

print('running inference on test dataset')
test_data_path = data_path / 'test_set'
inference_dataset=pa.load_dataset_dict(test_data_path)
NN_model.load_state_dict(torch.load(training_outputs_library_path / f'train_{hash_id}/trained_NN_model.pt'))

test_loss = pa.inference(NN_model,pa.flatten_hyperdict(hybrid_dict),inference_dataset)
print(test_loss)

torch.save(test_loss,training_outputs_library_path / f'train_{hash_id}/test_loss.pt')

# %%
