
import torch
import numpy as np
import uuid
import synthetic_data as sd
from pathlib import Path
import os

current_directory = Path(os.path.dirname(os.path.abspath(__file__)))

trait_NN=sd.synth_trait_network_test(hidden_dim=38, depth=3)



par_range = torch.tensor([[0.01,100],#mu_B
                                [1,500],#K_B
                                [0.20,0.995],#Y_B
                                [1e-7,0.3],#s_tr
                                [0.01,100],#r_B
                                [0.01,100],#d_B
                                [1e-6,10],#K_C_H
                                [1e-4,0.1],#vmax
                                [0.2,0.9],#f_m
                                [0.25,0.995],#Y_m
                                [1e-4,2],#m_max
                                [1e-4,1],#beta
                                [1e-4,1]])#tau


ini_val=[0.0000125,0.2499875,1,10,0,0]
                           
soilparams = [0.1,0.01,0.5,0,0]

trait_range=[5,200]

trait_state_dict_path = current_directory / 'data/synthetic_data/NNtrait_dict.pt'
print('initializing trait NN')

print('starting dataset generation')
X,dataset,filt,X_raw=sd.generate_datasetV2(par_range,trait_NN,7,ini_val,soilparams,2000000,str(trait_state_dict_path))

trait_data=dataset.traits
out_data=dataset.model_out[0]

sum_filt=np.sum(filt,axis=1)
n_constraints=len(filt)

X_less=X_raw[sum_filt==n_constraints-1]

print(len(X))

id=uuid.uuid1()#generate random id

#SAVE DATASET AS .PT FILES
out_dir=current_directory / 'data/data_generation/'
out_dir.mkdir(parents=True, exist_ok=True)
if len(X)>0 :
  torch.save(X,out_dir / f'/X_{id}.pt')
  torch.save(trait_data[0],out_dir / f'/trait_{id}.pt')
  torch.save(out_data,out_dir / f'/outputs_{id}.pt')
  print('parameter sets found ! saved dataset !')
  
if len(X)==0 :
  print('no parameter set found,saving parameters that fulfilled n-1 constraints and filter for inspection')
  torch.save(filt,out_dir / f'/filt_{id}')
  torch.save(X_less,out_dir / f'/X_not_filtered_{id}')
  

