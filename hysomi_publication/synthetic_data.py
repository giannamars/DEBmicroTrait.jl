
#%%
import torch
import torch.nn as nn
import torch.nn.functional as F
from torchdiffeq import odeint, odeint_adjoint
from torch.utils.data import Dataset, DataLoader
import numpy as np
import random
from scipy.stats import qmc
from hybrid_PBM import traitSOM_single_pool_BATCH_DEV_PROT_V3
#%%
torch.set_default_dtype(torch.float64)
torch.manual_seed(100)
random.seed(100)
torch.set_printoptions(sci_mode=False)
#%%
def flatten(A):
    rt = []
    for i in A:
        if isinstance(i,list): rt.extend(flatten(i))
        else: rt.append(i)
    return rt


#define a Dataset object to be used to store data for training
class sim_trait_process(Dataset):#dataset object to contain the synthetic data (trait-model)
    def __init__(self,traits,model_out,soil_param,ini_val,X):
        self.traits=traits,
        self.model_out=model_out,
        self.soil_param=soil_param,
        self.ini_val=ini_val,
        self.X=X,
    def __len__(self):
        return len(self.traits[0])
    
    def __getitem__(self, index):
        return self.traits[0][index],self.model_out[0][:,index],self.soil_param[0][index],self.ini_val[0][index],self.X[0][index]
    
class synth_trait_network(nn.Module):
    def __init__(self,ranges,input_dim=13, output_dim=38, hidden_dim=38, depth=3):
        super(synth_trait_network, self).__init__()
        self.depth = depth
        
        # Input layer
        self.input_layer = nn.Linear(input_dim, hidden_dim)
        
        # Hidden layers
        self.hidden_layers = nn.ModuleList([
            nn.Linear(hidden_dim,hidden_dim) for _ in range(self.depth - 1)])
        
        # Output layer
        self.output_layer = nn.Linear(hidden_dim, output_dim)

        self.sig_layer = nn.Sigmoid()

        self.range_mult=ranges[:,1]-ranges[:,0].detach()
        self.range_low=ranges[:,0].detach()

        #self.range_layer=enforce_ranges_V2(ranges)

    def forward(self, x):

        x = F.relu(self.input_layer(x))
        for layer in self.hidden_layers:
            x = F.relu(layer(x))
        x = self.output_layer(x)

        x = self.sig_layer(x)

        x = x*self.range_mult+self.range_low

        return x

def par_constraint(theta):#PARAMETER CONSTRAINTS
    mu,m_max,r_B,d_B=theta[0],theta[10],theta[4],theta[5]
    filter=[m_max<mu,mu<d_B,d_B<r_B]
    return not False in filter

def constrained_sampling(par_range,n_sample): #CONSTRAINED SAMPLING FOR SYNTHETIC DATA SELECTION

    l_bounds = par_range[:,0].detach().numpy()
    u_bounds = par_range[:,1].detach().numpy()
    sampler = qmc.LatinHypercube(d=len(par_range))
    sample = sampler.random(n=n_sample*7)
    sample_scaled = qmc.scale(sample, l_bounds, u_bounds)
    par_check=np.array(list(map(par_constraint,sample_scaled)))
    filtered_sample=sample_scaled[par_check]
    return filtered_sample[:n_sample]

def generate_datasetV2(ranges,trait_ranges,days,initial_values,soil_parameters,n,trait_state_dict_path): #CURRENTLY IN USE FUNCTION
    print('Hypercube sampling and parameter constraints...')

    X=torch.tensor(constrained_sampling(ranges,n))
    print('initializing model')

    ini=torch.tensor([initial_values]*len(X))
    soil_batch=torch.tensor([soil_parameters]*len(X))

    time_data = torch.linspace(0.0,days,days*24)

    model=traitSOM_single_pool_BATCH_DEV_PROT_V3()
    print(f'Running model for {n} parameter samples. This can take time')

    batch_ODE = odeint(lambda t,y: model(t,y,X,soil_batch),ini,time_data)

    print('done. Filtering for process constraints')
    sol_list=list(batch_ODE.detach().numpy().transpose())
    Y0,Y1,Y2,Y3,Y4,Y5=sol_list

    Btot=Y0+Y1
    Ctot=Y2+Y3
    whidth_B_A = np.sum(Y0/np.max(Y0,axis=1)[:,None]>0.5,axis=1)

    filt=np.array([np.sum(np.logical_and(Btot>0.05,Btot<3.125),axis=1)/len(time_data)==1,
            np.sum(np.logical_and(Ctot>1.5,Ctot<65),axis=1)/len(time_data)==1,
            np.max(Btot,axis=1)>0.50,
            Y1[:,-1]>0.5*np.max(Y0,axis=1),
            np.argmax(Y0,axis=1)<48,
            np.argmax(Y0,axis=1)>15,
            whidth_B_A<48,
            whidth_B_A>5]).transpose()
        
    sum_filt=np.sum(filt,axis=1)
    n_constraints=np.shape(filt)[1]

    X_filtered=X[sum_filt==n_constraints]
    batch_ODE_filtered=batch_ODE[:,sum_filt==n_constraints]
    soil_param_filtered=soil_batch[sum_filt==n_constraints]
    ini_filtered=ini[sum_filt==n_constraints]


    print(f'checked {n} ouputs for {np.shape(filt)[1]} process constraints. {len(X_filtered)} parameter sets satisfies all constraints')
    print('mapping filtered parameter sets to traits')

    trait_ranges_ready=torch.tensor([trait_ranges]*38)
    NN_trait=synth_trait_network(trait_ranges_ready)

    NN_trait.load_state_dict(torch.load(trait_state_dict_path,weights_only=True))

    trait_tensor=NN_trait(X_filtered)

    return X_filtered, sim_trait_process(trait_tensor,batch_ODE_filtered,soil_param_filtered,ini_filtered,X_filtered),filt,X

#%%
