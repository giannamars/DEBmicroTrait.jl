#%%
import torch
import torch.nn as nn
import torch.nn.functional as F

#%%
    
class MLP_V6(nn.Module): #PARAMETER ESTIMATION MODULE
    #OUTPUT 13 PARAMETERS FOR traitSOM_single_pool_BATCH_DEV_PROT_V3()
    def __init__(self,ranges,input_dim=38, output_dim=13, hidden_dim=64, depth=3,dropout=0):
        super(MLP_V6, self).__init__()
        self.depth = depth
        
        # Input layer
        self.input_layer = nn.Linear(input_dim, hidden_dim)
        
        # Hidden layers
        self.hidden_layers = nn.ModuleList([
            nn.Linear(hidden_dim,hidden_dim) for _ in range(self.depth - 1)])
        
        self.BN_layers = nn.ModuleList([
            nn.BatchNorm1d(hidden_dim) for _ in range(self.depth - 1)])
        
        #self.BN = nn.BatchNorm1d(hidden_dim)
        self.drop=nn.Dropout1d(p=dropout)
        
        # Output layer
        self.output_layer = nn.Linear(hidden_dim, output_dim)
        self.sig_layer = nn.Sigmoid()

        #Range transformation
        self.range_mult=ranges[:,1]-ranges[:,0].detach()
        self.range_low=ranges[:,0].detach()

    def forward(self, x):

        x = F.relu(self.input_layer(x))
        for i in range(len(self.hidden_layers)):
            x = self.drop(F.relu(self.BN_layers[i](self.hidden_layers[i](x))))
        x = self.output_layer(x)
        x = self.sig_layer(x)
        x = x*self.range_mult+self.range_low

        return x        

class MLP_V6_state(nn.Module): #PURE ML MODEL : OUTPUTS STATE VARIABLES DIRECTLY (USE FOR COMPARISON)
    def __init__(self,input_dim=38, output_dim=168*6, hidden_dim=64, depth=3,dropout=0):
        super(MLP_V6_state, self).__init__()
        self.depth = depth
        
        # Input layer
        self.input_layer = nn.Linear(input_dim, hidden_dim)
        
        # Hidden layers
        self.hidden_layers = nn.ModuleList([
            nn.Linear(hidden_dim,hidden_dim) for _ in range(self.depth - 1)])
        
        self.BN_layers = nn.ModuleList([
            nn.BatchNorm1d(hidden_dim) for _ in range(self.depth - 1)])
        
        #self.BN = nn.BatchNorm1d(hidden_dim)
        self.drop=nn.Dropout1d(p=dropout)
        
        # Output layer
        self.output_layer = nn.Linear(hidden_dim, output_dim)
        #self.sig_layer = nn.Sigmoid()

        #Range transformation
        #self.range_mult=ranges[:,1]-ranges[:,0].detach()
        #self.range_low=ranges[:,0].detach()

    def forward(self, x):

        x = F.relu(self.input_layer(x))
        for i in range(len(self.hidden_layers)):
            x = self.drop(F.relu(self.BN_layers[i](self.hidden_layers[i](x))))
        x = self.output_layer(x)
        #x = self.sig_layer(x)
        #x = x*self.range_mult+self.range_low
        x = x.view(x.size(0), 14, 6)

        return x
        

#%%
#################################################
####PROCESS BASED MODEL####
#################################################
        
class traitSOM_single_pool_BATCH_DEV_PROT_V3(nn.Module): #MICROBIAL EXPLICIT SOIL PBM
    def __init__(self):
        super(traitSOM_single_pool_BATCH_DEV_PROT_V3, self).__init__()

        self.batched_dot = torch.func.vmap(torch.dot)
    
    def forward(self, t, y,biopar,soilpar):
        
        mu = biopar[:,0]#define each parameters
        K_L = biopar[:,1]
        Y_mu = biopar[:,2]
        S_tr = biopar[:,3]
        r_react = biopar[:,4]
        r_deact = biopar[:,5]
        K_C_H = biopar[:,6]
        vmax = biopar[:,7]
        fm = biopar[:,8]
        Y_m = biopar[:,9]
        m_max = biopar[:,10]
        beta = biopar[:,11]
        tau= biopar[:,12]
            
        k_a = soilpar[:,0]
        k_d = soilpar[:,1]
        C_S_max = soilpar[:,2]
        I_C_L = soilpar[:,3]
        I_C_H = soilpar[:,4]
        
            
        B_a = y[:,0]#state variables
        B_i = y[:,1]
        C_l_L = y[:,2]
        C_H = y[:,3]
        C_s_L = y[:,4]
        CO2 = y[:,5]

        unit_vector=torch.ones(len(biopar),2).double()#unit tensor to be used with torch.dot, symplify the implementation

        dydt = torch.zeros_like(y)#ensure that output is in the right shape
        
        r_a=r_react*B_i*(1/(torch.exp((S_tr-C_l_L)/tau*S_tr)+1))#bacterial activation rates
        r_d=r_deact*B_a*(1-(1/(torch.exp((S_tr-C_l_L)/tau*S_tr)+1)))#bacterial deactivation rates
        
        r_m_a_B=m_max*B_a#maintenance rates
        r_m_a_M=((m_max*C_l_L*K_L)/(m_max+C_l_L*K_L))*B_a
        r_m_i_B=m_max*B_i*beta
        r_m_i_M=((m_max*C_l_L*K_L)/(m_max+C_l_L*K_L))*B_a*beta
        
        
        r_m=torch.vstack([r_m_a_M,r_m_i_M]).transpose(0,1).double()#maintenance vector (tensor)
        r_B=torch.vstack([r_m_a_B,r_m_i_B]).transpose(0,1).double()

        #MAIN EQUATIONS
        
        dydt[:,0] = mu*B_a*(K_L*C_l_L/(mu+K_L*C_l_L))+r_a-r_d-(1/Y_m)*(r_m_a_B-r_m_a_M)#active bacteria
            
        dydt[:,1] = r_d-r_a-(1/Y_m)*(r_m_i_B-r_m_i_M)#inactive bacteria
            
        dydt[:,2]=I_C_L+(1-fm)*(self.batched_dot(unit_vector,r_B-r_m))-\
            ((1/Y_mu)*mu*B_a*(K_L*C_l_L/(mu+K_L*C_l_L)))-self.batched_dot(unit_vector,r_m)+\
                vmax*B_a*(C_H/(K_C_H+C_H))-k_a*C_l_L*(C_S_max-C_s_L)+k_d*C_s_L#Low molecular weight carbon in solution
                
        dydt[:,3]=I_C_H+fm*(self.batched_dot(unit_vector,r_B-r_m))-\
            vmax*B_a*(C_H/(K_C_H+C_H))#High molecular weight carbon in solution
            
        dydt[:,4]=k_a*C_l_L*(C_S_max-C_s_L)-k_d*C_s_L#sorbed low molecular weight carbon
        
        dydt[:,5]=((1-Y_mu)/Y_mu)*mu*B_a*(K_L*C_l_L/(mu+K_L*C_l_L))+((1-Y_m)/Y_m)*self.batched_dot(unit_vector,r_B-r_m)+self.batched_dot(unit_vector,r_m)#microbial respiration (CO2)
        return dydt #return tensor shape = (time,batch_size,state variables)

