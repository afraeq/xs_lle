import numpy as np
from PolyLLE import PolyLLE

class Kamide (PolyLLE):
    
    def __init__ (self, label = 'pol', z_sol=0.98, 
                  shulz_kind = '1P', r_pol = np.arange(10,10**4), 
                  p = [0.999,0.995], 
                  teta1 = None, teta2 = None,
                  coarsen = True,
                  alpha = None, xs_fraction = 0.05,
                  equilibrium_method = 'optimize'):

        super().__init__ (label, z_sol,
                          shulz_kind, r_pol, p, teta1, teta2, 
                          coarsen, alpha, xs_fraction, equilibrium_method)
       
        self.n_comp = self.n_pol+1

    #########################
    
    def generate_pseudo_comp (self):
        
        self.r_pol = np.logspace(np.log10(self.r_pol[0]), 
                                 np.log10(self.r_pol[-1]), 
                                 num=len(self.r_pol), endpoint=True, 
                                 base=10.0)

        feed_exp_coarse = np.zeros_like(self.r_pol)
        xs_exp_coarse = np.zeros_like(self.r_pol)
        ins_exp_coarse = np.zeros_like(self.r_pol)
        
        for i in range(len(self.r_pol)-1):
            start = int(self.r_pol[i])
            end = int(self.r_pol[i+1])
            feed_exp_coarse[i] = np.mean(self.feed_exp[start:end])
            xs_exp_coarse[i] = np.mean(self.xs_exp[start:end])
            ins_exp_coarse[i] = np.mean(self.ins_exp[start:end])

        self.feed_exp = feed_exp_coarse
        self.xs_exp = xs_exp_coarse
        self.ins_exp = ins_exp_coarse
       
    #########################
    
    def equilibrium_equations ():
        pass
 
    #########################
           
    def gibbs_energy (self, vij_vetor):

        vij_vetor = vij_vetor
                
        nc = self.n_comp
        nf = self.nf
        vi = self.z
        n = self.r
        v_bar = self.r

        chi_ij = np.zeros((nc,nc))
        vij = np.zeros((nc,nf-1))
        phi_ij = np.zeros_like(vij)
        mi = np.zeros_like(vij)
        vi_nf = np.zeros(nc)
        phi_i_nf = np.zeros_like(vi_nf)
        mi_nf = np.zeros_like(vi_nf)

        for m in range(1,nc):
            chi_ij[0,m] = self.A
            chi_ij[m,0] = self.A

        for j in range (0,nf-1):
            for i in range (0,nc):
                vij[i,j] = vij_vetor[j+i*(nf-1)]
    
        dG = 0
    
        for j in range (0,nf-1):
            for i in range (0,nc):
                phi_ij[i,j] = vij[i,j]/sum(vij[:,j])
            for i in range (0,nc):
                sum1 = 0
                sum2 = 0
                sum3 = 0
                for k in range (0,nc):
                    if k==i:
                        continue
                    sum1 += (1-n[i]/n[k])*phi_ij[k,j]
                    sum2 += phi_ij[k,j]*chi_ij[i,k]
                    for l in range (0,k-1):
                        if l==i:
                            continue
                        sum3 += phi_ij[k,j]*phi_ij[l,j]*chi_ij[l,k]
                mi[i,j] = (np.log(phi_ij[i,j]) + 
                           sum1 + n[i]*((1-phi_ij[i,j])*sum2 - sum3))
                dG += (vij[i,j]/v_bar[i])*mi[i,j]
    
        for i in range(0,nc):
            vi_nf[i] = vi[i] - sum(vij[i,:])
    
        phi_i_nf = vi_nf/sum(vi_nf)
    
        for i in range (0,nc):
            sum1 = 0
            sum2 = 0
            sum3 = 0
            for k in range (0,nc):
                if k==i:
                    continue
                sum1 += (1-n[i]/n[k])*phi_i_nf[k]
                sum2 += phi_i_nf[k]*chi_ij[i,k]
                for l in range(0,k-1):
                    if l==i:
                        continue
                    sum3 += phi_i_nf[k]*phi_i_nf[l]*chi_ij[l,k]
            mi_nf[i] = (np.log(phi_i_nf[i]) + 
                        sum1 + n[i]*(((1-phi_i_nf[i])*sum2 - sum3)))
            dG += (vi_nf[i]/v_bar[i])*mi_nf[i]
                    
        return dG