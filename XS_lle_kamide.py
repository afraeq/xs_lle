#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Determination of Flory-Huggins Interaction Parameters 
for Polypropylene+Xylene Polydisperse Solutions from
Xylene Solubles (XS) Analysis Data

@author: afranio
         afrjr.weebly.com

Flory-Huggins model as provided by KAMIDE e DOBASHI, 2000 - 
Physical Chemistry of Polymer Solutions - Theoretical Background.

"""

import numpy as np

from XS_lle import XS_lle

class XS_lle_kamide (XS_lle):
    
    def __init__ (self, label, z_sol, teta1, teta2, alpha, mw, xs_fraction,
                  npol, rmin, rmax, equilibrium_method = 'optimize'):
       
        self.npol = npol
        self.rmin = rmin
        self.rmax = rmax

        super().__init__ (label, z_sol, teta1, teta2, alpha, mw, xs_fraction)
 
        self.equilibrium_method = equilibrium_method
      
        self.ncomp = npol+1

    #########################
    
    def generate_pseudo_comp (self):
        
        # solvent chain length
        self.r_sol = np.array([1])

        # polymer chain lengths
        self.r_pol = np.arange(1,40000)
        
        # uniform discretizing in logarithmic space
        
        nclass = self.npol
        
        rclass         = np.zeros(nclass)
        feed_exp_class = np.zeros(nclass)
        xs_exp_class   = np.zeros(nclass)
        ins_exp_class  = np.zeros(nclass)
        aux            = np.zeros(nclass)

        delta_log_r = (np.log(self.rmax)- np.log(self.rmin))/(nclass-1)

        for i in range(nclass):
            rclass[i] = np.floor(np.exp((i)*delta_log_r+np.log(self.rmin)))

        feed_exp_class[0] = sum(self.feed_exp[0:int(rclass[0])+1])
        xs_exp_class[0]   = sum(self.xs_exp[0:int(rclass[0])+1])
        ins_exp_class[0]  = sum(self.ins_exp[0:int(rclass[0])+1])

        for i in range(1,nclass):
            feed_exp_class[i] = sum(self.feed_exp[int(rclass[i-1]):
                                                  int(rclass[i])+1])
            xs_exp_class[i] = sum(self.xs_exp[int(rclass[i-1]):
                                              int(rclass[i])+1])
            ins_exp_class[i] = sum(self.ins_exp[int(rclass[i-1]):
                                                int(rclass[i])+1])

        feed_exp_class = feed_exp_class/sum(feed_exp_class)
        xs_exp_class = xs_exp_class/sum(xs_exp_class)
        ins_exp_class = ins_exp_class/sum(ins_exp_class)

        aux[:] = rclass

        aux[0] = np.floor(np.exp(((np.log(rclass[0])+np.log(1)/2))))

        for i in range (1,nclass):
            aux[i] =  np.floor(np.exp(((np.log(rclass[i])+
                                        np.log(rclass[i-1]))/2)))
            
        self.r_pol     = aux

        self.feed_exp  = feed_exp_class
        self.xs_exp    = xs_exp_class
        self.ins_exp   = ins_exp_class
       
    #########################
    
    def equilibrium_equations ():
        pass
 
    #########################
           
    def gibbs_energy (self, vij_vetor):
        
        nc = self.ncomp
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
            
        for m in range (1,nc):
            chi_ij[0,m], = self.chi
            chi_ij[m,0], = self.chi
    
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