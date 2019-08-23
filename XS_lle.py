#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Determination of Flory-Huggins Interaction Parameters 
for Polypropylene+Xylene Polydisperse Solutions from
Xylene Solubles (XS) Analysis Data

@author: afranio
         afrjr.weebly.com

This module was devoloped to perform calculations with data from three papers:

    MATOS, V., MATTOS NETO, A. G., PINTO, J. C., 2001
    “Method for Quantita- tive Evaluation of Kinetic Constants in
    Olefin Polymerizations
    I. Kinetic Study of a Conventional Ziegler-Natta Catalyst Used for
    Propylene Poly- merizations”,
    Journal of Applied Polymer Science, v. 79, pp. 2076–2018.

,

    MATOS, V., MATTOS NETO, A. G., PINTO, J. C., et al., 2002,
    “Method for Quantitative Evaluation of Kinetic Constants in
    Olefin Polymerizations.
    II. Kinetic Study of a High-Activity Ziegler–Natta Catalyst Used for
    Bulk Propylene Polymerizations”,
    Journal of Applied Polymer Science, v. 86, pp. 3226–3245.

and

    MATOS, V., MOREIRA, M., MATTOS NETO, A. G., et al., 2007,
    “Method for Quantitative Evaluation of Kinetic Constants in
    Olefin Polymerizations. III. Kinetic Study of the hiPP Synthesis”,
    Macromolecular Reaction Engineering, v. 1, pp. 137–159.

"""

import numpy as np
import scipy.optimize as opt
from abc import ABC, abstractmethod

class XS_lle (ABC):

    def __init__ (self, label, z_sol, teta1, teta2, alpha, mw, xs_fraction):

        # polymer label
        self.label = label

        # solvent extensive volume <---------------------------------------------------------------------------------------
        #estranho, ver linha 112
        self.z_sol = np.array([z_sol])

        # shulz-flory parameters
        # [0]: feed, [1]: solubles
        self.q1_exp = np.exp(-teta1)
        self.q2_exp = np.exp(-teta2)
        self.alpha_exp = alpha

        # feed mean molar mass
        self.mw = mw

        # xs fraction: polymer fraction extracted by xylene
        self.xs_fraction = xs_fraction

        # number of phases
        self.nf = 2
        
        # solvent chain length
        self.r_sol = np.array([1])

        # polymer chain lengths
        self.r_pol = np.arange(1,40000)
        
        # experimental feed distribution
        #NORMALIZACAO????????????????????????????????????????????????????????<------------------------------------------
        self.feed_exp = ((self.shulz_flory(self.r_pol,
                                           self.q1_exp[0],
                                           self.q2_exp[0],
                                           self.alpha_exp[0]))/
                         sum(self.shulz_flory(self.r_pol,
                                              self.q1_exp[0],
                                              self.q2_exp[0],
                                              self.alpha_exp[0])))
        
        
        #print(sum(self.shulz_flory(self.r_pol,self.q1_exp[0],self.q2_exp[0],self.alpha_exp[0])))
        # experimental xs distribution
        #nomeclatura não esta calara, na verdade xs exp é a curva de fracoes mássicas dos diversos tamanhos de polímeros que permaneceram na fase liquida apos o resfriamneto da solucao do teste xs, portanto esse vetor armazena a curva de distribuição de tamanhos dos polímeros que firaram na fase líquida ou seja e a curva do pp que e atatico
        self.xs_exp = ((self.shulz_flory(self.r_pol,
                                        self.q1_exp[1],
                                        self.q2_exp[1],
                                        self.alpha_exp[1]))/
                       sum(self.shulz_flory(self.r_pol,
                                            self.q1_exp[1],
                                            self.q2_exp[1],
                                            self.alpha_exp[1])))

        # experimental insolubles distribution
        #this expression came from a mass balance
        #look in comments in the end of file
        self.ins_exp = ((self.feed_exp-self.xs_fraction*self.xs_exp)/
                        (1-self.xs_fraction))
        for i in range(len(self.ins_exp)):
            if self.ins_exp[i]<0:
                self.ins_exp[i] = 0
                
        # generating pseudo components
        self.generate_pseudo_comp()

        # polymer extensive volumes
        #mass fraction given by the shultz flory model is equal to de volume fraction, in this special case as shown by afranio
        self.z_pol = (1-self.z_sol)*self.feed_exp

        # solvent + polymers chain lengths
        self.r = np.concatenate ([self.r_sol, self.r_pol])

        # solvent + polymers extensive volumes
        self.z = np.concatenate ([self.z_sol, self.z_pol])

        # number of components
        self.ncomp = len(self.r)

        # number of polymer pseudocomponents
        self.npol  = len(self.r_pol)    

    #########################

    def shulz_flory (self,r,q1,q2,alpha):
        #output-> mass fraction of a component with a chain size of r (mass fraction=(mass of polymer of size r in phase j)/(total polýmer mass of polymer in phase j))
        #input-> chain size (r)
        #input-> shultz flory parameter (q1 and q2) (probability of propagation of a polymer chain in an certain reactive site (site 1 or site 2) )
        #input-> relative activity of a catalitic site (alpha), in this case there are only two sites therfore de second site relative activity will be alpha-1
        return ((alpha*(1-q1)*(1-q1)*(q1**(r-1)))*r + 
                ((1-alpha)*(1-q2)*(1-q2)*(q2**(r-1)))*r)

    #########################

    def min_gibbs_energy (self):

        
        #porque esses limites?????????????????????<-------------------------------------------------------------------------------
        bounds = [(1e-10, self.z[0]-1e-10)]

        for m in range(1,self.ncomp):
            bounds.append((1e-10,self.z[m]-1e-10))


        # gibbs energy minimiztion for equilibrium calculations
        #differential evolution method is an genetic algorithm minimization strategy
        #for details look at: https://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.optimize.differential_evolution.html
        self.result_opt = opt.differential_evolution(self.gibbs_energy, 
                                                     bounds, popsize=40,
                                                     strategy = 'best2bin',
                                                     tol=1e-8, maxiter = 2000)

        zI = self.result_opt.x
        zII = self.z - zI

        self.phiI = zI/sum(zI)
        self.phiII = zII/sum(zII)

        return self.result_opt

    #########################
    
    def solve_equilibrium_equations (self, x0):

        self.result_eq = opt.root(self.equilibrium_equations, 
                                  x0, method='hybr')
        return self.result_eq

    #########################

    def estimation_objF (self, A):
        
        self.A = A
                
        if self.equilibrium_method == 'equations':
        
            if A>=0.51 and A<0.515:
                x0 = [-0.001, 0.0001, 5, 0.9]
            if A>=0.515 and A<0.52:
                x0 = [-0.03, 0.0003, 20, 0.5]
            if A>=0.52 and A<0.53:
                x0 = [-0.05, 0.0005, 35, 0.35]
            elif A>=0.53 and A<0.54:
                x0 = [-0.08, 0.0017, 90, 0.25]
            if A>=0.54 and A<0.55:
                x0 = [-0.12, 0.003, 150, 0.19]
            if A>=0.55 and A<0.56:
                x0 = [-0.15, 0.004, 225, 0.15]

            result = self.solve_equilibrium_equations(x0)
            
        elif self.equilibrium_method == 'optimize':  
            
            result = self.min_gibbs_energy()
                
        if result.success:
            if self.phiI[0] > self.phiII[0]:
                xs_model = self.phiI[1:]/sum(self.phiI[1:])
                ins_model = self.phiII[1:]/sum(self.phiII[1:])
            else:
                xs_model = self.phiII[1:]/sum(self.phiII[1:])
                ins_model = self.phiI[1:]/sum(self.phiI[1:])
            return (sum((xs_model-self.xs_exp)**2) + 
                    sum((ins_model-self.ins_exp)**2))
        else:
            return 1000

    #########################

    def estimation (self):

        bounds = [(0.51, 0.56)]
        result_estimation = opt.differential_evolution(self.estimation_objF, 
                                                       bounds)
        self.A = result_estimation.x[0]
        self.estimation_objF(self.A)
        return result_estimation

    #########################

    def plot_Experimental_Distributions (self,ax):

        ax.plot(self.r_pol, self.feed_exp, 
                 'k', label = 'Feed')
        ax.plot(self.r_pol, self.xs_exp, 
                 'b', label = 'Solubles, experimental')
        ax.plot(self.r_pol, self.ins_exp, 
                 'r', label = 'Insolubles, experimental')
        ax.set_title(self.label)
        ax.set_xlabel('Chain length, $r_i$')
        ax.set_ylabel('Mass fraction, $w_i$')
        ax.legend()

    #########################

    def plot_Calculated_Distributions (self,ax):

        if self.phiI[0] > self.phiII[0]:
            ax.plot(self.r_pol,self.phiI[1:]/sum(self.phiI[1:]),
                     '#00B0B2', label = 'Solubles, calculated')
            ax.plot(self.r_pol,self.phiII[1:]/sum(self.phiII[1:]),
                     '#FF6900',  label = 'Insolubles, calculated')
            ax.set_title(self.label+', A = '+str(self.A))
            ax.set_xlabel('Chain length, $r_i$')
            ax.set_ylabel('Mass fraction, $w_i$')
            ax.legend()
        else:
            ax.plot(self.r_pol,self.phiII[1:]/sum(self.phiII[1:]),
                     '#00B0B2',  label = 'Solubles, calculated')
            ax.plot(self.r_pol,self.phiI[1:]/sum(self.phiI[1:]),
                     '#FF6900',  label = 'Insolubles, calculated')
            ax.set_xlabel('Chain length, $r_i$')
            ax.set_ylabel('Mass fraction, $w_i$')
            ax.set_title(self.label+', A = '+str(self.A))
            ax.legend()

    #########################

    @abstractmethod
    def generate_pseudo_comp():
        pass

    #########################

    @abstractmethod
    def gibbs_energy ():
        pass
    
    #########################

    @abstractmethod
    def equilibrium_equations ():
        pass

#
#wi=(massa de polímero de tamanho i em determinada fase)/(massa total do polimero naquela fase)
#Xs=(massa de polímero obtida da fase líquida do teste xs)/(massa de polímero obtida da fase líquida do teste xs + massa de polímero obtida da fase solida do teste xs)
#Xs=(massa de polímero obtida da fase líquida do teste xs)/(massa total de polimero utilizada no teste xs)
#assim podemos dizer que Xs=mfl/(mfl+mfs)=mfl/mfa
#        
# wifs=fracao massica do polimero de tamanho i na fase solida obtida no teste xs
# wifl=fracao massica do polimero de tamanho i na fase liquida obtida no teste xs
# wifa=fracao masica do polimero de tamanho i na alimentacao (feed) 
# a distribuicao das fracoes massicas dos polimeros em cada fase sao obtidads da analise gpc         
# mfs=massa total dos polimeros ma fase solida obtida no teste xs 
# mfl=massa total dos polimeros ma fase liquida obtida no teste xs  
# mfa=massa total de polimero utilizada no teste xs 
#     
#Balanco de massa
#     
# wifx * mfs + wifl * mfl = wifa * mfa
#     
#mas
#   mfl= Xs * mfa
#     
# logo
# wifs * mfs = wifa * mfa - wifl * Xs * mfa 
#     
# wifs * mfs = (wifa - wifl * Xs) * mfa 
#     
# wifs = (wifa - wifl * Xs) * (mfa/mfs) 
#     
# wifs = (wifa - wifl * Xs) / (mfs/mfa) 
#    
# entretanto podemos dizer que:
#     
#    (mfs/mfa) + (mfl/mfa) = 1  
#    (mfs/mfa) + Xs = 1
#    (mfs/mfa) = 1 - Xs
#     
# Por fim obtemos    
#     
#   wifs = (wifa - wifl * Xs) / (1 - Xs)  
#     
#     Essa e a forma que esta no codigo
#     
#     
# 
#         
#    
#    
#    
#    
#    
#    
#   