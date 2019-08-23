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
from pyswarm import pso

class XS_lle (ABC):

    def __init__ (self, label, z_sol, teta1, teta2, alpha, mw, xs_fraction):

        # polymer label
        self.label = label

        # solvent extensive volume
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
        self.feed_exp = ((self.shulz_flory(self.r_pol,
                                           self.q1_exp[0],
                                           self.q2_exp[0],
                                           self.alpha_exp[0]))/
                         sum(self.shulz_flory(self.r_pol,
                                              self.q1_exp[0],
                                              self.q2_exp[0],
                                              self.alpha_exp[0])))

        # experimental xs distribution
        self.xs_exp = ((self.shulz_flory(self.r_pol,
                                        self.q1_exp[1],
                                        self.q2_exp[1],
                                        self.alpha_exp[1]))/
                       sum(self.shulz_flory(self.r_pol,
                                            self.q1_exp[1],
                                            self.q2_exp[1],
                                            self.alpha_exp[1])))

        # experimental insolubles distribution
        self.ins_exp = ((self.feed_exp-self.xs_fraction*self.xs_exp)/
                        (1-self.xs_fraction))
        for i in range(len(self.ins_exp)):
            if self.ins_exp[i]<0:
                self.ins_exp[i] = 0
                
        # generating pseudo components
        self.generate_pseudo_comp()

        # polymer extensive volumes
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
        return ((alpha*(1-q1)*(1-q1)*(q1**(r-1)))*r + 
                ((1-alpha)*(1-q2)*(1-q2)*(q2**(r-1)))*r)

    #########################

    def min_gibbs_energy (self):

        bounds = [(1e-10, self.z[0]-1e-10)]

        for m in range(1,self.ncomp):
            bounds.append((1e-10,self.z[m]-1e-10))
            
#        self.result_opt = opt.dual_annealing(self.gibbs_energy, 
#                                                     bounds)#, popsize=40,
#                                                     strategy = 'best2bin',
#                                                     tol=1e-8, maxiter = 2000)


#        print(self.result_opt.fun)
        
#        zI = self.result_opt.x

        lb = [1e-10 for _ in range(self.ncomp)]
        ub = [self.z[m]-1e-10 for m in range(self.ncomp)]

        xopt, fopt = pso(self.gibbs_energy, lb, ub, 
                         minfunc = 1e-16, minstep = 1e-14,
                         swarmsize = 100, maxiter = 200)        
        
#        print(fopt)
        zI = xopt    
        
        zII = self.z - zI

        self.phiI = zI/sum(zI)
        self.phiII = zII/sum(zII)
        
        # to mantain compatibility with
        # scipy.optimize usage
        class results(object):
            pass
        
        self.result_opt = results()
        self.result_opt.success = True
        self.result_opt.x = xopt

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
