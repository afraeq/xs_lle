import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
from abc import ABC, abstractmethod

class PolyLLE (ABC):

    def __init__ (self, label = 'pol1', z_sol = 0.98,
                  shulz_kind = '1P', r_pol = np.arange(10,10**4), 
                  p = [0.999,0.995],
                  teta1 = None, teta2 = None,
                  coarsen = True,
                  alpha = None, xs_fraction = 0.05,
                  equilibrium_method = 'equations'):

        # polymer label
        self.label = label

        # solvent extensive volume
        self.z_sol = np.array([z_sol])

        # shulz-flory parameters
        # [0]: feed, [1]: solubles

        if shulz_kind == '1P':
            self.p = p
        elif shulz_kind == '3P':
            self.q1_exp = np.exp(-teta1)
            self.q2_exp = np.exp(-teta2)
            self.alpha_exp = alpha

        # xs fraction: polymer fraction extracted by xylene
        self.xs_fraction = xs_fraction

        # number of phases
        self.nf = 2
        
        # solvent chain length
        self.r_sol = np.array([1])

        # polymer chain lengths

        self.r_pol = r_pol

        if shulz_kind == '1P':

            # experimental feed distribution

            dist = self.shulz_flory_1P(self.r_pol, self.p[0])

            self.feed_exp = dist/sum(dist)

            # experimental solubles distribution

            dist = self.shulz_flory_1P(self.r_pol, self.p[1])

            self.xs_exp = dist/sum(dist)

        elif shulz_kind == '3P':

            # experimental feed distribution

            dist = self.shulz_flory_3P(self.r_pol,
                                    self.q1_exp[0],
                                    self.q2_exp[0],
                                    self.alpha_exp[0])

            self.feed_exp = dist/sum(dist)

            # experimental solubles distribution

            dist = self.shulz_flory_3P(self.r_pol,
                                    self.q1_exp[1],
                                    self.q2_exp[1],
                                    self.alpha_exp[1])

            self.xs_exp = dist/sum(dist)

        # experimental insolubles distribution
        self.ins_exp = ((self.feed_exp-self.xs_fraction*self.xs_exp)/
                        (1-self.xs_fraction))
        for i in range(len(self.ins_exp)):
            if self.ins_exp[i]<0:
                self.ins_exp[i] = 0
                
        # generating pseudo components
        if coarsen:
            self.generate_pseudo_comp()

        # polymer extensive volumes
        self.z_pol = (1-self.z_sol)*self.feed_exp

        # solvent + polymers chain lengths
        self.r = np.concatenate ([self.r_sol, self.r_pol])

        # solvent + polymers extensive volumes
        self.z = np.concatenate ([self.z_sol, self.z_pol])

        # number of components
        self.n_comp = len(self.r)

        # number of polymer pseudocomponents
        self.n_pol  = len(self.r_pol)    

        # optimize or equations
        self.equilibrium_method = equilibrium_method

    #########################

    def shulz_flory_3P (self,r,q1,q2,alpha):
        return ((alpha*(1-q1)*(1-q1)*(q1**(r-1)))*r + 
                ((1-alpha)*(1-q2)*(1-q2)*(q2**(r-1)))*r)

    #########################

    def shulz_flory_1P (self, r, p):
        return r*((1-p)**2)*p**(r-1)

    #########################

    def min_gibbs_energy (self):

        lb = 1e-6*np.ones((self.n_comp,))
        ub = self.z - 1e-10

        bounds = [(lb[i], ub[i]) for i in range(self.n_comp)]

        self.result_opt = scipy.optimize.direct(self.gibbs_energy, bounds)

        zI = self.result_opt.x      
        zII = self.z - zI

        self.phiI = zI/sum(zI)
        self.phiII = zII/sum(zII)

        return self.result_opt

    #########################
    
    def solve_equilibrium_equations (self, x0):

        self.result_eq = scipy.optimize.root(self.equilibrium_equations, 
                                             x0, method='hybr')
        return self.result_eq

    #########################

    def estimation_objF (self, A):
        
        self.A = A
                
        if self.equilibrium_method == 'equations':
        
            if A<0.515:#A>=0.51 and 
                x0 = [-0.001, 0.0001, 5, 0.9]
            if A>=0.515 and A<0.52:
                x0 = [-0.03, 0.0003, 20, 0.5]
            if A>=0.52 and A<0.53:
                x0 = [-0.05, 0.0005, 35, 0.35]
            elif A>=0.53 and A<0.54:
                x0 = [-0.08, 0.0017, 90, 0.25]
            if A>=0.54 and A<0.55:
                x0 = [-0.12, 0.003, 150, 0.19]
            if A>=0.55:# and A<0.56:
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

        #bounds = [(0.51, 0.56)]
        #result_estimation = scipy.optimize.differential_evolution(
        #                    self.estimation_objF, bounds)
        self.result_estimation = \
            scipy.optimize.minimize_scalar(self.estimation_objF, 
                                           bounds=(0.51,0.56),
                                           method='bounded')
        self.A = self.result_estimation.x#[0]
        self.estimation_objF(self.A)
        return self.result_estimation

    #########################

    def plot_Experimental_Distributions (self,ax=None):

        if ax is None:
            ax=plt.gca()

        ax.plot(self.r_pol, self.feed_exp, 
                 'k',label = 'Feed')
        ax.plot(self.r_pol, self.xs_exp, 
                 'b', label = 'Solubles, experimental')
        ax.plot(self.r_pol, self.ins_exp, 
                 'r', label = 'Insolubles, experimental')
        ax.set_title(self.label)
        ax.set_xlabel('Chain length, $r_i$')
        ax.set_ylabel('Mass fraction, $w_i$')
        ax.legend()

    #########################

    def plot_Calculated_Distributions (self,ax = None):

        if ax is None:
            ax=plt.gca()

        ax.plot(self.r_pol, self.feed_exp, 
                 'k', label = 'Feed')

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




