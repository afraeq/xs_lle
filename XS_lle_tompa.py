#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Determination of Flory-Huggins Interaction Parameters 
for Polypropylene+Xylene Polydisperse Solutions from
Xylene Solubles (XS) Analysis Data

@author: afranio
         afrjr.weebly.com

Flory-Huggins model as provided by TOMPA, 1950 -
On The Precipitation Threshold of Solutions of Heterogeneous Polymers

Equilibrium equations based on HEIDEMANN et. al., 2006 -
An approach to expediting phase equilibrium calculations for
polydisperse polymers

"""

import numpy as np
from XS_lle import XS_lle

class XS_lle_tompa (XS_lle):

    def __init__ (self, label, z_sol, teta1, teta2, alpha, mw, xs_fraction,
                  equilibrium_method = 'equations'):

        super().__init__ (label, z_sol, teta1, teta2, alpha, mw, xs_fraction)
        
        self.equilibrium_method = equilibrium_method
        
        # quotients used in equilibrium equations
        self.quoc1 = np.zeros(self.npol)
        self.quoc2 = np.zeros(self.npol)

        for i in range(self.npol):
            self.quoc1[i] = ((self.r_pol[i] - self.r_pol[0])/
                             (self.r_pol[self.npol-1]-self.r_pol[0]))
            self.quoc2[i] = ((self.r_pol[self.npol-1] - self.r_pol[i])/
                             (self.r_pol[self.npol-1]-self.r_pol[0]))

        # equilibrium factors
        self.K   = np.zeros(self.ncomp)
        self.lnK = np.zeros(self.ncomp)

        # compositions of phases I and II
        self.phiI  = np.zeros(self.ncomp)
        self.phiII = np.zeros(self.ncomp)

        # mole numbers of phases I and II
        self.nI  = np.zeros(self.ncomp)
        self.nII = np.zeros(self.ncomp)
        
    #########################
    
    def generate_pseudo_comp (self):
        pass

    #########################

    def equilibrium_equations (self,x):

        self.lnK[0] = x[0]

        for i in range(1,self.ncomp):
            self.lnK[i] = self.quoc1[i-1] * x[2] + self.quoc2[i-1] * x[1]

        self.K = np.exp(self.lnK)

        for i in range(self.ncomp):
            self.phiI[i] = self.K[i]*self.z[i]/(self.K[i]*x[3]+(1-x[3]))
            self.phiII[i] = self.z[i]/(self.K[i]*x[3]+1-x[3])

        soma1 = 0
        soma2 = 0

        for i in range(self.npol):           
            soma1 += (((self.r_pol[i]-self.r[0])/
                          self.r_pol[i])*self.phiI[i+1])
            soma2 += (((self.r_pol[i]-self.r[0])/
                          self.r_pol[i])*self.phiII[i+1])

        y1 = (x[0] + 
              soma1 - 
              soma2 + 
              self.A*((self.r[0]-self.phiI[0])**2-
                      (self.r[0]-self.phiII[0])**2))
        y2 = ((self.r[0]/self.r[1])*x[1] - 
              x[0] + 
              2*self.A*(self.phiI[0]-self.phiII[0]))
        y3 = ((self.r[0]/self.r[self.ncomp-1])*x[2] - 
              x[0] + 
              2*self.A*(self.phiI[0]-self.phiII[0]))
        y4 = sum(self.phiI) - sum(self.phiII)

        return [y1,y2,y3,y4]

    #########################

    def gibbs_energy (self, x):

        for i in range(self.ncomp):
            self.phiI[i]  = x[i]/sum(x)
            self.phiII[i] = (self.z[i]-x[i])/sum(self.z-x)
            self.nI[i]    = x[i]/self.r[i]
            self.nII[i]   = self.z[i]/self.r[i] - self.nI[i]

        soma1 = 0
        soma2 = 0
        soma3 = 0
        soma4 = 0

        for i in range(self.npol):
            soma1 +=((self.r_pol[i]-self.r[0])/self.r_pol[i])*self.phiI[i+1]
            soma2+=((self.r_pol[i]-self.r[0])/self.r_pol[i])*self.phiII[i+1]
            soma3 += (self.r_pol[i]-self.r[0])*self.phiI[i+1]
            soma4 += (self.r_pol[i]-self.r[0])*self.phiII[i+1]

        delta_muI_0  = (np.log(self.phiI[0]) + 
                        soma1 + 
                        self.A*(self.r[0]-self.phiI[0])**2)
        delta_muII_0 = (np.log(self.phiII[0]) + 
                        soma2 + 
                        self.A*(self.r[0]-self.phiII[0])**2)

        dG = self.nI[0]*delta_muI_0 + self.nII[0]*delta_muII_0

        for i in range(1,self.ncomp):

            delta_muI  = (np.log(self.phiI[i]) + 
                          soma3 + 
                          self.A*self.r[i]*(self.phiI[0])**2)
            delta_muII = (np.log(self.phiII[i]) + 
                          soma4 + 
                          self.A*self.r[i]*(self.phiII[0])**2)

            dG += self.nI[i]*delta_muI + self.nII[i]*delta_muII

        return dG
