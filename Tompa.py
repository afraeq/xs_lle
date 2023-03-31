import numpy as np
from PolyLLE import PolyLLE

class Tompa (PolyLLE):

    def __init__ (self, label = '', z_sol=0.98,
                  shulz_kind = '1P', r_pol = np.arange(10,10**4), 
                  p = [0.999,0.995], 
                  teta1 = None, teta2 = None,
                  coarsen = False,
                  alpha = None, xs_fraction = 0.05,
                  equilibrium_method = 'equations'):

        super().__init__ (label, z_sol,
                          shulz_kind, r_pol, p, teta1, teta2, 
                          coarsen, alpha, xs_fraction, equilibrium_method)
                
        # quotients used in equilibrium equations
        self.quoc1 = np.zeros(self.n_pol)
        self.quoc2 = np.zeros(self.n_pol)

        for i in range(self.n_pol):
            self.quoc1[i] = ((self.r_pol[i] - self.r_pol[0])/
                             (self.r_pol[self.n_pol-1]-self.r_pol[0]))
            self.quoc2[i] = ((self.r_pol[self.n_pol-1] - self.r_pol[i])/
                             (self.r_pol[self.n_pol-1]-self.r_pol[0]))

        # equilibrium factors
        self.K   = np.zeros(self.n_comp)
        self.lnK = np.zeros(self.n_comp)

        # compositions of phases I and II
        self.phiI  = np.zeros(self.n_comp)
        self.phiII = np.zeros(self.n_comp)

        # mole numbers of phases I and II
        self.nI  = np.zeros(self.n_comp)
        self.nII = np.zeros(self.n_comp)
        
    #########################
    
    def generate_pseudo_comp (self):
        pass

    #########################

    def equilibrium_equations (self,x):

        self.lnK[0] = x[0]

        for i in range(1,self.n_comp):
            self.lnK[i] = self.quoc1[i-1] * x[2] + self.quoc2[i-1] * x[1]

        self.K = np.exp(self.lnK)

        for i in range(self.n_comp):
            self.phiI[i] = self.K[i]*self.z[i]/(self.K[i]*x[3]+(1-x[3]))
            self.phiII[i] = self.z[i]/(self.K[i]*x[3]+1-x[3])

        soma1 = 0
        soma2 = 0

        for i in range(self.n_pol):           
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
        y3 = ((self.r[0]/self.r[self.n_comp-1])*x[2] - 
              x[0] + 
              2*self.A*(self.phiI[0]-self.phiII[0]))
        y4 = sum(self.phiI) - sum(self.phiII)

        return [y1,y2,y3,y4]

    #########################

    def gibbs_energy (self, x):

        for i in range(self.n_comp):
            self.phiI[i]  = x[i]/sum(x)
            self.phiII[i] = (self.z[i]-x[i])/sum(self.z-x)
            self.nI[i]    = x[i]/self.r[i]
            self.nII[i]   = self.z[i]/self.r[i] - self.nI[i]

        soma1 = 0
        soma2 = 0
        soma3 = 0
        soma4 = 0

        for i in range(self.n_pol):
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

        for i in range(1,self.n_comp):

            delta_muI  = (np.log(self.phiI[i]) + 
                          soma3 + 
                          self.A*self.r[i]*(self.phiI[0])**2)
            delta_muII = (np.log(self.phiII[i]) + 
                          soma4 + 
                          self.A*self.r[i]*(self.phiII[0])**2)

            dG += self.nI[i]*delta_muI + self.nII[i]*delta_muII

        return dG
