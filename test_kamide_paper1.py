#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 20:53:15 2019

@author: afranio
"""

import numpy as np
import matplotlib.pyplot as plt

from XS_lle_kamide import XS_lle_kamide

nclass = 15
rmin = 10
rmax = 40000

####################
# EXPERIMENTAL DATA
####################

z_sol = 0.97738

labels = ['H1_9','H2_31','H3_8','H4_32','H5_29','H6_13','H8_11']

teta1 = 1e-4*np.array([[1.29, 0.000], [1.55, 0.0], [1.58, 1.163],
                       [1.63, 1.266], [1.37, 0.0], [1.39, 1.406],
                       [1.85, 0.0]])

teta2 = 1e-4*np.array([[6.08, 6.348], [5.6, 19.44], [6.57, 13.07],
                      [0.00, 10.52], [5.2, 12.38], [7.40, 16.20],
                      [6.30, 24.07]])

alpha = np.array([[0.868, 0.000],[0.237, 0.000],[0.549, 0.053],
                  [1.000, 0.145],[0.839, 0.000],[0.655, 0.083],
                  [0.353, 0.000]])

mw = np.array([582261, 242779, 349674, 514920, 540133, 435169, 246094])

xs = 1e-2*np.array([3.98, 4.58, 4.59, 3.43, 6.67, 6.48, 4.53])

# creating list of objects for each polymer
polymers = [XS_lle_kamide(labels[i],z_sol,teta1[i],teta2[i],
                          alpha[i],mw[i],xs[i], nclass, rmin, rmax) 
            for i in range(len(labels))]

#################################
# RUNNING TEST FLASH CALCULATION
#################################

A = 0.535

p = polymers[0]

p.estimation_objF(A)
p.plot_Experimental_Distributions(plt.gca())
p.plot_Calculated_Distributions(plt.gca())