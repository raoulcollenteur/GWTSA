# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 11:39:58 2014
this is a time series analysis model
@author: R.A. Collenteur
"""

#Import all the packages needed
#from datetime import datetime
import matplotlib.pyplot as plt
from GWTSA import *

# close the previous figures
plt.close('all')

bore = 'Test_Data/GW_Levels' #, 'B27B0238-002', 'B27B0081-001', 'B27C0002-001', 'B32F0002-001', 'B33A0113-001', 'B33C0140-001', 'B39E0117-001', 'B40B0304-001'
forcing = 'Test_Data/KNMI_Bilt'
 
Timeserie = Model(bore, forcing, rows=[1,8])

X0 = {'A': 2.6,'a': 10, 'n': 1.35,'Alpha': 1.0, 'S_cap': 0.0, 'K_sat':-3, 'Beta': 2.0, 'D': -3, 'f': 0.8} # initial parameters

Timeserie.solve('TFN4', X0, Opt = 0, Cor = 0)

Timeserie.plot_heads('r')

#TFN = 'TFN4' #TFN-model to use
#
#Bore1.Model(TFN,Par)
#Bore1.Plot_Heads('b')

#Timeserie.plot_forcings()

''' Optimal parameter set
# (A, a, n, d, Alpha, S_cap, K_sat, Beta, D)'''
Parameters = [400.0, 10.0, 1.35, 7.0, 1.0, 0.0, -3.0, 2.0, -3.0]



Timeserie.simulate(Parameters, 'TFN4')