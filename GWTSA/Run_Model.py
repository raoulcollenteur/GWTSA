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

bore = 'Test_Data/B27B0081-001' #, 'B27B0238-002', 'B27B0081-001', 'B27C0002-001', 'B32F0002-001', 'B33A0113-001', 'B33C0140-001', 'B39E0117-001', 'B40B0304-001'
forcing = 'Test_Data/KNMI_Bilt'
 
ts = Model(bore, forcing, rows=[5,8])

X0 = {'A': 2.0,'a': 20, 'n': 1.5,'Alpha': 0.5, 'S_cap': -0.50, 'K_sat':-2, 'Beta': 1.5, 'D': -3, 'f': 0.8} # initial parameters

ts.solve('TFN4', X0, method = 1, correlation = 0)

''' Optimal parameter set
# (A, a, n, d, Alpha, S_cap, K_sat, Beta, D)
Parameters = [2.6, 10.0, 1.35, 7.0, 1.0, 0.0, -3.0, 2.0, -3.0]'''



#Timeserie.simulate(Parameters, 'TFN4')
#
#X0 = np.array([[2.0, 8.0, 1.0, 7.0, 0.0, 0.0, -4.0, 1.0, -4.0],[2.6, 10.0, 1.35, 7.0, 1.0, 0.0, -3.0, 2.0, -3.0]])
#
#mont = Timeserie.monte_carlo('TFN4', X0, 1000)
#
#for i in range(9):
#    plt.figure()
#    plt.plot(Timeserie.Parameters[:,i], Timeserie.result, '.')

ts.simulate('TFN4', ts.parameters_opt)
ts.plot_heads()
