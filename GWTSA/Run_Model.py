# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 11:39:58 2014
this is a time series analysis model
@author: R.A. Collenteur
"""

#Import all the packages needed
#from datetime import datetime
from GWTSA import *

# close the previous figures
#plt.close('all')

bore = 'Test_Data/B33A0113-001.txt' # 'B27B0238-002', 'B27B0081-001', 'B27C0002-001', 'B32F0002-001', 'B33A0113-001', 'B33C0140-001', 'B39E0117-001', 'B40B0304-001'
forcing = 'Test_Data/KNMI_Bilt.txt'
TFN = 'TFN4'
 
ts = setup(bore, forcing, rows=[5,8], timestart=2000)

X0 = {'A': 3.0,'a': 2.0, 'n': 1.6,'Alpha': 1.0, 'S_cap': -2.0, 'K_sat':-3.0, 'Beta': 2, 'D': -3, 'f': -0.1} # initial parameters

ts.solve(TFN, X0, method=0)


''' Optimal parameter set
# (A, a, n, d, Alpha, S_cap, K_sat, Beta, D)
Parameters = [2.6, 10.0, 1.35, 17.0, 1.0, 0.0, -3.0, 2.0, -3.0]'''

#X0 = np.array([[100.0, 0.0, 0.8, 6.0, 1.7, -0.2],[1000, 1000.0, 2.0, 7.0, 2.3, 0.0]])
#mont = ts.monte_carlo('TFN2', X0, 10000)

#for i in range(len(X0[0])):
#    plt.figure()
#    plt.plot(ts.parameters[:,i], ts.result, '.')
#    plt.ylim(0,10)

#plt.figure()
#plt.scatter(ts.parameters[ts.result<1,0], ts.parameters[ts.result<1,1], c=ts.result[ts.result<1], cmap='RdYlGn')
#plt.xlabel('A')
#plt.ylabel('a')
#plt.colorbar()

ts.plot_results()
ts.plot_diagnostics()

#Check normal distribution of innovations
#plt.figure()
#plt.hist(ts.innovations, 50)

#
#
#Parameters = [3.0, 2.0, 1.6, 25.0, 2.0, 0.0, -3.0, 2.0, -3.0]
#ts.simulate(TFN,Parameters)
#ts.plot_results()