# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 11:39:58 2014
this is a time series analysis model
@author: R.A. Collenteur
"""

#Import all the packages needed
from GWTSA import *

#plt.close('all')

bores = ['Test_Data/B27B0238-002.txt', 'Test_Data/B27B0081-001.txt', 'Test_Data/B27C0002-001.txt', 'Test_Data/B32F0002-001.txt']#, 'Test_Data/B33A0113-001.txt', 'Test_Data/B33C0140-001.txt', 'Test_Data/B39E0117-001.txt', 'Test_Data/B40B0304-001.txt']

forcing = ['Test_Data/KNMI_Bilt.txt', 'Test_Data/KNMI_Bilt.txt','Test_Data/KNMI_Bilt.txt', 'Test_Data/KNMI_Bilt.txt','Test_Data/KNMI_Bilt.txt', 'Test_Data/KNMI_Bilt.txt','Test_Data/KNMI_Bilt.txt', 'Test_Data/KNMI_Bilt.txt']
       
ml = Model(bores,forcing)

TFN = ['TFN2','TFN2','TFN2','TFN2','TFN2','TFN2','TFN2','TFN2']
#TFN = ['TFN4','TFN4','TFN4','TFN4','TFN4','TFN4','TFN4','TFN4']
X0 = {'A': 2.6,'a': 2.0, 'n': 1.6,'Alpha': 150, 'S_cap': -2.0, 'K_sat':-3.0, 'Beta': 3.0, 'D': -3, 'f': -0.1} # initial parameters

ml.solve(TFN, X0, 0)

ml.plot()

#ml.bores_list[3].plot_results()
#ml.bores_list[3].plot_diagnostics()

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

#Check normal distribution of innovations
#plt.figure()
#plt.hist(ts.innovations, 50)

#Parameters = [3.0, 2.0, 1.6, 25.0, 2.0, 0.0, -3.0, 2.0, -3.0]
#ts.simulate(TFN,Parameters)
#ts.plot_results()