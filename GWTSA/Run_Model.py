# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 11:39:58 2014
this is a time series analysis model
@author: R.A. Collenteur
"""

#Import all the packages needed
from GWTSA import *

import glob
bores = glob.glob('Test_Data/Selectie/*.csv')

plt.close('all')

#bores = ['Test_Data/B33A0229001_1.csv']
#bores = ['Test_Data/B27B0081001_1.csv', 'Test_Data/B27B0081001_1.csv']#['Test_Data/B33A0229001_1.csv', 'Test_Data/B33C0141001_0.csv', 'Test_Data/B33A0114002_1.csv','Test_Data/B33A0227001_1.csv', 'Test_Data/B33A0229001_1.csv']

forcing = 'Test_Data/KNMI_Bilt.txt'
       
ml = Model(bores,forcing)

TFN = ['TFN4','TFN4','TFN2','TFN2','TFN2','TFN2','TFN2','TFN2','TFN2','TFN2','TFN2','TFN2','TFN2','TFN2','TFN2','TFN2']
X0 = {'A': 3.2,'a': 2.5, 'n': 1.6,'Alpha': 0, 'S_cap': -2.0, 'K_sat':-3.0, 'Beta': 3.0, 'D': -3, 'f': -0.1} # initial parameters

ml.solve(TFN, X0, 2)

ml.plot()

#X0 = ml.bores_list[0].parameters_opt
#ml.solve(TFN, X0, 2)
#ml.plot(0)

#ml.bores_list[0].plot_results()
#ml.bores_list[1].plot_results()
#ml.bores_list[0].plot_diagnostics()

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
