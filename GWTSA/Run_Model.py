# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 11:39:58 2014
this is a time series analysis model
@author: R.A. Collenteur
"""

#Import all the packages needed
from GWTSA import *
import glob

#bores = glob.glob('Test_Data/*.csv')

plt.close('all')

bores = ['Test_Data/B27C0049001_1.csv']
#bores = ['Test_Data/B27C0002001_1.csv', 'Test_Data/B27C0041001_1.csv', 'Test_Data/B27C0049001_1.csv', 'Test_Data/B27D0001001_1.csv', 'Test_Data/B33A0113001_1.csv', 'Test_Data/B33C0135001_1.csv', 'Test_Data/B33D0115001_1.csv', 'Test_Data/B33D0217001_1.csv', 'Test_Data/B40B0304001_1.csv']

forcing = 'Test_Data/KNMI_Bilt.txt'
       
ml = Model(bores,forcing, 3000)

#TFN = ['nonlinear', 'nonlinear', 'nonlinear', 'nonlinear', 'nonlinear', 'nonlinear', 'nonlinear', 'nonlinear', 'nonlinear', 'nonlinear']
TFN = ['nonlinear','linear','linear','linear','linear','linear','linear','linear','linear','linear']
#TFN = ['linear','nonlinear','linear','nonlinear','linear','linear','linear','linear','linear','linear']

X0 = Parameters()
#           (Name,  Value,  Vary,   Min,  Max,  Expr)
X0.add_many(('A',   3.1,    True,   None, None,  None),
           ('a',    2.5,    True,   None, None,  None),
           ('n',    1.6,    True,   None, None,  None),
           ('alpha',0.50,   True,   None, None,  None),
           ('scap', -1.0,   True,   None, None,  None),
           ('ksat', -2.0,   True,   None, None,  None),
           ('beta', 3.0,    True,   0.5,  3.0,   None),
           ('imax', -3.0,   False,  None, None,  None))
           #('f', -2.0,   True,   None, None,  None), )
           
ml.solve(TFN, X0, 0)

ml.plot(1)
ml.bores_list[0].plot_diagnostics()
ml.bores_list[0].plot_results()


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

