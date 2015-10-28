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

#plt.close('all')

bores = ['Test_Data/B27C0002001_1.csv']
#bores = ['Test_Data/B27C0002001_1.csv', 'Test_Data/B27C0041001_1.csv', 'Test_Data/B27C0049001_1.csv', 'Test_Data/B27D0001001_1.csv', 'Test_Data/B33A0113001_1.csv', 'Test_Data/B33C0135001_1.csv', 'Test_Data/B33D0115001_1.csv', 'Test_Data/B33D0217001_1.csv', 'Test_Data/B40B0304001_1.csv']

forcing = 'Test_Data/KNMI_Bilt.txt'
       
ml = Model(bores,forcing, 6000)

#TFN = ['nonlinear', 'nonlinear', 'nonlinear', 'nonlinear', 'nonlinear', 'nonlinear', 'nonlinear', 'nonlinear', 'nonlinear', 'nonlinear']
TFN = ['linear','nonlinear','nonlinear1','nonlinear2','linear','linear','linear','linear','linear','linear']
#TFN = ['linear','nonlinear','linear','nonlinear','linear','linear','linear','linear','linear','linear']

#print ml.bores_list[0].calcSumax()

X0 = Parameters()
#           (Name,  Value,  Vary,   Min,  Max,  Expr)
#X0.add_many(('A',   3.4,    True,   None, None,  None),
#           ('a',    2.5,    True,   None, None,  None),
#           ('n',    1.5,    True,   None, None,  None),
#           ('alpha',2.0,    True,   None, None,  None),
#           ('scap', -0.55,  False,   None, None,  None),
#           ('ksat', -2.0,   True,   None, None,  None),
#           ('beta', 2.0,    True,   None,  None,  None),
#           ('imax', 1e-3,   False,  None, None,  None))

X0.add_many(('A',   -1.4,    True,   None, None,  None),
           ('a',    2.5,    True,   None, None,  None),
           ('n',    1.5,    True,   None, None,  None),
           ('alpha',2.0,    True,   None, None,  None),
           ('f',    -2.0,   True,  None, None,  None))
           
#           
ml.solve(TFN, X0, 0)
latex_plot()
#ml.plot(1)
ml.bores_list[0].plot_diagnostics()
ml.bores_list[0].plot_results()


#Check normal distribution of innovations
#plt.figure()
#plt.hist(ts.innovations, 50)

R = ml.bores_list[0].result


