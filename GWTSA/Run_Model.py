# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 11:39:58 2014
this is a time series analysis model
@author: R.A. Collenteur
"""

#Import all the packages needed
from GWTSA import *

#import glob
#bores = glob.glob('Test_Data/*.csv')

#plt.close('all')

#bores = ['Test_Data/B33A0113001_1.csv','Test_Data/B33A0113001_1.csv','Test_Data/B33A0113001_1.csv','Test_Data/B33A0113001_1.csv']
#bores = [ 'Test_Data/B27D0001001_1.csv', 'Test_Data/B33D0217001_1.csv', 'Test_Data/B27C0049001_1.csv', 'Test_Data/B27C0041001_1.csv']#, 'Test_Data/B33A0113001_1.csv']
bores = ['Test_Data/B33C0135001_1.csv']


forcing = 'Test_Data/KNMI_Bilt.txt'
calibration = ['01-01-1960', '31-12-2004']
validation = ['01-01-2005', '31-12-2005']
       
ml = Model(bores, forcing, calibration, validation)

#TFN = ['linear','linear','linear','linear']
#TFN = ['preferential','preferential','preferential']
TFN = ['percolation','percolation','percolation','percolation']
#TFN = ['combination','combination','combination','combination']
#TFN = ['linear', 'preferential', 'percolation', 'combination']


#print ml.bores_list[0].calcSumax()

X1 = Parameters()
#           (Name,  Value,  Vary,   Min,  Max,  Expr)
X1.add_many(('A',   3.4,    True,   0.0,  None,  None),
           ('a',    2.0,    True,   None, None,  None),
           ('b',    0.0,    False,  0.0,  None,  None),
           ('n',    1.5,    True,   None, None,  None),
           ('alpha',2.0,    True,   0.0,  3.0,   None),
           ('Srmax', 0.28,  False,  None, None,  None),
           ('Kp', -2.25,    True,   None, -1.0,  None),
           ('beta', 1.0,    True,   None, None,  None),
           ('gamma', 1.0,   True,   None, None,  None),
           ('f',    1.0,    True,   None, 1.5,  None),
           ('imax', 1.5e-3, False,  None, None,  None))
#           (Name,  Value,  Vary,   Min,  Max,  Expr)
#X1.add_many(('A',   3.487,    False,   0.0,  None,  None),
#           ('a',    2.795,    False,   None, None,  None),
#           ('b',    0.0,    False,  0.0,  None,  None),
#           ('n',    1.403,    False,   None, None,  None),
#           ('alpha', 2.969,    True,   0.0,  3.0,   None),
#           ('Srmax', 0.28, False,  None, None,  None),
#           ('Kp', -1.721,    False,   None, -1.0,  None),
#           ('beta', 1.0,    False,   None, None,  None),
#           ('gamma', 27.3850,   False,   None, None,  None),
#           ('f',    1.0,    False,   None, 1.5,  None),
#           ('imax', 3e-3,   True,  None, None,  None))


X0 = Parameters()
#           (Name,  Value,  Vary,   Min,  Max,  Expr)
X0.add_many(('A',   3.4,    True,   None, None,  None),
           ('a',    2.5,    True,   None, None,  None),
           ('n',    1.5,    True,   None, None,  None),
           ('alpha',2.0,    True,   None, 3.2,  None),
           ('f',    1.0,    True,   None, None,  None))
                     
ml.solve(TFN, X1, 'leastsq')

latex_plot()
ml.plot()

ml.bores_list[0].plot_results(savefig=1)
ml.bores_list[1].plot_results(savefig=1)
ml.bores_list[2].plot_results(savefig=1)
ml.bores_list[3].plot_results(savefig=1)

#ml.bores_list[0].plot_diagnostics(savefig=0)
#ml.bores_list[1].plot_diagnostics()
#ml.bores_list[2].plot_diagnostics()

print np.mean(ml.bores_list[0].recharge.resample('A', how='sum'))
print np.mean(ml.bores_list[1].recharge.resample('A', how='sum'))
print np.mean(ml.bores_list[2].recharge.resample('A', how='sum'))
print np.mean(ml.bores_list[3].recharge.resample('A', how='sum'))


plt.plot(ml.bores_list[0].head_modeled[ml.bores_list[0]._index_observed],ml.bores_list[0].head_observed, 'bo')
plt.plot([18,21.5],[18,21.5])


