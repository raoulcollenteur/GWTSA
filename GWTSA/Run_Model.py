# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 11:39:58 2014
this is a time series analysis model
@author: R.A. Collenteur
"""

#Import all the packages needed
from GWTSA import *

import glob
bores = glob.glob('Test_Data/All/*.csv')

plt.close('all')

#bores = ['Test_Data/B33A0113001_1.csv','Test_Data/B33A0113001_1.csv','Test_Data/B33A0113001_1.csv','Test_Data/B33A0113001_1.csv']
#bores = [ 'Test_Data/B27D0001001_1.csv', 'Test_Data/B27D0001001_1.csv',  'Test_Data/B27D0001001_1.csv', 'Test_Data/B27D0001001_1.csv']
#bores = ['Test_Data/B27C0049001_1.csv','Test_Data/B27C0049001_1.csv','Test_Data/B27C0049001_1.csv','Test_Data/B27C0049001_1.csv']
#bores = [ 'Test_Data/B27D0001001_1.csv', 'Test_Data/B27C0049001_1.csv', 'Test_Data/B33A0113001_1.csv']
#bores = ['Test_Data/UH/B32C0058001_1.csv', 'Test_Data/UH/B32C0058001_1.csv'] 
#bores = [B33C0135001_1]

forcing = 'Test_Data/KNMI_Bilt.txt'
calibration = ['01-01-1960', '31-12-1990']
validation = ['01-01-1980', '31-12-2004']
       
ml = Model(bores, forcing, calibration, validation)
sl=ml

TFN = ['linear']*len(bores)
TFN1 = ['preferential']*len(bores)

#TFN = ['linear', 'preferential', 'percolation', 'combination']


#print ml.bores_list[0].calcSumax() = 0.27 for Deelen

X1 = Parameters()
#           (Name,  Value,  Vary,   Min,  Max,  Expr)
X1.add_many(('A',   2.4,    True,   0.0,  None,  None),
           ('a',    2.0,    True,   None, None,  None),
           #('mu',    150.0, True,  0.0,  None,  None),
           #('sig',   25.0,  True,  0.0,  None,  None),
           ('n',    1.5,    True,   None, None,  None),
           ('alpha',2.0,    True,   0.01,  3.0,   None),
           ('Srmax', 0.27,  False,  None, None,  None),
           ('Kp', -2.25,    True,   -4.0, -1.0,  None),
           ('beta', 3.0,    True,   0.0, None,  None),
           #('gamma', 3.0,   True,   0.0, None,  None),
           #('f',    1.0,    True,   0.0, 1.5,  None),
           ('imax', 1.5e-3, False,  None, None,  None))

X0 = Parameters()
#           (Name,  Value,  Vary,   Min,  Max,  Expr)
X0.add_many(('A',   3.4,    True,   None, None,  None),
           ('a',    2.5,    True,   None, None,  None),
           ('n',    1.5,    True,   None, None,  None),
           ('alpha',2.0,    True,   0.0, 3.2,  None),
           ('f',    1.0,    True,   0.0, 1.5,  None))
                     
ml.solve(TFN, X0, method='leastsq')
sl.solve(TFN1,X1, method='leastsq')

interface_plot()
cyan = [120/255.,196./255,1.]

#ml.plot(savefig=True)

ml.bores_list[0].plot_results(savefig=1)
ml.bores_list[1].plot_results(savefig=1)
ml.bores_list[2].plot_results(savefig=1)
#ml.bores_list[3].plot_results(savefig=1)

#ml.bores_list[0].plot_diagnostics(savefig=0)
#ml.bores_list[1].plot_diagnostics()
#ml.bores_list[2].plot_diagnostics()

print np.mean(ml.bores_list[0].recharge.resample('A', how='sum'))
print np.mean(ml.bores_list[1].recharge.resample('A', how='sum'))
print np.mean(ml.bores_list[2].recharge.resample('A', how='sum'))
#print np.mean(ml.bores_list[3].recharge.resample('A', how='sum'))

#plt.plot(ml.bores_list[0].head_modeled[ml.bores_list[0]._index_observed],ml.bores_list[0].head_observed, 'bo')
#plt.plot([18,21.5],[18,21.5])

#%% Plot the recharges
#
#plt.figure(figsize=(8.3,4))
#
#plt.subplot(411)
#plt.title('B27D00010')
#ml.bores_list[0].recharge.resample('A', how='sum').plot('bar', color = cyan)
#plt.subplot(412)
#plt.title('B27C00490')
#ml.bores_list[1].recharge.resample('A', how='sum').plot('bar', color = cyan)
#plt.subplot(413)
#plt.title('B33A01130')
#ax = ml.bores_list[2].recharge.resample('A', how='sum').plot('bar', color = cyan)
#
#ticks = ax.xaxis.get_ticklocs()
#ticklabels = ['1975', '1980', 1985, 1990, 1995, 2000, 2005]
#ax.xaxis.set_ticks(ticks[::5])
#ax.xaxis.set_ticklabels(ticklabels, rotation=0)
#plt.subplot(413)
#ml.plot()
#
#plt.savefig('Figures/recharge.eps', bbox_inches='tight') 
