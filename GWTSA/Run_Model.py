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

#plt.close('all')

#bores = ['Test_Data/B33A0113001_1.csv','Test_Data/B33A0113001_1.csv','Test_Data/B33A0113001_1.csv','Test_Data/B33A0113001_1.csv']
#bores = [ 'Test_Data/B27D0001001_1.csv', 'Test_Data/B27D0001001_1.csv',  'Test_Data/B27D0001001_1.csv', 'Test_Data/B27D0001001_1.csv']
#bores = ['Test_Data/B27C0049001_1.csv','Test_Data/B27C0049001_1.csv','Test_Data/B27C0049001_1.csv','Test_Data/B27C0049001_1.csv']
bores = [ 'Test_Data/B27D0001001_1.csv', 'Test_Data/B27C0049001_1.csv', 'Test_Data/B33A0113001_1.csv']
#bores = ['Test_Data/UH/B32C0058001_1.csv', 'Test_Data/UH/B32C0058001_1.csv'] 
bores = [ 'Test_Data/B27D0001001_1.csv']

forcing = 'Test_Data/KNMI_Bilt.txt'
calibration = ['01-01-1965', '31-12-1998']
validation = ['01-01-2000', '31-12-2004']
       
ml = Model(bores, forcing, calibration, validation)

TFN = ['linear']*len(bores)
TFN = ['preferential']*len(bores)
#TFN = ['percolation']*len(bores)
#TFN = ['combination']*len(bores)

#print ml.bores_list[0].calcSumax() #= 0.27 for Deelen

X0 = Parameters()
#           (Name,  Value,  Vary,   Min,  Max,  Expr)
X0.add_many(('A',   0.2,    True,   0.0,  None,  None),
            ('t_p',  2.1,    True,  0.0, 2.7,  None),
            ('n',    1.2,    True,   None, None,  None),
            ('alpha',1.0,    True,   0.01,  None,   None),
            # Recharge
            #('f',    1.0,    True,   0.0, 1.5,  None),
            #('mu',  150.0, True,  0.0,  None,  None),
            #('sig', 25.0,  True,  0.0,  None,  None),
            ('Srmax',0.27,  False,  None, None,  None),
            ('Imax', 1.5e-3, False,  None, None,  None),
            ('Kp',  -2.25,    True,   -4.0, -1.0,  None),
            #('Beta', 3.0,    True,   0.0, None,  None),
            ('gamma',3.0,   True,   0.0, None,  None),
            # Reclamation
            #('b',   3.0,    True,   None, None,  None),
            #('B',   -0.5,    True,   -0.5, 0.0,  None)
            # Linear Slope 
            #('slope', 0.0005,    True,   0.0001, None,  None),
            #('intercept', -1.5, -2.0,   0.50, None,  None)
            )
                     
ml.solve(X0, RM='percolation', method='leastsq')

interface_plot()
cyan = [120/255.,196./255,1.]

#ml.plot(modeled=0)

ml.bores_list[0].plot_results(savefig=1)
ml.bores_list[0].plot_diagnostics(savefig=0)

ml.bores_list[1].plot_results(savefig=1)
ml.bores_list[2].plot_results(savefig=1)

ml.bores_list[1].plot_diagnostics()
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
#plt.subplot(311)
#plt.title('B27D00010')
#ax = ml.bores_list[0].recharge.resample('A', how='sum').plot('bar', color = cyan)
#plt.subplot(312)
#plt.title('B27C00490')
#ml.bores_list[1].recharge.resample('A', how='sum').plot('bar', color = cyan)
#plt.subplot(313)
#plt.title('B33A01130')
#ax = ml.bores_list[2].recharge.resample('A', how='sum').plot('bar', color = cyan)
#ticks = ax.xaxis.get_ticklocs()
#ticklabels = ['1975', '1980', 1985, 1990, 1995, 2000, 2005]
#ax.xaxis.set_ticks(ticks[::5])
#ax.xaxis.set_ticklabels(ticklabels, rotation=0)
#plt.subplot(413)
#ml.plot()
#plt.savefig('Figures/recharge.eps', bbox_inches='tight') 

#%% Plot individual daily recharges
latex_plot()
plt.figure(figsize=(8.3,4))
plt.title('B27D00010 - Linear Recharge')
RU = ml.bores_list[0].recharge[ml.bores_list[0].recharge.index > '1997-01-01']
plt.bar(R.index, R.values, color = 'k', lw=1)

plt.fill_between(R.index, 0.0, 0.1, where=(R.index.month > 10), facecolor=cyan, alpha=0.5, label='winter', lw=0)
plt.fill_between(R.index, 0.0, 0.1, where=(R.index.month <4), facecolor=cyan, alpha=0.5, lw=0)
plt.ylim(0,0.035)
plt.ylabel('Recharge [m/d]')
plt.xlabel('Time [years]')
p = plt.Rectangle((0, 0), 1, 1, fc=cyan)
plt.legend([p], ['November-March'])
plt.savefig('Figures/recharge_detail.eps', bbox_inches='tight') 

ticks = ax.xaxis.get_ticklocs()
ax.xaxis.set_ticks(ticks[::365])
ax.xaxis.set_ticklabels(ticklabels, rotation=0)
plt.savefig('Figures/recharge.eps', bbox_inches='tight') 

