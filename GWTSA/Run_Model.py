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
bores = [ 'Test_Data/B27D0001001_1.csv', 'Test_Data/B27C0049001_1.csv', 'Test_Data/B33A0113001_1.csv']
#bores = ['Test_Data/UH/B32C0058001_1.csv', 'Test_Data/UH/B32C0058001_1.csv'] 
#bores = [ 'Test_Data/B33A0113001_1.csv']

forcing = 'Test_Data/KNMI_Bilt.txt'
calibration = ['01-01-1974', '31-12-1994']
#calibration = ['01-01-1974', '31-12-2003']
validation = ['01-01-1995', '31-12-2004']
#validation = ['01-01-2004', '31-12-2004']
       
ml = Model(bores, forcing, calibration, validation)

trend = ['linear']*len(bores)
RM = ['preferential']*len(bores)
#TFN = ['percolation']*len(bores)
#TFN = ['combination']*len(bores)

#print ml.bores_list[0].calcSumax() #= 0.27 for Deelen

X0 = Parameters()
#           (Name,  Value,  Vary,   Min,  Max,  Expr)
X0.add_many(('A',  -1.0,    True,   None, 0.0,  None),
            ('a',  2.0,    True,  0.0, None,  None),
            ('n',    2.0,    True,   None, None,  None),
            #('A1',   0.4,    True,   0.0,  None,  None),
            #('t_p1',  2.4,    True,  0.0, 3.0,  None),
            #('n1',    1.5,    True,   None, None,  None),
            ('alpha',1.0,    True,   0.01,  3.0,   None),
            # Recharge
            #('f',    1.0,    True,   0.0, 1.5,  None),
            #('mu',  150.0, True,  0.0,  None,  None),
            #('sig', 25.0,  True,  0.0,  None,  None),
            ('Srmax',0.14,  False,  None, None,  None),
            ('Imax', 1.5e-3, False,  None, None,  None),
            ('Kp',  -2.25,    True,   -4.0, -1.0,  None),
            ('Beta', 1.0,    True,   0.0, None,  None),
            ('gamma',3.0,   True,   0.0, None,  None),
            # Reclamation
            #('b',   3.0,    True,   None, None,  None),
            #('B',   -0.5,    True,   -0.5, 0.0,  None)
            # Linear Slope 
            #('slope', 0.00015, False, 0.00001, 0.00025,  None),
            #('intercept', -1.10, False, None, None,  None)
            )

ml.solve(X0, RM= 'combination', method='leastsq')

latex_plot()
cyan = [120/255.,196./255,1.]

ml.plot(modeled=1, savefig=1)

ml.bores_list[0].plot_results(savefig=1)
ml.bores_list[0].plot_diagnostics(savefig=1)
#ml.bores_list[0].recharge_uncertainty()

ml.bores_list[1].plot_results(savefig=1)
ml.bores_list[1].plot_diagnostics(savefig=1)
#ml.bores_list[1].recharge_uncertainty()

ml.bores_list[2].plot_results(savefig=1)
ml.bores_list[2].plot_diagnostics(savefig=1)
#ml.bores_list[2].recharge_uncertainty()

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
#latex_plot()
#plt.figure(figsize=(8.3,4))
#plt.title('B27D00010 - Linear Recharge')
#R = ml.bores_list[0].recharge[ml.bores_list[0].recharge.index > '1997-01-01']
#plt.bar(R.index, R.values, color = 'k', lw=1)

#plt.fill_between(R.index, 0.0, 0.1, where=(R.index.month > 10), facecolor=cyan, alpha=0.5, label='winter', lw=0)
#plt.fill_between(R.index, 0.0, 0.1, where=(R.index.month <4), facecolor=cyan, alpha=0.5, lw=0)
#plt.ylim(0,0.035)
#plt.ylabel('Recharge [m/d]')
#plt.xlabel('Time [years]')
#p = plt.Rectangle((0, 0), 1, 1, fc=cyan)
#plt.legend([p], ['November-March'])
#plt.savefig('Figures/recharge_detail.eps', bbox_inches='tight') 
#
#ticks = ax.xaxis.get_ticklocs()
#ax.xaxis.set_ticks(ticks[::365])
#ax.xaxis.set_ticklabels(ticklabels, rotation=0)
#plt.savefig('Figures/recharge.eps', bbox_inches='tight') 

R = ml.bores_list[0].recharge.resample('A', how='sum')
plt.bar(R.index, R.values, lw=1, color=[120/255.,196./255,1.], width = -250)