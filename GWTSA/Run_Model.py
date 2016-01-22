# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 11:39:58 2014
this is a time series analysis model
@author: R.A. Collenteur
"""

#Import all the packages needed
from GWTSA import *

#import glob
#bores = glob.glob('Test_Data/All/*.csv')

plt.close('all')

#bores = ['Test_Data/B33A0113001_1.csv','Test_Data/B33A0113001_1.csv','Test_Data/B33A0113001_1.csv','Test_Data/B33A0113001_1.csv']
#bores = [ 'Test_Data/B27D0001001_1.csv', 'Test_Data/B27D0001001_1.csv',  'Test_Data/B27D0001001_1.csv', 'Test_Data/B27D0001001_1.csv']
#bores = ['Test_Data/B27C0049001_1.csv','Test_Data/B27C0049001_1.csv','Test_Data/B27C0049001_1.csv','Test_Data/B27C0049001_1.csv']
bores = [ 'Test_Data/B27D0001001_1.csv', 'Test_Data/B27C0049001_1.csv', 'Test_Data/B33A0113001_1.csv']
#bores = ['Test_Data/UH/B32C0058001_1.csv', 'Test_Data/UH/B32C0058001_1.csv'] 
bores = [ 'Test_Data/B27D0001001_1.csv']

forcing = 'Test_Data/KNMI_Bilt.txt'
calibration = ['01-01-1974', '31-12-1994']
#calibration = ['01-01-1974', '31-12-2003']
validation = ['01-01-1995', '31-12-2004']
#validation = ['01-01-2004', '31-12-2004']

ml = Model(bores, forcing, calibration, validation)

#trend = ['linear']*len(bores)
#RM = ['preferential']*len(bores)
#TFN = ['percolation']*len(bores)
#TFN = ['combination']*len(bores)

#print ml.bores_list[0].calcSumax() #= 0.35 for Deelen

X0 = Parameters()
#           (Name,  Value,  Vary,   Min,  Max,  Expr)
X0.add_many(('A',  3,    True,   1e-6, None,  None),
            ('a',  200.,    True,  0.0, None,  None),
            ('n',    1.0,    True,   None, None,  None),
            ('A1',   3,    True,   1e-6,  None,  None),
            ('a1',  200,    True,  0.0, None,  None),
            ('n1',    1.0,    True,   None, None,  None),
            ('alpha',100.0,    True,   0.0,  1000.0,   None),
            # Recharge
            #('f',    1.0,    True,   0.0, 1.5,  None),
            #('mu',  150.0, True,  0.0,  None,  None),
            #('sig', 25.0,  True,  0.0,  None,  None),
            ('Srmax',0.305,  False,  None, None,  None),
            ('Imax', 1.5e-3, False,  None, None,  None),
            ('Kp',  0.005,    True,   1e-4, 0.1,  None),
            ('Beta', 2.0,    True,   0.0, None,  None),
            ('gamma',3.0,   True,   0.0, None,  None),
            # Reclamation
            #('b',   3.0,    True,   None, None,  None),
            #('B',   -0.5,    True,   -0.5, 0.0,  None)
            # Linear Slope 
            #('slope', 0.00015, False, 0.00001, 0.00025,  None),
            #('intercept', -1.10, False, None, None,  None)
            )

ml.solve(X0, RM='preferential', method='leastsq')

latex_plot()
cyan = [120/255.,196./255,1.]

ml.plot(modeled=1, savefig=1)

ml.bores_list[0].plot_results(savefig=1)
ml.bores_list[0].plot_diagnostics(savefig=1)


ml.bores_list[1].plot_results(savefig=1)
ml.bores_list[1].plot_diagnostics(savefig=1)


ml.bores_list[2].plot_results(savefig=1)
ml.bores_list[2].plot_diagnostics(savefig=1)


print np.mean(ml.bores_list[0].recharge.resample('A', how='sum'))
print np.mean(ml.bores_list[1].recharge.resample('A', how='sum'))
print np.mean(ml.bores_list[2].recharge.resample('A', how='sum'))
#print np.mean(ml.bores_list[3].recharge.resample('A', how='sum'))

#plt.plot(ml.bores_list[0].head_modeled[ml.bores_list[0]._index_observed],ml.bores_list[0].head_observed, 'bo')
#plt.plot([18,21.5],[18,21.5])

#%% Plot the recharges uncertainty
#
plt.figure(figsize=(8.3,4))

ax1 = plt.subplot(311)
plt.title('B27D00010')
ml.bores_list[0].recharge_uncertainty()
ax1.xaxis.set_ticklabels([])
ax2 = plt.subplot(312, sharey=ax1, sharex=ax1)
plt.title('B27C00490')
ml.bores_list[1].recharge_uncertainty()
ax2.xaxis.set_ticklabels([])     
plt.ylabel('recharge [m/y]')
plt.subplot(313, sharey=ax1)
plt.title('B33A01130')
ml.bores_list[2].recharge_uncertainty()

plt.xlabel('Time [years]')
plt.savefig('Figures/recharge_boxplot.eps', bbox_inches='tight') 

#%% Plot individual daily recharges
latex_plot()
plt.figure(figsize=(8.3,2.0))
plt.title('B27D00010 - Combination')
R = ml.bores_list[0].recharge[ml.bores_list[0].recharge.index > '1997-01-01']
plt.bar(R.index, R.values, color = 'k', lw=1)

plt.fill_between(R.index, 0.0, 0.1, where=(R.index.month > 10), facecolor=cyan, alpha=0.5, label='winter', lw=0)
plt.fill_between(R.index, 0.0, 0.1, where=(R.index.month <4), facecolor=cyan, alpha=0.5, lw=0)
plt.ylim(-0.005,0.05)
plt.ylabel('Recharge [m/d]')
plt.xlabel('Time [years]')
p = plt.Rectangle((0, 0), 1, 1, fc=cyan)
plt.legend([p], ['November-March'])
plt.savefig('Figures/recharge_detail.eps', bbox_inches='tight') 

#%% Plot the residual series
latex_plot()
plt.figure(figsize=(8.3,2.0))
plt.plot(md.num2date(ml.bores_list[0]._time_axis), ml.bores_list[0].residuals, 'k', label='B27D00010')
plt.plot(md.num2date(ml.bores_list[1]._time_axis), ml.bores_list[1].residuals, color=cyan, label='B27C00490')
plt.plot(md.num2date(ml.bores_list[2]._time_axis), ml.bores_list[2].residuals, 'gray', label='B33A01130')
plt.xlabel('Time [years]')
plt.ylabel('Error [m]')
plt.legend(loc=(0,1), ncol=3, frameon=False, handlelength=3)
plt.axhline(0.0, linestyle='--', color='dimgray')
plt.savefig('Figures/residuals.eps', bbox_inches='tight') 

#%%
R = ml.bores_list[0].df
R = R.T
R.hist()


#%%
#X0 = Parameters()
##           (Name,  Value,  Vary,   Min,  Max,  Expr)
#X0.add_many(('A',  -1.0,    True,   None, 0.0,  None),
#            ('a',  2.0,    True,  0.0, None,  None),
#            ('n',    2.0,    True,   None, None,  None),
#            #('A1',   0.4,    True,   0.0,  None,  None),
#            #('t_p1',  2.4,    True,  0.0, 3.0,  None),
#            #('n1',    1.5,    True,   None, None,  None),
#            ('alpha',1.0,    True,   0.01,  3.0,   None),
#            # Recharge
#            #('f',    1.0,    True,   0.0, 1.5,  None),
#            #('mu',  150.0, True,  0.0,  None,  None),
#            #('sig', 25.0,  True,  0.0,  None,  None),
#            ('Srmax',0.14,  False,  None, None,  None),
#            ('Imax', 1.5e-3, False,  None, None,  None),
#            #('Kp',  -2.25,    True,   -4.0, -1.0,  None),
#            ('Beta', 1.0,    True,   0.0, None,  None),
#            #('gamma',3.0,   True,   0.0, None,  None),
#            # Reclamation
#            #('b',   3.0,    True,   None, None,  None),
#            #('B',   -0.5,    True,   -0.5, 0.0,  None)
#            # Linear Slope 
#            #('slope', 0.00015, False, 0.00001, 0.00025,  None),
#            #('intercept', -1.10, False, None, None,  None)
#            )

plt.figure(figsize=(8.3,2.0))
plt.plot(md.num2date(ml.bores_list[0]._time_axis), ml.bores_list[0].residuals, 'k', label='residuals')
plt.plot(md.num2date(ml.bores_list[0]._time_axis[0:-1]), ml.bores_list[0].innovations, color=cyan, label='innovations')
plt.legend(loc=(0,1), ncol=3, frameon=False, handlelength=3)
plt.ylabel('Error [m]')
plt.xlabel('Time [Years]')   
plt.savefig('Figures/Error_Model.eps', bbox_inches='tight') 