# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 11:39:58 2014
this is a time series analysis model
@author: R.A. Collenteur
"""

# Import all the packages needed
from GWTSA import *

#import glob
#bores = glob.glob('Test_Data/All/*.csv')

plt.close('all')

bores = ['Test_Data/B33A0113001_1.csv', 'Test_Data/B33A0113001_1.csv', 'Test_Data/B33A0113001_1.csv', 'Test_Data/B33A0113001_1.csv']
#bores = [ 'Test_Data/B27D0001001_1.csv', 'Test_Data/B27D0001001_1.csv',  'Test_Data/B27D0001001_1.csv',  'Test_Data/B27D0001001_1.csv']
#bores = ['Test_Data/B27C0049001_1.csv', 'Test_Data/B27C0049001_1.csv', 'Test_Data/B27C0049001_1.csv', 'Test_Data/B27C0049001_1.csv']
#bores = [ 'Test_Data/B27D0001001_1.csv', 'Test_Data/B27C0049001_1.csv', 'Test_Data/B33A0113001_1.csv']
#bores = ['Test_Data/B33A0113001_1.csv']

forcing = 'Test_Data/KNMI_Apeldoorn.txt'
calibration = ['01-01-1974', '31-12-1994']
#calibration = ['01-01-1974', '31-12-2003']
validation = ['01-01-1995', '31-12-2004']
#validation = ['01-01-2004', '31-12-2004']

ml = Model(bores, forcing, calibration, validation, discharge='Test_Data/Discharge_Tongeren.csv', Cf=[1.0,1.0])

RM = ['combination']*len(bores)
RM = ['linear', 'preferential', 'percolation', 'combination']

#print ml.bores_list[0].calcSumax() #= 0.305, 0.226, 0.262

X0 = Parameters()
#           (Name,  Value,  Vary,   Min,  Max,  Expr)
X0.add_many(('A',  3.,    True,   1e-6, None,  None),
            ('a',  200.,    True,  0.0, None,  None),
            ('n',    1.0,    False,   None, None,  None),
            #('A1',   3,    True,   1e-6,  None,  None),
            #('a1',  200,    True,  0.0, None,  None),
            #('n1',    1.0,    True,   None, None,  None),
            ('alpha', 100.0,    True,   0.0,  None,   None),
            # Recharge
            ('f',    1.0,    True,   0.0, 1.5,  None),
            ('mu',  15.0, True,  0.0,  None,  None),
            ('sig', 25.0,  True,  0.0,  None,  None),
            ('Srmax',0.262,  False,  None, None,  None),
            ('Imax', 1.5e-3, False,  None, None,  None),
            ('Kp',  0.005,    True,   1e-4, 0.1,  None),
            ('Beta', 3.0,    True,   0.0, None,  None),
            ('gamma',3.0,   True,   0.0, None,  None),
            # Reclamation
            ('b',   4.0,    True,   0.0, None,  None),
            ('B',   -0.50,    True,   -1.0, 0.0,  None),
            # Pumping well
            #('u',   6.55,    True,   None, None,  None),
            #('v',   -0.0193,    True,   -10.0, 0.0,  None)
            # Linear Slope 
            #('slope', -1.0e-4, False, -0.1, 0.00025,  None),
            #('intercept', 0.50, False, None, None,  None)
            )

ml.solve(X0, RM=RM, trend='reclamation', method='leastsq')

latex_plot()
cyan = [120/255., 196./255, 1.]

#ml.plot(modeled=1, savefig=1)

ml.bores_list[0].plot_results(savefig=1)
ml.bores_list[0].plot_diagnostics(savefig=1)

plt.figure(figsize=(8.3,2.0))
ax1 = plt.subplot(111)
plt.title('%s %s' %(ml.bores_list[0].bore, ml.bores_list[0]._RM))
ml.bores_list[0].recharge_uncertainty()
plt.ylabel('recharge [m/y]')
plt.ylim(-0.4,1.2)
plt.xlabel('Time [years]')
plt.savefig('Figures/recharge_uncertainty_%s.eps' % (ml.bores_list[0]._RM), bbox_inches='tight') 

ml.bores_list[1].plot_results(savefig=1)
ml.bores_list[1].plot_diagnostics(savefig=1)

ml.bores_list[3].plot_results(savefig=1)
ml.bores_list[2].plot_diagnostics(savefig=1)

print np.mean(ml.bores_list[0].recharge.resample('A', how='sum'))
print np.mean(ml.bores_list[1].recharge.resample('A', how='sum'))
print np.mean(ml.bores_list[2].recharge.resample('A', how='sum'))
#print np.mean(ml.bores_list[3].recharge.resample('A', how='sum'))

# %% Plot the recharges uncertainty
plt.figure(figsize=(8.3,2.0))
ax1 = plt.subplot(111)
plt.title('%s %s' %(ml.bores_list[0].bore, ml.bores_list[0]._RM))
ml.bores_list[0].recharge_uncertainty()
plt.ylabel('recharge [m/y]')
plt.ylim(-0.4,1.2)
plt.xlabel('Time [years]')
plt.savefig('Figures/recharge_uncertainty_%s.eps' % (ml.bores_list[0]._RM), bbox_inches='tight') 

# %% Plot individual daily recharges

for i in range(len(RM)):
    plt.figure(figsize=(8.3, 2.0))
    plt.title('%s - %s' %(ml.bores_list[i].bore, RM[i]))
    R = ml.bores_list[i].recharge[ml.bores_list[i].recharge.index > '1997-01-01']
    plt.bar(R.index, R.values, color='k', lw=1)
    plt.fill_between(R.index, 0.0, 0.1, where=(R.index.month > 10),
                     facecolor=cyan, alpha=0.5, label='winter', lw=0)
    plt.fill_between(R.index, 0.0, 0.1, where=(R.index.month <4),
                     facecolor=cyan, alpha=0.5, lw=0)
    plt.ylim(-0.005, 0.05)
    plt.ylabel('Recharge [m/d]')
    plt.xlabel('Time [years]')
    p = plt.Rectangle((0, 0), 1, 1, fc=cyan)
    plt.legend([p], ['November-March'])
    plt.savefig('Figures/recharge_detail_%s.eps' %RM[i], bbox_inches='tight')

#%% Plot the residual series
presentation_plot()
plt.figure(figsize=(8.3,2.0))
plt.plot(md.num2date(ml.bores_list[0]._time_axis), ml.bores_list[0].residuals, 'r', label='linear')
plt.plot(md.num2date(ml.bores_list[1]._time_axis), ml.bores_list[1].residuals, color=cyan, label='preferential')
plt.plot(md.num2date(ml.bores_list[2]._time_axis), ml.bores_list[2].residuals, 'gray', label='percolation')
plt.plot(md.num2date(ml.bores_list[3]._time_axis), ml.bores_list[3].residuals, 'k--', markersize=2, label='combination')
plt.xlabel('Time [years]')
plt.ylabel('Error [m]')
plt.legend(loc=(0,1), ncol=4, frameon=False, handlelength=3)
plt.axhline(0.0, linestyle='--', color='dimgray')
plt.savefig('Figures/residuals.eps', bbox_inches='tight') 

plt.figure(figsize=(8.3,2.0))
plt.plot(md.num2date(np.arange(ml.bores_list[0]._time_begin, ml.bores_list[0]._time_end+1)), ml.bores_list[0].trend, 'r', label='linear')
plt.plot(md.num2date(np.arange(ml.bores_list[1]._time_begin, ml.bores_list[1]._time_end+1)), ml.bores_list[1].trend, color=cyan, label='preferential')
plt.plot(md.num2date(np.arange(ml.bores_list[2]._time_begin, ml.bores_list[2]._time_end+1)), ml.bores_list[2].trend, 'gray', label='percolation')
plt.plot(md.num2date(np.arange(ml.bores_list[3]._time_begin, ml.bores_list[3]._time_end+1)), ml.bores_list[3].trend, 'k--', markersize=2, label='combination')
plt.xlabel('Time [years]')
plt.ylabel('Trend [m]')
plt.legend(loc=(0,1), ncol=4, frameon=False, handlelength=3)
plt.axhline(0.0, linestyle='-', color='dimgray')
#plt.axvline('01-01-1967', linestyle=':', color='dimgray')
plt.savefig('Figures/trend.eps', bbox_inches='tight') 
#%% Plot histograms of the recharge monte carlo runs per year.
#R = ml.bores_list[0].df
#R = R
##R.hist(bins=20)
#plt.figure()
#plt.bar(R.index, R.quantile(q=0.5, axis=1), yerr=[(R.quantile(q=0.5, axis=1)-R.quantile(q=0.025, axis=1)),(R.quantile(q=0.975, axis=1)-R.quantile(q=0.5, axis=1))], color=[120/255.,196./255,1.], width = -250, lw=0, error_kw={'ecolor': 'dimgray', 'ewidth': '5'})
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

#%%Plot all boreholes with manual control
plt.figure(figsize=(8.3,2.0))
plt.plot(ml.bores_list[0]._time_steps, 'k')
plt.axvline(466, color='k', linestyle='--')
plt.ylabel('Time step [days]')
plt.xlabel('Groundwater observations')
plt.savefig('Figures/Timesteps.eps', bbox_inches='tight') 

plt.figure()
plt.plot(md.num2date(np.arange(ml.bores_list[0]._time_begin, ml.bores_list[0]._time_end+1)), ml.bores_list[0].head_modeled, 'r', label='Linear')
plt.plot(md.num2date(np.arange(ml.bores_list[1]._time_begin, ml.bores_list[1]._time_end+1)), ml.bores_list[1].head_modeled, color=cyan, label='Preferential')
plt.plot(md.num2date(np.arange(ml.bores_list[2]._time_begin, ml.bores_list[2]._time_end+1)), ml.bores_list[2].head_modeled, 'gray', label='Percolation')
plt.plot(md.num2date(np.arange(ml.bores_list[3]._time_begin, ml.bores_list[3]._time_end+1)), ml.bores_list[3].head_modeled, 'k-', label='Combination')
plt.plot(md.num2date((ml.bores_list[0]._time_axis)), ml.bores_list[0].head_observed, 'k.', label='%s Observed' %ml.bores_list[0].bore,markersize=2)
plt.axvline(ml.bores_list[0]._date_calibration, color='k', linestyle='--')
plt.xlabel('Time [years]')
plt.ylabel('Groundwater level [m]')
plt.legend(loc=(0,1), ncol=5, frameon=False, handlelength=2, )
plt.savefig('Figures/%s_All.eps' %ml.bores_list[0].bore, bbox_inches='tight') 