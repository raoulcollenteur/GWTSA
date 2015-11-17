# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 09:25:59 2015

@author: Raoul

In this file all the different plots are made for the report
"""

cyan = [120/255.,196./255,1.]

from GWTSA import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
from pandas import read_csv

plt.close('all')
latex_plot()

#%% Plot the impulse response function
X0 = Parameters()
X0.add_many(('A',   3.0,    True,   None, None,  None),
           ('a',    2.48,    True,   None, None,  None),
           ('b',    0,    True,   None, None,  None),
           ('n',    1.5,    True,   None, None,  None))
x1 = TFN_Model.IRF(X0)

X0 = Parameters()
X0.add_many(('A',   3.0,    True,   None, None,  None),
           ('a',    2.48,    True,   None, None,  None),
           ('b',    0.0,    True,   None, None,  None),
           ('n',    3.0,    True,   None, None,  None))
x2 = TFN_Model.IRF(X0)

X0 = Parameters()
X0.add_many(('A',   2.7,    True,   None, None,  None),
           ('a',    2.0,    True,   None, None,  None),
           ('b',    0.0,    True,   None, None,  None),
           ('n',    1.5,    True,   None, None,  None))
x3 = TFN_Model.IRF(X0)

X0 = Parameters()
X0.add_many(('A',   3.0,    True,   None, None, None),
           ('a',    2.48,    True,   None, None,  None),
           ('b',    2.5,    True,   None, None,  None),
           ('n',    1.0,    True,   None, None,  None))

x4 = TFN_Model.IRF(X0)


plt.figure(figsize=(4.15,3))    
plt.plot(x1, '-', color=cyan, label='A=1000, a=300, n=1.5')
plt.plot(x2, '--', color=cyan, label='A=1000, a=300, n=3.0')
plt.plot(x3, 'k-', label='A=500, a=100, n=1.5')
plt.plot(x4, 'k--', label='A=1000, a=300, n=1.0')

plt.xlim(0,1500)
plt.legend()
plt.xlabel('Time [Days]')
plt.ylabel('Response [-]')
plt.savefig('Figures/impulse_response.eps', bbox_inches='tight')

#%% Plot the recharge coefficients

S = np.arange(0,1.,0.01)
Beta = [0.5, 1.0, 3.0, 10.0, 15.0]
Cr = np.ones((len(Beta),100))
Scap=1.
for i in range(len(Beta)):
    Cr[i] = (S/Scap)**Beta[i]
plt.figure(figsize=(4.15,3))    
plt.plot(Cr[0],S, '-', color='k', label=r'$\beta$=0.5')
plt.plot(Cr[1],S, '-', color=cyan, label=r'$\beta$=1.0')
plt.plot(Cr[2],S, '--', color='k', label=r'$\beta$=3.0')
plt.plot(Cr[3],S, '--', color=cyan, label=r'$\beta$=10.0')
plt.plot(Cr[4],S, '-', color='k', label=r'$\beta$=15.0')
#plt.plot((1-Cr).T, S)
plt.ylabel(r'$S_r/S_{rmax}$ [-]'), plt.xlabel(r'$C_r$ [-]')
plt.yticks([0.2,0.4,0.6,0.8,1.0],[0.2,0.4,0.6,0.8,1.0], rotation=0)
plt.legend(loc=0)


plt.savefig('Figures/beta_parameter.eps', bbox_inches='tight')


S = np.arange(0,1.,0.01)
Beta = np.arange(0, 1.5, 0.05)
Cr = np.ones((len(Beta),100))
for i in range(len(Beta)):
    Cr[i] = (1-(1-(S/Scap))**Beta[i])
plt.figure()
plt.plot(Cr.T, S)
#plt.plot(1-Cr, S)
plt.ylabel('St'), plt.xlabel('Cr,1-Cr')
plt.legend(['Recharge', 'UZ'])


#r[i] = 1. / (1+np.exp((-S/Scap + 0.5)/Beta[i])) # The soil retention curve

#%% Plot the precipitation and evapotranspiration

data = read_csv('Test_Data/KNMI_Deelen.txt', skipinitialspace=True, skiprows=11, delimiter=',', parse_dates=['YYYYMMDD'], index_col=['YYYYMMDD'])

data.RH = data.RH/10.
data.EV24 = data.EV24/10.

P_avg = data.RH.resample('A', how='sum', kind='YYYYMMDD')
E_avg = data.EV24.resample('A', how='sum', kind='YYYYMMDD')
R_avg = P_avg - E_avg

P_mean = P_avg.mean()
E_mean = E_avg.mean()
R_mean = R_avg.mean()

plt.figure(figsize=(8.3,3))

plt.subplot(211)
P_avg[P_avg.index>1986].plot(kind='bar', color=cyan, label='Precipitation')
plt.ylabel(r'$P$ [mm/year]')
plt.legend(loc=4)

plt.subplot(212)
E_avg[E_avg.index>1986].plot(kind='bar', color='lightcoral', label='Evapotranspiration')
plt.ylabel(r'$E_p$ [mm/year]')
plt.xlabel('Time [years]')
plt.legend(loc=4)
plt.savefig('Figures/climate_deelen.eps', bbox_inches='tight')

#%% Plot the boreholes used in this study

