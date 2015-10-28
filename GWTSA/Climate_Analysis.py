# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 11:22:33 2015

@author: Raoul
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
from pandas import read_csv

plt.close('all')

# Import the Climate data
data = read_csv('Test_Data/KNMI_Deelen.txt', skipinitialspace=True, skiprows=11, delimiter=',', parse_dates=['YYYYMMDD'], index_col=['YYYYMMDD'])

data.RH = data.RH/10.
data.EV24 = data.EV24/10.

P_avg = data.RH.resample('A', how='sum', kind='YYYYMMDD')
E_avg = data.EV24.resample('A', how='sum', kind='YYYYMMDD')
R_avg = P_avg - E_avg

P_mean = P_avg.mean()
E_mean = E_avg.mean()
R_mean = R_avg.mean()

#plt.bar(md.num2date(data[:,0]),data[:,1])

P_mean - E_mean




plt.figure(figsize=(15,9))

plt.subplot(211)
P_avg[P_avg.index>1986].plot(kind='bar', color='lightskyblue', label='Precipitation')
plt.ylabel('annual precipitation [mm]')
plt.legend()

plt.subplot(212)
E_avg[E_avg.index>1986].plot(kind='bar', color='lightcoral', label='Evapotranspiration')
plt.ylabel('annual evapotranspiration [mm]')
plt.xlabel('Time [years]')
plt.legend()
plt.savefig('climate_deelen.eps', format='eps', bbox_inches='tight')

PW_avg = data.RH.resample('M', how='sum', kind='YYYYMMDD')
PW_avg[(PW_avg.index>1992) & (PW_avg.index<1996)].plot(kind='bar', color='lightskyblue', label='Precipitation')


# Mass Curve Technique to estimate S_cap
P = data.RH.resample('A', how='cumsum', kind='YYYYMMDD')
E = data.EV24.resample('A', how='cumsum', kind='YYYYMMDD')
R = P-E
R.plot()


#%% INFLUENCE 
S = np.arange(0,1.,0.01)
Beta = np.arange(0, 15, 0.5)
Cr = np.ones((len(Beta),100))
Scap=1.
for i in range(len(Beta)):
    Cr[i] = (S/Scap)**Beta[i]
plt.figure()    
plt.plot(Cr.T,S)
#plt.plot((1-Cr).T, S)
plt.ylabel('St'), plt.xlabel('Cr')
plt.legend(['Recharge', 'UZ'])


S = np.arange(0,1.,0.01)
Beta = np.arange(0, 1.0, 0.01)
Cr = np.ones((len(Beta),100))
for i in range(len(Beta)):
    Cr[i] = (1-(1-(S/Scap))**Beta[i])
plt.figure()
plt.plot(Cr.T, S)
#plt.plot(1-Cr, S)
plt.ylabel('St'), plt.xlabel('Cr,1-Cr')
plt.legend(['Recharge', 'UZ'])







