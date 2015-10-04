# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 11:22:33 2015

@author: Raoul
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
from pandas import read_csv



# Import the Climate data
data = read_csv('Test_Data/KNMI_Bilt.txt', skipinitialspace=True, skiprows=11, delimiter=',', parse_dates=['YYYYMMDD'], index_col=['YYYYMMDD'])

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
P_avg[P_avg.index>1960].plot(kind='bar', color='lightskyblue', label='Precipitation')
plt.ylabel('annual precipitation [mm]')
plt.legend()

plt.subplot(212)
E_avg[E_avg.index>1960].plot(kind='bar', color='lightcoral', label='Evapotranspiration')
plt.ylabel('annual evapotranspiration [mm]')
plt.xlabel('Time [years]')
plt.legend()
