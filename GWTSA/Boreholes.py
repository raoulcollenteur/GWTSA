# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 15:26:29 2015

@author: Raoul
"""

from pandas import read_csv
import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.stats import linregress


plt.close('all')

bores = glob.glob('Test_Data/*.csv')

bore = [] #Store de different boreholes data
peak = [] #Store the 1995 peak moment


for i in range(len(bores)):
    if i == 0:
        bore.append(read_csv('%s' %bores[0], parse_dates=True, index_col=2, skiprows=17, skipinitialspace=True))
    else:
        bore.append(read_csv('%s' %bores[i], parse_dates=True, index_col=2, skiprows=15, skipinitialspace=True))
    bore[i].rename(columns={'Stand (cm t.o.v. MV)': 'h'}, inplace=True)
    bore[i].h = bore[i].h/100.
    bore[i].h.plot()
    peak.append(bore[i].h[bore[i].index > '1995-01-01 00:00:00'].argmax())
    bores[i] = bores[i][-17:-8]
    
Thickness = [10.35, 34.02, 29.24, 6.34, 49.09, 30.95, 33.06, 23.02, 47.10] #Estimated based on ground level minus highest groundwater level 

#plt.figure()
#plt.plot(peak, Thickness, 'o')

days = []
for i in range(len(peak)):
    days.append(peak[i] - min(peak))
    days[i] = days[i].days

plt.figure()
#plt.xkcd()
plt.plot(Thickness, days, 'k+', markersize=10, markeredgewidth=3)
plt.ylabel('Relative Delay Time [Days]')
plt.xlabel('Thickness Unsaturated Zone [Meters]')
plt.ylim(-5,300)
plt.xlim(0,55)
for i in range(len(bores)):
    plt.annotate(bores[i],(Thickness[i],days[i]))
plt.savefig('delay.eps', format='eps', bbox_inches='tight')

slope, intercept, r_value, p_value, std_er = linregress(days, Thickness)
y = intercept + slope * (np.linspace(0,300))
plt.plot(y,np.linspace(0,300), 'k--')

W = bore[0].resample('2W')