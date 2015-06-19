# -*- coding: utf-8 -*-
"""
@author: R.A. Collenteur
"""

#Import all the packages needed
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import gammainc
from GWTSA import *
from Unsat_Zone import percolation
from datetime import date

plt.close('all')

Bore = 'Test_Data/GW_Levels' # For Time Data to make realistic test series
forcing = 'Test_Data/KNMI_Bilt'

# (A, a, n, d, Alpha, S_cap, K_sat, Beta, D)
Parameters = [2.6, 10.0, 1.35, 7.0, 1.0, 0.0, -3.0, 2.0, -3.0]
TFN = 'TFN4' #TFN-model to use

Bore1 = Model(Bore, forcing, rows=[1,8])

# Unpack all the parameters that should be calibrated
A = 10**Parameters[0]
a = Parameters[1]
n = Parameters[2]
d = Parameters[3]
alpha = 10**Parameters[4]
S_cap = Parameters[5]
K_sat = Parameters[6]
Beta = Parameters [7]
Imax = Parameters[8]
dt= 1 

#Recharge model 
recharge = percolation(Bore1._time_model, Bore1.precipitation, Bore1.evaporation, S_cap, K_sat, Beta, Imax, dt, solver=0)[0]

# Set the value for the timestep to calculate the innovations
Fs = A * gammainc(n,Bore1._time_model/a) # Step response function based on pearsonIII
Fb = Fs[1:] - Fs[0:-1] #block reponse function
H = d + np.convolve(recharge,Fb)

np.random.seed(22)
r = 0.00 * np.random.standard_normal(len(Bore1._time_model)) # Residuals
NM = np.exp(-Bore1._time_model/alpha) #Noise Model
e = np.convolve(r,NM)
Ho = H[0:len(Bore1._time_model)] + e[0:len(Bore1._time_model)]

plt.plot(Ho)
plt.plot(H[0:len(Bore1._time_model)])
plt.plot(e[0:len(Bore1._time_model)])

Ho = Ho[Bore1._time_observed]

# Make a datetime 
Time = md.num2date(Bore1._time_axis)

# Create a list with string of YY-MM-DD format
for i in range(len(Time)):
    Time[i] = date.strftime(Time[i], format = '%Y%m%d')
    
X = zip(Time,Ho)   

np.savetxt('Test_Data/GW_Levels.txt', X, header='[2.6, 10.0, 1.35, 7.0, 1.0, 0.0, -3.0, 2.0, -3.0]', fmt = '%s', delimiter = ',', newline = '\r\n')

Hm = Bore1.simulate('TFN4', Parameters)
Bore1.plot_heads(newfig = 0)