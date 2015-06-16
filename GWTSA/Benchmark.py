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
from TFN_Model import TFN4
from datetime import date


print 'Import of packages succesfull!'
#plt.close('all')

Bore = 'Test_Data/B27B0081-001' # For Time Data to make realistic test series
forcing = 'Test_Data/KNMI_Bilt'

# (A, a, n, d, Alpha, S_cap, K_sat, Beta, D)
Parameters = [400.0, 10.0, 1.35, 7.0, 10.0, 0.0, -3.0, 2.0, -3.0]
TFN = 'TFN4' #TFN-model to use

Bore1 = Model(Bore, forcing)

# Unpack all the parameters that should be calibrated
A = Parameters[0]
a = Parameters[1]
n = Parameters[2]
d = Parameters[3]
Alpha = Parameters[4]
S_cap = Parameters[5]
K_sat = Parameters[6]
Beta = Parameters [7]
D = Parameters[8]
dt= 1 

#Recharge model 
R = percolation(Bore1.Time_Model, Bore1.P, Bore1.E, S_cap, K_sat, Beta, D , dt, solver = 1)[0]

#plt.plot((Bore1.P-Bore1.E))
#plt.plot(R, 'r')
#print sum(R[2000:]), sum(Bore1.P[2000:]) #/ 10000.0 - Bore1.E[2000:] / 10000.0)

# Set the value for the timestep to calculate the innovations
Fs = A * gammainc(n,Bore1.Time_Model/a) # Step response function based on pearsonIII
Fb = Fs[1:] - Fs[0:-1] #block reponse function
H = d + np.convolve(R,Fb)

np.random.seed(22)
r = 0.005 * np.random.standard_normal(len(Bore1.Time_Model)) # Residuals
NM = np.exp(-Bore1.Time_Model/Alpha) #Noise Model
e = np.convolve(r,NM)
Ho = H[0:len(Bore1.Time_Model)] + e[0:len(Bore1.Time_Model)]

plt.plot(Ho)
plt.plot(H[0:len(Bore1.Time_Model)])
plt.plot(e[0:len(Bore1.Time_Model)])

Ho = Ho[Bore1.Time_Observed]

# Make a datetime 
Time = md.num2date(Bore1.Time_Axis)

# Create a list with string of YY-MM-DD format
for i in range(len(Time)):
    Time[i] = date.strftime(Time[i], format = '%Y%m%d')
    
X = zip(Time,Ho)   

np.savetxt('testserie.txt', X, header='[400.0, 10.0, 1.35, 7.0, 50.0, 0.0, -3.0, 2.0, -3.0]', fmt = '%s', delimiter = ',', newline = '\r\n')
    
Hm = Bore1.test(Parameters, 'TFN4')
plt.plot(Hm,'r')

#R1 = percolation(Bore1.Time_Model, Bore1.P, Bore1.E, S_cap, K_sat, Beta, D , dt)
#R2 = percolation(Bore1.Time_Model, Bore1.P, Bore1.E, S_cap, K_sat, Beta, D , dt)
#R3 = percolation(Bore1.Time_Model, Bore1.P, Bore1.E, S_cap, K_sat, Beta, D , dt)
#plt.figure()
#plt.plot(R)
#plt.plot(R1)
#plt.plot(R2)
#plt.plot(R3)
#plt.legend(['R', 'R1', 'R2'])

'''
#For plotting the time series imported and created above, comment out.
#fig = figure(figsize=(40,10));
#subplot(311)
#bar(t,P)
#subplot(312)
#bar(t,E,color='red', lw=0)
#subplot(313)
#bar(t,R,color='green', lw=0)

#Create a benchmark groundwater time series
A=1.0; a=400.0; n=1.4; alpha=50; d=1.0; f=0.7
R = P - f * E
Fs = A * a * gammainc(n,t/a) # Step response function based on pearsonIII
Fb = Fs[1:] - Fs[0:-1] # block reponse function      
H = np.convolve(R,Fb)
H = H[0:tmax]
np.random.seed(22)
r = 0.05 * np.random.standard_normal(tmax+1) # Residuals
NM = np.exp(-t/alpha) #Noise Model
e = np.convolve(r,NM)
Ho = d + H + e[0:tmax]
H += d
np.savetxt('testserie.txt',[ClimateData[:,0],Ho],header='Head, # A=10; a=0.5; n=2; alpha=2; d=1')

print np.std(e)
print np.std(Ho)
print np.std(e)/np.std(Ho)

Here we try to find the values used above to see if our benchmark problem works. We use the formulas used in other applications, see the import of those function at the top of this script. 

#Initial_Parameters = [2,250,1,np.mean(Ho),50,0.5] #par = [A,a,n,d,alpha] initial parameters
#InputData = [t,ClimateData,Ho,t,dt,0] 
##data = [timesteps, Recharge, tmax, Observed GWL, Timesteps]
##Parameters = fmin(SWSI,Initial_Parameters,args=(InputData,))
#Parameters = cma.fmin(SWSI2,Initial_Parameters,2,args=(InputData,))[0]
#InputData = [t,ClimateData,Ho,t,dt,0] #use parameters for entire dataset
#[Hm, v, R] = TFN2(Parameters,InputData); #model the GWL

fig = plt.figure(figsize=(20,8),fontsize=15)
ax1 = plt.gca()
plt.plot(Ho[0:tmax]); 
plt.plot(H[0:tmax],'b-');
plt.legend(['Modeled Head','Observed Head'],loc=4)
plt.xlabel('Time [T]')
plt.ylabel('Groundwater Head [L]')


fig = plt.figure(figsize=(20,8))
ax1 = plt.gca()
plt.plot(r[0:500], 'r'); 
plt.plot(e[0:500],'b-');
plt.legend(['Innovations','Residuals'],loc=4)
plt.xlabel('Time [T]', size=20)
plt.ylabel('Groundwater Head Error [L]', size=20)
plt.title('Residuals and Innovations', size=20)
fig.set_fontsize(50)

Print the output parameters and compare to the real values
parname = ['A','a','n','d','alpha','f','B','b','m']
print 'the sum of least squares of the innovations is', SWSI2(Parameters,InputData)
#print 'the sum of least squares of the innovations is', LI(par,data)

for i in range(0,len(Parameters)):
    print parname[i], '=', Parameters[i]


plt.figure()
Fs = Parameters[0] * Parameters[1] * gammainc(Parameters[2],t/Parameters[1])
Fb = Fs[1:] - Fs[0:-1]
plt.plot(t[0:-1],Fb)

'''


