# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 15:37:24 2016

@author: Raoul
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from scipy.signal import fftconvolve
from GWTSA import *

cyan = [120/255., 196./255, 1.]
presentation_plot()

t = np.arange(0,1000)
y = np.exp(-t/150.)
plt.figure(figsize=(3,2))
plt.plot(t,y, color=cyan, lw=2)
plt.xlabel('Time [days]')
plt.ylabel('Response [-]')
plt.yticks(np.arange(0.00, 1.1, 0.5))
plt.xticks(np.arange(0.00, 1001, 500))
plt.savefig('Figures/exponential_decay.eps', bbox_inches='tight') 

t = np.arange(0,300)
y = norm.pdf(t, 150,50)
plt.figure(figsize=(3,2))
ax = plt.plot(t,y, color=cyan, lw=2)
plt.xlabel('Time [days]')
plt.ylabel('Response [-]')
plt.yticks(np.arange(0.00, 0.009, 0.004))
plt.xticks(np.arange(0.00, 301, 100))

plt.savefig('Figures/gaussian.eps', bbox_inches='tight') 

#%% Convolution
A = 0.003
a = 200.0
n = 1.5
t = np.arange(1.,1000.)
Fs = A * t**n * (t/a)**-n * gammainc(n, t/a) # Step response function based on pearsonIII
Fb = np.append(0, Fs[1:] - Fs[0:-1]) #block reponse function 
presentation_plot()

plt.figure(figsize=(2,1.5))    
plt.plot(Fb, '-', color=cyan, lw=2)

z = [0,0,1,1,0,0]
x = [0,1,1,2,2,3]
plt.plot(x,z)

plt.xlim(0,1000)
plt.legend()
plt.xlabel('Time [T]')
plt.ylabel('Response [-]')
plt.xticks([])
plt.yticks([])
plt.savefig('Figures/impulse_response.eps', bbox_inches='tight')

z = [0,0,1,1,0,0,0,0]
x = [0,1,1,2,2,3,4,5]
plt.figure(figsize=(2,1.5))
plt.plot(x,z, 'k', lw=2)
plt.xlim(0,5)
plt.ylim(0.0,1.1)
plt.xlabel('Time [T]')
plt.ylabel('Recharge [L/T]')
plt.xticks([1,2],['0', '$t+\Delta t$'])
plt.yticks([])
plt.savefig('Figures/block.eps', bbox_inches='tight')

R = [0,0,1,0,0,2,0,0,1,1,1,0,0]
x = 0.001*np.ones((1000,13))
R = R*x
R = R.flatten(order='F')
y = fftconvolve(R,Fb)

mu = 120.12
sig= 79.23      
Fs4 = norm.cdf(range(500), mu, sig) #/norm.pdf(mu,mu,sig)
Fb4 = Fs4[1:] - Fs4[0:-1] #block reponse function
plt.figure(figsize=(10,5))
plt.plot(Fb4, '-', color='r', lw=2)
plt.plot(Fb4, '-', color=cyan, lw=2)
plt.legend(['Linear Model','Combination Model'])
plt.xlabel('Time [Days]')
plt.ylabel('Response [-]')
plt.yticks([0,0.005])
plt.xticks([0,500])
plt.savefig('Figures/percolation_response.eps', bbox_inches='tight')

plt.figure()
plt.subplot(311)
plt.plot(R, color='k')
plt.ylim(0.0001,0.0021)
plt.xticks([])
plt.yticks([])
plt.ylabel('Recharge [L/T]')

plt.subplot(312)
plt.plot(0.001+y, color=cyan)
plt.ylim(0.0,0.03)
plt.xticks([])
plt.yticks([])

plt.ylabel('Head [L]')

plt.subplot(313)
plt.plot(0.007+y, color='k')
plt.ylim(0.0,0.03)
plt.xticks([])
plt.yticks([])
plt.ylabel('Head [L]')
plt.xlabel('Time [T]')
plt.savefig('Figures/convolution.eps', bbox_inches='tight')