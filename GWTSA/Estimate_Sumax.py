# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 19:04:50 2015

This function calculates the value for Sumax (The storage size of the root zone) based on the precipitation deficits. 

This function was originally written by T. Euser from the TU Delft (2015). The theory is further explained here:
Euser, T. et al. (2015) estimating 

@author: Raoul
"""

#Import all the packages needed for the calculation
import numpy as np
from pandas import Series, read_csv
import matplotlib.pyplot as plt
from scipy.stats import gumbel_r

def calcSumax(P,E, Cr=0.35, imax=0.001, T=20, EaMethod='gradual'):
    
    #%% Interception Reservoir
    Si = np.zeros(len(P))
    Pe = np.zeros(len(P))
    Ei = np.zeros(len(P))
    Ep = np.zeros(len(P))
    
    for t in range(len(P)-1):    
        Si[t+1] = Si[t]+P[t+1]                      # Fill interception bucket with new rain
        Pe[t+1] = np.max(((Si[t+1]-imax), 0.0))     # Calculate effective precipitation
        Si[t+1] = Si[t+1] - Pe[t+1]                 # Update interception state
        Ei[t+1] = np.min((Si[t+1], E[t+1]))         # Evaporation from interception
        Si[t+1] = Si[t+1] - Ei[t+1]                 # Update interception state
        Ep[t+1] = E[t+1] - Ei[t+1]                  # Update potential evapotranspiration   
    
    #%% Estimate Actual Evaporation
    
    # Averag actual evaporation [L/T]
    Ea = np.min((sum((1.-Cr)*Pe), sum(Ep)))/len(P)  # 1-Cr because Q = Cr * Pe
    #Get seasonality back in there
    EaAvgC = np.min((sum((1.-Cr)*Pe), sum(Ep))) / len(P)  # 1-Cr because Q = Cr * Pe
    
    #Get seasonality back in there
    if EaMethod == 'gradual':
        EaAvg = np.ones(len(Ep)) * EaAvgC           # get timeseries with constant average Ea
        a = np.min((Ep,EaAvg), axis=0)
        A = sum(EaAvg - a)
        B = sum(Ep - a)
        Ea = A / B * (Ep - a) + a                   # Calculate Ea with similar pattern as Ep

    if  EaMethod == 'flatline':
        EaAvg = np.ones(len(Ep)) * EaAvgC           # get timeseries with constant average Ea
        a = np.min((Ep,EaAvg))
        EaMiss = EaAvg - a
        EpMiss = Ep - a
        x = EaAvgC / 8.
        i = 1
        err = np.max(Ep)
        while abs(err) > (0.01 * x):
            x = x + err/(len(Ep)*0.8)
            addEa = np.min((x, EpMiss))
            A = sum(EaMiss)
            B = sum(addEa)
            err = A - B;
            
            i = i+1;
            print i
            if EpMiss == addEa:
                break
        
    soildeficit = np.zeros(len(P))   
    for t in range(len(P)-1):
        soildeficit[t+1] = np.min(((soildeficit[t] + Pe[t] - Ea[t]), 0.0))     
    
    soildeficit = Series(soildeficit, index=P.index)
    Sumax = np.sqrt(soildeficit.resample('A', how='min') ** 2. ) # Work with positive values

    Sumax.sort()
    mu, sigma = gumbel_r.fit(Sumax)
    y = gumbel_r.pdf(Sumax, mu, sigma)
    #plt.plot(Sumax,y)
    #Sumax.hist()
    #plt.axvline(gumbel_r.isf(1./T, loc=mu, scale=sigma)) #Causes trouble for multiple return periods

    return gumbel_r.isf(1./T, loc=mu, scale=sigma), Ea  

#%% Import data`
data = read_csv('Test_Data/KNMI_Bilt.txt', skipinitialspace=True, skiprows=11, delimiter=',', parse_dates=['YYYYMMDD'], index_col=['YYYYMMDD'])

P = data.RH[data.index>'1960-01-01 00:00:00']/10000.
E = data.EV24[data.index>'1960-01-01 00:00:00']/10000.
Cr = 0.40
imax = 0.001
T = np.array([20, 30])
EaMethod = 'gradual'

sumax, Ea = calcSumax(P,E, Cr, imax, T, EaMethod)
