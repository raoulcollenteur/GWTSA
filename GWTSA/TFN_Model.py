# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 14:34:13 2014
@author: rcollenteur
"""

#Import all the packages needed for the simulation
import numpy as np
from scipy.special import gammainc
from Unsat_Zone import percolation

'''
TFN1 defines what is in this research the most simple transfer function noise (TFN) model. The model is linear and the recharge is calculated by extracting the evaporation from the precipitation
'''

def TFN1(Parameters,InputData, solver = 0):
    
    #Unpack all the parameters that should be calibrated
    A = Parameters[0]
    a = Parameters[1]
    n = Parameters[2]
    d = Parameters[3]
    alpha = Parameters[4]
    
    #unpack all the data that is needed for the simulation
    t = InputData[1]
    P = InputData[2]
    E = InputData[3]
    Ho = InputData[4]
    to = InputData[5] 
    dt = InputData[6]
    tstart = InputData[7]
    istart = np.where(to > tstart)[0][0]
    
    #Recharge model
    R = P - E  
    
    # Set the value for the timestep to calculate the innovations
    Fs = A * gammainc(n,t/a) # Step response function based on pearsonIII
    Fb = Fs[1:] - Fs[0:-1] #block reponse function
    Hm = d + np.convolve(R,Fb)
    r = Ho - Hm[to] #Calculate the residuals at the timesteps with an observation
    v = r[1:] - ( r[0:-1] * np.exp(-dt/alpha) )
    v = v[istart:len(v):3] #give back the innovations for every xth time step
    return [Hm, v]

'''
TFN2 defines a TFN model that deals with the recharge a little more, adding a parameter 'f' that determines what part of the evaporation is extracted from the recharge. 
'''

def TFN2(Parameters,InputData):
    
    #Unpack all the parameters that should be calibrated
    A = Parameters[0]
    a = Parameters[1]
    n = Parameters[2]
    d = Parameters[3]
    alpha = Parameters[4]
    f = Parameters[5]
    
    #unpack all the data that is needed for the simulation
    t = InputData[1]
    P = InputData[2]
    E = InputData[3]
    Ho = InputData[4]
    to = InputData[5] 
    dt = InputData[6]
    tstart = InputData[7]
    istart = np.where(to > tstart)[0][0]
    
    #Recharge model
    R = P - f * E     
    
    # Set the value for the timestep to calculate the innovations
    Fs = A * gammainc(n,t/a) # Step response function based on pearsonIII
    Fb = Fs[1:] - Fs[0:-1] #block reponse function
    Hm = d + np.convolve(R,Fb)
    r = Ho - Hm[to] #Calculate the residuals at the timesteps with an observation
    v = r[1:] - ( r[0:-1] * np.exp(-dt/alpha) )
    v = v[istart:len(v):1] #give back the innovations for every xth time step
    return [Hm, v]

'''
TFN3 defines a TFN model that deals with the recharge a little more, adding a parameter 'f' that determines what part of the evaporation is extracted from the recharge. 
'''
    
    
def TFN3(Parameters,InputData):
    
    #Unpack all the parameters that should be calibrated
    A = Parameters[0]
    a = Parameters[1]
    n = Parameters[2]
    d = Parameters[3]
    alpha = Parameters[4]
    f = Parameters[5]
    B = Parameters[6]
    b = Parameters[7]
    m = Parameters[8]
    
    #unpack all the data that is needed for the simulation
    t = InputData[1]
    P = InputData[2]
    E = InputData[3]
    Ho = InputData[4]
    to = InputData[5] 
    dt = InputData[6]
    tstart = InputData[7]
    istart = np.where(to > tstart)[0][0]
    
    #Precipitation response
    Fs = A * gammainc(n,t/a) # Step response function based on pearsonIII
    Fb = Fs[1:] - Fs[0:-1] #block reponse function
    
    #Evaporation response
    Fs = B * gammainc(m,t/b) # Step response function based on pearsonIII
    Fe = Fs[1:] - Fs[0:-1] #block reponse function
    
    Hm = d + np.convolve(P,Fb) - f * np.convolve(E,Fe)
    r = Ho - Hm[to] #Calculate the residuals at the timesteps with an observation
    v = r[1:] - r[0:-1] * np.exp(-dt/alpha)
    v = v[istart:len(v):1] #give back the innovations for every xth time step
    return [Hm, v]    
    
def TFN4(Parameters,InputData, solver = 0):
    # Unpack all the parameters that should be calibrated
    A = 10**Parameters[0]
    a = Parameters[1]
    n = Parameters[2]
    d = Parameters[3]
    Alpha = 10**Parameters[4]
    S_cap = Parameters[5]
    K_sat = Parameters[6]
    Beta = Parameters [7]
    Imax = Parameters[8]
    
    # unpack all the data that is needed for the simulation
    Time_Model = InputData[1]
    P = InputData[2]
    E = InputData[3]
    Ho = InputData[4]
    to = InputData[5] 
    Time_Step = InputData[6]
    tstart = InputData[7]
    istart = np.where(to > tstart)[0][0]
    dt= 1 
  
    #Recharge model
    R = percolation(Time_Model, P, E, S_cap, K_sat, Beta, Imax , dt, solver = 0)[0]
    
    # Set the value for the timestep to calculate the innovations
    Fs = A * gammainc(n,Time_Model/a) # Step response function based on pearsonIII
    Fb = Fs[1:] - Fs[0:-1] #block reponse function
    Hm = d + np.convolve(R,Fb)
    r = Ho - Hm[to] #Calculate the residuals at the timesteps with an observation
    v = r[1:] - ( r[0:-1] * np.exp(-Time_Step/Alpha) )
    v = v[istart:len(v):3] #give back the innovations for every xth time step
    return [Hm, v]


    