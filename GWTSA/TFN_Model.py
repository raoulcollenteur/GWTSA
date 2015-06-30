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

def TFN1(parameters,InputData, solver = 0):
    
    #Unpack all the parameters that should be calibrated
    A = parameters[0]
    a = parameters[1]
    n = parameters[2]
    d = parameters[3]
    alpha = parameters[4]
    
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
    
def TFN2(parameters,InputData, solver=0):
    
    #Unpack all the parameters that should be calibrated
    A = parameters[0]
    a = parameters[1]
    n = parameters[2]
    d = parameters[3]
    alpha = 10**parameters[4]
    f = 10**parameters[5]
    
    #unpack all the data that is needed for the simulation
    time_model = InputData[1]
    precipitation = InputData[2]
    evaporation = InputData[3]
    head_observed = InputData[4]
    time_observed = InputData[5] 
    time_step = InputData[6]
    
    #Recharge model
    recharge = precipitation - f * evaporation     
    
        # Set the value for the timestep to calculate the innovations
    Fs = A * gammainc(n, time_model/a) # Step response function based on pearsonIII
    Fb = Fs[1:] - Fs[0:-1] #block reponse function
    head_modeled = d + np.convolve(recharge,Fb)
    residuals = head_observed - head_modeled[time_observed]
    innovations = residuals[1:] - ( residuals[0:-1] * np.exp(-time_step/alpha) )
    
    return [head_modeled, innovations, recharge, residuals]

'''
TFN3 defines a TFN model that deals with the recharge a little more, adding a parameter 'f' that determines what part of the evaporation is extracted from the recharge. 
'''
    
    
def TFN3(parameters,InputData):
    
    #Unpack all the parameters that should be calibrated
    A = parameters[0]
    a = parameters[1]
    n = parameters[2]
    d = parameters[3]
    alpha = parameters[4]
    f = parameters[5]
    B = parameters[6]
    b = parameters[7]
    m = parameters[8]
    
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
    
''' TFN4 has the non-linear unsaturated zone model to calculate the recharge  
'''    
    
def TFN4(parameters,InputData, solver = 0):
    # Unpack all the parameters that should be calibrated
    A = 10**parameters[0]
    a = 10**parameters[1]
    n = parameters[2]
    d = parameters[3]
    alpha = 10**parameters[4]
    S_cap = parameters[5]
    K_sat = parameters[6]
    Beta = parameters [7]
    Imax = -3 #parameters[8]

    # unpack all the data that is needed for the simulation
    time_model = InputData[1]
    P = InputData[2]
    E = InputData[3]
    head_observed = InputData[4]
    time_observed = InputData[5] 
    time_step = InputData[6]
    dt= 1 
  
    #Recharge model
    recharge = percolation(time_model, P, E, S_cap, K_sat, Beta, Imax , dt, solver = 0)[0]
    
    Fs = A * gammainc(n, time_model/a) # Step response function based on pearsonIII
    Fb = Fs[1:] - Fs[0:-1] #block reponse function
    head_modeled = d + np.convolve(recharge,Fb)
    residuals = head_observed - head_modeled[time_observed]
    innovations = residuals[1:] - ( residuals[0:-1] * np.exp(-time_step/alpha) )
    
    return [head_modeled, innovations, recharge, residuals]
    