# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 14:34:13 2014
@author: rcollenteur
"""

#Import all the packages needed for the simulation
import numpy as np
from scipy.special import gammainc
from Unsat_Zone import percolation
from scipy.signal import fftconvolve


'''
TFN2 defines a TFN model that deals with the recharge a little more, adding a parameter 'f' that determines what part of the evaporation is extracted from the recharge. 
'''
    
def linear(parameters, InputData, solver=0):
    
    #Unpack all the parameters that should be calibrated
    A = 10.0** parameters['A'].value
    a = 10.0** parameters['a'].value
    n = parameters['n'].value
    d = parameters['d'].value
    f = 10.0** parameters['f'].value
    
    #unpack all the data that is needed for the simulation
    time_model = np.arange(0,10000)
    P = InputData[1]
    E = InputData[2]
    
    #Recharge model
    recharge = P - f * E     
    
        # Set the value for the timestep to calculate the innovations
    Fs = A * gammainc(n, time_model/a) # Step response function based on pearsonIII
    Fb = Fs[1:] - Fs[0:-1] #block reponse function
    head_modeled = d + fftconvolve(recharge,Fb)  
    return [head_modeled, recharge]

   
''' TFN4 has the non-linear unsaturated zone model to calculate the recharge  
'''    
    
def nonlinear(parameters,InputData, solver):
    # Unpack all the parameters that should be calibrated    
    A = 10.0**parameters['A'].value
    a = 10.0**parameters['a'].value
    n = parameters['n'].value
    d = parameters['d'].value
    S_cap = 10.0**parameters['scap'].value
    K_sat = 10.0**parameters['ksat'].value
    Beta = parameters['beta'].value       
    Imax = 10.0** parameters['imax'].value

    # unpack all the data that is needed for the simulation
    time_model = InputData[0]
    P = InputData[1]
    E = InputData[2]
    dt= 1 
    
    #Recharge model
    recharge = percolation(time_model, P, E, S_cap, K_sat, Beta, Imax , dt, solver)[0]
    time_model = np.arange(0,10000)
    Fs = A * gammainc(n, time_model/a) # Step response function based on pearsonIII
    Fb = Fs[1:] - Fs[0:-1] #block reponse function
    head_modeled = d + fftconvolve(recharge,Fb)    
    return [head_modeled, recharge]
    

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
    