# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 14:34:13 2014
@author: rcollenteur
"""

#Import all the packages needed for the simulation
import numpy as np
from scipy.special import gammainc, gamma
from Unsat_Zone import perc, pref, comb
from scipy.signal import fftconvolve

def IRF(parameters):
    # Unpack all the parameters that should be calibrated
    A = 10.0** parameters['A'].value
    a = 10.0** parameters['a'].value
    n = parameters['n'].value
    b = parameters['b'].value
    
    time_model = np.arange(0,10000)
    Fs = A * gammainc(n, time_model/a) # Step response function based on pearsonIII
    Fs = np.append(np.zeros(b), Fs)
    return Fs[1:] - Fs[0:-1] #block reponse function
    #Fs = np.cumsum((A*time_model**(n-1)*np.exp(-time_model/a))/(((n-1)*a)**(n-1)*np.exp(1-n)))
    #Fs[time_model>=b] = np.cumsum(A * (time_model[time_model>=b]-b)**(n-1)*np.exp(-(time_model[time_model>=b]-b)/a)/gamma(n))


def linear(parameters, InputData):
    # Unpack all the parameters that should be calibrated  
    d = parameters['d'].value
    f = parameters['f'].value
    
    #unpack all the data that is needed for the simulation
    P = InputData[1]
    E = InputData[2]
    
    #Recharge model
    recharge = P - f * E     
    
    Fb = IRF(parameters) #block reponse function
    head_modeled = d + fftconvolve(recharge,Fb)  
    return [head_modeled, recharge]

  
    
def preferential(parameters, InputData):
    # Unpack all the parameters that should be calibrated    
    d = parameters['d'].value
    Srmax = parameters['Srmax'].value
    Beta = parameters['beta'].value
    Imax = parameters['imax'].value

    # unpack all the data that is needed for the simulation
    time_model = InputData[0]
    P = InputData[1]
    E = InputData[2]
    solver = InputData[3]    
    dt= 1 
    
    #Recharge model
    recharge = pref(time_model, P, E, Srmax, Beta, Imax , dt, solver)[0]
    Fb = IRF(parameters) #block reponse function
    head_modeled = d + fftconvolve(recharge,Fb)    
    return [head_modeled, recharge]    
    
def percolation(parameters, InputData):
    # Unpack all the parameters that should be calibrated    
    d = parameters['d'].value
    Srmax = parameters['Srmax'].value
    Kp = 10.0**parameters['Kp'].value
    Gamma = parameters['gamma'].value         
    Imax = parameters['imax'].value

    # unpack all the data that is needed for the simulation
    time_model = InputData[0]
    P = InputData[1]
    E = InputData[2]
    solver = InputData[3]
    dt= 1 
    
    #Recharge model
    recharge = perc(time_model, P, E, Srmax, Kp, Gamma, Imax , dt, solver)[0]
    
    Fb = IRF(parameters) #block reponse function
    head_modeled = d + fftconvolve(recharge,Fb)    
    return [head_modeled, recharge]
    
    
def combination(parameters, InputData):
    # Unpack all the parameters that should be calibrated    
    d = parameters['d'].value
    Srmax = parameters['Srmax'].value
    Kp = 10.0**parameters['Kp'].value
    Beta = parameters['beta'].value
    Gamma = parameters['gamma'].value       
    Imax = parameters['imax'].value

    # unpack all the data that is needed for the simulation
    time_model = InputData[0]
    P = InputData[1]
    E = InputData[2]
    solver = InputData[3]    
    dt= 1 
    
    #Recharge model
    recharge = comb(time_model, P, E, Srmax, Kp, Beta, Gamma, Imax , dt, solver)[0]

    Fb = IRF(parameters) #block reponse function
    head_modeled = d + fftconvolve(recharge,Fb)    
    return [head_modeled, recharge]         