# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 14:34:13 2014
@author: rcollenteur
"""

#Import all the packages needed for the simulation
import numpy as np
from scipy.special import gammainc
from scipy.stats import norm
from Unsat_Zone import perc, pref, comb
from scipy.signal import fftconvolve

#%% Define the basic model that construct the model from different options in Impulse response and recharge simulation models. 
def construct_model(parameters, InputData):
    d = parameters['d'].value    
    time_model = InputData[0]
    
    ImpulseResponse = eval(InputData[4])
    RechargeModel = eval(InputData[5])
    trend = InputData[6]
    
    Fb = ImpulseResponse(parameters) #block reponse function
    recharge = RechargeModel(parameters, InputData) #Simulate the recharge series
    
    head_modeled = (d + fftconvolve(recharge,Fb))[time_model] 
    if trend != None:    
        trend = eval(trend)(parameters, InputData)
        head_modeled += trend
        return [1, head_modeled, recharge, trend]
    else:
        return [0, head_modeled, recharge, trend]

#%% Define the different impulse Response Function
   
def IRF(parameters):
    # Unpack all the parameters that should be calibrated
    A = parameters['A'].value
    n = parameters['n'].value
    t_p = 10**parameters['t_p'].value
    t = np.arange(1.,10000)
    Fs = A * t**n * (t*(n-1)/t_p)**-n * gammainc(n, (t*(n-1)/t_p)) 
    return np.append(0, Fs[1:] - Fs[0:-1]) #block reponse function    

def IRF1(parameters):
    # Unpack all the parameters that should be calibrated
    A = 10.0** parameters['A'].value
    a = 10.0** parameters['a'].value
    n = parameters['n'].value
    time_model = np.arange(0,10000)
    Fs = A * (a**n) * gammainc(n, time_model/a) # Step response function based on pearsonIII
    return Fs[1:] - Fs[0:-1] #block reponse function

def IRF2(parameters, recharge):    
    mu = parameters['mu'].value
    sig = parameters['sig'].value
    
    #Percolation impulse response    
    Fs = norm.cdf(range(1000), mu, sig) #/norm.pdf(mu,mu,sig)
    Fb = Fs[1:] - Fs[0:-1] #block reponse function
    return fftconvolve(recharge,Fb)[0:len(recharge)]
    
#%% Define the different recharge models

def linear(parameters, InputData):
    # Unpack all the parameters that should be calibrated  
    f = parameters['f'].value    
    #unpack all the data that is needed for the simulation
    P = InputData[1]
    E = InputData[2]    
    #Recharge model
    recharge = P - f * E     
    return recharge

def preferential(parameters, InputData):
    # Unpack all the parameters that should be calibrated    
    Srmax = parameters['Srmax'].value
    Beta = parameters['Beta'].value
    Imax = parameters['Imax'].value
    # unpack all the data that is needed for the simulation
    time_model = InputData[0]
    P = InputData[1]
    E = InputData[2]
    solver = InputData[3]    
    dt= 1   
    #Recharge model
    recharge = pref(time_model, P, E, Srmax, Beta, Imax , dt, solver)[0]
    return recharge   
    
def percolation(parameters, InputData):
    # Unpack all the parameters that should be calibrated    
    Srmax = parameters['Srmax'].value
    Kp = 10.0**parameters['Kp'].value
    Gamma = parameters['gamma'].value         
    Imax = parameters['Imax'].value
    # unpack all the data that is needed for the simulation
    time_model = InputData[0]
    P = InputData[1]
    E = InputData[2]
    solver = InputData[3]
    dt= 1 
    #Recharge model
    recharge = perc(time_model, P, E, Srmax, Kp, Gamma, Imax , dt, solver)[0]
    return recharge
    
def combination(parameters, InputData):
    # Unpack all the parameters that should be calibrated    
    Srmax = parameters['Srmax'].value
    Kp = 10.0**parameters['Kp'].value
    Beta = parameters['Beta'].value
    Gamma = parameters['gamma'].value       
    Imax = parameters['Imax'].value
    # unpack all the data that is needed for the simulation
    time_model = InputData[0]
    P = InputData[1]
    E = InputData[2]
    solver = InputData[3]    
    dt= 1 
    #Recharge model
    Rs, Rf= comb(time_model, P, E, Srmax, Kp, Beta, Gamma, Imax , dt, solver)[0:1]  
    recharge = Rs + Rf
    return recharge      
    
#%% Define Alternative model additions (E.g. Wells, linear slope, reclamation)    
    
def linear_slope(parameters, InputData):
    slope = parameters['slope'].value
    intercept = parameters['intercept'].value
    time_model = InputData[0]    
    head_modeled = slope * time_model + intercept
    head_modeled[head_modeled<0.0] = 0.0
    return head_modeled
          
    
def reclamation(parameters, InputData):
    B = parameters['B'].value
    b = 10**parameters['b'].value
    time_model = InputData[0]
    t_start = parameters['t_start'].value
    Fb = B*(1.-np.exp(-(time_model-t_start)/b))
    Fb[Fb>0.0]=0.0
    return Fb    
        