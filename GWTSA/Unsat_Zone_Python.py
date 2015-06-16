# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 18:20:35 2015

This python script is used to solve the differential equations that are used for the unsaturated zone module. The output of the function is the soil state (how full it is) at each timestep, and the groundwater recharge R.

-------------------------------------Two models can be used:
- Piston Flow
dS/dt = (P-I) - K_sat * (S/S_cap)**Beta - E_p * min(1, S/0.5S_cap)

- Percolation
dS/ Dt = P[t+1] * (1 - (S[t] / S_cap)**Beta)- E_p * min(1, S/0.5S_cap)

-------------------------------------Numerical info:
The Soil module is is solved with an implicit euler and Newton Raphson iteration. The initial guesstimate for the the NR-iteration is provided by an Explicit Euler solution of the above differential equation. 

-------------------------------------To Do:
    - Solve the Zero Division error by applying a bisection method
    - Solve when S[t+1] > S_Cap or S[t+1] < 0.00
    - Compile a C-file of this functin for quick calculation using Cython
    - Built in more external / internal checks for water balance

-------------------------------------References:     
- Kavetski, D., Kuczera, G. & Franks, S.W. [2006]. Calibration of conceptuyal hydrological models revisited: 1. Overcoming numerical artefacts. Journal of Hydrology, 320, p. 173-186.  

@author: Raoul Collenteur
"""


import numpy as np

    
''' -----------------------------------------------
In this section the percolation model is defined. 
dS/dt = (P-I) - K_sat * (S/S_cap)**Beta - E_p * min(1, S/0.5S_cap) 
----------------------------------------------- '''  

def percolation(Time_Model, P, E, S_cap = 1.0, K_sat = -1.5, Beta = 2.0, D = 0.0001 , dt = 1):
    """ This is the percolation function docstring 
    """    
    K_sat = 10.0**K_sat
    S_cap = 10.0**S_cap
    D = 10.0**D
    error = 1.0e-5
    # Create an empty array to store the soil state in
    S = np.zeros(len(Time_Model)/dt) 
    S[0] = 0.5 * S_cap   #Set the initial system state
    R = np.zeros(len(Time_Model)/dt)
    R[0] = 0.0    
    # Calculate the interception, minimum of P, E or D
    D1 = np.ones(len(P)) * D # Create an array for element-wise comparison
    I = np.amin([P,D1,E], axis = 0)
    
    Ei = E - I # Update the amount of potential evaporation

    for t in np.arange(0,len(Time_Model)-1,dt):
        assert S[t] <= S_cap, 'Error, soil is saturated and overland flow will occur S[t], K_sat, S_cap, Beta, D^= %r,' % [S[t], K_sat, S_cap, Beta, D, t]
        Last_S = S[t]
        iteration = 0
        bisection = 1
        #Use explicit Euler scheme to find an initial estimate for the newton raphson-method
        
        S[t+1] = np.max([0.0, S[t]+ dt * ( (P[t]-I[t]) - K_sat * (S[t] / S_cap)**Beta - Ei[t] * np.min([1.0, (S[t] / (0.5 * S_cap))]))])
       
        #Start the while loop for the newton-Raphson iteration   
        while (abs(Last_S - S[t+1]) > error):
            if iteration > 100:
                break #Check if the number of iterations is not too high  
            iteration += 1
            Last_S = S[t+1]      
       
            g = Last_S - S[t] - dt *( (P[t]-I[t]) - K_sat * (Last_S / S_cap)**Beta - Ei[t] * np.min([1, (Last_S / (0.5 * S_cap))]) )       
            # Derivative depends on the state of the system
            if Last_S > (0.5 * S_cap):
                g_derivative = 1.0 - dt * ( -Beta * K_sat * (Last_S / S_cap)**(Beta-1))
            else:        
                g_derivative = 1.0 - dt * ( -Beta * K_sat * (Last_S / S_cap)**(Beta-1) - Ei[t] * (0.5 * S_cap) ) 
            
            # Check if there is no zero-division error            
            if np.isnan(g / g_derivative):
                bisection = 0
                break
            # if there is no zero-division error                
            else: # use newton raphson
                S[t+1] = Last_S - g / g_derivative                
        
        if bisection == 0:
            iteration = 0
            a = S[t]
            b = S[t+1]
            c = a + b / 2.0      
            
            while ((b - a) / 2.0) > 1e-6:
                if iteration > 100:
                    print 'iteration in bisection method exceeded 100', iteration                        
                    break   
                iteration += 1 #increase the number of iterations by 1

                if (c - S[t] - dt *( (P[t]-I[t]) - K_sat * (c / S_cap)**Beta - Ei[t] * np.min([1, (c / (0.5 * S_cap))]) )) == 0:
                    return c # Return the current value if it is correct
                elif (a - S[t] - dt *( (P[t]-I[t]) - K_sat * (a / S_cap)**Beta - Ei[t] * np.min([1, (a / (0.5 * S_cap))]) )) * (c - S[t] - dt *( (P[t]-I) - K_sat * (c / S_cap)**Beta - Ei[t] * np.min([1.0, (c / (0.5 * S_cap))]) )) > 0.0 :
                    b = c
                else : 
                    a = c
                    
                c = a + b / 2.0 
                
            S[t+1] = c    
                
        assert ~np.isnan(S[t+1]), 'NaN-value calculated for soil state'       
        
        S[t+1] = np.min([S_cap, np.max([0.0,S[t+1]])]) #Make sure the solution is larger then 0.0 and smaller than S_cap

        R[t+1] = K_sat * dt / 2 * ((S[t] + S[t+1]) / S_cap) ** Beta #This can be written to be outside the loop for quicker calculation
    
    return R, S


''' -----------------------------------------------
In this section the piston flow model is defined.
dS/ Dt = P[t+1] * (1 - (S[t] / S_cap)**Beta)- E_p * min(1, S/0.5S_cap)
----------------------------------------------- '''


def piston_flow(Time_Model, P, E, S_cap = 1.0, Beta = 2.0, D = 0.001 , dt = 1):

    e = 1.0e-5 
    
    # Create an empty array to store the soil state in
    S = np.zeros(len(Time_Model)/dt) 
    S[0] = 0.2 * S_cap   #Set the initial system state
    R = np.zeros(len(Time_Model)/dt)
    R[0] = 0
    
    for t in np.arange(0,len(Time_Model)-1,dt):
        assert S[t] <= S_cap, 'Error, soil is saturated and overland flow will occur'
        Last_S = S[t]
        iteration = 0
        
        # Calculate the interception, minimum of P, E or D
        I = np.min([P[t],D,E[t]])        
        E[t] = E[t] - I # Update the amount of potential evaporation
        assert E[t] >= 0.0, 'Error,A negative value in the Interception bucket is detected'
        
        
        #Use explicit Euler scheme to find an initial estimate for the newton raphson-method
        
        S[t+1] = S[t]+ dt * ( (P[t]-I) * (1 - (S[t] / S_cap)**Beta) - E[t] * np.min([1, (S[t] / (0.5 * S_cap))]))
       
    #Start the while loop for the newton-Raphson iteration   
        while (abs(Last_S - S[t+1]) > e):
            #print iteration
            if iteration > 100:
                #Check if the number of iterations is not too high
                break   
            iteration += 1 #increase the number of iterations by 1
            Last_S = S[t+1]       
       
            g = Last_S - S[t] - dt *( (P[t]-I) * ( 1 - (Last_S / S_cap)**Beta) - E[t] * np.min([1, (Last_S / (0.5 * S_cap))]) )       
       
       
            # Derivative depends on the state of the system
            if Last_S > (0.5 * S_cap):
                g_derivative = 1.0 - dt * ( -(P[t]-I) * Beta  * (Last_S / S_cap)**(Beta-1))
            else:        
                g_derivative = 1.0 - dt * ( -(P[t]-I) * Beta * (Last_S / S_cap)**(Beta-1) - E[t] * (0.5 * S_cap) ) 
        
            S[t+1] = Last_S - (g / g_derivative)
       
        R[t+1] = (P[t+1] - I) * dt / 2 * ( (S[t]+S[t+1]) / S_cap) ** Beta      
       
    return R