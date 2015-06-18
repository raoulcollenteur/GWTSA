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

 -----------------------------------------------
In this section the percolation model is defined. 
dS/dt = (P-I) - K_sat * (S/S_cap)**Beta - E_p * min(1, S/0.5S_cap) 
----------------------------------------------- """


# cython: profile=True
# filename: Unsat_Zone.pyx

import numpy as np
cimport numpy as np


# Define some C-function for more efficient computatation
cdef inline double c_max(double a, double b): return a if a >= b else b
cdef inline double c_min(double a, double b): return a if a <= b else b


def percolation(np.ndarray[np.int_t] Time_Model, np.ndarray[np.float_t, ndim=1] P, np.ndarray[np.float_t, ndim=1] Ep, double S_cap = 1.0, double K_sat = -1.5, double Beta = 2.0, double Imax = -3, int dt = 1, int solver = 0):
    
    cdef int t, iteration, bisection, n
    cdef double error, Last_S, g, g_derivative, a, b, c
    
    n = len(Time_Model) / dt

    K_sat = 10.0**K_sat
    S_cap = 10.0**S_cap
    Imax = 10.0**Imax
    error = 1.0e-5
    
    # Create an empty array to store the soil state in
    cdef np.ndarray[np.float_t] S = np.zeros(n)
    S[0] = 0.5 * S_cap   #Set the initial system state
    cdef np.ndarray[np.float_t] Si = np.zeros(n)
    Si[0] = 0.0
    cdef np.ndarray[np.float_t] Pe = np.zeros(n)
    Pe[0] = 0.0    
    cdef np.ndarray[np.float_t] Ei = np.zeros(n)
    Ei[0] = 0.0  
    
    for t in range(n-1):
        Si[t+1] = Si[t] + P[t+1]                # Fill interception bucket with new rain
        Pe[t+1] = c_max(0.0, Si[t+1] - Imax)    # Calculate effective precipitation
        Si[t+1] = Si[t+1] - Pe[t+1]             
        Ei[t+1] = c_min(Si[t+1], Ep[t+1])       # Evaporation from interception
        Si[t+1] = Si[t+1] - Ei[t+1]             # Update interception state
        Ep[t+1] = Ep[t+1] - Ei[t+1]             # Update potential evapotranspiration    
    
        Last_S = S[t]
        iteration = 0
        bisection = 1
        #Use explicit Euler scheme to find an initial estimate for the newton raphson-method
        
        S[t+1] = c_max(0.0, S[t] + dt * ( Pe[t] - K_sat * (S[t] / S_cap)**Beta - Ep[t] * c_min(1.0, (S[t] / (0.5 * S_cap)) )))
        
        if solver == 1:
            #Start the while loop for the newton-Raphson iteration   
            while abs(Last_S - S[t+1]) > error:
                if iteration > 100:
                    break #Check if the number of iterations is not too high  
                iteration += 1
                Last_S = S[t+1]      
           
                g = Last_S - S[t] - dt *( Pe[t] - K_sat * (Last_S / S_cap)**Beta - Ep[t] * c_min(1, (Last_S / (0.5 * S_cap))) )       
                # Derivative depends on the state of the system
                if Last_S > (0.5 * S_cap):
                    g_derivative = 1.0 - dt * ( -Beta * K_sat * (Last_S / S_cap)**(Beta-1))
                else:        
                    g_derivative = 1.0 - dt * ( -Beta * K_sat * (Last_S / S_cap)**(Beta-1) - Ep[t] * (0.5 * S_cap) ) 
                
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
                
                while ((b - a)/2.0) > error:
                    if iteration > 100:
                        print 'iteration in bisection method exceeded 100', iteration                        
                        break   
                    iteration += 1 #increase the number of iterations by 1
    
                    if (c - S[t] - dt *( Pe[t] - K_sat * (c / S_cap)**Beta - Ep[t] * c_min(1, (c / (0.5 * S_cap))) )) == 0.0:
                        return c # Return the current value if it is correct
                    elif (a - S[t] - dt *( Pe[t] - K_sat * (a / S_cap)**Beta - Ep[t] * c_min(1.0, (a / (0.5 * S_cap))) )) * (c - S[t] - dt *( Pe[t] - K_sat * (c / S_cap)**Beta - Ep[t] * c_min(1.0, (c / (0.5 * S_cap))) )) > 0.0 :
                        b = c
                    else : 
                        a = c
                        
                    c = a + b / 2.0 
                    
                S[t+1] = c    
                    
            assert ~np.isnan(S[t+1]), 'NaN-value calculated for soil state'       
        
        S[t+1] = c_min(S_cap, c_max(0.0,S[t+1])) #Make sure the solution is larger then 0.0 and smaller than S_cap
        
    cdef np.ndarray[np.float_t] R = np.append(0.0, K_sat * dt / 2 * ((S[:-1] + S[1:]) / S_cap) ** Beta )  
    
    return R, S