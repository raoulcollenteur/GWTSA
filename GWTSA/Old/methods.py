# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 12:26:48 2016

@author: Raoul
"""
import numpy as np
from scipy.special import gammainc
from scipy.stats import norm
from scipy.signal import fftconvolve
from Unsat_Zone import perc, pref, comb


class Methods(object):
    def __init__(self):
        pass
# %% Define de impulse response function (IR) here

    def IRF(self, parameters):
        # Unpack all the parameters that should be calibrated
        A = parameters['A'].value
        a = 10**parameters['a'].value
        n = parameters['n'].value
        t = np.arange(1.,10000)
        Fs = A * t**n * (t/a)**-n * gammainc(n, t/a) # Step response function based on pearsonIII
        return np.append(0, Fs[1:] - Fs[0:-1]) #block reponse function
        
    def IRF1(self, parameters):
        """The IRF1 function is a modified version of the general IRF, 
        where the T_peak is given as an parameter instead of 'a'. 
        """
        # Unpack all the parameters that should be calibrated
        A = parameters['A'].value
        n = parameters['n'].value
        t_p = 10**parameters['t_p'].value
        t = np.arange(1.,10000)
        Fs = A * t**n * (t*(n-1)/t_p)**-n * gammainc(n, (t*(n-1)/t_p)) 
        return np.append(0, Fs[1:] - Fs[0:-1]) #block reponse function    

    def IRF2(self, parameters, recharge):
        """The IRF2 function can be used to simulate the effect of the
        percolation zone. Followed by an exponential decay function.
        """
        mu = parameters['mu'].value
        sig = parameters['sig'].value

        # Percolation impulse response    
        Fs = norm.cdf(range(1000), mu, sig)  # /norm.pdf(mu,mu,sig)
        Fb = Fs[1:] - Fs[0:-1]  # block reponse function
        return fftconvolve(recharge, Fb)[0:len(recharge)]

    def IRF3(self, parameters):
        A1 = parameters['A1'].value
        a1 = parameters['a1'].value
        n1 = parameters['n1'].value
        t = np.arange(1.0, 10000)
        Fs = A1 * t**n1 * (t/a1)**-n1 * gammainc(n1, t/a1)  # Step response function based on pearsonIII
        return np.append(0, Fs[1:] - Fs[0:-1])  # block reponse function

# %% Define the Recharge Models (RM) here.

    def linear(self, parameters):
        f = parameters['f'].value
        P = np.array(self.data.P)
        E = np.array(self.data.E)
        recharge = P - f * E
        return recharge

    def preferential(self, parameters):
        Srmax = parameters['Srmax'].value
        Beta = parameters['Beta'].value
        Imax = parameters['Imax'].value
        t = np.array(self._t)
        P = np.array(self.data.P)
        E = np.array(self.data.E)
        solver = 1
        dt = 1
        recharge = pref(t, P, E, Srmax, Beta, Imax, dt, solver)[0]
        return recharge

    def percolation(self, parameters):
        Srmax = parameters['Srmax'].value
        Kp = parameters['Kp'].value
        Gamma = parameters['Gamma'].value
        Imax = parameters['Imax'].value
        t = np.array(self._t)
        P = np.array(self.data.P)
        E = np.array(self.data.E)
        solver = 1
        dt = 1
        recharge = perc(t, P, E, Srmax, Kp, Gamma, Imax, dt, solver)[0]
        return recharge

    def combination(self, parameters):
        Srmax = parameters['Srmax'].value
        Kp = parameters['Kp'].value
        Beta = parameters['Beta'].value
        Gamma = parameters['Gamma'].value
        Imax = parameters['Imax'].value
        t = np.array(self._t)
        P = np.array(self.data.P)
        E = np.array(self.data.E)
        solver = 1
        dt = 1
        Rs, Rf = comb(t, P, E, Srmax, Kp, Beta, Gamma, Imax, dt, solver)[0:2]
        recharge = Rs + Rf
        return recharge
        
    

# %% Define Alternative model additions (E.g. Wells, linear slope, reclamation)

    def linear_slope(self, parameters):
        slope = parameters['slope'].value
        intercept = parameters['intercept'].value
        t = np.array(self._t)
        trend = slope * t + intercept
        trend[trend > 0.0] = 0.0
        trend[-25*365:] = trend[-15*365]  # Alternatively, let the rise last for a couple of years
        return trend

    def reclamation(self, parameters):
        B = parameters['B'].value
        b = 10**parameters['b'].value
        t = np.array(self._t)
        t_start = parameters['t_start'].value
        Fb = B*(1.-np.exp(-(t-t_start)/b))
        Fb[Fb > 0.0] = 0.0
        return Fb

    def well(self, parameters):
        B = parameters['B'].value
        b = parameters['b'].value
        t = np.array(self._t)
        Fi = Fi = B/t[1:] * np.exp(-b/t[1:])
        Fs = np.cumsum(Fi)
        Fb = np.append(0, Fs[1:] - Fs[0:-1])
        drawdown = fftconvolve(self.data.discharge, Fb)[t]
        return drawdown
