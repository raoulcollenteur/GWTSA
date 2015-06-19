# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 10:08:36 2015

@author: Raoul
"""
import numpy as np
import matplotlib.dates as md
import matplotlib.pyplot as plt
from scipy.optimize import fmin
from scipy.special import gammainc
import TFN_Model
import cma


class Model:

    def __init__(self, bore, forcing, rows=[5, 8], timestart=1000):
        
        """
        Prepares the time series model with the forcings and observed heads
        
        Parameters
        ----------
        bore: string
            name of the txt file containing the observed head values. 
        forcing: string
            name of the txt file containing the forcing values. 
        rows: list
            list with the number of rows to skip when reading the txt files
            E.g. [5,8] #skip 5 lines in bore and 8 lines in forcing txt file
        timestart: int
            warmup time for model in days    
        
        Returns
        -------
           - Precipitation values (self.P)
           - Potential Evaporation values (self.E)
           - Name of borehole (self.bore)
        
        See Also
        --------
        See the documentation of GWTSA for a description and examples of how the 
        input files for the observed heads and the forcings should look.
        
        """
        
        Ho = np.genfromtxt(('./%s.txt' % (bore)), delimiter=',', skiprows=rows[0], usecols=[0, 1], converters={0: md.strpdate2num('%Y%m%d')});
        Ho = Ho[Ho[:,1] > -999] #Select only real values
        self._time_begin = Ho[0,0] #In time number, not date
        self._time_end = Ho[-1,0] #In time number, not date
        self._time_steps = Ho[1:,0] - Ho[0:-1,0] #Calculate the timesteps
        self._time_axis = Ho[:,0] # Real Time in time number format for plots
        self._time_observed = np.append(0, [np.cumsum(self._time_steps)]).astype(int) #Timesteps with observed heads
        self.__time_model = np.arange(0,sum(self._time_steps,1), dtype=int)     
        self._time_start = timestart; #warmup time of the model [Days]
        self.head_observed = Ho[:,1]
               
        # Import the precipitation and evaporation data and calculate the recharge
        ClimateData = np.genfromtxt('./%s.txt' % forcing , delimiter=',', skiprows=rows[1], converters={1: md.datestr2num}, usecols=[1, 2, 3]);
        ClimateData = ClimateData[:,:][(ClimateData[:,0] >= self._time_begin) & (ClimateData[:,0] <= self._time_end)]
        self.precipitation = ClimateData[:,1] / 10000.0 # divided by 10000 to make it in meters
        self.precipitation[self.precipitation < 0.0] = 0.0
        self.evaporation = ClimateData[:,2] / 10000.0 # divided by 10000 to make it in meters
        self.evaporation[self.evaporation < 0.0] = 0.0

        self.bore = bore
        
    def __repr__(self):
        return 'Time Series Data of Bore: %s' %self.bore
        
    def solve(self, TFN, X0, method=0, correlation=1):
        
        """ 
        Solves the time series model
        
        parameters
        ----------     
        TFN: string 
            define the name of the Transfer function noise model you want to use (see TFN_Model.py)
        X0: dictionary with initial parameter guesses
            E.g. X0 = {'A': 20,'a': 10, 'n': 1.5,'Alpha':40}. The order does not matter as solve() 
            constructs an array with the initial guesses based on the names.
        method: integer 0 or 1
            0 = fmin method, 1 = cma.fmin method (much slower, but does not get stuck in local opt.)
        correlation: 0 or 1
            return correlation matrix (self.correlation_matrix)
            
        Returns
        -------
            - an array of the optimum parameters set (self.Parameters)
            - array with modeled heads (self.head_modeled)
            - array with innivations (self.innovations)
        
        See Also
        --------
        Help(TFN_Model) can give you information on how the parameter dictionary should look, as
        that depends on the Transfer Functions Noise Model that is used. 
        
        """
        # Define the TFN Model and the Initial Parameters
        if TFN == 'TFN1':
            initial_parameters = [X0['A'],X0['a'],X0['n'],np.mean(self.head_observed),X0['Alpha']]
        elif TFN == 'TFN2':
            initial_parameters = [X0['A'],X0['a'],X0['n'],np.mean(self.head_observed),X0['Alpha'], X0['f']]
        elif TFN == 'TFN3':
            initial_parameters = [X0['A'],X0['a'],X0['n'],np.mean(self.head_observed),X0['Alpha']]
        elif TFN == 'TFN4':
            initial_parameters = [X0['A'],X0['a'],X0['n'],np.mean(self.head_observed),X0['Alpha'],X0['S_cap'],X0['K_sat'], X0['Beta'], X0['D']]           
        else:
            print 'Error: TFN model does not exist, chose another one or check TFN_Model.py!'
        
        self.TFN = getattr(TFN_Model, TFN)    
    
        InputData = [self.TFN, self._time_model, self.precipitation, self.evaporation, self.head_observed,self._time_observed, self._time_steps, self._time_start]
        
        if method == 0:
            self.Parameters_opt = fmin(self.swsi, initial_parameters, args= (InputData,), maxiter= 1000 )      
        elif method == 1: 
            res = cma.fmin(self.swsi, initial_parameters, 0.5, args=(InputData,), options={'ftarget': 1e-5})
            self.Parameters_opt = res[0]

        #Print the output parameters and compare to the real values'''
        self.Parameter_Names = ['A','a','n','d','alpha','S_cap', 'K_sat','Beta','D']
        #print 'Explained Variance Percentage is:', self.explained_variance(self.Parameters,InputData)
        for i in range(0,len(self.Parameters_opt)):
            print self.Parameter_Names[i], '=', self.Parameters_opt[i]  
        
        if Cor == 0 and method == 1:
            self.correlation_matrix = res[-2].correlation_matrix()
            print 'The Correlation Matrix:'
            print self.correlation_matrix
        
    def monte_carlo(self, TFN, X0, n=1000):
        self.TFN = getattr(TFN_Model, TFN)    
        InputData = [self.TFN, self._time_model, self.P, self.E, self.head_observed, self._time_observed, self._time_steps, self._time_start]
        
        x = np.random.uniform(0.0, 1.0, (n,len(X0[0])))
        self.Parameters = x * (X0[1] - X0[0]) + X0[0]
        self.result = []
        for i in range(n):
            self.result = np.append(self.result, self.swsi(self.Parameters[i], InputData))
        self.Parameter_opt = self.Parameters[self.result.argmin()]
        print 'SWSI is:', self.result[self.result.argmin()]
        return self.result, self.Parameters    
        
        
    
    '''This section contains the objective functions and diagnostic tests.'''

# swsi constitutes the adapted version of the Sum of weighted squared innovations (swsi) developed by asmuth et al. (2005). For large values of alpha and small timesteps the numerator approaches zero. Therefore, Peterson et al. (2014) adapted the swsi function, making the numerator a natural log and changing the product operator to a summations.
   
    def swsi(self, Parameters, InputData):
        TFN = InputData[0]
        innovation = TFN(Parameters, InputData, solver=0)[1]
        N = len(innovation)
        alpha = 10**Parameters[4]
        dt = InputData[6][-N:]
        x = np.exp(sum(np.log(1 - np.exp(-2.0  / alpha * dt)))*(1.0/N))
        swsi = sum( (x / (1 - np.exp(-2.0 / alpha * dt ))) * innovation**2)
        #print swsi2 #Print innovation to explore source of nan-values?!
        return swsi 

# The Explained Variance Percentage (explained_variance) method compares variance of the
# observed values to the variance of the residuals, that is observed values
# minus the modelled values. (Asmuth et al, 2012)

    def explained_variance(self, Parameters, InputData):
        TFN = InputData[0]
        to = InputData[5] #timesteps of observation
        Ho = InputData[4]
        tstart = InputData[7]
        istart = np.where(to > tstart)[0][0]
        Hm = TFN(Parameters, InputData, method = 0)[0]
        Hm = Hm[to[istart:]]
        Ho = Ho[istart:]
        E = np.array((np.var(Ho)**2 - np.var(Hm-Ho)**2)/np.var(Ho)**2*100)
        return E

# the test function is sometimes helpfull when you have parameter set and want to see what the modeled heads will result in. This function is experimental. You need to run a model first..    

    def simulate(self, Xt, TFN):
        self.TFN = getattr(TFN_Model, TFN)   
        InputData = [self.TFN, self._time_model, self.P, self.E, self.head_observed, self.Time_Observed, self.Time_Steps, self.Time_Start] 
        self.head_modeled, self.Innovations = self.TFN(Xt, InputData); #model the GWL
              
        ''' In this section the functions are defined that relate to the plotting of different forcings and results. Each function starts with plot_function to be able to quickly find these modules. '''        

    def plot_heads(self,color='r',observed=0, modeled=0, newfig=0):
        assert modeled == 0 or observed == 0, 'No heads are chosen to be plotted'
        if newfig == 0:
            plt.figure()
        if observed == 0:
            plt.plot(md.num2date(self.Time_Axis), self.head_observed, 'k.')
        if modeled == 0:
            plt.plot(md.num2date(np.arange(self.Time_Begin, self.Time_End+1)), 
                 self.head_modeled[self._time_model], '-', color=color)
        plt.legend(['Observed Head','Modeled Head'], loc=1)
        plt.xlabel('Time [T]', size=20)
        plt.ylabel('Groundwater Head [L]', size=20)
        plt.title('%s' % (self.bore))
        
    def plot_innovations(self):
        plt.plot(self.Innovations)
        plt.title('Innovations of the time series')
        plt.xlabel('Innovations [-]')             
        
    def plot_forcings(self):
        plt.figure()
        plt.bar(md.num2date( self.ClimateData[:,0]), self.ClimateData[:,1], color='b', lw=0)
        plt.ylabel('P [0.1 mm/D]', size=20)
        plt.ylim(0,max(self.ClimateData[:,1]))        
        ax1 = plt.gca()
        ax2 = ax1.twinx() # To create a second axis on the right
        plt.bar(md.num2date( self.ClimateData[:,0]), self.ClimateData[:,2], color='red', lw=0)
        ax2.set_ylim((0,400)[::-1])
        plt.ylabel('E [0.1 mm/D]',size=20)
        plt.xlabel('Time [Years]',size=20)
        plt.legend('Precipitation','Potential Evapotranspiration')
        plt.title('Forcings',size=20)
        
    def plot_impulseResponseFunction(self):
        Fs = self.Parameters[0] * self.Parameters[1] * gammainc(self.Parameters[2], self._time_model/self.Parameters[1])
        Fb = Fs[1:] - Fs[0:-1]
        plt.plot(self._time_model[0:-1],Fb)
