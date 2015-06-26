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

#%%
class Model:

    def __init__(self, bore, forcing, rows=[5, 8], timestart=2000):
        
        """
        Prepares the time series model with the forcings and observed heads
        
        parameters
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
        self._time_model = np.arange(0,sum(self._time_steps,1), dtype=int)     
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
 
#%% 
       
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
            - an array of the optimum parameters set (self.parameters)
            - array with modeled heads (self.head_modeled)
            - array with innivations (self.innovations)
        
        See Also
        --------
        Help(TFN_Model) can give you information on how the parameter dictionary should look, as
        that depends on the Transfer Functions Noise Model that is used. 
        
        """
        # Define the TFN Model and the Initial parameters
        if TFN == 'TFN1':
            initial_parameters = [X0['A'],X0['a'],X0['n'],np.mean(self.head_observed),X0['Alpha']]
        elif TFN == 'TFN2':
            initial_parameters = [X0['A'],X0['a'],X0['n'],np.mean(self.head_observed),X0['Alpha'], X0['f']]
            self.Parameter_Names = ['A','a','n','d','alpha','f']
        elif TFN == 'TFN3':
            initial_parameters = [X0['A'],X0['a'],X0['n'],np.mean(self.head_observed),X0['Alpha']]
        elif TFN == 'TFN4':
            initial_parameters = [X0['A'],X0['a'],X0['n'],np.mean(self.head_observed),X0['Alpha'],X0['S_cap'],X0['K_sat'], X0['Beta'], X0['D']]
            self.Parameter_Names = ['A','a','n','d','alpha','S_cap', 'K_sat','Beta','D']
        else:
            print 'Error: TFN model does not exist, chose another one or check TFN_Model.py!'
        
        self.TFN = getattr(TFN_Model, TFN)
        self._TFN = TFN
    
        InputData = [self.TFN, self._time_model, self.precipitation, self.evaporation, self.head_observed,self._time_observed, self._time_steps, self._time_start]
        
        if method == 0:
            self.parameters_opt= fmin(self.swsi, initial_parameters, args= (InputData,), maxiter= 10000)      
        elif method == 1: 
            res = cma.fmin(self.swsi, initial_parameters, 0.5, args=(InputData,), options={'ftarget': 1e-5})
            self.parameters_opt = res[0]

        if correlation == 0 and method == 1:
            self.correlation_matrix = res[-2].correlation_matrix()
            print 'The Correlation Matrix:'
            print self.correlation_matrix
        
    def monte_carlo(self, TFN, X0, n=1000):
        """
        Runs a monte carlo analysis
        
        parameters
        ----------     
        TFN: string 
            define the name of the Transfer function noise model you want to use (see TFN_Model.py)
        X0: array
            An two row array with the lower and upper bounds of the parameters            
            E.g. X0 = np.array([[2.0, 8.0, 1.0],[2.6, 10.0]])
        n: integer
            number of runs that are performed.
            
        Returns
        -------
            - an array of the optimum parameters set (self.parameters_opt)
            - an array with the results of the objective function
            - an array with al the used parameters
        
        See Also
        --------
        
        """
        self.TFN = getattr(TFN_Model, TFN)    
        InputData = [self.TFN, self._time_model, self.precipitation, self.evaporation, self.head_observed,self._time_observed, self._time_steps, self._time_start]
        
        rnd = np.random.uniform(0.0, 1.0, (n,len(X0[0])))
        self.parameters = rnd * (X0[1] - X0[0]) + X0[0]
        self.result = []
        
        for i in range(n):
            self.result = np.append(self.result, self.swsi(self.parameters[i], InputData))
        
        self.parameters_opt = self.parameters[self.result.argmin()]
        print 'SWSI is:', self.result[self.result.argmin()]
        return self.result, self.parameters            
    
    def simulate(self, TFN, parameters):
        """
        Simulate the groundwater levels with a certain parameter set
        """
        self.TFN = getattr(TFN_Model, TFN)   
        InputData = [self.TFN, self._time_model, self.precipitation, self.evaporation, self.head_observed,self._time_observed, self._time_steps, self._time_start]
        
        self.head_modeled, self.innovations, self.recharge, self.residuals = self.TFN(parameters, InputData);   
        
        # Calculate the Explained Variance Percentage
        self.explained_variance = (np.var(self.head_observed)**2 - np.var(self.head_modeled[self._time_observed]-self.head_observed)**2)/np.var(self.head_observed)**2*100
        
        self.swsi_value = self.swsi(self.parameters_opt, InputData)

    '''This section contains the objective functions and diagnostic tests.'''

#%% swsi constitutes the adapted version of the Sum of weighted squared innovations (swsi) developed by asmuth et al. (2005). For large values of alpha and small timesteps the numerator approaches zero. Therefore, Peterson et al. (2014) adapted the swsi function, making the numerator a natural log and changing the product operator to a summations.
   
    def swsi(self, parameters, InputData):
        TFN = InputData[0]
        spinup = np.where(InputData[5]  > InputData[7])[0][0] # Where 
        innovation = TFN(parameters, InputData, solver=0)[1][spinup:-1:3]
        N = len(innovation)
        alpha = 10**parameters[4]
        dt = InputData[6][-N:]
        x = np.exp(sum(np.log(1 - np.exp(-2.0  / alpha * dt)))*(1.0/N))
        swsi = sum( (x / (1 - np.exp(-2.0 / alpha * dt ))) * innovation**2)
        return swsi 
              
#%% In this section the functions are defined that relate to the plotting of different forcings and results. Each function starts with plot_function to be able to quickly find these modules.        

    def plot_heads(self,color='r',observed=0, modeled=0, newfig=0):
        assert modeled == 0 or observed == 0, 'No heads are chosen to be plotted'
        if newfig == 0:
            plt.figure()
        if observed == 0:
            plt.plot(md.num2date(self._time_axis), self.head_observed, 'k.')
        if modeled == 0:
            plt.plot(md.num2date(np.arange(self._time_begin, self._time_end+1)), 
                 self.head_modeled[self._time_model], '-', color=color)
        plt.legend(['Observed Head','Modeled Head'], loc=0)
        plt.xlabel('Time [T]', size=20)
        plt.ylabel('Groundwater Head [L]', size=20)
        plt.title('%s' % (self.bore))
        
    def plot_results(self):
        
        plt.figure(figsize=(15,9))
        plt.suptitle('GWTSA %s Model Results' %self._TFN, fontsize=16, fontweight='bold')
        gs = plt.GridSpec(3, 4, wspace=0.2)

        # Plot the recharge
        ax1 = plt.subplot(gs[0,:-1])
        plt.bar(md.num2date(np.arange(self._time_begin, self._time_end+1)), self.recharge, lw=0)
        plt.ylabel('Recharge [m]')
        plt.legend(['Recharge'])
        ax1.xaxis.set_visible(False)
        plt.title('%s' % (self.bore))   
        
        # Plot the Groundwater levels
        ax2 = plt.subplot(gs[1,:-1])
        plt.plot(md.num2date(self._time_axis), self.head_observed, 'k.')
        plt.plot(md.num2date(np.arange(self._time_begin, self._time_end+1)), 
                 self.head_modeled[self._time_model], '-')
        ax2.xaxis.set_visible(False)         
        plt.legend(['Observed Head','Modeled Head'], loc=0)
        plt.ylabel('Groundwater head [m]')
              
        # Plot the residuals and innovations  
        ax3 = plt.subplot(gs[2,:-1])      
        plt.plot(md.num2date(self._time_axis), self.residuals, 'b')
        plt.plot(md.num2date(self._time_axis[0:-1]), self.innovations, 'orange')
        plt.legend(['residuals','innovations'], loc=0)
        plt.ylabel('Error [m]')
        plt.xlabel('Time')                         
        
        # Plot the Impulse Response Function
        ax4 = plt.subplot(gs[0,-1])    
        Fs = self.parameters_opt[0] * gammainc(self.parameters_opt[2], self._time_model/self.parameters_opt[1])
        Fb = Fs[1:] - Fs[0:-1]
        plt.plot(self._time_model[0:-1],Fb)
        plt.xlim(0,np.where(np.cumsum(Fb)>0.999*sum(Fb))[0][0]) # cut off plot after 99.9% of the response
        plt.title('Impulse Response')
        
        # Plot the Model Parameters (Experimental)
        ax5 = plt.subplot(gs[1,-1])
        ax5.xaxis.set_visible(False)
        ax5.yaxis.set_visible(False)
        text = np.vstack((self.Parameter_Names,self.parameters_opt)).T
        colLabels=("Parameter", "Value")
        ax5.table(cellText=text, colLabels=colLabels, loc='center') 
    
        # Table of the numerical diagnostics.
        ax6 = plt.subplot(gs[2,-1])   
        ax6.xaxis.set_visible(False)
        ax6.yaxis.set_visible(False)        
        plt.text(0.1,0.33,'SWSI: %s' %self.swsi_value, fontsize=12 )
        plt.text(0.1,0.66,'Explained variance: %.2f %s' %( self.explained_variance, '%'), fontsize=12 )        
        
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
        plt.figure()
        Fs = self.parameters_opt[0] * gammainc(self.parameters_opt[2], self._time_model/self.parameters_opt[1])
        Fb = Fs[1:] - Fs[0:-1]
        plt.plot(self._time_model[0:-1],Fb)
