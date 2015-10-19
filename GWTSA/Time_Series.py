# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 10:45:17 2015

@author: Raoul
"""
import numpy as np
import matplotlib.dates as md
import matplotlib.pyplot as plt
from scipy.optimize import fmin, leastsq
from scipy.special import gammainc
import TFN_Model
import cma
from statsmodels.tsa.stattools import acf
from lmfit import minimize, Parameters, Parameter, report_fit

#%% Model Class
def latex_plot():
    #from matplotlib import rcParams
    params = {'backend': 'ps',
              #'text.latex.preamble': ['\usepackage{amsmath}','\usepackage[utf8]{inputenc}'],
              #'text.latex.unicode': True,
              'axes.labelsize': 10, 
              'axes.titlesize': 10,
              'font.size': 10, 
              'font.family': 'serif',
              #'font.serif': 'Bookman',
              'legend.fontsize': 10,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              #'text.usetex': 0,
              #'text.dvipnghack' : True,
              'figure.figsize': [8.29,5],
              'figure.dpi': 300
    }
    return plt.rcParams.update(params)

class Model:
    def __init__(self, bores, forcing, timestart=2000):
        self.bores_list = []            # Create an empty list for all the bores instances
        self.bores_number = len(bores)  # Determine how many boreholes are entered
        
        for i in range(self.bores_number):
            self.bores_list.append(TimeSeries(bores[i], forcing, timestart=timestart))
            
    def add_bore(self, bores, forcing):
        self.bores_list.append(TimeSeries(bores, forcing))
        self.bores_number += 1          # Increase number of boreholes
        
    def delete_bore(self, boreid):      #Not functional yet!!!
        self.bores_list.remove(boreid)
        self.bores_number -= 1
    
    def solve(self, TFN, X0, method=0):       
        for i in range(self.bores_number):
            self.bores_list[i].solve(TFN[i], X0, method)
    
    def plot(self, modeled=1):
        fig2 = plt.figure('Boreholes', figsize=(15,9))
        colors=plt.cm.nipy_spectral(np.linspace(0,1,self.bores_number))
        ax = fig2.add_subplot(111)
        ax.set_position([0.05,0.1,0.8,0.8])
        for i in range(self.bores_number):
            ax.plot(md.num2date(self.bores_list[i]._time_axis), self.bores_list[i].head_observed,'.', c=colors[i])
            if modeled == 1:            
                ax.plot(md.num2date(np.arange(self.bores_list[i]._time_begin, self.bores_list[i]._time_end+1)), self.bores_list[i].head_modeled[self.bores_list[i]._time_model], c=colors[i], label='%s model' %self.bores_list[i].bore) 
        plt.ylabel('Groundwater head [m]')
        plt.xlabel('Time [Years]')
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1.014))
        fig2.savefig('boreholes.eps', format='eps', bbox_inches='tight')       
    
#%% Time series class    
            
class TimeSeries:
    def __init__(self, bore, forcing, timestart, rows=[50, 8], cl_unit=10000.0, gw_unit=100.0):
        
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
        
        Ho = np.genfromtxt(('./%s' % (bore)), delimiter=',', skiprows=rows[0], usecols=[2, 5], converters={2: md.strpdate2num('%d-%m-%Y')});
        Ho = Ho[Ho[:,1] > -999] #Select only real values
        self._time_begin = Ho[0,0] #In time number, not date
        self._time_end = Ho[-1,0] #In time number, not date
        self._time_steps = Ho[1:,0] - Ho[0:-1,0] #Calculate the timesteps
        self._time_axis = Ho[:,0] # Real Time in time number format for plots
        self._time_observed = np.append(0, [np.cumsum(self._time_steps)]).astype(int) #Timesteps with observed heads
        self._time_model = np.arange(0,sum(self._time_steps,1), dtype=int)     
        self._time_start = timestart; #warmup time of the model [Days]
        self._spinup = np.where(self._time_observed  > self._time_start)[0][0]

        self.head_observed = Ho[:,1] / gw_unit
               
        # Import the precipitation and evaporation data and calculate the recharge
        ClimateData = np.genfromtxt('./%s' % forcing , delimiter=',', skip_header=rows[1], converters={1: md.datestr2num}, usecols=[1, 2, 3]);
        ClimateData = ClimateData[:,:][(ClimateData[:,0] >= self._time_begin) & (ClimateData[:,0] <= self._time_end)]
        self.precipitation = ClimateData[:,1] / cl_unit 
        self.precipitation[self.precipitation < 0.0] = 0.0
        self.evaporation = ClimateData[:,2] / cl_unit 
        self.evaporation[self.evaporation < 0.0] = 0.0
        self.bore = bore[-17:-8]
        print 'Model setup for bore id %s completed succesfully.' %self.bore

    def __repr__(self):
        return 'Time Series Data of Bore: %s' %self.bore
        
    def solve(self, TFN, X0, method=0):
        
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
            - an array of the optimum parameters set (self.parameters_optimized)
            - array with modeled heads (self.head_modeled)
            - array with innivations (self.innovations)
        
        See Also
        --------
        Help(TFN_Model) can give you information on how the parameter dictionary should look, as
        that depends on the Transfer Functions Noise Model that is used. 
        
        """
        # Define the TFN Model
        self.TFN = getattr(TFN_Model, TFN)
        self._TFN = TFN
    
        InputData = [self._time_model, self.precipitation, self.evaporation]
        
        #self.parameters = initial_parameters
        
        if method == 0:
            X0.add('d', value=np.mean(self.head_observed))
            self.result = minimize(self.objective_function, X0, args=(InputData,), method='leastsq')
            self.parameters_optimized = self.result.params.valuesdict()
            if self.result.success: 
                print 'Optimization completed succesfully!'
                print(report_fit(self.result))
        elif method == 1:  # Support depreciated when lmfit is 
            res = cma.fmin(self.swsi, X0, 2.0, args=(InputData,), options={'ftarget': 1e-5})
            self.parameters_optimized = res[0]
            self.parameters = np.loadtxt('outcmaesxrecentbest.dat', skiprows=2)[:,5:]
            self.correlation_matrix = res[-2].correlation_matrix()
            
        self.simulate(InputData, self._TFN, self.result.params)    
    
    def simulate(self, InputData, TFN, parameters, solver=1):
        """
        Simulate the groundwater levels with a certain parameter set
        """        
        self.swsi_value = self.swsi(parameters, InputData ) #Also calculates a new model? Is this still usefull>?
        
        i = self._time_observed > self._time_start
        self.explained_variance = (np.var(self.head_observed[i])**2 - np.var(self.head_modeled[self._time_observed[i]]-self.head_observed[i])**2)/np.var(self.head_observed[i])**2*100. # Calculate the Explained Variance Percentage
                
        self.rmse = np.sqrt(sum((self.head_modeled[self._time_observed[i]]-self.head_observed[i])**2) / len(self.head_observed[i]))
        self.avg_dev = sum(self.head_modeled[self._time_observed[i]]-self.head_observed[i]) / len(self.head_observed[i])
        
    '''This section contains the objective functions and diagnostic tests.'''

#%% swsi constitutes the adapted version of the Sum of weighted squared innovations (swsi) developed by asmuth et al. (2005). For large values of alpha and small timesteps the numerator approaches zero. Therefore, Peterson et al. (2014) adapted the swsi function, making the numerator a natural log and changing the product operator to a summation.
   
    def swsi(self, parameters, InputData, solver=1):
        alpha = 10.0 ** parameters['alpha'].value
        self.head_modeled, self.recharge = self.TFN(parameters, InputData, solver=solver)
        self.residuals = self.head_observed - self.head_modeled[self._time_observed]
        self.innovations = self.residuals[1:] - (self.residuals[0:-1] * np.exp(-self._time_steps/alpha))
        innovations = self.innovations[self._spinup::1]  # Select only innovations after the spinup period
        N = len(innovations)                        # Number of innovations
        dt = self._time_steps[-N:]                  # 
        numerator = np.exp((1.0/N) * sum(np.log(1 - np.exp(-2.0 / alpha * dt))))
        self.swsi_value = np.sqrt(sum( (numerator / (1 - np.exp(-2.0 / alpha * dt ))) * innovations**2))
        return self.swsi_value 
#       
    def objective_function(self, parameters, InputData, solver=1):
        #self.parameters = np.vstack((self.parameters, parameters))
        alpha = 10.**parameters['alpha'].value
        self.head_modeled, self.recharge = self.TFN(parameters, InputData, solver=solver)
        self.residuals = self.head_observed - self.head_modeled[self._time_observed]
        self.innovations = self.residuals[1:] - (self.residuals[0:-1] * np.exp(-self._time_steps/alpha))
        innovations = self.innovations[self._spinup::1]  # Select only innovations after the spinup period
        N = len(innovations)                        # Number of innovations
        dt = self._time_steps[-N:]                  # 
        numerator = np.exp((1.0/N) * sum(np.log(1 - np.exp(-2.0 / alpha * dt))))
        innovations= (numerator / (1 - np.exp(-2.0 / alpha * dt ))) * innovations**2.
        return innovations        
              
#%% In this section the functions are defined that relate to the plotting of different forcings and results. Each function starts with plot_function to be able to quickly find these modules.        

    def plot_heads(self,color='r',observed=0, modeled=0, newfig=0):
        assert modeled == 0 or observed == 0, 'No heads are chosen to be plotted'
        if newfig == 0:
            plt.figure()
        if observed == 0:
            plt.plot(md.num2date(self._time_axis), self.head_observed, 'k')
        if modeled == 0:
            plt.plot(md.num2date(np.arange(self._time_begin, self._time_end+1)), 
                 self.head_modeled[self._time_model], '-', color=color)
        plt.legend(['Observed Head','Modeled Head'], loc=0)
        plt.xlabel('Time [T]', size=20)
        plt.ylabel('Groundwater Head [L]', size=20)
        plt.title('%s' % (self.bore))
        
    def plot_results(self, ylim=[]):
        
        plt.figure('%s_%s' %(self.bore,self._TFN), figsize=(15,9))
        plt.suptitle('GWTSA %s Model Results' %self._TFN, fontsize=16, fontweight='bold')
        gs = plt.GridSpec(3, 4, wspace=0.2)

        # Plot the recharge
        ax1 = plt.subplot(gs[0,:-1])
        plt.bar(md.num2date(np.arange(self._time_begin, self._time_end+1)), self.recharge, lw=0)
        plt.ylabel('Recharge [m/d]')
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
        plt.axvline(md.num2date(self._time_begin + self._time_start), c='grey', linestyle='--', label='Spinup period')
        #plt.text(md.num2date(self._time_begin + self._time_start), 0.0, 'Spinup period2')  # Not displaying label yet?!
        plt.ylim(min(self.head_observed), max(self.head_observed))
        if ylim == '':
            plt.ylim(ylim[0],ylim[1])
              
        # Plot the residuals and innovations  
        ax3 = plt.subplot(gs[2,:-1])      
        plt.plot(md.num2date(self._time_axis), self.residuals, 'b')
        plt.plot(md.num2date(self._time_axis[0:-1]), self.innovations, 'orange')
        plt.legend(['residuals','innovations'], loc=0)
        plt.ylabel('Error [m]')
        plt.xlabel('Time [Years]')                         
        
        # Plot the Impulse Response Function
        ax4 = plt.subplot(gs[0,-1])    
        A = 10**self.parameters_optimized['A']
        a = 10**self.parameters_optimized['a']
        n = self.parameters_optimized['n']
        Fs = A * gammainc(n, self._time_model/a)
        Fb = Fs[1:] - Fs[0:-1]
        plt.plot(self._time_model[0:-1],Fb)
        plt.xlim(0,np.where(np.cumsum(Fb)>0.99*sum(Fb))[0][0]) # cut off plot after 99.0% of the response
        plt.title('Impulse Response')
        
        # Plot the Model Parameters (Experimental)
        ax5 = plt.subplot(gs[1,-1])
        ax5.xaxis.set_visible(False)
        ax5.yaxis.set_visible(False)
        text = np.vstack((self.parameters_optimized.keys(),self.parameters_optimized.values())).T
        colLabels=("Parameter", "Value")
        ax5.table(cellText=text, colLabels=colLabels, loc='center') 
    
        # Table of the numerical diagnostics.
        ax6 = plt.subplot(gs[2,-1])   
        ax6.xaxis.set_visible(False)
        ax6.yaxis.set_visible(False)        
        plt.text(0.05, 0.80, 'SWSI: %.2f meter' %self.swsi_value, fontsize=12 )
        plt.text(0.05, 0.60, 'Explained variance: %.2f %s' %( self.explained_variance, '%'), fontsize=12 )
        plt.text(0.05, 0.40, 'RMSE: %.2f meter' %self.rmse, fontsize=12)        
        plt.text(0.05, 0.20, 'Average deviation: %.2f meter' %self.avg_dev, fontsize=12)
        
        plt.savefig('%s_%s.eps' %(self.bore,self._TFN), format='eps', bbox_inches='tight')

    def plot_diagnostics(self):
        plt.figure('Diagnostics %s' %self._TFN, figsize=(15,9))
        gs = plt.GridSpec(4,4, wspace=0.4, hspace=0.4)
        plt.suptitle('GWTSA Parameter diagnostics', fontsize=16, fontweight='bold')
        
        # Plot the parameter evolutions
        ax1 = plt.subplot(gs[0:2,0:2])
        #plt.plot(self.parameters)
        plt.title('parameters evolution')
        plt.ylabel('parameter value')
        plt.xlabel('iteration')
        plt.legend(self.result.var_names)
        
        try:
            self.result.covar
            ax2 = plt.subplot(gs[0:2,2:4])
            plt.pcolormesh(self.result.covar, cmap='coolwarm', alpha=0.6)
            plt.xticks(np.arange(0.5,len(self.result.covar)+0.5,1),self.result.var_names, rotation=40, ha='right')
            plt.yticks(np.arange(0.5,len(self.result.covar)+0.5,1),self.result.var_names, rotation=40, ha='right')
            plt.colorbar()
            plt.title('covariance matrix')
            
            for (i, j), z in np.ndenumerate(self.result.covar):
                plt.text(j+0.5, i+0.5, '{:0.2f}'.format(z), ha='center', va='center', color='darkslategrey')
                
        except:
            pass
        
        ax3 = plt.subplot(gs[2:4,0:2])        
        x = acf(self.innovations, nlags=3650)
        plt.title('Autocorrelation graph of the innovations')
        plt.plot(x)
        plt.xlabel('Time [Days]')        
        
        ax4 = plt.subplot(gs[2:4,2:4])  
        ax4.xaxis.set_visible(False)
        ax4.yaxis.set_visible(False)
        text = np.vstack((self.parameters_optimized.keys(),self.parameters_optimized.values())).T
        colLabels=("Parameter", "Value")
        ax4.table(cellText=text, colLabels=colLabels, loc='center')                   
        
    def plot_forcings(self):
        plt.figure()
        plt.bar(md.num2date( self.ClimateData[:,0]), self.ClimateData[:,1], color='b', lw=0)
        plt.ylabel('P [m/d]', size=20)
        plt.ylim(0,max(self.ClimateData[:,1]))        
        ax1 = plt.gca()
        ax2 = ax1.twinx() # To create a second axis on the right
        plt.bar(md.num2date( self.ClimateData[:,0]), self.ClimateData[:,2], color='red', lw=0)
        ax2.set_ylim((0,400)[::-1])
        plt.ylabel('E [m/d]',size=20)
        plt.xlabel('Time [Years]',size=20)
        plt.legend('Precipitation','Potential Evapotranspiration')
        plt.title('Forcings',size=20)
#%%



