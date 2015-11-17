# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 10:45:17 2015

@author: Raoul

Notes To self:
--------------------------------
- Eventually, the package has to deal with multiple model part, wells, transfer functions, error models. This can easily (?) be implemented by using the "Model" approach of lmfit.

"""
import numpy as np
import matplotlib.dates as md
import matplotlib.pyplot as plt
import TFN_Model
from statsmodels.tsa.stattools import acf
from lmfit import minimize, Parameters, Parameter, report_fit, fit_report
from pandas import Series



#%% Model Class

class Model:
    def __init__(self, bores, forcing, calibration, validation):
        self.bores_list = []            # Create an empty list for all the bores instances
        self.bores_number = len(bores)  # Determine how many boreholes are entered
        
        for i in range(self.bores_number):
            self.bores_list.append(TimeSeries(bores[i], forcing, calibration, validation))
            
    def add_bore(self, bores, forcing):
        self.bores_list.append(TimeSeries(bores, forcing))
        self.bores_number += 1          # Increase number of boreholes
        
    def delete_bore(self, boreid):      #Not functional yet!!!
        self.bores_list.remove(boreid)
        self.bores_number -= 1
    
    def solve(self, TFN, X0, method=0, solver=1):       
        for i in range(self.bores_number):
            self.bores_list[i].solve(TFN[i], X0, method, solver)
    
    def plot(self, modeled=1, savefig=False):
        fig2 = plt.figure('Boreholes', figsize=(8.3,5))
        colors=plt.cm.nipy_spectral(np.linspace(0,1,self.bores_number))
        colors[1] = [0.47058823529411764, 0.7686274509803922, 1.0, 1.0]
        for i in range(self.bores_number):
            plt.plot(md.num2date(self.bores_list[i]._time_axis), self.bores_list[i].head_observed,'.', c=colors[i])
            if modeled == 1:            
                plt.plot(md.num2date(np.arange(self.bores_list[i]._time_begin, self.bores_list[i]._time_end+1)), self.bores_list[i].head_modeled, c=colors[i], label='%s, %s' %(self.bores_list[i].bore, self.bores_list[i]._TFN), linestyle='-') 
        plt.ylabel('Groundwater head [m]')
        plt.xlabel('Time [Years]')
        plt.legend(loc=(0,1), ncol=3, frameon=False, handlelength=3)
        plt.axvline(self.bores_list[0]._date_calibration, color='k', linestyle='--')

        if savefig:
            fig2.savefig('Figures/boreholes.eps', format='eps', bbox_inches='tight')       
    
#%% Time series class    
            
class TimeSeries:
    def __init__(self, bore, forcing, calibration, validation, rows=[50, 8], cl_unit=10000.0, gw_unit=100.0):
        
        """
        Prepares the time series model with the forcings and observed heads
        
        parameters
        ----------
        bore: string
            name of the csv file containing the observed head values. 
        forcing: string
            name of the csv file containing the forcing values. 
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
        
        Ho = np.genfromtxt(('./%s' % (bore)), delimiter=',', skip_header=rows[0], usecols=[2, 5], converters={2: md.strpdate2num('%d-%m-%Y')})
        Ho = Ho[Ho[:,1] > -999] #Select only real values
        Ho = Ho[(Ho[:,0] > md.datestr2num(calibration[0])) & (Ho[:,0] < md.datestr2num(validation[1]))] # Select all observed GWL in the calibration and validation period
        
        self.head_observed = Ho[:,1] / gw_unit # Observed groundwater levels [L]
        self._time_begin = Ho[0,0] # In time number format
        self._time_end = Ho[-1,0] # In time number format      
        self._time_steps = Ho[1:,0] - Ho[0:-1,0] # Calculate the timesteps [Days]  
        self._time_axis = Ho[:,0] # Real Time in time number format for plots

        # Import the precipitation and evaporation data and calculate the recharge
        ClimateData = np.genfromtxt('./%s' % forcing , delimiter=',', skip_header=rows[1], converters={1: md.datestr2num}, usecols=[1, 2, 3]);
        ClimateData = ClimateData[(ClimateData[:,0] >= md.datestr2num('01-01-1960')) & (ClimateData[:,0] <= self._time_end)]
        
        self.precipitation = ClimateData[:,1] / cl_unit 
        self.precipitation[self.precipitation < 0.0] = 0.0
        self.evaporation = ClimateData[:,2] / cl_unit 
        self.evaporation[self.evaporation < 0.0] = 0.0
        
        self._time_climate = md.num2date(ClimateData[:,0]) # Time for climate series in real time [dd-mm-yyyy]
        
        # Times with reference to model start time
        self._time_model = np.arange(0,(ClimateData[-1,0]-ClimateData[0,0])+1, dtype=int) # Model time [Days]
        self._time_spinup = self._time_begin - ClimateData[0,0]


        self._index_observed = np.append(0, [np.cumsum(self._time_steps)-1]).astype(int) # Indexes with observed heads        
        self._index_calibration = self._index_observed[1:] < (md.datestr2num(calibration[1])-self._time_begin)
        self._index_validation = self._index_observed[1:] > (md.datestr2num(calibration[1])-self._time_begin)
        self._date_calibration = calibration[1]
        
        # Save a name for the borehole
        self.bore = bore[-17:-8]
        print 'Model setup for bore id %s completed succesfully.' %self.bore

    def __repr__(self):
        return 'Time Series Data of Bore: %s' %self.bore
        
    def solve(self, TFN, X0, method='leastsq', solver=1):
        
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
            0 = Levenberg-Marquardt method, 1 =  Nelder-Mead method
        solver: integer 0 or 1
            0 = explicit euler, 1 = implicit euler
            
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
    
        InputData = [self._time_model, self.precipitation, self.evaporation, solver]
        
        if method == 'leastsq':
            X0.add('d', value=np.mean(self.head_observed), vary=True) 
            self.result = minimize(self.objective_function, X0, args=(InputData,), method='leastsq')
            self.parameters_optimized = self.result.params.valuesdict()
            if self.result.success: 
                print 'Optimization completed succesfully!'
                print(report_fit(self.result))
                #np.savetxt('fit_report_%s_%s.txt' %(self.bore, self._TFN),(fit_report(self.result),), fmt='%str')
        else:
            X0.add('d', value=np.mean(self.head_observed))
            self.result = minimize(self.objective_function, X0, args=(InputData,), method=method)
            self.parameters_optimized = self.result.params.valuesdict()
         
        # Calculate statistics for both the calibration and validation period 
        self.SWSI_Cal = self.swsi(self._index_calibration)
        self.SWSI_Val = self.swsi(self._index_validation)
        
        self.RMSE_Cal = self.rmse(self._index_calibration)
        self.RMSE_Val = self.rmse(self._index_validation)
        
        self.AVGDEV_Cal = self.avg_deviation(self._index_calibration)
        self.AVGDEV_Val = self.avg_deviation(self._index_validation)      
        
        self.EXPVAR_Cal = self.explained_variance(self._index_calibration)
        self.EXPVAR_Val = self.explained_variance(self._index_validation)  
    
        
    '''This section contains the objective functions and diagnostic tests.'''

#%% swsi constitutes the adapted version of the Sum of weighted squared innovations (swsi) developed by asmuth et al. (2005). For large values of alpha and small timesteps the numerator approaches zero. Therefore, Peterson et al. (2014) adapted the swsi function, making the numerator a natural log and changing the product operator to a summation.
      
    def objective_function(self, parameters, InputData):
        alpha = 10.0**parameters['alpha'].value
        self.head_modeled, recharge = self.TFN(parameters, InputData)
        self.recharge = Series(recharge, index=self._time_climate)
        self.recharge = self.recharge[self.recharge.index > md.num2date(self._time_begin)]
        self.head_modeled = self.head_modeled[self._time_spinup:self._time_model[-1]+1] #Select entire period
        self.residuals = self.head_observed - self.head_modeled[self._index_observed] 
        self.innovations = self.residuals[1:] - (self.residuals[0:-1] * np.exp(-self._time_steps/alpha))
        
        # Select the period for which to calibrate
        innovations = self.innovations[self._index_calibration]
        dt = self._time_steps[self._index_calibration]    
        N = len(innovations) # Number of innovations
        
        # Weighing the innovations for optimization
        numerator = np.exp((1.0/N) * sum(np.log(1 - np.exp(-2.0 / alpha * dt))))
        innovations = (numerator / (1 - np.exp(-2.0 / alpha * dt ))) * innovations**2.
        return innovations      

    def swsi(self, period):
        alpha = 10.0**self.parameters_optimized['alpha']
        # Select the period to calculate statistics
        innovations = self.innovations[period]
        dt = self._time_steps[period]    
        N = len(innovations) # Number of innovations
        
        numerator = np.exp((1.0/N) * sum(np.log(1 - np.exp(-2.0 / alpha * dt))))
        return np.sqrt(sum( (numerator / (1 - np.exp(-2.0 / alpha * dt ))) * innovations**2)) 

    def rmse(self, period):
        return np.sqrt(sum((self.head_modeled[self._index_observed][period]-self.head_observed[period])**2) / len(self.head_observed[period]))
      
    def avg_deviation(self, period):
        return sum(self.head_modeled[self._index_observed][period]-self.head_observed[period]) / len(self.head_observed[period])

    def explained_variance(self, period):
        return (np.var(self.head_observed[period])**2 - np.var(self.residuals[period])**2)/np.var(self.head_observed[period])**2*100. 
             
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
        
    def plot_results(self, savefig=True):      
        plt.figure('%s_%s' %(self.bore,self._TFN))
        plt.suptitle('GWTSA %s Model Results for bore %s' %(self._TFN, self.bore), fontweight='bold')
        gs = plt.GridSpec(3, 4, wspace=0.2)

        # Plot the recharge
        ax1 = plt.subplot(gs[0,:-1])
        self.recharge.resample('A', how='sum').plot('bar')
        plt.ylabel('Recharge [m/year]')
        plt.text(0.1,0.03, 'average %.2f m/y' %np.mean(self.recharge.resample('A', how='sum')), backgroundcolor='w')
        ax1.xaxis.set_visible(False)
        
       # Plot the Groundwater levels
        ax2 = plt.subplot(gs[1,:-1])
        plt.plot(md.num2date(self._time_axis), self.head_observed, 'k.')
        plt.plot(md.num2date(np.arange(self._time_begin, self._time_end+1)), 
                 self.head_modeled, '-')
        ax2.xaxis.set_visible(False)         
        plt.legend(['Observed Head','Modeled Head'], loc=2)
        plt.ylabel('Gwl [m]')
        plt.ylim(min(self.head_observed), max(self.head_observed))
        plt.axvline(self._date_calibration, color='k', linestyle='--')
        ymin, ymax = plt.ylim()
        plt.text(self._date_calibration,ymin+0.1, 'validation', verticalalignment = 'bottom')
              
        # Plot the residuals and innovations  
        ax3 = plt.subplot(gs[2,:-1])      
        plt.plot(md.num2date(self._time_axis), self.residuals, 'b')
        plt.plot(md.num2date(self._time_axis[0:-1]), self.innovations, 'orange')
        plt.legend(['residuals','innovations'], loc=2)
        plt.ylabel('Error [m]')
        plt.xlabel('Time [Years]')                        
        
        # Plot the Impulse Response Function
        ax4 = plt.subplot(gs[0,-1])    
        IRF = getattr(TFN_Model, 'IRF')
        Fb = IRF(self.result.params)
        plt.plot(Fb)
        plt.xticks(range(0,10000,500))
        plt.xlim(0,np.where(np.cumsum(Fb)>0.99*sum(Fb))[0][0]) # cut off plot after 99.0% of the response
        plt.text(5,0.0,'Peak Time: %i' %Fb.argmax(), verticalalignment='bottom')
        plt.title('Impulse Response')
        
        # Plot the Model Parameters (Experimental)
        ax5 = plt.subplot(gs[1,-1])
        ax5.xaxis.set_visible(False)
        ax5.yaxis.set_visible(False)
        text = np.vstack((self.parameters_optimized.keys(),[round(float(i), 3) for i in self.parameters_optimized.values()])).T
        colLabels=("Parameter", "Value")
        ytable = plt.table(cellText=text, colLabels=colLabels, loc='center')
        ytable.scale(1,0.5)
    
        # Table of the numerical diagnostics.
        ax6 = plt.subplot(gs[2,-1])   
        ax6.xaxis.set_visible(False)
        ax6.yaxis.set_visible(False)        
        plt.text(0.05, 0.84, 'SWSI: %.2f / %.2f meter' %(self.SWSI_Cal, self.SWSI_Val))
        plt.text(0.05, 0.68, 'Expl. var: %.2f / %.2f pr' %(self.EXPVAR_Cal, self.EXPVAR_Val))
        plt.text(0.05, 0.52, 'RMSE: %.2f / %.2f meter' %(self.RMSE_Cal, self.RMSE_Val))      
        plt.text(0.05, 0.36, 'Avg dev: %.2f / %.2f meter' %(self.AVGDEV_Cal, self.AVGDEV_Val))
        plt.text(0.05, 0.20, 'AIC: %.2f' %self.result.aic)
        plt.text(0.05, 0.04, 'BIC: %.2f' %self.result.bic)
        
        if savefig:
            plt.savefig('Figures/%s_%s.eps' %(self.bore,self._TFN), format='eps', bbox_inches='tight')

    def plot_diagnostics(self, savefig=True):
        plt.figure('Diagnostics %s_%s' %(self.bore,self._TFN))
        plt.suptitle('GWTSA %s Model Diagnostics for bore %s' %(self._TFN, self.bore), fontweight='bold')
        gs = plt.GridSpec(4,4, wspace=0.4, hspace=0.4)
        
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
        text = np.vstack((self.parameters_optimized.keys(),[round(float(i), 3) for i in self.parameters_optimized.values()])).T
        colLabels=("Parameter", "Value")
        ax4.table(cellText=text, colLabels=colLabels, loc='center', fontsize=4)   

        if savefig:
            plt.savefig('%s_%s_diagnostics.eps' %(self.bore,self._TFN), format='eps', bbox_inches='tight')                
        
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
#%% Estiamte Sumax
    def calcSumax(self, Cr=0.40, imax=0.001, T=20, EaMethod='gradual'):
        from scipy.stats import gumbel_r
        # Interception Reservoir
        n = len(self.evaporation)
        Si = np.zeros(n)
        Pe = np.zeros(n)
        Ei = np.zeros(n)
        Ep = np.zeros(n)
        
        for t in range(n-1):    
            Si[t+1] = Si[t]+self.precipitation[t+1]                      # Fill interception bucket with new rain
            Pe[t+1] = np.max(((Si[t+1]-imax), 0.0))     # Calculate effective precipitation
            Si[t+1] = Si[t+1] - Pe[t+1]                 # Update interception state
            Ei[t+1] = np.min((Si[t+1], self.evaporation[t+1]))         # Evaporation from interception
            Si[t+1] = Si[t+1] - Ei[t+1]                 # Update interception state
            Ep[t+1] = self.evaporation[t+1] - Ei[t+1]                  # Update potential evapotranspiration   
        
        # Estimate Actual Evaporation
        
        # Averag actual evaporation [L/T]
        Ea = np.min((sum((1.-Cr)*Pe), sum(Ep))) / n  # 1-Cr because Q = Cr * Pe
        #Get seasonality back in there
        EaAvgC = np.min((sum((1.-Cr)*Pe), sum(Ep))) / n  # 1-Cr because Q = Cr * Pe
        
        #Get seasonality back in there
        if EaMethod == 'gradual':
            EaAvg = np.ones(len(Ep)) * EaAvgC           # get timeseries with constant average Ea
            a = np.min((Ep,EaAvg), axis=0)
            A = sum(EaAvg - a)
            B = sum(Ep - a)
            Ea = A / B * (Ep - a) + a                   # Calculate Ea with similar pattern as Ep
            
        soildeficit = np.zeros(n)   
        for t in range(n-1):
            soildeficit[t+1] = np.min(((soildeficit[t] + Pe[t] - Ea[t]), 0.0))     
        
        soildeficit = Series(soildeficit, index=self._time_climate)
        Sumax = np.sqrt(soildeficit.resample('A', how='min') ** 2. ) # Work with positive values
    
        Sumax.sort()
        mu, sigma = gumbel_r.fit(Sumax)
        y = gumbel_r.pdf(Sumax, mu, sigma)
        plt.plot(Sumax,y)
        Sumax.hist()
        plt.axvline(gumbel_r.isf(1./T, loc=mu, scale=sigma)) #Causes trouble for multiple return periods
    
        return gumbel_r.isf(1./T, loc=mu, scale=sigma) 

#%% Make a nicely looking latex plot! Just run this function before plotting :)

def latex_plot():
    #from matplotlib import rcParams
    params = {'backend': 'ps',
              #'text.latex.preamble': ['\usepackage{amsmath}','\usepackage[utf8]{inputenc}'],
              #'text.latex.unicode': True,
              'axes.labelsize': 8, 
              'axes.titlesize': 8,
              'font.size': 6, 
              'font.family': 'serif',
              #'font.serif': 'Bookman',
              'legend.fontsize': 8,
              'xtick.labelsize': 6,
              'ytick.labelsize': 6,
              #'text.usetex': 0,
              #'text.dvipnghack' : True,
              'figure.figsize': [8.29,5],
              'figure.dpi': 150
    }
    return plt.rcParams.update(params)

