# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 10:45:17 2015

@author: Raoul

Notes To self:
--------------------------------


"""
import numpy as np
import matplotlib.dates as md
import matplotlib.pyplot as plt
from TFN_Model import *
from statsmodels.tsa.stattools import acf
from lmfit import minimize, Parameters, Parameter, report_fit, fit_report
from pandas import Series, DataFrame, read_csv

cyan = [120/255.,196./255,1.]
#%% Model Class

class Model:    
    """ The Model class creates instances for all observations well.
    Parameters
    ----------
    bores, forcing: string
        The filename and path provided as a string
    calibration, validation: list of strings
        list that contains the start and end date 
        e.g. ['01-01-1974', '31-12-1994']
    discharge: string, optional
        keyword-argument, possibility to add discharge from a well
    Cf: list, optional
        Cropfactor, list of one or two values. With two values a linearly 
        changing cropfactor is applied over the entire period.
    
    Returns
    -------
    Object M that contains a list (M.bores_list) with an instance for each 
    borehole entered.
    
    Example
    -------
    M = Model('bores.csv', 'forcing.csv', ['01-01-1974', '31-12-1994'], 
               ['01-01-1974', '31-12-2003'], Discharge.csv', [1.2,1.0])
    
    See Also
    --------
    See the documentation of GWTSA for a description and examples of how the 
    input files for the observed heads and the forcings should look.
    """
    
    def __init__(self, bores, forcing, calibration, validation, discharge=None, 
                 Cf=[1.0]):
        self.bores_list = []            # Create an empty list for all the bores instances
        self.bores_number = len(bores)  # Determine how many boreholes are entered
        
        for i in range(self.bores_number):
            self.bores_list.append(TimeSeries(bores[i], forcing, calibration, 
                                   validation, discharge, Cf=Cf))
            
    def add_bore(self, bores, forcing, calibration, validation, discharge=None, 
                 Cf=[1.0]):
        """ Add a single borehole to the bores_list. Input is similar to that 
        of the Model class.
        """
        self.bores_list.append(TimeSeries(bores, forcing, calibration, 
                               validation, discharge, Cf))
        self.bores_number += 1          # Increase number of boreholes
        
    def delete_bore(self, boreid):      #Not functional yet!!!
        self.bores_list.pop(boreid)
        self.bores_number -= 1
    
    def solve(self, X0, IR='IRF', RM='linear', trend=None, method='leastsq', 
              solver=1):
        """ The solve function fits a model on each of the boreholes in the 
        bores_list.
        Parameters
        ----------
        X0: Parameters object
            Parameters object from the LMFIT package. See help(Parameters) 
            and the GWTSA manual on how to construct this input
        IR: string, optional
            The Impulse Response (IR) function to apply. e.g. 'IRF'.
        RM: list, optional
            list of strings with the Recharge Model (RM) to use.
        trend: string, optional
            string-name of the trend to use. options are: 'reclamation', 
            'linear_slope', 'well' or a another self-defined function.
        method: string, optional
            optimization method used by LMFIT. standard is 'leastsq'.
        solver: int 0 or 1
            choose how to solve the recharge model. 0 for explicit Euler and
            1 for implicit Euler.
        
        Returns
        -------
        An optimized model that is stored in the instance object for each 
        borehole.

        """
        for i in range(self.bores_number):
            self.bores_list[i].solve(X0, IR, RM[i], trend, method, solver)
    
    def plot(self, modeled=1, savefig=False):
        """ Plot the observed and modelled groundwater levels of all boreholes 
        in the bores_list. keyword-argument: 'modeled', True or False to plot
        the modelled gwl's. Can also be used when no model is fitted. 'savefig'
        is used to save the figure or not.
        """
        fig2 = plt.figure('Boreholes')
        colors=plt.cm.nipy_spectral(np.linspace(0,1,self.bores_number))
        if self.bores_number>1: colors[1] = [0.47058823529411764, 0.7686274509803922, 1.0, 1.0] #Add cyan TUDelft color
        for i in range(self.bores_number):
            plt.plot(md.num2date(self.bores_list[i]._time_axis), self.bores_list[i].head_observed,'.', markersize=2, c=colors[i], label='%s observed' %self.bores_list[i].bore)
            if modeled == 1:            
                plt.plot(md.num2date(np.arange(self.bores_list[i]._time_begin, self.bores_list[i]._time_end+1)), self.bores_list[i].head_modeled, c=colors[i], label='%s, %s' %(self.bores_list[i].bore, self.bores_list[i]._RM), linestyle='-') 
        plt.ylabel('Groundwater head [m]')
        plt.xlabel('Time [Years]')
        plt.legend(loc=(0,1), ncol=3, frameon=False, handlelength=3)
        plt.axvline(self.bores_list[0]._date_calibration, color='k', linestyle='--')

        if savefig:
            fig2.savefig('Figures/boreholes.eps', bbox_inches='tight')
    
#%% Time series class    
            
class TimeSeries:
    def __init__(self, bore, forcing, calibration, validation, Discharge=None, rows=[50, 8], cl_unit=10000.0, gw_unit=100.0, Cf=[1.0]):     
        """ The Model class creates instances for all observations well.
        Parameters
        ----------
        bores, forcing: string
            The filename and path provided as a string
        calibration, validation: list of strings
            list that contains the start and end date 
            e.g. ['01-01-1974', '31-12-1994']
        discharge: string, optional
            keyword-argument, possibility to add discharge from a well
        Cf: list, optional
            Cropfactor, list of one or two values. With two values a linearly 
            changing cropfactor is applied over the entire period.
        
        Returns
        -------
        Object M that contains a list (M.bores_list) with an instance for each 
        borehole entered.
        
        Example
        -------
        M = Model('bores.csv', 'forcing.csv', ['01-01-1974', '31-12-1994'], 
                   ['01-01-1974', '31-12-2003'], Discharge.csv', [1.2,1.0])
        
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
        ClimateData = ClimateData[(ClimateData[:,0] <= self._time_end) & (ClimateData[:,0] > -999) & (ClimateData[:,2] >-999)] # Get climate series untill Time_end
        
        self.precipitation = ClimateData[:,1] / cl_unit *1.0
        self.precipitation[self.precipitation < 0.0] = 0.0
        if len(Cf)==2:
            self.Cf = np.linspace(Cf[0],Cf[1], len(ClimateData[:,2]))
        else:
            self.Cf = np.ones(len(ClimateData[:,2]))* Cf
        self.evaporation = self.Cf * ClimateData[:,2] / cl_unit # 0.65 is the cropfactor
        self.evaporation[self.evaporation < 0.0] = 0.0
        
        self._time_climate = md.num2date(ClimateData[:,0]) # Time for climate series in real time [dd-mm-yyyy]
        
        # Times with reference to model start time
        self._time_model = np.arange(0,(ClimateData[-1,0]-ClimateData[0,0])+1, dtype=int) # Model time [Days]
        self._time_spinup = self._time_begin - ClimateData[0,0]

        self._index_observed = np.append(0, [np.cumsum(self._time_steps)-1]).astype(int) # Indexes with observed heads        
        self._index_calibration = self._index_observed < (md.datestr2num(calibration[1])-self._time_begin)
        self._index_validation = self._index_observed > (md.datestr2num(calibration[1])-self._time_begin)
        self._date_calibration = calibration[1]
        
        # Figure out the discharge per day and for the correct period
        if Discharge != None:         
            data = read_csv('%s' %Discharge, delimiter=';', index_col=0, parse_dates=True)
            data = data.resample('1d', fill_method='ffill') # Downsample discharge to daily time step if necessary
            data = data[(data.index > md.num2date(ClimateData[0,0])) & (data.index < md.num2date(self._time_end))]
            i = md.date2num(data.index[0])-ClimateData[0,0]
            self.Discharge = np.append(np.zeros(i), data.values)
        
        # Save a name for the borehole
        self.bore = bore[-17:-8]
        print 'Model setup for bore id %s completed succesfully.' %self.bore

    def __repr__(self):
        return 'Time Series Data of Bore: %s' %self.bore
        
    def solve(self, X0, IR='IRF', RM='linear', trend=None, method='leastsq', solver=1):
        """ The solve function fits a model.
        Parameters
        ----------
        X0: Parameters object
            Parameters object from the LMFIT package. See help(Parameters) 
            and the GWTSA manual on how to construct this input
        IR: string, optional
            The Impulse Response (IR) function to apply. e.g. 'IRF'.
        RM: list, optional
            list of strings with the Recharge Model (RM) to use.
        trend: string, optional
            string-name of the trend to use. options are: 'reclamation', 
            'linear_slope', 'well' or a another self-defined function.
        method: string, optional
            optimization method used by LMFIT. standard is 'leastsq'.
        solver: int 0 or 1
            choose how to solve the recharge model. 0 for explicit Euler and
            1 for implicit Euler.
        
        Returns
        -------
        An optimized model that is stored in the instance object for each 
        borehole. See plot_result and plot_diagnostics for more info on
        showing the results.
        """        

        # Define the TFN Model
        self.TFN = eval(IR)     # Save the impulse response function
        self.RM = eval(RM)      # Save the recharge calculation function
        self._TFN = IR          # Save the name of the impulse response function
        self._RM = RM           # Save the name of the recharge model

        InputData = [self._time_model, self.precipitation, self.evaporation,
                     solver, IR, RM, trend]
        
        if trend == 'well':
            
            InputData.append(self.Discharge)
        
        if method == 'leastsq':
            X0.add('d', value=np.mean(self.head_observed), vary=True)
            if trend=='reclamation': X0.add('t_start', value=(md.datestr2num('01-01-1975') - md.date2num(self._time_climate[0])), vary=False)
            self.result = minimize(self.objective_function, X0, args=(InputData,), method='leastsq', scale_covar=True)
            self.parameters_optimized = self.result.params.valuesdict()
            if self.result.success: 
                print 'Optimization completed succesfully!'
                print(report_fit(self.result))
                np.savetxt('Figures/fit_report_%s_%s.txt' % (self.bore, self._TFN), (fit_report(self.result),), fmt='%str')
        else:
            X0.add('d', value=np.mean(self.head_observed))
            self.result = minimize(self.objective_function, X0,
                                   args=(InputData,), method=method)
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
        alpha = parameters['alpha'].value
        output = construct_model(parameters, InputData)
        #output = combi_model(parameters, InputData)        
        self.head_modeled = output[1][self._time_spinup:self._time_model[-1]+1] #Select entire period        
        self.recharge = Series(output[2], index=self._time_climate)
        self.recharge = self.recharge[self.recharge.index > md.num2date(self._time_begin)]
        if output[0]==1:
            self.trend = output[3][self._time_spinup:self._time_model[-1]+1]
        self.residuals = self.head_observed - self.head_modeled[self._index_observed] 
        self.innovations = self.residuals[1:] - (self.residuals[0:-1] * np.exp(-self._time_steps/alpha))
        
        # Select the period for which to calibrate
        innovations = self.innovations[self._index_calibration[:-1]]
        dt = self._time_steps[self._index_calibration[:-1]]    
        N = len(innovations) # Number of innovations
        
        # Weighing the innovations for optimization
        numerator = np.exp((1.0/N) * sum(np.log(1 - np.exp(-2.0 / alpha * dt))))
        innovations = (numerator / (1 - np.exp(-2.0 / alpha * dt ))) * innovations**2.
        return innovations      

    def swsi(self, period):
        alpha = self.parameters_optimized['alpha']
        # Select the period to calculate statistics
        innovations = self.innovations[period[:-1]]
        dt = self._time_steps[period[:-1]]    
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
        
    def plot_results(self, savefig=True):      
        plt.figure('%s_%s_%s' %(self.bore,self._TFN,self._RM))
        plt.suptitle('GWTSA %s Model Results for bore %s' %(self._TFN, self.bore), fontweight='bold')
        gs = plt.GridSpec(3, 4, wspace=0.2)
 
       # Plot the Groundwater levels
        ax2 = plt.subplot(gs[1,:-1])
        plt.plot(md.num2date(self._time_axis), self.head_observed, 'k.', markersize=2, label='observed head')
        plt.plot(md.num2date(np.arange(self._time_begin, self._time_end+1)), 
                 self.head_modeled, '-', color=cyan, label='modeled head')
        try:
            plt.plot(md.num2date(np.arange(self._time_begin, self._time_end+1)),self.trend+np.mean(self.head_observed))  
        except:
            pass
        #std = np.std(self.residuals)
        #plt.fill_between(md.num2date(np.arange(self._time_begin, self._time_end+1)),self.head_modeled-2*std, self.head_modeled+2*std, color='gray', alpha=0.3, label='95% confidence interval')
        ax2.xaxis.set_visible(False)         
        plt.legend(loc=(0,1), ncol=3, frameon=False, handlelength=3)
        plt.ylabel('Head [m]')
        #plt.ylim(min(self.head_modeled-2.2*std), max(self.head_modeled+2.2*std))
        plt.axvline(self._date_calibration, color='gray', linestyle=':')
        ymin, ymax = plt.ylim()
        #plt.text(self._date_calibration+2,ymin+0.1, 'Validation', verticalalignment = 'bottom')
           
        # Plot the recharge
        ax1 = plt.subplot(gs[0,:-1])
        self.recharge.resample('A', how='sum').plot('bar', color=cyan, linewidth=0)
        ax1.grid(0)
        #ax1.bar(R.index.to_pydatetime(), R.values, lw=1, color=[120/255.,196./255,1.], width=-250.)
        plt.ylabel('Recharge [m/year]')
        plt.text(0.1,0.1, 'average %.2f m/y' %np.mean(self.recharge.resample('A', how='sum')), backgroundcolor='w', alpha=0.5)
        ax1.xaxis.set_visible(False)           
           
        # Plot the residuals and innovations  
        ax3 = plt.subplot(gs[2,:-1], sharex=ax2)      
        plt.plot(md.num2date(self._time_axis), self.residuals, 'k', label='residuals')
        plt.plot(md.num2date(self._time_axis[0:-1]), self.innovations, color=cyan, label='innovations')
        plt.legend(loc=(0,1), ncol=3, frameon=False, handlelength=3)
        plt.ylabel('Error [m]')
        plt.xlabel('Time [Years]')                        
        
        # Plot the Impulse Response Function
        ax4 = plt.subplot(gs[0,-1])    
        Fb = self.TFN(self.result.params)
        plt.plot(Fb, 'k')
        plt.xticks(range(0,10000,1000))
        plt.xlim(0,np.where(np.cumsum(Fb)>0.99*sum(Fb))[0][0]) # cut off plot after 99.0% of the response
        plt.ylim(0.0)
        plt.text(10,0.1,'Peak Time: %i' %Fb.argmax(), verticalalignment='bottom')
        plt.title('Impulse Response')
        
        # Plot the Model Parameters (Experimental)
        ax5 = plt.subplot(gs[1,-1])
        ax5.xaxis.set_visible(False)
        ax5.yaxis.set_visible(False)
        text = np.vstack((self.parameters_optimized.keys(),[round(float(i), 4) for i in self.parameters_optimized.values()])).T
        colLabels=("Parameter", "Value")
        ytable = ax5.table(cellText=text, colLabels=colLabels, loc='center')
        ytable.scale(1,0.6)

        # Table of the numerical diagnostics.
        ax6 = plt.subplot(gs[2,-1])   
        ax6.xaxis.set_visible(False)
        ax6.yaxis.set_visible(False)        
        plt.text(0.05, 0.84, 'SWSI: %.3f / %.3f meter' %(self.SWSI_Cal, self.SWSI_Val))
        plt.text(0.05, 0.68, r'Expl. var: %.2f / %.2f %%' %(self.EXPVAR_Cal, self.EXPVAR_Val))
        plt.text(0.05, 0.52, 'RMSE: %.3f / %.3f meter' %(self.RMSE_Cal, self.RMSE_Val))      
        plt.text(0.05, 0.36, 'Avg dev: %.3f / %.3f meter' %(self.AVGDEV_Cal, self.AVGDEV_Val))
        plt.text(0.05, 0.20, 'AIC: %.2f' %self.result.aic)
        plt.text(0.05, 0.04, 'BIC: %.2f' %self.result.bic)
        
        if savefig:
            plt.savefig('Figures/%s_%s_%s.eps' %(self.bore,self._TFN,self._RM), bbox_inches='tight')

    def plot_diagnostics(self, savefig=True):
        plt.figure('Diagnostics %s_%s' %(self.bore,self._TFN))
        plt.suptitle('GWTSA %s Model Diagnostics for bore %s' %(self._TFN, self.bore), fontweight='bold')
        gs = plt.GridSpec(4,4, wspace=0.4, hspace=0.6)
        
        # Plot the parameter evolutions
        ax1 = plt.subplot(gs[0:2,0:2])
        Hm = self.head_modeled[self._index_observed]
        Hm = np.sort(Hm)
        FOE = (1-np.arange(0.0,len(Hm),1) / len(Hm))*100
        plt.plot(FOE, Hm, 'x', color='k')

        Ho = self.head_observed
        Ho = np.sort(Ho)
        FOE = (1-np.arange(0.0,len(Ho),1) / len(Ho))*100
        plt.plot(FOE, Ho, 'k-')

        plt.title('Frequency of Exceedence graph')
        plt.ylabel('Head [m]')
        plt.xlabel('frequency of exceedence ($\%$)')
        plt.legend(['Modeled', 'Observed'])
        
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
            pass #This should give some message that the covariance matrix is unavailable
        
        ax3 = plt.subplot(gs[2:4,0:2])        
        x = acf(self.innovations, nlags=730)
        plt.title('Autocorrelation graph of the innovations')
        plt.stem(x, color='k', markerfmt=' ' )
        plt.xlabel('Time [Days]')        
        
        ax4 = plt.subplot(gs[2:4,2:4])  
        ax4.xaxis.set_visible(False)
        ax4.yaxis.set_visible(False)
        text = np.vstack((self.parameters_optimized.keys(),[round(float(i), 3) for i in self.parameters_optimized.values()])).T
        colLabels=("Parameter", "Value")
        ax4.table(cellText=text, colLabels=colLabels, loc='center', fontsize=4)   

        if savefig:
            plt.savefig('Figures/%s_%s_diagnostics.eps' %(self.bore,self._TFN), format='eps', bbox_inches='tight')                

#%% Plot the recharge
    def recharge_uncertainty(self, n=1000, fig=1):
        keys = self.result.var_names
        par = [self.result.params.valuesdict()[x] for x in keys]
        par = np.random.multivariate_normal(par, self.result.covar, 1000)
        InputData = [self._time_model, self.precipitation, self.evaporation, 1]
        df = DataFrame(index=self._time_climate)
        for i in range(n):
            parameters = self.result.params
            for j in range(len(keys)): parameters['%s' %keys[j]].value = par[i,j]
            df['r%s'%i] = self.RM(parameters, InputData)
        df = df[df.index > md.num2date(self._time_begin)]
        df = df.resample('A', how='sum')
        self.df = df
        df = df[(df>-5.0) & (df<5.0)]
#        self.recharge_std = df.std(1)
#        self.recharge_mean = df.mean(1)
         
        if fig==1:
            return plt.bar(df.index, df.quantile(q=0.5, axis=1), yerr=[(df.quantile(q=0.5, axis=1)-df.quantile(q=0.025, axis=1)),(df.quantile(q=0.975, axis=1)-df.quantile(q=0.5, axis=1))], color=[120/255.,196./255,1.], width = -250, lw=0, error_kw={'ecolor': 'dimgray', 'ewidth': '5'})

#%% Estimate Sumax     
    def calcSumax(self, Cr=0.40, imax=0.001, T=20, fc=1.0, EaMethod='gradual'):
        from scipy.stats import gumbel_r
        # Interception Reservoir
        n = len(self.evaporation)
        Si = np.zeros(n)
        Pe = np.zeros(n)
        Ei = np.zeros(n)
        Ep = np.zeros(n)
        evaporation = fc * self.evaporation
        for t in range(n-1):    
            Si[t+1] = Si[t]+self.precipitation[t+1]                      # Fill interception bucket with new rain
            Pe[t+1] = np.max(((Si[t+1]-imax), 0.0))     # Calculate effective precipitation
            Si[t+1] = Si[t+1] - Pe[t+1]                 # Update interception state
            Ei[t+1] = np.min((Si[t+1], evaporation[t+1]))         # Evaporation from interception
            Si[t+1] = Si[t+1] - Ei[t+1]                 # Update interception state
            Ep[t+1] = evaporation[t+1] - Ei[t+1]                  # Update potential evapotranspiration   
        
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
        
        self.soildeficit = Series(soildeficit, index=self._time_climate)
        Sumax = np.sqrt(self.soildeficit.resample('A', how='min') ** 2. ) # Work with positive values
    
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
              'legend.scatterpoints': 3,
              'xtick.labelsize': 6,
              'xtick.color': 'gray',     
              'ytick.labelsize': 6,
              'ytick.color': 'gray',              
              #'text.usetex': 0,
              #'text.dvipnghack' : True,
              'figure.figsize': [8.29,5],
              'figure.dpi': 300,
              'figure.facecolor' : 'white',
              'axes.edgecolor': 'lightgray',
              'axes.facecolor': 'white',
              'axes.labelcolor': 'dimgray',
    }
    return plt.rcParams.update(params)

def presentation_plot():
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
              'legend.scatterpoints': 3,
              'xtick.labelsize': 10,
              'xtick.color': 'k',     
              'ytick.labelsize': 10,
              'ytick.color': 'k',              
              #'text.usetex': 0,
              #'text.dvipnghack' : True,
              'figure.figsize': [8.29,5],
              'figure.dpi': 300,
              'figure.facecolor' : 'white',
              'axes.edgecolor': 'k',
              'axes.facecolor': 'white',
              'axes.labelcolor': 'k',
    }
    return plt.rcParams.update(params)

def interface_plot():
    #from matplotlib import rcParams
    params = {'backend': 'ps',
              #'text.latex.preamble': ['\usepackage{amsmath}','\usepackage[utf8]{inputenc}'],
              #'text.latex.unicode': True,
              'axes.labelsize': 12, 
              'axes.titlesize': 12,
              'font.size': 10, 
              'font.family': 'serif',
              #'font.serif': 'Bookman',
              'legend.fontsize': 12,
              'xtick.labelsize': 10,
              'ytick.labelsize': 10,
              #'text.usetex': 0,
              #'text.dvipnghack' : True,
              'figure.figsize': [16.29,10],
              'figure.dpi': 300,
              'figure.facecolor' : 'white'
    }
    return plt.rcParams.update(params)
