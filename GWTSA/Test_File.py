# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 18:58:37 2016

@author: Raoul
"""
import numpy as np
from pandas import read_csv, to_datetime
import matplotlib.pyplot as plt
import matplotlib.dates as md
#from TFN_Model import *
from statsmodels.tsa.stattools import acf
from lmfit import minimize, Parameters, Parameter, report_fit, fit_report

from scipy.special import gammainc
from scipy.stats import norm
from scipy.integrate import quad
from scipy.signal import fftconvolve
from Unsat_Zone import perc, pref, comb

from methods import *

#%% Time series class    
            
class TimeSeries(Methods):
    def __init__(self, bore, forcing, Discharge=None, rows=[50, 8], 
                 cl_unit=10000.0, gw_unit=100.0):
        # Load the precipitation and evaporation data and check values
        self.data = read_csv('%s' %forcing, skipinitialspace=True, skiprows=12, 
                             delimiter=',', parse_dates=[0], index_col=['Date'],                                                         usecols=[1,2,3], names=['Date','P','E'])
        self.data['E'] = self.data['E'] / (cl_unit * 1.0)   # Get the correct units 
        self.data['E'][self.data['E'] < 0.0] = 0.0          # No negative values
        self.data['P'] = self.data['P'] / (cl_unit * 1.0)   # Get the correct units  
        self.data['P'][self.data['P'] < 0.0] = 0.0          # No negative values
        
        # Load observed heads and check values
        parse = lambda x: md.datetime.datetime.strptime(x, '%d-%m-%Y')
        self.data['Ho'] = read_csv('%s' %bore, index_col=0, skiprows=50,
                                   usecols=[2, 5], parse_dates=True,
                                   date_parser=parse)  
        self.data['Ho'] = self.data['Ho'] / gw_unit # Get the correct units 

        # Load well discharge data if available
        if Discharge != None:   
            D = read_csv('%s' %Discharge, delimiter=';', index_col=0, 
                         parse_dates=True)
            self.data['Discharge'] = D.resample('1d', fill_method='ffill',
                                                base=0.0) # Downsample discharge to daily time step if necessary
            self.data['Discharge'].replace('nan', 0.0, inplace=True)
        
        self.data = self.data[(self.data.P>-999) & (self.data.E>-999)] # Only model where P and E are available
        
        self.data['t'] = range(len(self.data.P))                                    # Index
        self._index_observed = self.data.t[self.data.Ho > -999].tolist()   # Index list with observed heads
        self._dt = np.array(self._index_observed[1:]) - np.array(self._index_observed[0:-1]) # Time Steps
        self._t = self.data.t.tolist()

        # Save a name for the borehole
        self.bore = bore[-17:-8]
        print 'Model setup for bore id %s completed succesfully.' % self.bore

    def __repr__(self):
        return 'Time Series Data of Bore: %s' % self.bore
#%%
    def solve(self, X0, IR='IRF', RM='linear', trend=None, Cf=[1.0], 
              method='leastsq', period = ['01-01-1900', '01-01-2000', '01-01-2100'], 
              solver=1):
        # Define the TFN Model
        #self.IR = eval(IR) # Save the impulse response function
        self.IR = getattr(self, IR)
        self.RM = getattr(self, RM)  # Save the recharge calculation function
        self._TFN = IR      # Save the name of the impulse response function
        self._RM = RM       # Save the name of the recharge model

        self._index_calibration = self.data.t[(self.data.index > period[0]) & 
        (self.data.index < period[1]) & (self.data.Ho > -999)].tolist()
        self._index_validation = self.data.t[(self.data.index > period[1]) & 
        (self.data.index < period[2]) & (self.data.Ho > -999)].tolist()
        
        
        if method == 'leastsq':
            X0.add('d', value=self.data.Ho.mean(), vary=True)
            self.result = minimize(self.objective_function, X0, 
                                   method='leastsq', scale_covar=True)
            self.parameters_optimized = self.result.params.valuesdict()
            if self.result.success: 
                print 'Optimization completed succesfully!'
                print(report_fit(self.result))
                np.savetxt('Figures/fit_report_%s_%s.txt' %(self.bore, self._TFN),
                           (fit_report(self.result),), fmt='%str')
        else:
            X0.add('d', value=np.mean(self.head_observed))
            self.result = minimize(self.objective_function, X0, args=(InputData,), 
                                   method=method)
            self.parameters_optimized = self.result.params.valuesdict()

        # Calculate statistics for both the calibration and validation period 
        self.result.SWSI_Cal = sum(self.swsi(self._index_calibration))
        self.result.SWSI_Val = sum(self.swsi(self._index_validation))
        
        self.result.RMSE_Cal = sum(self.rmse(self._index_calibration))
        self.result.RMSE_Val = sum(self.rmse(self._index_validation))
        
        self.result.AVGDEV_Cal = self.avg_deviation(self._index_calibration)
        self.result.AVGDEV_Val = self.avg_deviation(self._index_validation)      
        
        self.result.EXPVAR_Cal = self.explained_variance(self._index_calibration)
        self.result.EXPVAR_Val = self.explained_variance(self._index_validation)  

# %% swsi constitutes the adapted version of the Sum of weighted squared innovations (swsi) developed by asmuth et al. (2005). For large values of alpha and small timesteps the numerator approaches zero. Therefore, Peterson et al. (2014) adapted the swsi function, making the numerator a natural log and changing the product operator to a summation.
      
    def objective_function(self, parameters):       
        self.alpha = parameters['alpha'].value
        self.construct_model(parameters)
        innovations = self.swsi(self._index_calibration)
        return innovations
        
    def construct_model(self, parameters):
        d = parameters['d'].value       
        self.data['recharge'] = self.RM(parameters) #Simulate the recharge series          
        self.Fb = self.IR(parameters) #block reponse function
        self.data['Hm'] = (d + fftconvolve(self.data.recharge, self.Fb))[self._t] 
        self.data['residuals'] = self.data.Ho - self.data.Hm
        self.data['innovations'] = self.data.residuals[self._index_observed][1:]
        - (self.data.residuals[self._index_observed][0:-1] * np.exp(-self._dt /
        self.alpha))

    def swsi(self, period):
        # Select the period to calculate statistics
        innovations = self.data.innovations[period[1:]]
        dt = np.array(period[1:]) - np.array(period[0:-1])    
        N = len(innovations) # Number of innovations
        numerator = np.exp((1.0/N) * sum(np.log(1 - np.exp(-2.0 / self.alpha * dt))))
        return (numerator / (1 - np.exp(-2.0 / self.alpha * dt ))) * innovations**2

    def rmse(self, period):
        return (self.data.Hm[period] - 
        self.data.Ho[period]**2) / len(self.data.Ho[period])

    def avg_deviation(self, period):
        return sum(self.data.Hm[period] - 
        self.data.Ho[period]) / len(self.data.Ho[period])

    def explained_variance(self, period):
        return (np.var(self.data.Ho[period])**2 - np.var(self.data.residuals[period])**2) / np.var(self.data.Ho[period])**2*100.


    

        
#%%
bores = 'Test_Data/B33A0113001_1.csv'
forcing = 'Test_Data/KNMI_Bilt.txt'

ml = TimeSeries(bores, forcing)

X0 = Parameters()
#           (Name,  Value,  Vary,   Min,  Max,  Expr)
X0.add_many(('A',  3,    True,   1e-6, None,  None),
            ('a',  200.,    True,  0.0, None,  None),
            ('n',    1.0,    False,   None, None,  None),
            ('alpha', 100.0,    True,   0.0,  1000.0,   None),
            # Recharge
            ('f',    1.0,    True,   0.0, 1.5,  None),
            ('Srmax',0.305,  False,  None, None,  None),
            ('Imax', 1.5e-3, False,  None, None,  None),
            #('Kp',  0.005,    True,   1e-4, 0.1,  None),
            ('Beta', 2.0,    True,   0.0, None,  None)
            )

ml.solve(X0, RM='linear',  method='leastsq')