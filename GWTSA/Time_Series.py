# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 18:58:37 2016

@author: Raoul
"""
import numpy as np
from pandas import read_csv, to_datetime
import matplotlib.pyplot as plt
import matplotlib.dates as md
from statsmodels.tsa.stattools import acf
from lmfit import minimize, Parameters, Parameter, report_fit, fit_report
from methods import *

cyan = [120/255.,196./255,1.]

#%% Time series class

class TimeSeries(Methods):
    def __init__(self, bore, forcing, discharge=None, rows=[50, 5],
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
        self.data = self.data[(self.data.P>-999) & (self.data.E>-999)] # Only model where P and E are available
        # Load well discharge data if available
        if discharge != None:
            D = read_csv('%s' %discharge, delimiter=';', index_col=0,
                         parse_dates=True)
            self.data['discharge'] = D.resample('1d', fill_method='ffill',
                                                base=0.0) # Downsample discharge to daily time step if necessary
            self.data['discharge'].replace('nan', 0.0, inplace=True)

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
    def solve(self, X0, IR='IRF', RM='linear', TR=None, Cf=[1.0],
              method='leastsq', period = ['01-01-1900', '01-01-2000', '01-01-2100'],
              solver=1):
        # Define the TFN Model
        #self.IR = eval(IR) # Save the impulse response function
        self.IR = getattr(self, IR)
        self.RM = getattr(self, RM)  # Save the recharge calculation function
        if TR != None:
            self.TR = getattr(self, TR)
        self._TFN = IR      # Save the name of the impulse response function
        self._RM = RM       # Save the name of the recharge model

        self._index_calibration = self.data.t[(self.data.index > period[0]) &
        (self.data.index < period[1]) & (self.data.Ho > -999)].tolist()
        self._index_validation = self.data.t[(self.data.index > period[1]) &
        (self.data.index < period[2]) & (self.data.Ho > -999)].tolist()
        self.period = period

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
        self.result.SWSI_Cal = np.sqrt(sum(self.swsi(self._index_calibration)))
        self.result.SWSI_Val = np.sqrt(sum(self.swsi(self._index_validation)))

        self.result.RMSE_Cal = np.sqrt(sum(self.rmse(self._index_calibration)))
        self.result.RMSE_Val = np.sqrt(sum(self.rmse(self._index_validation)))

        self.result.AVGDEV_Cal = self.avg_deviation(self._index_calibration)
        self.result.AVGDEV_Val = self.avg_deviation(self._index_validation)

        self.result.EXPVAR_Cal = self.explained_variance(self._index_calibration)
        self.result.EXPVAR_Val = self.explained_variance(self._index_validation)

# %% swsi constitutes the adapted version of the Sum of weighted squared innovations (swsi) developed by asmuth et al. (2005). For large values of alpha and small timesteps the numerator approaches zero. Therefore, Peterson et al. (2014) adapted the swsi function, making the numerator a natural log and changing the product operator to a summation.

    def objective_function(self, parameters):
        self.alpha = 10**parameters['alpha'].value
        self.construct_model(parameters)
        innovations = self.swsi(self._index_calibration)
        self.test = innovations
        return np.array(innovations)

    def construct_model(self, parameters):
        d = parameters['d'].value
        self.data['recharge'] = self.RM(parameters) #Simulate the recharge series
        self.Fb = self.IR(parameters) #block reponse function
        self.data['Hm'] = (d + fftconvolve(self.data.recharge, self.Fb))[self._t]

        try:
            trend = self.TR(parameters)
            head_modeled += trend
        except:
            pass

        self.data['residuals'] = self.data.Ho - self.data.Hm
        self.data['innovations'] = self.data.residuals[self._index_observed][1:] - (self.data.residuals[self._index_observed][0:-1] * np.exp(-self._dt *1.0 / self.alpha))

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

# %% Different plotting functions go here
    def plot_results(self, savefig=True):
        plt.figure('%s_%s_%s' %(self.bore,self._TFN,self._RM))
        plt.suptitle('GWTSA %s Model Results for bore %s' %(self._TFN, self.bore), fontweight='bold')
        gs = plt.GridSpec(3, 4, wspace=0.2)

       # Plot the Groundwater levels
        ax2 = plt.subplot(gs[1,:-1])
        self.data.Ho[self._index_observed].plot(marker='.', color='k', markersize=2, label='observed head')
        self.data.Hm.plot(color=cyan, label='modeled head')
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
        plt.axvline(self.period[1], color='gray', linestyle=':', label='calibration/validation')
        plt.ylim(self.data.Ho.min(), self.data.Ho.max())

        # Plot the recharge
        ax1 = plt.subplot(gs[0,:-1], sharex=ax2)
        x = self.data.recharge.resample('A', how='sum')
        plt.bar(x.index, x.values, color=cyan, edgecolor='gray', width=250)
        ax1.grid(0)
        plt.ylabel('Recharge [m/year]')
        plt.text(0.1,0.1, 'average %.2f m/y' %np.mean(self.data.recharge.resample('A', how='sum')), backgroundcolor='w', alpha=0.5)
        ax1.xaxis.set_visible(False)

        # Plot the residuals and innovations
        ax3 = plt.subplot(gs[2,:-1], sharex=ax2)
        self.data.residuals[self._index_observed].plot(color='k')
        self.data.innovations[self._index_observed].plot(color=cyan, label='innovations')
        plt.legend(loc=(0,1), ncol=3, frameon=False, handlelength=3)
        plt.ylabel('Error [m]')
        plt.xlabel('Time [Years]')

        # Plot the Impulse Response Function
        ax4 = plt.subplot(gs[0,-1])
        plt.plot(self.Fb, 'k')
        plt.xticks(range(0,10000,1000))
        plt.xlim(0,np.where(np.cumsum(self.Fb)>0.99*sum(self.Fb))[0][0]) # cut off plot after 99.0% of the response
        plt.ylim(0.0)
        plt.text(10,0.1,'Peak Time: %i' %self.Fb.argmax(), verticalalignment='bottom')
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
        plt.text(0.05, 0.84, 'SWSI: %.3f / %.3f meter' %(self.result.SWSI_Cal, self.result.SWSI_Val))
        plt.text(0.05, 0.68, r'Expl. var: %.2f / %.2f %%' %(self.result.EXPVAR_Cal, self.result.EXPVAR_Val))
        plt.text(0.05, 0.52, 'RMSE: %.3f / %.3f meter' %(self.result.RMSE_Cal, self.result.RMSE_Val))
        plt.text(0.05, 0.36, 'Avg dev: %.3f / %.3f meter' %(self.result.AVGDEV_Cal, self.result.AVGDEV_Val))
        plt.text(0.05, 0.20, 'AIC: %.2f' %self.result.aic)
        plt.text(0.05, 0.04, 'BIC: %.2f' %self.result.bic)

        if savefig:
            plt.savefig('Figures/%s_%s_%s.eps' %(self.bore,self._TFN,self._RM), bbox_inches='tight')

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
# %%


