# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 10:08:36 2015

@author: Raoul
"""
import numpy as np
import matplotlib.dates as md
import matplotlib.pyplot as plt
from scipy.optimize import fmin, fmin_l_bfgs_b
from Objective_Function import SWSI2,EVP
from scipy.special import gammainc
import TFN_Model
import cma


class TimeSeries:


    '''Here we initiate the Time Series Analysis for a specific borehole. We import the data from the text-file and get all the time data necessary for the analysis and plotting later. '''

    def __init__(self, Bore):
        
        Ho = np.genfromtxt(('./Data/%s.txt' % (Bore)), delimiter=',', skiprows=5, usecols=[0, 1], converters={0: md.strpdate2num('%Y%m%d')});
        Ho = Ho[Ho[:,1] > -999] #Select only real values
        self.Time_Begin = Ho[0,0] #In time number, not date
        self.Time_End = Ho[-1,0] #In time number, not date
        self.Time_Steps = Ho[1:,0] - Ho[0:-1,0] #Calculate the timesteps
        self.Time_Axis = Ho[:,0] # Real Time in time number format for plots
        self.Time_Observed = np.append(0, [np.cumsum(self.Time_Steps)]).astype(int)                 #Timesteps with observed heads
        self.Time_Model = np.arange(0,sum(self.Time_Steps,1), dtype=int)     
        self.Head_Observed = Ho[:,1] #Save Observed Heads column only
        self.Time_Start = 2000; #warmup time of the model [Days] not to consider in the calibration. This parameter is hardcoded into the program.
               
        # Import the precipitation and evaporation data and calculate the recharge
        ClimateData = np.genfromtxt('./Data/KNMI_Bilt.txt', delimiter=',', skiprows=8, converters={1: md.datestr2num}, usecols=[1,2,3]);
        self.ClimateData = ClimateData[:,:][(ClimateData[:,0] >= self.Time_Begin) & (ClimateData[:,0] <= self.Time_End)]

        self.bore = Bore


    '''In this part the actual model is run and the parameters are optimized. This function takes an array of initial parameters and the TFN-model that has been chosen for the Time Series Analysis. '''
        
    def Model(self, TFN,Par):
        
        # Define the TFN Model and the Initial Parameters
        if TFN == 'TFN1':
            self.Initial_Parameters = [Par['A'],Par['a'],Par['n'],np.mean(self.Head_Observed),Par['Alpha']]
        elif TFN == 'TFN2':
            self.Initial_Parameters = [Par['A'],Par['a'],Par['n'],np.mean(self.Head_Observed),Par['Alpha']]
        elif TFN == 'TFN3':
            self.Initial_Parameters = [Par['A'],Par['a'],Par['n'],np.mean(self.Head_Observed),Par['Alpha']]
        elif TFN == 'TFN4':
            self.Initial_Parameters = [Par['A'],Par['a'],Par['n'],np.mean(self.Head_Observed),Par['Alpha'],Par['S_cap'],Par['Beta']]
            mybounds = [(0,10),(0,10),(0,10),(0,None),(0,10),(0,10),(0,10)]    
           
        else:
            print 'Error: TFN model does not exist, chose another one or check TFN_Model.py!'
        
        self.TFN = getattr(TFN_Model, TFN)          
        
        InputData = [self.TFN, self.Time_Model, self.ClimateData, self.Head_Observed, self.Time_Observed, self.Time_Steps, self.Time_Start] 
        #self.Parameters = fmin(SWSI2,self.Initial_Parameters,args=(InputData,))     
        #self.Parameters = fmin_l_bfgs_b(SWSI2, x0=self.Initial_Parameters, args=(InputData,), bounds=mybounds, approx_grad=True)  
        self.Parameters = cma.fmin(SWSI2, self.Initial_Parameters, 2, args=(InputData,))[0] #options={'boundary_handling': 'BoundTransform','bounds': [np.zeros(10), np.ones(10)] }
        [self.Head_Modeled, self.Innovations] = self.TFN(self.Parameters, InputData); #model the GWL

        #Print the output parameters and compare to the real values'''
        self.Parameter_Names = ['A','a','n','d','alpha','S_cap','Beta']
        print 'Explained Variance Percentage is:', EVP(self.Parameters,InputData)
        for i in range(0,len(self.Parameters)):
            print self.Parameter_Names[i], '=', self.Parameters[i]  



    '''This function is used for easy plotting of the observed and the modeled head. the Modeled head is drawn as a line whereas the observed head (often fewer data points) is drawn in dots. Both have the same color, making it easy to draw more boreholes/models in one figure. '''                     
           
    def Plot_Heads(self,color):
        plt.plot(md.num2date(self.Time_Axis), self.Head_Observed,'k.')
        plt.plot(md.num2date(np.arange(self.Time_Begin,self.Time_End+1)), self.Head_Modeled[self.Time_Model],'-',color=color)
        plt.legend(['Observed Head','Modeled Head'],loc=4)
        plt.xlabel('Time [T]',size=20)
        plt.ylabel('Groundwater Head [L]',size=20)
        plt.title('%s' % (self.bore))
        
    def Plot_Innovations(self):
        plt.plot(self.Innovations)
        plt.title('Innovations of the time series')
        plt.xlabel('Innovations [-]')
                
        
    def Plot_Forcings(self):
        plt.figure()
        plt.bar(md.num2date( self.ClimateData[:,0]),self.ClimateData[:,1],color='b',lw=0)
        plt.ylabel('P [0.1 mm/D]',size=20)
        plt.ylim(0,max(self.ClimateData[:,1]))        
        ax1 = plt.gca()
        ax2 = ax1.twinx() # To create a second axis on the right
        plt.bar(md.num2date( self.ClimateData[:,0]), self.ClimateData[:,2], color='red', lw=0)
        ax2.set_ylim((0,400)[::-1])
        plt.ylabel('E [0.1 mm/D]',size=20)
        plt.xlabel('Time [Years]',size=20)
        plt.legend('Precipitation','Potential Evapotranspiration')
        plt.title('Forcings',size=20)
        
    def Plot_ImpulseResponseFunction(self):
        Fs = self.Parameters[0] * self.Parameters[1] * gammainc(self.Parameters[2], self.Time_Model/self.Parameters[1])
        Fb = Fs[1:] - Fs[0:-1]
        plt.plot(self.Time_Model[0:-1],Fb)
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        