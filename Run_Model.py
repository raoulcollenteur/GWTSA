# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 11:39:58 2014
this is a time series analysis model
@author: R.A. Collenteur
"""

#Import all the packages needed
#from datetime import datetime
from Time_Series import *

# close the previous figures
plt.close('all')

Bore = 'B27C0002-001' #, 'B27B0238-002', 'B27B0081-001', 'B27C0002-001', 'B32F0002-001', 'B33A0113-001', 'B33C0140-001', 'B39E0117-001', 'B40B0304-001'

Par = {'A': 1,'a': 20, 'n': 1.35,'Alpha':7,'K_sat':1.18,'Beta': 1.0,'S_cap': 10.0} # initial parameters
TFN = 'TFN4' #TFN-model to use

Bore1 = TimeSeries(Bore)
Bore1.Model(TFN,Par)

Bore1.Plot_Heads('r')

#TFN = 'TFN4' #TFN-model to use
#
#Bore1.Model(TFN,Par)
#Bore1.Plot_Heads('b')

#Bore1.Plot_Forcings()
