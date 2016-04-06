# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 11:39:58 2014
this is a time series analysis model
@author: R.A. Collenteur
"""

# Import all the packages needed
from GWTSA_new import *

plt.close('all')

bores = 'Data/B27D0001001_1.csv'
forcing = 'Data/KNMI_Bilt.txt'

ml = TimeSeries(bores, forcing, discharge = None)

X0 = Parameters()
#           (Name,  Value,  Vary,   Min,  Max,  Expr)
X0.add_many(('A',  3,    True,   1e-6, None,  None),
            ('a',  3.0,    True,  0.0, None,  None),
            ('n',    1.0,    True,   None, None,  None),
            ('alpha', 2.0, True,   0.0,  3.0,   None),
            # Recharge
            ('f',    1.0,    True,   0.0, 1.5,  None),
            ('Srmax', 0.305,  False,  None, None,  None),
            ('Imax', 1.5e-3, False,  None, None,  None),
            ('Kp',  0.005,    True,   1e-4, 0.1,  None),
            #('Beta', 2.0,    True,   0.0, None,  None),
            ('Gamma', 3.0,    True,   0.0, None,  None)
            )

ml.solve(X0, RM='percolation', period = ['01-01-1975', '31-12-2003', '31-12-2004'], method='leastsq')
interface_plot()
ml.plot_results()