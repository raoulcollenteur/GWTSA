# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:45:57 2015

@author: Raoul
"""

__name__ ='GWTSA_new'
__description__='Software to perform time series analysis on groundwater levels'
__author__='Raoul Collenteur'
__author_email__='info@raoulcollenteur.nl'
__url__='https://github.com/raoulcollenteur/GWTSA'
__license__='GNU General Public License'
__version__='1.0.0.dev1'


# Explicitly tell what modules should be imported
from time_series import *

# Cythonized Unsat_Zone Module has to be imported seperately as it is an different file extension
from distutils.core import Extension
ext_modules = Extension('Unsat_Zone', ['Unsat_Zone.so', 'Unsat_Zone.c'])

