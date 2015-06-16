# -*- coding: utf-8 -*-
"""
Created on Sat Feb 28 23:15:23 2015

Run this script in the terminal (on a Mac) by typing (without brackets): "python setup.py build_ext --inplace"
    - Make sure you are in the directory where this setup.py file is placed!
@author: Raoul
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np
from Cython.Build import cythonize

ext  =  [Extension( "Unsat_Zone", sources=["Unsat_Zone.pyx"] )]

setup(
    name = 'Unsat_Zone',
    cmdclass={'build_ext' : build_ext}, 
    include_dirs = [np.get_include()],    
    ext_modules = cythonize("Unsat_Zone.pyx"),
)
