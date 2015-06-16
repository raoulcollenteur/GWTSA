# -*- coding: utf-8 -*-

"""Setup script to install the GWTSA module in python
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'DESCRIPTION.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='GWTSA',
    # https://packaging.python.org/en/latest/single_source_version.html
    version='1.0.0.dev1',
    long_description=long_description,
    description='Software for time series analysis of groundwater levels',
    url='https://github.com/raoulcollenteur/GWTSA',
    author='Raoul Collenteur',
    author_email='info@raoulcollenteur.nl',
    license='GNU General Public License',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers :: Geohydrologists :: Hydrologists',
        'Topic :: Software Development :: Hydrology :: ',
        'License :: GNU General Public License',
        'Programming Language :: Python :: 2.7',
    ],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['cma', 'numpy', 'matplotlib', 'math', 'scipy'],

)