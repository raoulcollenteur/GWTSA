# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 11:22:18 2016

@author: Raoul
"""

from Test_File import method    
class method(object):
    def plus(self):
        print self.x + self.y 
    def multiply(self):
        print self.x * self.y       
        
class Main(method):
    def __init__(self, x, y):
        self.x = x
        self.y = y

  
        
ml = Main(10, 5)