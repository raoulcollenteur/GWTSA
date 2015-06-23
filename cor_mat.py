## -*- coding: utf-8 -*-
#"""
#Created on Tue Jun 23 12:29:14 2015
#
#@author: Raoul
#"""
#
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#import numpy as np
#

#
#surf = ax.plot_surface(a,b,y, rstride=1, cstride=1)
#ax.set_zlim(-1.01, 1.01)
#
#fig.colorbar(surf, shrink=0.5, aspect=5)
#
#plt.show()

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.gca(projection='3d')

def func(a,b):
    y = a + a*b**2
    return y
    
a = np.random.uniform(0,1,10)    
b = np.random.uniform(0,1,10)



a, b = np.meshgrid(a, b)
y = func(a,b)

surf = ax.plot_surface(a,b,y, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)


ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

