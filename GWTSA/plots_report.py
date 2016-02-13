# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 09:25:59 2015

@author: Raoul

In this file all the different plots are made for the report
"""

cyan = [120/255.,196./255,1.]

from GWTSA import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
from pandas import read_csv
from scipy.special import gammainc, gamma

plt.close('all')
latex_plot()

#%% Plot the impulse response function
X0 = Parameters()
X0.add_many(('A',   0.3,    True,   None, None,  None),
           ('a',    200,    True,   None, None,  None),
           ('b',    0.0,    True,   None, None,  None),
           ('n',    1.5,    True,   None, None,  None))
x1 = TFN_Model.IRF(X0)

X0 = Parameters()
X0.add_many(('A',   2e-4,    True,   None, None,  None),
           ('a',    200,    True,   None, None,  None),
           ('b',    0.0,    True,   None, None,  None),
           ('n',    3.0,    True,   None, None,  None))
x2 = TFN_Model.IRF(X0)

X0 = Parameters()
X0.add_many(('A',   0.3,    True,   None, None,  None),
           ('a',    300,    True,   None, None,  None),
           ('b',    0.0,    True,   None, None,  None),
           ('n',    1.5,    True,   None, None,  None))
x3 = TFN_Model.IRF(X0)

X0 = Parameters()
X0.add_many(('A',   3.0,    True,   None, None, None),
           ('a',    200,    True,   None, None,  None),
           ('b',    2.5,    True,   None, None,  None),
           ('n',    1.0,    True,   None, None,  None))

x4 = TFN_Model.IRF(X0)

plt.figure(figsize=(4.15,3))    
plt.plot(x1, '-', color=cyan, label='A=0.3, a=200, n=1.5')
plt.plot(x2, '--', color=cyan, label='A=2e-4, a=200, n=3.0')
plt.plot(x3, 'k-', label='A=0.3, a=300, n=1.5')
plt.plot(x4, 'k--', label='A=3.0, a=200, n=1.0')

plt.xlim(0,1500)
plt.ylim(0,3.5)
plt.legend()
plt.xlabel('Time [Days]')
plt.ylabel('Response [-]')
plt.savefig('Figures/impulse_response.eps', bbox_inches='tight')

#%% Plot the recharge coefficients

S = np.arange(0,1.,0.01)
Beta = [0.5, 1.0, 3.0, 10.0, 15.0]
Cr = np.ones((len(Beta),100))
Scap=1.
for i in range(len(Beta)):
    Cr[i] = (S/Scap)**Beta[i]
plt.figure(figsize=(4.15,3))    
plt.plot(Cr[0],S, '-', color='k', label=r'$\beta$=0.5')
plt.plot(Cr[1],S, '-', color=cyan, label=r'$\beta$=1.0')
plt.plot(Cr[2],S, '--', color='k', label=r'$\beta$=3.0')
plt.plot(Cr[3],S, '--', color=cyan, label=r'$\beta$=10.0')
plt.plot(Cr[4],S, '-', color='k', label=r'$\beta$=15.0')
#plt.plot((1-Cr).T, S)
plt.ylabel(r'$S_r/S_{rmax}$ [-]'), plt.xlabel(r'$C_r$ [-]')
plt.yticks([0.2,0.4,0.6,0.8,1.0],[0.2,0.4,0.6,0.8,1.0], rotation=0)
plt.legend(loc=0)


plt.savefig('Figures/beta_parameter.eps', bbox_inches='tight')


S = np.arange(0,1.,0.01)
Beta = np.arange(0, 1.5, 0.05)
Cr = np.ones((len(Beta),100))
for i in range(len(Beta)):
    Cr[i] = (1-(1-(S/Scap))**Beta[i])
plt.figure()
plt.plot(Cr.T, S)
#plt.plot(1-Cr, S)
plt.ylabel('St'), plt.xlabel('Cr,1-Cr')
plt.legend(['Recharge', 'UZ'])


#r[i] = 1. / (1+np.exp((-S/Scap + 0.5)/Beta[i])) # The soil retention curve

#%% Plot the precipitation and evapotranspiration

data = read_csv('Test_Data/KNMI_Bilt.txt', skipinitialspace=True, skiprows=11, delimiter=',', parse_dates=['YYYYMMDD'], index_col=['YYYYMMDD'])

data.RH = data.RH/10000.
data.EV24 = data.EV24/10000.

P_avg = data.RH.resample('A', how='sum', kind='YYYYMMDD')
P_avg = P_avg[(P_avg.index>1958) & (P_avg.index<2015)]
E_avg = data.EV24.resample('A', how='sum', kind='YYYYMMDD')
E_avg = E_avg[(E_avg.index>1958) & (E_avg.index<2015)]
R_avg = P_avg - E_avg

P_mean = P_avg.mean()
E_mean = E_avg.mean()
R_mean = R_avg.mean()

presentation_plot()
plt.figure(figsize=(8.3,4))

plt.subplot(211)
P = data[(data.index > '1997-01-01') & (data.index < '2005-01-01')]
plt.bar(P.index, P.RH, color='k', label='Precipitation')

#ax = P_avg.plot(kind='bar', color=cyan, label='Precipitation')
plt.ylabel(r'$P$ [m/d]')
plt.ylim(-0.005,0.05)
plt.legend(loc='best')
#ax.grid(False)

plt.subplot(212)
plt.bar(P.index, P.EV24,color='k', label='Evaporation')
#ax1 = E_avg.plot(kind='bar', color='lightcoral', label='Evapotranspiration')
plt.ylabel(r'$E_p$ [m/d]')
plt.xlabel('Time [years]')
plt.ylim(-0.005,0.05)
#ax1.grid(False)
plt.legend(loc='best')
plt.savefig('Figures/climate_deelen.eps', bbox_inches='tight')

from pandas.stats.moments import rolling_mean

b = bore[1]
b['E']=rolling_mean(data.EV24, 3650)
b['P']=rolling_mean(data.RH, 3650)
b['R']=b.P-b.E
b.corr()
b[b.index<'2005-09-23 00:00:00'].corr()
#%% Plot the boreholes used in this study

#bores = glob.glob('Test_Data/*.csv')
bores = [ 'Test_Data/B27D0001001_1.csv', 'Test_Data/B27C0049001_1.csv', 'Test_Data/B33A0113001_1.csv']

bore = [] #Store de different boreholes data
peak = [] #Store the 1995 peak moment


parse = lambda x: md.datetime.datetime.strptime(x, '%d-%m-%Y')

for i in range(len(bores)):
    bore.append(read_csv('%s' %bores[i], parse_dates=True, index_col=2, skiprows=16, skipinitialspace=True, usecols=[2,5], date_parser=parse))
    bore[i].rename(columns={'Stand (cm t.o.v. NAP)': 'h'}, inplace=True)
    bore[i].h = bore[i].h/100.
    #bore[i].h = bore[i].h-bore[i].h.mean()
    #bore[i].h.plot()
    bore[i] = bore[i][(bore[i].index > '1958-01-01 00:00:00') & (bore[i].index < '2005-01-01 00:00:00')]
    peak.append(bore[i].h[(bore[i].index > '1975-01-01 00:00:00') & (bore[i].index < '1980-01-01 00:00:00')].argmax())
    bores[i] = bores[i][-17:-8]
    
plt.figure(figsize=(8.3,2))
plt.plot(bore[0].index, bore[0].h, linestyle='-', color='k', label='B27D00010')
plt.plot(bore[1].index, bore[1].h, linestyle='--', color='k', label='B27C00490')
plt.plot(bore[2].index, bore[2].h, linestyle='-', color=cyan, label='B33A01130')
plt.xlabel('Time [Years]')
plt.ylabel('GWL [m]')
plt.legend(loc=(0,1), ncol=3, frameon=False, handlelength=3)
plt.savefig('Figures/boreholes_plot.eps', bbox_inches='tight') 
 
#%% 
from pandas.stats.moments import rolling_mean 
 
b = bore[0]
dt=5
b['E']=rolling_mean(data.EV24, dt*365)
b['P']=rolling_mean(data.RH, dt*365)
b['R']=b.P-b.E
print b.corr().iloc[0]
label = []
period = range(1970,2006,1)
cor =[0,0,0,0]
for i in range(len(period)-dt):
    label = np.append(label,'%s-%s' %(period[i], period[i+dt]))
    cor = np.vstack((cor, b[(b.index>'%s' %period[i]) & (b.index<'%s' %period[i+dt])].corr().iloc[0]))

latex_plot()
plt.figure(figsize=(4.15,3.0))
plt.plot(period[0:-dt],cor[1:,1], 'k-o', label='Evaporation')
plt.plot(period[0:-dt],cor[1:,2], 'k', label='Precipitation')
plt.plot(period[0:-dt],cor[1:,3], 'k--', label='Recharge')
plt.xticks(period[0::5],label[0::5], rotation=0)
plt.ylabel('Correlation')
plt.xlabel('Period [years]')
plt.xlim(1970,2000) 
plt.legend(loc=0)
plt.savefig('Figures/correlation.eps', bbox_inches='tight')

latex_plot()
plt.figure(figsize=(4.15,3.0))
plt.subplot(211)
ax = (b.R[b.index>'1967-07-20 00:00:00']*1000).plot(color='k')
ax.set_xticklabels([])
plt.ylabel('Recharge [mm/d]')
plt.subplot(212)
b.h[b.index>'1967-07-20 00:00:00'].plot(color='k')
plt.ylabel('Head [m]')
plt.xlabel('Time [years]')
plt.savefig('Figures/mov_avg_10yr_recharge.eps', bbox_inches='tight') 

plt.plot(b.R, b.h)
 
#%% Plot of the thickness of the unsaturated zone against the T_peak

Thickness = [10.35, 34.02, 29.24, 6.34, 49.09, 30.95, 33.06, 23.02, 47.10] #Estimated based on ground level minus highest groundwater level 

days = []
for i in range(len(peak)):
    days.append(peak[i] - min(peak))
    days[i] = days[i].days

plt.figure()
#plt.xkcd()
plt.plot(Thickness, days, 'k+', markersize=10, markeredgewidth=3)
plt.ylabel('Relative Delay Time [Days]')
plt.xlabel('Thickness Unsaturated Zone [Meters]')
plt.ylim(-5,300)
plt.xlim(0,55)
for i in range(len(bores)):
    plt.annotate(bores[i],(Thickness[i],days[i]))
plt.savefig('delay.eps', format='eps', bbox_inches='tight')

slope, intercept, r_value, p_value, std_er = linregress(days, Thickness)
y = intercept + slope * (np.linspace(0,300))
plt.plot(y,np.linspace(0,300), 'k--')
#
#a = bore[0].resample('2W')
#b = bore[1].resample('2W')
#c = bore[2].resample('2W')
#
#a = a.interpolate(method='cubic')
#a = a[(a.index > '1975-01-15 00:00:00') & (a.index < '2004-01-15 00:00:00')]
#b = b.interpolate(method='cubic')
#b = b[(b.index > '1975-01-15 00:00:00') & (b.index < '2004-01-15 00:00:00')]
#c = c.interpolate(method='cubic')
#c = c[(c.index > '1975-01-15 00:00:00') & (c.index < '2004-01-15 00:00:00')]
#
#import scipy
#af = scipy.fft(a)
#bf = scipy.fft(c)
#d = scipy.ifft(af * scipy.conj(bf))
#
#time_shift = np.argmax(abs(d))

#%% Plot the Impulse, Step and block response

# First start with the input parameters
A = 10.
a = 5.
n = 1.5

# Create the impulse response curve
t1 = np.arange(0,30.,0.1)
Fi = A*1./a**n*t1**(n-1) * np.exp(-t1/a)/gamma(n)
Fi = A * t1**(n-1)*np.exp(-t1/a)

# Create the step response curve
t = np.arange(0,30.,1)
Fs = A * gammainc(n,t/a) # Step response function based on pearsonIII
#Fs = np.cumsum(Fi)

# Create the block response curve
Fb = Fs[1:] - Fs[0:-1]
Fb = np.append(0, Fb) #This is only done for drawing the graph as you normally Convolute 

# Make the plots
plt.subplot(321)
plt.bar([0],[1],width=0.1, color='gray')
plt.xlim(-1,8)
plt.ylim(0,1.2)
plt.yticks([1],['$\infty$'])
plt.ylabel('N [L]')

plt.subplot(322)
plt.plot(t1,Fi, 'k')
plt.ylabel(r'$\theta$ [-]')

plt.subplot(323)
plt.bar([0],[1],width=8, color='gray')
plt.ylim(0,1.2)
plt.xlim(-1,8)
plt.xticks([0, 8],[0, '$\infty$'])
plt.ylabel('N [L]')

plt.subplot(324)
plt.plot(np.arange(len(Fs)),Fs, 'k')
plt.ylabel(r'$\theta_s$ [-]')

plt.subplot(325)
plt.bar([0],[1],width=1, color='gray')
plt.ylim(0,1.2)
plt.xlim(-1,8)
plt.ylabel('N [L]')
plt.xlabel(' Time [T]')

plt.subplot(326)
plt.plot(np.arange(len(Fb)),Fb, 'k')
plt.ylabel(r'$\theta_b$ [-]')
plt.xlabel(' Time [T]')
  
plt.savefig('Figures/responses.eps', bbox_inches='tight')   

#%% 

X0 = Parameters()
X0.add_many(('A',   0.20,    True,   None, None,  None),
           ('t_p',    2.0,    True,   None, None,  None),
           ('n',    1.5,    True,   None, None,  None))
x1 = TFN_Model.IRF2(X0)

X0 = Parameters()
X0.add_many(('A',   0.0003,    True,   None, None,  None),
           ('t_p',    2.5,    True,   None, None,  None),
           ('n',    2.8,    True,   None, None,  None))
x2 = TFN_Model.IRF2(X0)


plt.figure(figsize=(4.15,3))    
plt.plot(x1, '-', color='k', label='A=0.2, $T_{Peak}$=100, n=1.5')
plt.plot(x2, '-', color=cyan, label='A=3e-4, $T_{Peak}$=316, n=2.8')


plt.annotate(r'$\Delta T_{peak}$', xy=(100.0, 1.0), xytext=(316.0, 1.0), arrowprops={'arrowstyle': '<|-|>', 'color': 'k'}, va='center', fontsize=12)
plt.annotate(r'$T_{peak}$', xy=(0.0, 1.2), xytext=(100.0, 1.2), arrowprops={'arrowstyle': '<|-|>', 'color': 'k'}, va='center', fontsize=12)
plt.xlim(0,500)
plt.ylim(0)
plt.legend(loc=0)
plt.xlabel('Time [Days]')
plt.ylabel('Response [-]')
plt.savefig('Figures/IRF_Peak.eps', bbox_inches='tight')


#%% moving averages plots
from pandas.stats.api import ols

df = bore[0]
df1 = bore[1]
df2 = bore[2]

df['date'] = df.index.to_datetime() 
df['date_delta'] = (df['date'] - df['date'].min())  / np.timedelta64(1,'D')
res = ols(y=df.h, x=df.date_delta)

print res
print res.summary_as_matrix
plt.figure(figsize=(8.3,2))
res.y.plot(color='k', label='B27D00010')
res.y_fitted.plot(color='k', label='linear_fit')


df1['date'] = df1.index.to_datetime() 
df1['date_delta'] = (df1['date'] - df1['date'].min())  / np.timedelta64(1,'D')
res = ols(y=df1.h, x=df1.date_delta)

print res
print res.summary_as_matrix
res.y.plot(color='gray', label='B27C00490')
res.y_fitted.plot(color='gray', label='linear_fit')


df2['date'] = df2.index.to_datetime() 
df2['date_delta'] = (df2['date'] - df2['date'].min())  / np.timedelta64(1,'D')
res = ols(y=df2.h, x=df2.date_delta)
print res
print res.summary_as_matrix

res.y.plot(color=cyan, label='B33A01130')
res.y_fitted.plot(color=cyan, label='linear_fit')

plt.xlabel('Time [Years]')
plt.ylabel('GWL [m]')
plt.legend(loc=(0,1), ncol=3, frameon=False, handlelength=3)

plt.savefig('Figures/linear_trend.eps', bbox_inches='tight')

#%% Plot the discharges

data = read_csv('Test_Data/Discharge_Tongeren.csv', delimiter=';', index_col=0, parse_dates=True)

ax = data.plot(kind='bar', color='k', figsize=(8.3,2.0))
plt.ylabel('Discharge M $m^3$/year')
plt.xlabel('Year')
label = ax.get_xticklabels()
xticks = ax.get_xticks()
ax.set_xticks(xticks[0::5])
ax.set_xticklabels(range(1958,2012,5))
plt.legend(['Discharge B27C0049'])
plt.savefig('Figures/Discharge_Tongeren.eps', bbox_inches='tight')