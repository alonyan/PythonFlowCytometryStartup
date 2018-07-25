# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 15:17:22 2016

 Load standards, fit N gaussians for APC channel to locate centers, 
 gate +-1 sigma for each one, Find mean/gmean/whatever, Fit hill curve
 Load samples, gate, interpolate. return table of values
@author: Alonyan
"""
#%%
import numpy as np
import sys
sys.path.append("/Users/oyleryaniva/Documents/Python/Modules")
from FlowAnalysis import LoadFlowSamples as LFS
import matplotlib.pyplot as plt
import os
os.chdir('/Users/oyleryaniva/Documents/Python/example')

#%% Load standards 

datadir = '/Users/oyleryaniva/Documents/Python/example/FCS files for CBA'
pattern = 'Calibrations_Tube_'
Standards = LFS.FlowData(datadir, pattern)

#%%  Make histogram of one of the standards, remove debrees, and plot
xChannel = 'APC-A' #Channel for bead identifier
NSamp = 0 #shouldn't matter too much, as long as it's in range
datax=Standards.samples[NSamp].data[xChannel] #data
hist, edges = np.histogram(datax, bins=200)
hist = np.double(hist)
centers = (edges[1:]+edges[:-1])/2
hist[centers<1]=0
hist = hist/sum(hist)
#% use cluster analysis to find initial guess and fit to individual gaussians
from scipy import optimize
from Utilities import kmeans

def Gaussians(x, mu, sigma, Amp):# general function of N gaussians
    if len(Amp)==len(sigma)==len(mu):
        return np.array([(Amp/sigma)*np.exp((-(x0-mu)**2)/(2*sigma**2))/np.sqrt(2*np.pi) for x0 in x]).sum(axis=1)
    else:
        print('must have 3n fitting parameters')
def Gau6(x, *p): return Gaussians(x, np.array(p[0:6]), np.array(p[6:12]), np.array(p[12:18]))#stupid function for 6 gaussians


pvar = np.ones((2,2))#some initialization
while max(sum(pvar))>1e-4:#check that the fit is good
    mu = kmeans(datax[datax>1],6).mu #automated guessing
    beta00 = mu+[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, .02, .02, .02, .02, .02, .02]
    p, pvar = optimize.curve_fit(Gau6, centers, hist,p0=beta00)    

import pylab as pl
fig = pl.figure(num=None, figsize=(5, 5), dpi=80, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(111)
ax1.set_position([0.15,0.15,1,0.6])
plt.plot(centers, hist,centers, Gau6(centers,*p))

ax1.set_xlabel(xChannel)
a = np.outer(np.arange(1,10),10**np.arange(1,2)).T.reshape((1,-1)).squeeze()
ticks = np.append(-a[::-1],0)
ticks = np.append(-100,ticks)
a = np.outer(np.arange(1,10),10**np.arange(1,6)).T.reshape((1,-1)).squeeze()
ticks = np.append(ticks,a[:])
emptvec = ['','','','','','','','']
tickmarks = ['-0.1']+emptvec+['']+['0']+emptvec+['']+['0.1']+emptvec+['1']+emptvec+['10']+emptvec+['100']+emptvec
#ax1.set_xticks(np.arcsinh(np.array([-100,-90,-80,-70,-60,-50,-40,-30,-20, -10, 0,10,20,30,40,50,60,70,80,90,100, 1000, 10000, 100000])/150.))
ax1.set_xticks(np.arcsinh(ticks/150.))
ax1.set_xticklabels(tickmarks)
ax1.set_xlim(left=-0.7, right = 7.5)
#%reshape and sort
p = p.reshape((3,6))

p[2][:]= p[2][p[0][:].argsort()]
p[1][:]= p[1][p[0][:].argsort()]
p[0][:]= p[0][p[0][:].argsort()]

#%% make gates based on mean+-2std
from FlowCytometryTools import IntervalGate
gates = []
for i in np.arange(0,6):
    gates = gates + [IntervalGate((p[0][i]-2*p[1][i], p[0][i]+2*p[1][i]), ('APC-A'), region='in')]

#%% Make tables for MFIs for all the different beads and assign names using robust geomeans
import pandas as pd
MFItableCals = Standards.gate(gates[0]).geoMeans()['PE-A']
for i in np.arange(1,6):
    MFItableCals = pd.concat([MFItableCals, Standards.gate(gates[i]).geoMeans()['PE-A']], axis=1, ignore_index=True)
MFItableCals.columns = ['IL12-p70', 'TNF', 'IFNg', 'MCP-1','IL10', 'IL6']
MFItableCals
#%% cool, now just put in the concentrations and plot the different standards
gpl = np.append(5000/2**np.arange(0,10),0) #pg/l
Mw = np.array([70, 18, 18, 8.5, 20, 22]) #kDa

ConcTableCals = pd.DataFrame(gpl/Mw[0]) #pM
for i in np.arange(1,6):
    ConcTableCals = pd.concat([ConcTableCals, pd.DataFrame(gpl/Mw[i])], axis=1, ignore_index=True)
ConcTableCals.columns = ['IL12-p70', 'TNF', 'IFNg', 'MCP-1','IL10', 'IL6']
ConcTableCals.index = MFItableCals.index

fig = plt.figure(num=None, figsize=(5, 4), dpi=100, facecolor='w', edgecolor='k')  

for param in MFItableCals.columns:
    plt.semilogx(ConcTableCals[param], MFItableCals[param],'-o')
plt.legend(MFItableCals.columns, bbox_to_anchor=(0.4, 1),frameon=False)
    
#%% Fit all to hill functions to get extrapolation functions and make interpolators out of the fitted functions
from scipy import interpolate

def Hill1(x, amp, EC50, BG): #Hill function
    return BG+amp*1/(1+(EC50/x))

FitParam = {}
FitFunc = {}
InterpFun = {}
x = np.append(0, 10**np.linspace(-1, 3, 10000))
for param in MFItableCals.columns:
    print(param)
    p, pvar = optimize.curve_fit(Hill1,ConcTableCals[param],MFItableCals[param],p0=[100000, 20, 100],bounds=([0, 0, 0], [ np.inf, np.inf, np.inf]))
    FitParam[param] = p 
    FitFunc[param] = Hill1(x,*p)
    InterpFun[param] = interpolate.interp1d(FitFunc[param],x, kind='nearest', fill_value='extrapolate')

#%% If you really feel the need to plot something right now, please go ahead
param = MFItableCals.columns[5]
plt.semilogy(MFItableCals[param],ConcTableCals[param],'-o', MFItableCals[param],InterpFun[param](MFItableCals[param]),'-d')
plt.xlabel('gMGI')
plt.ylabel(param+'[pM]')
plt.legend(['Measured', 'Postdicted'],bbox_to_anchor=(1, 0.6),frameon=False)

#%% Now, if the postdicted concentrations are close enough,  load samples
pattern = 'R'
Samples = LFS.FlowData(datadir, pattern)
#%% Make tables for MFIs for all the different beads and assign names using robust geomeans
import pandas as pd
MFItableSamples = Samples.gate(gates[0]).geoMeans()['PE-A']
for i in np.arange(1,6):
    MFItableSamples = pd.concat([MFItableSamples, Samples.gate(gates[i]).geoMeans()['PE-A']], axis=1, ignore_index=True)
MFItableSamples.columns = ['IL12-p70', 'TNF', 'IFNg', 'MCP-1','IL10', 'IL6']


#Calculate concentrations
ConcTableSamples = pd.DataFrame() #pM
for param in MFItableCals.columns:
    ConcTableSamples = pd.concat([ConcTableSamples, pd.DataFrame(InterpFun[param](MFItableSamples[param]))], axis=1, ignore_index=True)
ConcTableSamples.columns = MFItableSamples.columns
ConcTableSamples.index = MFItableSamples.index
#%% some weird plot...
import pylab as pl

fig = pl.figure(num=None, figsize=(5.1, 6.6), dpi=80, facecolor='w', edgecolor='k')  
ax1 = fig.add_subplot(1,1,1)
heatmap = plt.pcolor(np.arcsinh(ConcTableSamples), cmap=plt.cm.inferno)
ax1.set_yticks(np.arange(ConcTableSamples.shape[0]) + 0.5, minor=False)
ax1.set_xticks(np.arange(ConcTableSamples.shape[1]) + 0.5, minor=False)
ax1.set_xticklabels(ConcTableSamples.columns.get_values())
ax1.set_yticklabels(ConcTableSamples.index.get_values())



#%%
import pandas as pd
store = pd.HDFStore(datadir+'ProcessCBA.h5')

store.put('ConcTableSamples',ConcTableSamples,format='t')
store.put('MFItableSamples',MFItableSamples,format='t')
store.put('ConcTableCals',ConcTableCals,format='t')
store.put('MFItableCals',MFItableCals,format='t')

store.close()

#%%
import pandas as pd

store2 = pd.HDFStore(datadir+'ProcessCBA.h5')

ConcTableSamples = store2.get('ConcTableSamples')
MFItableSamples = store2.get('MFItableSamples')
ConcTableCals = store2.get('ConcTableCals')
MFItableCals = store2.get('MFItableCals')