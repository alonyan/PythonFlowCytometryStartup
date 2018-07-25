# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 10:20:26 2016

@author: Alonyan
"""
#%%
import numpy as np
import sys
sys.path.append("/Users/oyleryaniva/Documents/Python/Modules")
from FlowAnalysis import LoadFlowSamples as LFS
import matplotlib.pyplot as plt
#%%
datadir = '/Users/oyleryaniva/Documents/Python/Comps'
pattern = 'Compensation'
Comps = LFS.FlowData(datadir, pattern)
#Comps = Comps.asinhtform(Comps.fluorescent_channel_names,150)

#%%
Comps.samples[0].view_interactively()
#%% general gates
from FlowCytometryTools import PolyGate

gate1 = PolyGate([(3.030e+04, 6.309e+04), (3.266e+04, 2.436e+04), (4.636e+04, 1.281e+04), (1.531e+05, 1.893e+04), (2.609e+05, 1.147e+05), (2.481e+05, 2.248e+05), (1.834e+05, 2.445e+05), (9.597e+04, 2.105e+05)], ('FSC-A', 'SSC-A'), region='in', name='gate1')
gate2 = PolyGate([(2.102e+04, 7.858e+04), (1.157e+05, 8.865e+04), (1.477e+05, 7.507e+04), (1.414e+05, 6.456e+04), (1.884e+04, 6.368e+04), (1.884e+04, 6.368e+04)], ('SSC-H', 'SSC-W'), region='in', name='gate2')
Comps = Comps.gate(gate1).gate(gate2)
#Comps1 = Comps1.gate(gate1)

#%% Gating individual samples
from FlowCytometryTools import ThresholdGate

gatePB = ThresholdGate(4.024e+00, ('Pacific Blue-A'), region='above')
Comps.samples[3] = Comps.samples[3].gate(gatePB)
#%% calculate compensation matrix
           
IdxMatch = [0, 3, 4, 1] #make diagonal, index of channel in samples
IndBlank = 2
        
meansMatrix = Comps.aSinhmeans().as_matrix()[IdxMatch].T
BGMatrix = np.tile(Comps.aSinhmeans().as_matrix()[IndBlank],(np.shape(meansMatrix)[1],1)).T
Hu = meansMatrix-BGMatrix
        
SpilloverM = Hu/Hu.diagonal()
SpilloverM = np.matrix(SpilloverM)
from numpy.linalg import inv
        
CompMat = inv(SpilloverM)
        
import pandas as pd

SpilloverM = pd.DataFrame(SpilloverM)
SpilloverM.columns = Comps.fluorescent_channel_names
SpilloverM.index = Comps.fluorescent_channel_names
CompMat = pd.DataFrame(CompMat)
CompMat.columns = Comps.fluorescent_channel_names
CompMat.index = Comps.fluorescent_channel_names
        
#%% Applying compensation matrix
Comps1  = Comps.compensate(CompMat)

#%%
        
sample = Comps.samples[1]
sampleT = Comps1.samples[1]
import pylab as pl

fig = pl.figure(num=None, figsize=(5.1, 6.6), dpi=80, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(111)
ax1.set_position([0.15,0.6,0.5,0.4])
#sampleT.plot(['Alexa Fluor 700-A', 'APC-A'], kind='scatter', color='green',  alpha=0.6, label='Comped')

#sample.plot(['Alexa Fluor 700-A', 'APC-A'], kind='scatter', color='gray',  alpha=0.6, label='Comped')
ax1.scatter(sample.data[['Alexa Fluor 700-A']], sample.data[['APC-A']],c='gray',alpha=0.4,linewidth=0, label='NoComp')
ax1.scatter(sampleT.data[['Alexa Fluor 700-A']], sampleT.data[['APC-A']],c='green',alpha=0.4,linewidth=0, label='Comped')

ax1.set_xlabel('Alexa 700')
ax1.set_ylabel('APC')

#ax1.legend(loc='best')
ax1.grid(True)
      #  ax1.set_ylim([-1000,50000])
      #  ax1.set_xlim([-1000,50000])


#%%
def num2str(num, precision): 
    return "%0.*f" % (precision, num)
#%%
sample = Comps.samples[1]
sampleT = Comps1.samples[1]

fig = pl.figure(num=None, figsize=(5, 5), dpi=80, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(111)
ax1.set_position([0.15,0.15,0.6,0.6])

xChannel = 'Alexa Fluor 700-A'
yChannel = 'APC-A'

datax=sample.data[xChannel] #data
datay=sample.data[yChannel]

from Utilities import colorcode
z, idx = colorcode(datax, datay)
ax1.scatter(datax[idx], datay[idx],c=z[idx]**10,alpha=0.3, s=5,linewidth=0, label='NoComp', cmap='viridis')


datax=sampleT.data[xChannel] #data
datay=sampleT.data[yChannel]


z, idx = colorcode(datax, datay)

ax1.scatter(datax[idx], datay[idx],c=z[idx]**10,alpha=0.3, s=5,linewidth=0, label='Comp', cmap='plasma')

ax1.set_xlabel(xChannel)
ax1.set_ylabel(yChannel)
ax1.set_xticks(np.arcsinh(np.array([-1000,-100,-90,-80,-70,-60,-50,-40,-30,-20, -10, -1, 0,1,10,20,30,40,50,60,70,80,90,100, 1000, 10000, 100000])/150.))
ax1.set_xticklabels(['-1','-.1','','','','','','','','','','0','','','','','','','','','','','','.1','1','10','100'])
ax1.set_yticks(np.arcsinh(np.array([-1000,-100,-90,-80,-70,-60,-50,-40,-30,-20, -10, -1, 0,1,10,20,30,40,50,60,70,80,90,100, 1000, 10000, 100000])/150.))
ax1.set_yticklabels(['-1','-.1','','','','','','','','','','0','','','','','','','','','','','','.1','1','10','100'])

#ax1.legend(loc='best')
#ax1.grid(True)
#   ax1.set_ylim([-1000,15000])
ax1.set_xlim([np.arcsinh(-200/150),np.arcsinh(270000/150)])
ax1.set_ylim([np.arcsinh(-200/150),np.arcsinh(270000/150)])   


#fig.savefig('/Users/Alonyan/Desktop/Comps.pdf')

#%%

sample = Comps1.samples[1]

xChannel = 'Pacific Blue-A'
yChannel = 'APC-A'


datax=sample.data[xChannel] #data
datay=sample.data[yChannel]

    
fig = pl.figure(num=None, figsize=(5, 5), dpi=80, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(111)
ax1.set_position([0.15,0.15,0.6,0.6])

from Utilities import colorcode
z, idx = colorcode(datax, datay)
ax1.scatter(datax[idx], datay[idx],c=z[idx]**10,alpha=0.3, s=5,linewidth=0, label='NoComp', cmap='plasma')#scatter colorcoaded 

ax1.set_xlabel(xChannel)
ax1.set_ylabel(yChannel)
ax1.set_xticks(np.arcsinh(np.array([-1000,-100,-90,-80,-70,-60,-50,-40,-30,-20, -10, -1, 0,1,10,20,30,40,50,60,70,80,90,100, 1000, 10000, 100000])/150.))
ax1.set_xticklabels(['-1','-.1','','','','','','','','','','0','','','','','','','','','','','','.1','1','10','100'])
ax1.set_yticks(np.arcsinh(np.array([-1000,-100,-90,-80,-70,-60,-50,-40,-30,-20, -10, -1, 0,1,10,20,30,40,50,60,70,80,90,100, 1000, 10000, 100000])/150.))
ax1.set_yticklabels(['-1','-.1','','','','','','','','','','0','','','','','','','','','','','','.1','1','10','100'])

#ax1.legend(loc='best')
ax1.set_xlim([np.arcsinh(-200/150),np.arcsinh(300000/150)])
ax1.set_ylim([np.arcsinh(-200/150),np.arcsinh(300000/150)])   


#fig.savefig('/Users/Alonyan/Desktop/Comps.pdf')


#%%

