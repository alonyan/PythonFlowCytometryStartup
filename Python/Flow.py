# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 18:53:58 2016

@author: Alonyan
"""
#%%
import FlowCytometryTools


datadir = '/Users/Alonyan/GoogleDrive/Experiments/Th1Th2'

datafile = '/Users/Alonyan/GoogleDrive/Experiments/Th1Th2/Expt_Tube.001.fcs'

import sys
sys.path.append("/Users/Alonyan/Documents/Python/FlowAnalysis")
#%%
from FlowCytometryTools import FCMeasurement

sample = FCMeasurement(ID='Test Sample', datafile=datafile)

print sample.channel_names


print sample.data[['PE-Texas Red-A', 'Pacific Blue-A']][:10]

print sample['PE-Texas Red-A'].describe()


#%%
import os
samples = []
pattern = 'Tube'

for file in os.listdir(datadir):
    if file.endswith(".fcs") and pattern in file:
        print(file)
        samples.append(FCMeasurement(ID=file, datafile=datadir+"/"+ file))
        
        
#%%
        from LoadFlowSamples_v1 import Samples
        datadir = '/Users/Alonyan/GoogleDrive/Experiments/Th1Th2'
        pattern = 'Tube'
        FS = Samples(datadir, pattern)
#%%
tsample = sample.transform('hlog', channels=['PE-Texas Red-A', 'Pacific Blue-A'], b=500.0)
#%%
import pylab as pl
#%%
pl.figure()
tsample.plot('Pacific Blue-A', bins=150)

#%%

pl.figure(num=None, figsize=(5.1, 6.6), dpi=80, facecolor='w', edgecolor='k')

tsample.plot(['PE-Texas Red-A', 'Pacific Blue-A'], cmap=pl.cm.gist_rainbow, bins=150)

#%%
import wx

tsample.view_interactively()

#%%
from FlowCytometryTools import PolyGate

gate1 = PolyGate([(3.864e+03, 4.313e+03), (-2.091e+02, 3.875e+03), (-2.813e+01, 7.716e+03), (5.040e+03, 9.256e+03), (7.778e+03, 6.499e+03), (7.778e+03, 6.499e+03), (7.778e+03, 6.499e+03)], ('PE-Texas Red-A', 'Pacific Blue-A'), region='in', name='gate1')

gated_tsample = tsample.gate(gate1)

pl.figure(num=None, figsize=(5.1, 6.6), dpi=80, facecolor='w', edgecolor='k')


tsample.plot(['PE-Texas Red-A', 'Pacific Blue-A'], kind='scatter', color='gray',  alpha=0.6, label='Original')

gated_tsample.plot(['PE-Texas Red-A', 'Pacific Blue-A'], kind='scatter', color='green',  alpha=0.6, label='gated' ,gates=[gate1])
pl.legend(loc='best')
pl.grid(True)



#%%






#%%
import FlowCytometryTools


datadir = '/Users/Alonyan/GoogleDrive/Experiments/Th1Th2'

datafile = '/Users/Alonyan/GoogleDrive/Experiments/Th1Th2/Expt_Tube.001.fcs'

import sys
sys.path.append("/Users/Alonyan/Documents/Python/FlowAnalysis")

#%%
import LoadFlowSamples_v1 as LFS

from LoadFlowSamples_v1 import Samples
datadir = '/Users/Alonyan/GoogleDrive/Experiments/Th1Th2'
pattern = 'Tube'
FS = Samples(datadir, pattern)


#%%
import wx

GatedSamps[0].view_interactively()
#%%
from FlowCytometryTools import PolyGate

gate1 = PolyGate([(3.864e+03, 4.313e+03), (-2.091e+02, 3.875e+03), (-2.813e+01, 7.716e+03), (5.040e+03, 9.256e+03), (7.778e+03, 6.499e+03), (7.778e+03, 6.499e+03), (7.778e+03, 6.499e+03)], ('PE-Texas Red-A', 'Pacific Blue-A'), region='in', name='gate1')


GatedSamps = LFS.gate(FS.samples,gate1)
