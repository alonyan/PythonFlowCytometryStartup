# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 10:20:26 2016

@author: Alonyan
"""
import numpy as np
import FlowCytometryTools
import sys
sys.path.append("/Users/Alonyan/Documents/Python/FlowAnalysis")
import LoadFlowSamples_v1 as LFS
import flowstats as fs
#%% load samples by pattern
        #from LoadFlowSamples_v1 import Samples
        datadir = '/Users/Alonyan/Documents/Python/Comps'
        pattern = 'Compensation'
        CompensationSamples = np.array(LFS.FlowData(datadir, pattern).samples)
        Channels = np.array(LFS.FlowData(datadir, pattern).fluorescent_channel_names)
        Names = np.array(LFS.FlowData(datadir, pattern).IDlist)
#% hlog transform
        CompensationSamples = LFS.transform(CompensationSamples, 'hlog', Channels, 500)
#%% calculate compensation matrix
           
        IdxMatch = [0, 3, 4, 1] #make diagonal, index of channel in samples
        IndBlank = 2
        
        meansMatrix = fs.means(CompensationSamples[IdxMatch],Channels).as_matrix()
        BGMatrix = np.tile(fs.means(CompensationSamples[[IndBlank]],Channels).as_matrix(),(np.shape(meansMatrix)[1],1)).T
        Hu = meansMatrix-BGMatrix
        
        SpilloverM = Hu/Hu.diagonal()
        SpilloverM = np.matrix(SpilloverM)
        from numpy.linalg import inv
        
        CompMat = inv(SpilloverM)
        
        import pandas as pd

        SpilloverM = pd.DataFrame(SpilloverM)
        SpilloverM.columns = Channels
        SpilloverM.index = Channels
        CompMat = pd.DataFrame(CompMat)
        CompMat.columns = Channels
        CompMat.index = Channels
        
#%% Applying compensation matrix
        NewSamps = LFS.compansate(CompensationSamples, CompMat, Channels)
#%%
        
        sample = CompensationSamples[1]
        sampleT = NewSamps[1]
        import pylab as pl

        fig = pl.figure(num=None, figsize=(5.1, 6.6), dpi=80, facecolor='w', edgecolor='k')
        ax1 = fig.add_subplot(111)
        ax1.set_position([0.15,0.6,0.3,0.22])
        #sampleT.plot(['Alexa Fluor 700-A', 'APC-A'], kind='scatter', color='green',  alpha=0.6, label='Comped')

        #sample.plot(['Alexa Fluor 700-A', 'APC-A'], kind='scatter', color='gray',  alpha=0.6, label='Comped')
        ax1.scatter(sample.data[['Alexa Fluor 700-A']], sample.data[['APC-A']],c='green',alpha=0.6,linewidth=0, label='NoComp')
        ax1.scatter(sampleT.data[['Alexa Fluor 700-A']], sampleT.data[['APC-A']],c='gray',alpha=0.6,linewidth=0, label='Comped')

        ax1.set_xlabel('Alexa 700')
        ax1.set_ylabel('APC')
        
        #ax1.legend(loc='best')
        ax1.grid(True)
        ax1.set_ylim([-1000,15000])
        ax1.set_xlim([-1000,15000])


#%%
        counts, xedges, yedges, Image = plt.hist2d(datax,datay,bins=50)

#%%
        from scipy import interpolate

        fig = pl.figure(num=None, figsize=(5.1, 6.6), dpi=80, facecolor='w', edgecolor='k')
        ax1 = fig.add_subplot(111)
        ax1.set_position([0.15,0.6,0.9,0.9])

        datax=sample.data['Alexa Fluor 700-A'] #data
        datay=sample.data['APC-A']
        H, xedges, yedges = np.histogram2d(datax,datay, bins=10)
        xedges, yedges = np.meshgrid((xedges[:-1]+xedges[1:])/2, (yedges[:-1]+yedges[1:])/2)
        f = interpolate.interp2d(xedges, yedges, H, kind='linear')
        z = f(datax[0],datay[0])
        for i in np.arange(1,len(datax)):
            z = np.append(z,f(datax[i],datay[i]))   
        z=(z-min(z))/(max(z)-min(z))     
        idx = z.argsort()
        ax1.scatter(datax[idx], datay[idx],c=z[idx],alpha=0.6,linewidth=0, label='NoComp')


        datax=sampleT.data['Alexa Fluor 700-A'] #data
        datay=sampleT.data['APC-A']
        H, xedges, yedges = np.histogram2d(datax,datay, bins=10)
        xedges, yedges = np.meshgrid((xedges[:-1]+xedges[1:])/2, (yedges[:-1]+yedges[1:])/2)
        %f1 = interpolate.interp2d(xedges, yedges, H, kind='linear')
        %z = f1(datax[0],datay[0])
        %for i in np.arange(1,len(datax)):
        %    z = np.append(z,f(datax[i],datay[i])) 
        %z=(z-min(z))/(max(z)-min(z))
        %idx = z.argsort()
        ax1.scatter(datax[idx], datay[idx],c=z[idx],alpha=0.6,linewidth=0, label='NoComp')



        ax1.set_xlabel('Alexa 700')
        ax1.set_ylabel('APC')
        ax1.legend(loc='best')
        ax1.grid(True)
        ax1.set_ylim([-1000,15000])
        ax1.set_xlim([-1000,15000])
       




