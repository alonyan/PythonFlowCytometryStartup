# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 16:35:20 2016

@author: Alonyan
"""

import numpy as np
import pandas as pd

          
def means(samples,ChanList):
        means = samples[0][ChanList].mean()      
        IDlist = [samples[0].ID[0:-4]]
        for i in np.arange(1,len(samples)):
            means = pd.concat([means, samples[i][ChanList].mean()], axis=1, ignore_index=True)
            IDlist = IDlist + [samples[i].ID[0:-4]]
        means.columns = IDlist
        return means
        
def medians(samples,ChanList):
        medians = samples[0][ChanList].median()      
        IDlist = [samples[0].ID[0:-4]]
        for i in np.arange(1,len(samples)):
            medians = pd.concat([medians, samples[i][ChanList].median()], axis=1, ignore_index=True)
            IDlist = IDlist + [samples[i].ID[0:-4]]
        medians.columns = IDlist
        return medians
        
def stds(samples,ChanList):
        stds = samples[0][ChanList].std()      
        IDlist = [samples[0].ID[0:-4]]
        for i in np.arange(1,len(samples)):
            stds = pd.concat([stds, samples[i][ChanList].std()], axis=1, ignore_index=True)
            IDlist = IDlist + [samples[i].ID[0:-4]]
        stds.columns = IDlist
        return stds
        
def counts(samples,ChanList):
        counts = samples[0][ChanList].count()      
        IDlist = [samples[0].ID[0:-4]]
        for i in np.arange(1,len(samples)):
            counts = pd.concat([counts, samples[i][ChanList].count()], axis=1, ignore_index=True)
            IDlist = IDlist + [samples[i].ID[0:-4]]
        counts.columns = IDlist
        return counts
        

          