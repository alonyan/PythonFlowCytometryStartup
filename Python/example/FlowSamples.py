# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 15:22:47 2016

@author: Alonyan
"""

from FlowCytometryTools import FCMeasurement
import os
import numpy as np
import pandas as pd


class Samples:
    def __init__(self, datadir, pattern):
        self.datadir = datadir
        self.pattern = pattern
        self.samples = []
        self.samples = self._load()
                
        blankInd = []
        if 'FSC-A' in self.samples[1].channel_names:
            blankInd.append(self.samples[1].channel_names.index('FSC-A'))
        if 'FSC-H' in self.samples[1].channel_names:
            blankInd.append(self.samples[1].channel_names.index('FSC-H'))
        if 'FSC-W' in self.samples[1].channel_names:
            blankInd.append(self.samples[1].channel_names.index('FSC-W'))
        if 'SSC-A' in self.samples[1].channel_names:
            blankInd.append(self.samples[1].channel_names.index('SSC-A'))
        if 'SSC-H' in self.samples[1].channel_names:
            blankInd.append(self.samples[1].channel_names.index('SSC-H'))
        if 'SSC-W' in self.samples[1].channel_names:
            blankInd.append(self.samples[1].channel_names.index('SSC-W'))
        if 'Time' in self.samples[1].channel_names:
            blankInd.append(self.samples[1].channel_names.index('Time'))
            
        ChanInd = list(set(np.arange(0,len(self.samples[1].channel_names)))-set(blankInd))
        ChanList = list(np.array(self.samples[1].channel_names)[ChanInd])
        self.fluorescent_channel_names = ChanList
        self.channel_names = self.samples[1].channel_names
        print ChanList
        self.means = self.samples[0][ChanList].mean()
        self.medians = self.samples[0][ChanList].median()
        self.stds = self.samples[0][ChanList].std()
        self.counts = self.samples[0][ChanList].count()        
        self.IDlist = [self.samples[0].ID[0:-4]]
        for i in np.arange(1,len(self.samples)):
            self.means = pd.concat([self.means, self.samples[i][ChanList].mean()], axis=1, ignore_index=True)
            self.medians = pd.concat([self.medians, self.samples[i][ChanList].median()], axis=1, ignore_index=True)
            self.stds = pd.concat([self.stds, self.samples[i][ChanList].std()], axis=1, ignore_index=True)
            self.counts = pd.concat([self.counts, self.samples[i][ChanList].count()], axis=1, ignore_index=True)

            self.IDlist = self.IDlist + [self.samples[i].ID[0:-4]]
        self.means.columns = self.IDlist
        self.medians.columns = self.IDlist
        self.stds.columns = self.IDlist
        self.counts.columns = self.IDlist
        

        
    def _load(self):
        for file in os.listdir(self.datadir):
            if file.endswith(".fcs") and self.pattern in file:
                print(file)
                self.samples.append(FCMeasurement(ID=file, datafile=self.datadir+"/"+ file))
        return self.samples
    
